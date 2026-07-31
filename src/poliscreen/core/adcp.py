"""AutoDock CrankPep (ADCP): acoplamiento específico de péptidos.

Vina trata el ligando como un árbol de torsiones independientes y su muestreo se degrada rápido con
la flexibilidad: sobre saFtsZ, un péptido de 5 residuos (23 rotables) ya no termina en dos minutos y
uno de 10 no converge. ADCP modela el péptido con rotámeros y muestrea conformación y posición a la
vez, que es la única forma razonable de acoplar algo tan flexible.

No sustituye a Vina para moléculas pequeñas, y por debajo de MIN_RESIDUOS su muestreador no aplica
(ver esa constante). Coste sobre 8HTB, octapéptido, seis hilos: 35 s (250.000 pasos × 10 réplicas) y
208 s (1.000.000 × 20), con la energía mejorando de forma monótona (-13,2 → -17,4 → -19,8 kcal/mol),
señal de que el muestreo converge.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess  # noqa: F401  (shutil se usa para copiar la diana junto al trabajo)
from pathlib import Path
from typing import Optional, Sequence

CARPETA = Path.home() / "poliscreen_tools" / "adfrsuite"
MAX_RESIDUOS = 20
# Límite inferior del muestreador, comprobado sobre 8HTB: con 3-4 residuos aborta al construir la
# tabla de movimientos ("Maybe too short chains?") sin escribir nada; desde 5 da poses en segundos.
# El crankshaft necesita esqueleto entre sus dos puntos de giro: es límite del algoritmo, no de la
# instalación. Por debajo de 5, Vina sigue siendo practicable.
MIN_RESIDUOS = 5


def _raiz() -> Optional[Path]:
    """Carpeta de la suite ADFR instalada.

    Se admiten varias instalaciones y se elige la más reciente que tenga el ejecutable: probar
    versiones distintas es habitual (no todas funcionan), y mirar un nombre fijo obligaría a
    exportar una variable en cada sesión.
    """
    env = os.environ.get("POLISCREEN_ADCP")
    if env and Path(env).exists():
        p = Path(env)
        return p if p.is_dir() else p.parent.parent

    cand = []
    base = CARPETA.parent                       # ~/poliscreen_tools
    if base.exists():
        for d in base.glob("adfrsuite*"):
            cand.append(d)
            cand.extend(d.glob("ADFRsuite-*"))
    validos = [d for d in cand if (d / "bin" / "adcp").exists()]
    if validos:
        return max(validos, key=lambda d: d.stat().st_mtime)

    exe = shutil.which("adcp")
    return Path(exe).parent.parent if exe else None


def bin_dir() -> Optional[Path]:
    r = _raiz()
    b = (r / "bin") if r else None
    return b if b and (b / "adcp").exists() else None


def available() -> bool:
    return bin_dir() is not None


def _entorno() -> dict:
    """PATH y LD_LIBRARY_PATH que ADCP necesita.

    Los binarios de la suite requieren libgomp (OpenMP), que no siempre esta en el sistema pero si
    en el entorno de conda; se anade su lib para no depender de una instalación con privilegios.
    """
    env = dict(os.environ)
    b = bin_dir()
    if b:
        env["PATH"] = str(b) + os.pathsep + env.get("PATH", "")
    libs = [p for p in (Path(os.environ.get("CONDA_PREFIX", "")) / "lib",) if p.exists()]
    if libs:
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            [str(p) for p in libs] + ([env["LD_LIBRARY_PATH"]] if env.get("LD_LIBRARY_PATH") else []))
    return env


def _archivo_rama() -> Optional[Path]:
    """Tabla de probabilidades de Ramachandran que el motor de muestreo espera en su directorio
    de trabajo. Viene con la suite, pero el binario no la busca en su propia instalación."""
    r = _raiz()
    if not r:
        return None
    cand = r / "CCSBpckgs" / "ADCP" / "ramaprob.data"
    if cand.exists():
        return cand
    return next(iter(r.rglob("ramaprob.data")), None)


class AdcpError(RuntimeError):
    pass


def prepare_target(receptor_pdb, box, out_dir, nombre: str = "diana", timeout: int = 900) -> Optional[Path]:
    """Prepara la diana con AGFR y devuelve el archivo .trg que consume ADCP.

    ADCP no acepta un PDB directamente: necesita los mapas de afinidad precalculados en una caja,
    que es lo que genera AGFR. Se reutiliza la misma caja de busqueda que el resto de PoliScreen,
    de modo que ambos motores exploran exactamente la misma region.
    """
    b = bin_dir()
    if b is None:
        raise AdcpError("ADCP no esta instalado. Ejecuta scripts/get_adcp.sh.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    env = _entorno()

    # AGFR parte de un pdbqt de receptor; se genera con la utilidad de la propia suite para que la
    # asignación de tipos de átomo sea la que ADCP espera.
    recq = out_dir / f"{nombre}.pdbqt"
    if not recq.exists():
        r = subprocess.run([str(b / "prepare_receptor"), "-r", str(receptor_pdb), "-o", str(recq)],
                           capture_output=True, text=True, env=env, timeout=timeout)
        if not recq.exists():
            raise AdcpError(f"prepare_receptor fallo: {(r.stderr or r.stdout or '')[-300:]}")

    destino = out_dir / nombre
    cmd = [str(b / "agfr"), "-r", str(recq), "-o", str(destino),
           "-b", f"user {box.cx} {box.cy} {box.cz} {box.sx} {box.sy} {box.sz}"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    trg = next(iter(sorted(out_dir.glob(f"{nombre}*.trg"))), None)
    if trg is None:
        raise AdcpError(f"AGFR no genero la diana: {(r.stderr or r.stdout or '')[-300:]}")
    return trg


# Tabla final de ADCP: rango, energía, RMSD y más columnas. Se exigen tres numeros para no confundir
# las líneas sueltas del preambulo del muestreador, que también empiezan por un entero y un decimal.
_RE_ENERGIA = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s", re.M)
_RE_RANK = re.compile(r"_ranked_(\d+)\.pdb$")


def dock_peptide(secuencia: str, target_trg, out_dir, nombre: Optional[str] = None,
                 n_pasos: int = 250000, n_repeticiones: int = 10, semilla: int = 42,
                 nucleos: int = 1, helice: bool = False, ciclico: bool = False,
                 cistina: bool = False, timeout: int = 3600) -> dict:
    """Acopla un péptido a partir de su SECUENCIA. Devuelve {'energía', 'poses', 'salida'}.

    ADCP genera la conformación durante el acoplamiento, así que no hace falta pasarle una
    estructura 3D previa. En la secuencia, MAYUSCULA indica que ese residuo parte de conformación
    helicoidal y minuscula de ovillo; útil cuando se sabe que el péptido es helicoidal.

    ciclico: cierre cabeza-cola por el esqueleto. Lo aplica el propio muestreador, que restringe la
    conformación al ciclo en vez de acoplar el péptido lineal; pasarlo importa porque un ciclo
    tiene muchos menos grados de libertad y un modo de unión distinto al de su análogo abierto.
    cistina: cierre por puente disulfuro entre dos cisteinas.
    """
    b = bin_dir()
    if b is None:
        raise AdcpError("ADCP no esta instalado.")
    seq = "".join(ch for ch in secuencia if ch.isalpha())
    if not MIN_RESIDUOS <= len(seq) <= MAX_RESIDUOS:
        raise AdcpError(f"ADCP acopla de {MIN_RESIDUOS} a {MAX_RESIDUOS} residuos; se pidieron "
                        f"{len(seq)}. Por debajo de {MIN_RESIDUOS} usa Vina.")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    nombre = nombre or seq

    # ADCP compone rutas relativas a su directorio de trabajo (./tmp_<trabajo>/), así que se ejecuta
    # DENTRO de out_dir con nombres sin ruta; con rutas absolutas no halla sus intermedios y aborta.
    trg = Path(target_trg)
    local_trg = out_dir / trg.name
    if not local_trg.exists():
        shutil.copy(trg, local_trg)

    # El muestreador lee ramaprob.data (Ramachandran) del directorio de trabajo. El lanzador de la
    # suite intenta copiarlo, pero se repone aquí por si no ocurre: sin él el binario aborta.
    rama = _archivo_rama()
    if rama and not (out_dir / rama.name).exists():
        shutil.copy(rama, out_dir / rama.name)

    cmd = [str(b / "adcp"), "-t", trg.name,
           "-s", seq.upper() if helice else seq.lower(),
           "-N", str(n_repeticiones), "-n", str(n_pasos),
           "-o", nombre, "-c", str(nucleos), "-S", str(semilla)]
    if ciclico:
        cmd.append("-cyc")
    if cistina:
        cmd.append("-cys")
    r = subprocess.run(cmd, capture_output=True, text=True, env=_entorno(),
                       cwd=str(out_dir), timeout=timeout)
    salida = (r.stdout or "") + (r.stderr or "")

    # ADCP escribe dos cosas: <nombre>_<i>.pdb, trayectorias completas de cada réplica (cientos de
    # modelos), y <nombre>_ranked_<k>.pdb, las poses finales ya ordenadas. Solo estas son resultados.
    poses = sorted((p for p in out_dir.glob(f"{nombre}_ranked_*.pdb") if _tiene_atomos(p)),
                   key=lambda p: int(_RE_RANK.search(p.name).group(1)))
    energias = [float(m.group(2)) for m in _RE_ENERGIA.finditer(salida)]
    energias.sort()
    return {"secuencia": seq, "energia": energias[0] if energias else None,
            "energias": energias[:len(poses)] or None,
            "poses": [str(p) for p in poses], "salida": salida[-1500:],
            "ok": bool(poses)}


def _tiene_atomos(pdb: Path) -> bool:
    try:
        return any(l.startswith(("ATOM", "HETATM")) for l in pdb.read_text(errors="ignore").splitlines())
    except Exception:
        return False


def diagnostico(receptor_pdb, box, timeout: int = 900) -> tuple:
    """(funciona, mensaje). Comprueba de extremo a extremo que ADCP produce poses utilizables.

    Se prueba con un pentapeptido, no con uno más corto: por debajo de MIN_RESIDUOS el muestreador
    aborta siempre, de modo que una prueba más corta diagnosticaría como averiada una instalación
    correcta.
    """
    import tempfile
    if not available():
        return False, "ADCP no esta instalado. Ejecuta scripts/get_adcp.sh."
    tmp = Path(tempfile.mkdtemp())
    try:
        trg = prepare_target(receptor_pdb, box, tmp, nombre="diag")
    except Exception as e:
        return False, f"AGFR no pudo preparar la diana: {str(e)[:160]}"
    try:
        r = dock_peptide("KWKLF", trg, tmp, nombre="diag_pep", n_pasos=50000, n_repeticiones=2,
                         timeout=timeout)
    except Exception as e:
        return False, f"ADCP no pudo acoplar: {str(e)[:160]}"
    if r["ok"]:
        return True, (f"ADCP operativo: {len(r['poses'])} poses, mejor energia "
                      f"{r['energia']:.1f} kcal/mol.")
    if "too short chains" in r["salida"]:
        return False, (f"El muestreador rechazo la secuencia de prueba por corta; ADCP necesita al "
                       f"menos {MIN_RESIDUOS} residuos.")
    return False, "ADCP no genero poses con coordenadas."


def dock_batch(secuencias: Sequence[str], target_trg, out_dir, on_progress=None, **kw) -> list:
    """Acopla varias secuencias contra la misma diana. Secuencial: ADCP ya paraleliza internamente."""
    filas = []
    for i, seq in enumerate(secuencias, 1):
        try:
            filas.append(dock_peptide(seq, target_trg, out_dir, **kw))
        except Exception as e:
            filas.append({"secuencia": seq, "energia": None, "poses": [],
                          "salida": str(e)[:300], "ok": False})
        if on_progress:
            on_progress(i, len(secuencias), seq)
    return filas


# Un enlace amida mide 1,33 A. Se admite holgura para la imprecision del muestreo, pero no tanta
# como para aceptar dos extremos meramente proximos.
CIERRE_MAX = 1.8


def _cerrar_anillo(pdb) -> Optional[float]:
    """Distancia C(último)-N(primero); si es la de un enlace, añade un CONECT.

    ADCP no escribe el enlace de cierre, así que el visor —que deduce enlaces por distancia— dibuja
    el péptido abierto aunque el anillo esté formado; el CONECT lo muestra como el ciclo que es. Solo
    se añade si la distancia lo justifica: forzarlo entre átomos separados dibujaría una barra larga
    y engañaría sobre la geometría real.
    """
    pdb = Path(pdb)
    n_serial = c_serial = None
    n_xyz = c_xyz = None
    n_res = []
    try:
        lineas = pdb.read_text(errors="ignore").splitlines()
    except Exception:
        return None
    for l in lineas:
        if l.startswith("ATOM"):
            n_res.append(int(l[22:26]))
    if not n_res:
        return None
    a, b = min(n_res), max(n_res)
    for l in lineas:
        if not l.startswith("ATOM"):
            continue
        nombre, resi = l[12:16].strip(), int(l[22:26])
        try:
            serial = int(l[6:11])
        except ValueError:
            continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if resi == a and nombre == "N":
            n_serial, n_xyz = serial, xyz
        elif resi == b and nombre == "C":
            c_serial, c_xyz = serial, xyz
    if not (n_xyz and c_xyz):
        return None
    d = sum((x - y) ** 2 for x, y in zip(n_xyz, c_xyz)) ** 0.5
    if d <= CIERRE_MAX and n_serial and c_serial:
        cuerpo = [l for l in lineas if not l.strip().startswith(("END", "CONECT"))]
        cuerpo.append(f"CONECT{c_serial:>5}{n_serial:>5}")
        cuerpo.append(f"CONECT{n_serial:>5}{c_serial:>5}")
        cuerpo.append("END")
        pdb.write_text("\n".join(cuerpo) + "\n")
    return d


def dock_sitios(targets, peptidos: dict, out_dir, receptor_por_sitio: dict,
                n_poses: Optional[int] = None, on_progress=None, **kw) -> tuple:
    """Acopla péptidos con ADCP; devuelve (filas, errores) con la forma de docking.dock_batch.

    Nombra las poses con el convenio del resto de la aplicación para que la fusión de complejos, PLIP
    y la puntuación las traten como las de cualquier motor. La energía de ADCP no es comparable con
    la de Vina, así que cada fila declara su motor.

    targets: [(ruta_receptor, id_sitio, caja)]. péptidos: {nombre: (secuencia, ciclico)}.
    """
    out_dir = Path(out_dir)
    poses_dir = out_dir / "poses"
    poses_dir.mkdir(parents=True, exist_ok=True)
    trabajo = out_dir / "adcp"
    trabajo.mkdir(parents=True, exist_ok=True)

    filas, errores = [], []
    total = len(targets) * max(1, len(peptidos))
    hecho = 0
    for rpath, sid, box in targets:
        try:
            trg = prepare_target(receptor_por_sitio.get(sid, rpath), box, trabajo / sid, nombre=sid)
        except Exception as e:
            errores.append((sid, f"AGFR no pudo preparar la diana: {str(e)[:160]}"))
            continue
        for nombre, (seq, ciclico) in peptidos.items():
            hecho += 1
            base = f"docking_{sid}_compounds_a_{nombre}"
            try:
                r = dock_peptide(seq, trg, trabajo / sid / nombre, nombre=nombre,
                                 ciclico=ciclico, **kw)
            except Exception as e:
                errores.append((base, str(e)[:200]))
                if on_progress:
                    on_progress(hecho, total, base, str(e)[:200])
                continue
            if not r["ok"]:
                errores.append((base, "ADCP no genero poses"))
                if on_progress:
                    on_progress(hecho, total, base, "sin poses")
                continue

            # El número de poses que se reportan es una decision del usuario y es independiente del
            # esfuerzo de muestreo: en ADCP, -N son las REPLICAS de la busqueda, y recortarlas para
            # obtener menos poses empeoraria el resultado en vez de solo acortar la lista. Se
            # muestrea con las replicas configuradas y se entregan las n mejores ya ordenadas.
            salida_poses = r["poses"][:n_poses] if n_poses else r["poses"]
            energias = (r.get("energias") or [r["energia"]] * len(r["poses"]))[:len(salida_poses)]
            aviso = None
            abiertas = 0
            for i, (p, e) in enumerate(zip(salida_poses, energias), 1):
                destino = poses_dir / f"{base}-model{i}.pdb"
                destino.write_bytes(Path(p).read_bytes())
                if ciclico:
                    # El enlace cabeza-cola se añade solo si la geometría lo justifica: con cierre
                    # real el visor dibuja el anillo; con los extremos separados no se falsea un
                    # enlace inexistente y la pose se cuenta para avisar. ADCP restringe al ciclo
                    # pero no siempre cierra los anillos pequeños.
                    d = _cerrar_anillo(destino)
                    if d is None or d > CIERRE_MAX:
                        abiertas += 1
                    if i == 1 and d is not None and d > CIERRE_MAX:
                        aviso = (f"se pidio ciclado y la pose sale con los extremos a {d:.2f} A "
                                 f"(un enlace amida mide 1,33 A): trata la geometria con reserva")
                filas.append({"receptor": sid, "pose_name": f"{base}-model{i}",
                              "compound_name": nombre, "docking_score": e, "motor": "adcp"})
            if ciclico and abiertas == len(salida_poses):
                aviso = (f"ADCP no cerro el anillo en ninguna de las {abiertas} poses; los ciclos de "
                         f"pocos residuos quedan tensos y sus extremos no llegan a enlazarse. "
                         f"La estructura para descriptores si es ciclica; la pose acoplada no.")
            if aviso:
                errores.append((base, aviso))
            if on_progress:
                on_progress(hecho, total, base, None)
    return filas, errores
