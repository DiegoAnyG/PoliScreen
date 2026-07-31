"""Pruebas del nucleo: química, geometría y enrutado. No requieren binarios externos.

Cada prueba corresponde a un fallo que se detecto en uso real; el comentario dice cual, porque una
prueba sin ese contexto acaba borrandose cuando estorba.
"""
from pathlib import Path

import pytest

from poliscreen.core import docking as dk
from poliscreen.core import peptides as pp
from poliscreen.core import session as ss


# --------------------------------------------------------------- rutas escritas a mano
@pytest.mark.parametrize("entrada,esperado", [
    (r"\\wsl.localhost\Ubuntu-24.04\home\u\proy", "/home/u/proy"),
    (r"\\wsl$\Ubuntu\home\u\proy", "/home/u/proy"),
    ("/home/u/proy", "/home/u/proy"),
    ('"/home/u/proy"', "/home/u/proy"),
])
def test_rutas_de_windows_se_traducen(entrada, esperado):
    """Una ruta de Windows creaba UNA carpeta con barras invertidas en el nombre."""
    assert str(ss.normalizar_ruta(entrada)[0]) == esperado


def test_ruta_vacia_da_un_destino_valido():
    p, _ = ss.normalizar_ruta("")
    assert p.is_absolute()


# --------------------------------------------------------------- química de péptidos
def test_ciclacion_es_cabeza_cola_con_cadenas_laterales_reactivas():
    """Buscar el COOH por su forma cerraba por el aspartato o el glutamato."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors as rmd
    RDLogger.DisableLog("rdApp.*")
    for seq in ("SDCEFGQ", "ELQGRAK", "KWKLF", "ICIWDDS"):
        lin = Chem.MolFromSmiles(pp.to_smiles(seq))
        cic = Chem.MolFromSmiles(pp.to_smiles(seq, ciclico=True))
        macro = max((len(r) for r in cic.GetRingInfo().AtomRings() if len(r) > 8), default=0)
        assert macro == 3 * len(seq), f"{seq}: anillo de {macro}, esperado {3 * len(seq)}"
        perdida = Descriptors.MolWt(lin) - Descriptors.MolWt(cic)
        assert abs(perdida - 18.02) < 0.1, f"{seq}: perdio {perdida:.2f}, esperado una agua"
        assert rmd.CalcNumRings(cic) == rmd.CalcNumRings(lin) + 1


def test_los_extremos_se_modifican_solos():
    """Acetilar tocaba la lisina y amidar el aspartato."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    AMINA = Chem.MolFromSmarts("[NX3;H2;!$(N-C=O)]")
    ACIDO = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    def cuenta(smi, patron):
        return len(Chem.MolFromSmiles(smi).GetSubstructMatches(patron))

    assert cuenta(pp.to_smiles("KKKGG"), AMINA) - cuenta(pp.to_smiles("KKKGG", n_acetil=True), AMINA) == 1
    assert cuenta(pp.to_smiles("DDEGG"), ACIDO) - cuenta(pp.to_smiles("DDEGG", c_amida=True), ACIDO) == 1


def test_la_carga_refleja_la_quimica_de_los_extremos():
    """La tabla ignoraba acetilación y ciclación al calcular la carga."""
    base = pp.propiedades("KWKLF")["carga_neta"]
    assert pp.propiedades("KWKLF", n_acetil=True)["carga_neta"] == base - 1
    assert pp.propiedades("KWKLF", c_amida=True)["carga_neta"] == base + 1
    assert pp.propiedades("KWKLF", ciclico=True)["carga_neta"] == base


def test_los_nombres_distinguen_las_variantes():
    """Dos moléculas distintas compartian nombre en las tablas y en los archivos de poses."""
    nombres = {pp.etiqueta("KWKLF", **kw) for kw in
               ({}, {"n_acetil": True}, {"c_amida": True},
                {"n_acetil": True, "c_amida": True}, {"ciclico": True})}
    assert len(nombres) == 5


def test_smiles_falla_en_vez_de_devolver_el_lineal():
    """Devolver el péptido abierto con nombre de ciclico ocultaba el error."""
    assert pp.to_smiles("") is None


# --------------------------------------------------------------- reparto de recursos
def test_el_paralelismo_baja_con_ligandos_flexibles(tmp_path):
    """Repartir solo por nucleos agoto la memoria y el sistema mato el proceso."""
    caja = dk.Box(0, 0, 0, 24, 24, 24)
    assert dk.coste_memoria_gb(caja, 25) > dk.coste_memoria_gb(caja, 5)
    grande = dk.Box(0, 0, 0, 60, 60, 60)
    assert dk.coste_memoria_gb(grande, 5) > dk.coste_memoria_gb(caja, 5)
    assert dk.paralelismo_seguro([caja], []) >= 1


def test_torsdof_se_lee_del_pdbqt(tmp_path):
    f = tmp_path / "l.pdbqt"
    f.write_text("ATOM      1  C   LIG A   1       0.0   0.0   0.0\nTORSDOF 22\n")
    assert dk.torsdof(f) == 22
    assert dk.torsdof(tmp_path / "no_existe.pdbqt") == 0


# --------------------------------------------------------------- reconocer péptidos
def test_una_molecula_pequena_no_se_toma_por_peptido(tmp_path):
    """El reconocimiento decide el motor: un falso positivo mandaria un ligando al motor erroneo."""
    f = tmp_path / "lig.pdb"
    f.write_text("HETATM    1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\nEND\n")
    assert pp.secuencia_de_estructura(f) is None
    assert pp.secuencia_de_estructura(tmp_path / "x.sdf") is None


def test_se_reconoce_un_peptido_con_esqueleto_completo(tmp_path):
    f = tmp_path / "pep.pdb"
    filas = []
    n = 1
    for i, (res, ) in enumerate([("ALA",), ("GLY",), ("LYS",), ("TRP",), ("PHE",)], start=1):
        for at in ("N", "CA", "C"):
            filas.append("ATOM  %5d  %-3s %3s B%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                         % (n, at, res, i, i * 3.0, 0.0, 0.0))
            n += 1
    f.write_text("\n".join(filas) + "\nEND\n")
    seq, ciclo = pp.secuencia_de_estructura(f)
    assert seq == "AGKWF"
    assert ciclo is False


# --------------------------------------------------------------- exportación
def test_el_paquete_se_arma_en_memoria(tmp_path):
    """La exportación dejaba copias sueltas en la carpeta de resultados."""
    (tmp_path / "ranking.csv").write_text("compound,best_dock\na,-8.0\n")
    datos, incluidos = ss.paquete(tmp_path, ["resultados_csv"], methods_text="# m")
    assert datos[:2] == b"PK"
    assert any("ranking.csv" in i for i in incluidos)
    assert not list(tmp_path.glob("export_*"))


def test_el_catalogo_marca_lo_recomendado(tmp_path):
    (tmp_path / "ranking.csv").write_text("a,b\n1,2\n")
    cat = ss.catalogo(tmp_path)
    assert cat["resultados_csv"]["hay"] and cat["resultados_csv"]["motivo"]
    assert not cat["poses_zip"]["hay"]


# --------------------------------------------------------------- residuos modificados
def test_no_se_duplican_atomos_al_preparar(tmp_path):
    """Un residuo modificado pedido como cofactor se anadia sobre el que ya estaba en la cadena."""
    from poliscreen.core import receptor as rc
    f = tmp_path / "r.pdb"
    f.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
        "END\n")
    assert rc._quitar_superpuestos(f) == 1
    assert rc._quitar_superpuestos(f) == 0        # idempotente


def test_conect_solo_cuando_el_anillo_cierra(tmp_path):
    """El visor dibujaba abierto un péptido ciclico; se cierra solo si la geometría lo justifica."""
    from poliscreen.core import adcp

    def pose(dist):
        return (f"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                f"ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C\n"
                f"ATOM      3  C   ALA A   2       {dist:.3f}   0.000   0.000  1.00  0.00           C\n"
                f"END\n")
    cerrada = tmp_path / "cerrada.pdb"
    cerrada.write_text(pose(1.5))                 # N(res1)-C(res2) a 1.5 A
    d = adcp._cerrar_anillo(cerrada)
    assert d is not None and "CONECT" in cerrada.read_text()

    abierta = tmp_path / "abierta.pdb"
    abierta.write_text(pose(4.0))                 # a 4 A: no se enlaza
    adcp._cerrar_anillo(abierta)
    assert "CONECT" not in abierta.read_text()


def test_se_explica_por_que_no_hay_rmsd(tmp_path):
    """'no calculable' no decia nada; la causa habitual es comparar moléculas distintas."""
    from poliscreen.core import validation as vl
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    a.write_text("ATOM      1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\n"
                 "ATOM      2  P1  LIG A   1       1.0   0.0   0.0  1.00  0.00           P\n")
    b.write_text("ATOM      1  C1  LIG A   1       0.0   0.0   0.0  1.00  0.00           C\n")
    assert "no son la misma molécula" in vl._motivo_sin_rmsd(a, b)
    assert vl._motivo_sin_rmsd(a, a) == "no calculable"


def test_un_peptido_largo_llega_a_estructura_3d(tmp_path):
    """Dos tridecapeptidos de una biblioteca de diez desaparecian sin mensaje: fallaba el encaje."""
    from poliscreen.core import ligands as lg
    seqs = ["DHITYAVHVQIRW", "WMHSPRFKIVVKW"]
    smis = [pp.to_smiles(s, c_amida=True) for s in seqs]
    assert all(smis)
    made = lg.materialize(smis, tmp_path, names=seqs)
    assert len(made) == len(seqs), "un ligando se perdio al generar la estructura 3D"


def test_se_prefiere_el_envoltorio_de_gnina(tmp_path, monkeypatch):
    """El binario suelto de gnina no es autocontenido: elegirlo obligaba a exportar una variable."""
    herramientas = tmp_path / "poliscreen_tools"
    herramientas.mkdir()
    for nombre in ("gnina", "gnina-run"):
        f = herramientas / nombre
        f.write_text("#!/bin/sh\n")
        f.chmod(0o755)
    monkeypatch.delenv("POLISCREEN_GNINA", raising=False)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    assert dk.gnina_exe().endswith("gnina-run")


def test_los_modificados_se_detectan_sin_listas_propias(tmp_path):
    """La detección usa la tabla de PDBFixer; sin estructura valida devuelve lista vacia, no falla."""
    from poliscreen.core import receptor as rc
    f = tmp_path / "vacio.pdb"
    f.write_text("END\n")
    assert rc.modified_residues(f) == []
