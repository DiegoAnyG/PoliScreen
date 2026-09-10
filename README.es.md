![PoliScreen](docs/banner.png)

# PoliScreen (Español)

[Read in English](README.md)

Cribado virtual reproducible que cierra el ciclo **diseño → filtro de sintetizabilidad
→ acoplamiento molecular (docking) → evaluación de calidad de interacción → ADMET**,
con una función de puntuación objetiva por cavidad, métrica de confianza ortogonal y
optimización multi-objetivo mediante el Frente de Pareto.

> **v1.1.0** — Paisaje interactivo de Pareto, polígonos geométricos de huella de
> interacción, túneles de transporte (CAVER / CaverDock) y distribución reproducible en contenedor.

---

## Inicio Rápido (Quick Start)

La vía de distribución recomendada es **Docker** (o Linux/WSL). Debido a que el docking
y el cálculo de interacciones involucran operaciones de punto flotante en múltiples
librerías de C/Fortran, ejecutarlo dentro del contenedor oficial garantiza que los
resultados numéricos sean exactamente idénticos en cualquier computadora.

### Windows (Sin necesidad de terminal)

1. Instala e inicia [Docker Desktop](https://www.docker.com/products/docker-desktop).
2. Descarga los archivos `PoliScreen-Docker.bat` y `PoliScreen.ico` desde la [sección de Releases](https://github.com/DiegoAnyG/PoliScreen/releases/latest).
3. Haz doble clic en `PoliScreen-Docker.bat`. El lanzador descargará la imagen oficial,
   creará automáticamente un acceso directo en tu escritorio con el logotipo de PoliScreen
   y abrirá la aplicación en tu navegador (`http://localhost:8501`).

### Linux / WSL2 (Un solo comando)

```bash
docker run --rm -it --init -p 127.0.0.1:8501:8501 \
  -v "$PWD/proyectos:/data" ghcr.io/diegoanyg/poliscreen:latest
# Luego abre http://localhost:8501 en tu navegador
```

Para instalación nativa en entornos de desarrollo Conda: **[docs/INSTALL.md](docs/INSTALL.md)**.
Requisitos del sistema: **[docs/REQUIREMENTS.md](docs/REQUIREMENTS.md)**.

---

## Flujo de Trabajo en 4 Pasos

PoliScreen cuenta con una interfaz web modular disponible en **Español e Inglés**
(Configuración → Idioma):

1. **Receptores** — Descarga una estructura del PDB (ej. `8HTB`) o sube tu propio archivo `.pdb`.
   Permite inspeccionar cadenas, conservar cofactores y extraer el ligando cristalográfico como control.
   *(Consejo: haz clic en el botón **Ejemplo (8HTB)** para cargar el caso de estudio de inmediato).*
2. **Ligandos** — Sube tus moléculas preparadas (`.sdf`, `.mol2`, `.smi`), construye análogos mediante
   reacciones químicas (con filtro de viabilidad sintética), filtra fármacos aprobados de ChEMBL
   o genera péptidos.
   *(Consejo: haz clic en **Ligandos de ejemplo (8HTB)** para cargar 3 inhibidores activos conocidos y 3 señuelos).*
3. **Ejecutar (Run)** — Ajusta la caja de búsqueda alrededor de la cavidad e inicia el acoplamiento
   (AutoDock Vina para moléculas pequeñas, ADCP para péptidos) y opcionalmente el cálculo de túneles (CAVER).
4. **Resultados** — Analiza la clasificación a través de herramientas avanzadas:
   - **Frente de Pareto Interactivo**: Paisaje multi-objetivo que equilibra afinidad,
     calidad de interacción y confianza, mostrando la estructura química 2D al posar el cursor sobre cada punto.
   - **Leaderboard TOP de Interacciones**: Polígonos geométricos donde los vértices y aristas
     representan los residuos clave de contacto y el tipo de enlace molecular.
   - **Perfil ADMET**: Estimación de propiedades farmacocinéticas y diagramas de radar.
   - **Túneles de Transporte**: Barrera de activación ($E_a$) y perfil energético hacia el sitio activo.

---

## Principios Científicos de PoliScreen

La mayoría de programas de cribado clasifican compuestos basándose únicamente en la *afinidad*
y el *número de contactos*. Esto premia moléculas promiscuas: una molécula grande que toca muchos
residuos irrelevantes puede superar a un inhibidor que interactúa con precisión catalítica.
PoliScreen se fundamenta en cuatro pilares:

1. **Puntuación objetiva por cavidad, no por similitud al control:** La calidad de interacción
   suma cada contacto ponderado por su **tipo de enlace** (puente salino > puente de hidrógeno >
   interacción pi > hidrofóbico) y por el **rol del residuo** (catalítico, secundario, cavidad o externo).
   Se normaliza contra la huella cristalográfica del ligando de referencia en su pose real.

2. **Métrica de confianza ortogonal al puntaje:** La `confianza` (0 a 1) es la media geométrica
   de la convergencia de poses, concordancia entre afinidad e interacciones, y consenso Vina-red neuronal
   (si se activa gnina). Mide **qué tan fiable es el resultado**, no su magnitud.

3. **Sintetizabilidad desde el diseño:** Los análogos se filtran por viabilidad química real
   (regioselectividad, clasificación de grupos hidroxilo, impedimento estérico), asegurando que lo
   que se acopla sea sintetizable en el laboratorio.

4. **Reproducibilidad explícita:** Semilla fija, un hilo por tarea de acoplamiento y exportación
   de un archivo de Métodos detallado con cada parámetro y versión. Dentro del contenedor Docker,
   dos corridas con la misma configuración generan exactamente los mismos resultados numéricos.

---

## Motores Computacionales Integrados

- **AutoDock Vina 1.2.5**: Motor de acoplamiento para moléculas pequeñas.
- **AutoDock CrankPep (ADCP)**: Motor de acoplamiento para péptidos (5 a 20 residuos).
- **PLIP**: Perfilador de interacciones macromoleculares (puentes de H, puentes salinos, pi-stacking, contactos hidrofóbicos).
- **fpocket**: Detección geométrica de cavidades y bolsillos de unión.
- **CAVER 3.0.2**: Detección de túneles de transporte (incluido en la imagen Docker).
- **CaverDock**: Simulación de trayectorias de transporte y cálculo de barrera energética ($E_a$).
- **gnina**: Re-puntuación opcional basada en redes neuronales profundas (acelerado por GPU).
- **ADMET-AI**: Predicción de propiedades fisicoquímicas y perfiles farmacocinéticos.

---

## Limitaciones

- **Receptor rígido**: No modela flexibilidad en cadenas laterales durante el acoplamiento estándar.
- **Acoplamiento no covalente**: Vina y motores compatibles no modelan formación de enlaces covalentes.
- **Predicciones ADMET estimadas**: Diseñadas para priorización de compuestos, no sustituyen pruebas experimentales.

---

## Cita

Consulta **[CITATION.cff](CITATION.cff)** para citar PoliScreen en publicaciones científicas.
Recuerda citar también las herramientas subyacentes que PoliScreen invoca; la bibliografía completa
se encuentra en el panel *Cómo citar* de la aplicación y en el archivo de Métodos exportado.

## Licencia

**GNU GPL v3 o posterior** (ver [LICENSE](LICENSE)). Los programas científicos que PoliScreen
ejecuta (Vina, ADCP, PLIP, RDKit, Open Babel, fpocket, gnina, CAVER) son software independiente
que conservan sus propias licencias de distribución.

## Autor

**Diego Cesar Anaya Guerrero**

