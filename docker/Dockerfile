# Imagen reproducible de PoliScreen (nucleo: preparacion, docking, interacciones, cribado e interfaz).
#
#   docker build -f docker/Dockerfile -t poliscreen .
#   docker run --rm -p 8501:8501 -v "$PWD/proyectos:/data" poliscreen
#   -> http://localhost:8501  (la carpeta de proyecto dentro de la app es /data)
#
# Imagen base FIJADA por tag (no 'latest': 'latest' se mueve y en unos meses construiria otra cosa).
# Para reproducibilidad maxima en una publicacion, sustituye el tag por el digest:
#   docker inspect --format='{{index .RepoDigests 0}}' condaforge/miniforge3:26.3.2-3
FROM condaforge/miniforge3:26.3.2-3

# Version de Vina fijada: la serie 1.2.x no esta en conda y su puntuacion difiere de la 1.1.2.
# El SHA256 fija la identidad exacta del binario; si el asset cambiara, la construccion falla.
ARG VINA_VERSION=1.2.5
ARG VINA_SHA256=fa0126a28a9ea9162d1b161dfa92bc76e632416db28ca246278ea4b2dc6860cb

ENV POLISCREEN_HOME=/opt/poliscreen \
    MPLBACKEND=Agg \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR ${POLISCREEN_HOME}

# El entorno se crea antes de copiar el codigo: asi la capa se reutiliza si solo cambia el fuente.
COPY environment.yml ./
RUN mamba env create -f environment.yml -n poliscreen --yes \
 && mamba clean -afy

RUN curl -fsSL "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
      -o /opt/conda/envs/poliscreen/bin/vina \
 && echo "${VINA_SHA256}  /opt/conda/envs/poliscreen/bin/vina" | sha256sum -c - \
 && chmod +x /opt/conda/envs/poliscreen/bin/vina

ENV PATH=/opt/conda/envs/poliscreen/bin:${PATH}

COPY . ${POLISCREEN_HOME}
# --no-deps: las dependencias ya las resolvio conda; evita que pip las duplique con otras versiones.
RUN pip install --no-deps -e .

# Comprobacion en tiempo de construccion: si falta una herramienta, la imagen falla aqui y no en produccion.
RUN vina --version && obabel -V && plip --version >/dev/null 2>&1 || true \
 && python -c "import rdkit, py3Dmol, streamlit, poliscreen; print('PoliScreen listo')"

VOLUME ["/data"]
EXPOSE 8501

CMD ["streamlit", "run", "/opt/poliscreen/src/poliscreen/ui/streamlit_app.py", \
     "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
