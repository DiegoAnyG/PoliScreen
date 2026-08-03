# Reproducible PoliScreen image.
#
#   docker build -f docker/Dockerfile -t poliscreen .
#   docker run --rm -p 8501:8501 -v "$PWD/proyectos:/data" poliscreen
#   -> http://localhost:8501  (projects are created under /data, the mounted volume)
#
# Base image pinned by tag: 'latest' moves and would build something else in a few months.
# For a publication, replace the tag with the digest:
#   docker inspect --format='{{index .RepoDigests 0}}' condaforge/miniforge3:26.3.2-3
FROM condaforge/miniforge3:26.3.2-3

# Vina pinned: 1.2.x is not on conda and its scoring function differs from 1.1.2.
# The SHA256 pins the exact binary: if the asset changed, the build fails.
ARG VINA_VERSION=1.2.5
ARG VINA_SHA256=fa0126a28a9ea9162d1b161dfa92bc76e632416db28ca246278ea4b2dc6860cb

ENV POLISCREEN_HOME=/opt/poliscreen \
    POLISCREEN_PROJECTS=/data \
    MPLBACKEND=Agg \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR ${POLISCREEN_HOME}

COPY environment.yml ./
RUN mamba env create -f environment.yml -n poliscreen --yes \
 && mamba clean -afy

RUN curl -fsSL "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
      -o /opt/conda/envs/poliscreen/bin/vina \
 && echo "${VINA_SHA256}  /opt/conda/envs/poliscreen/bin/vina" | sha256sum -c - \
 && chmod +x /opt/conda/envs/poliscreen/bin/vina

ENV PATH=/opt/conda/envs/poliscreen/bin:${PATH}

COPY . ${POLISCREEN_HOME}
RUN pip install --no-deps -e .

RUN vina --version && obabel -V && plip --version >/dev/null 2>&1 || true \
 && python -c "import rdkit, py3Dmol, streamlit, poliscreen; print('PoliScreen ready')"

VOLUME ["/data"]
EXPOSE 8501

CMD ["streamlit", "run", "/opt/poliscreen/src/poliscreen/ui/streamlit_app.py", \
     "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
