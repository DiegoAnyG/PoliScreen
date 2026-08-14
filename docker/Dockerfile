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
 && apt-get install -y --no-install-recommends curl ca-certificates git \
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

# Everything below is ordered by how often it changes, because Docker throws away every layer
# after the first one that does. The source used to be copied in here, before the engines, so
# editing a single Python file invalidated the ADMET layer and downloaded torch again -- a minute
# of work costing the better part of an hour, every time. What nothing in the source can affect
# now comes first, the scripts next, and the source last.

# The ADMET engine, on by default: without it the screening cycle still runs, but the
# reaction-based analogue builder and the ADMET report are disabled, which is the gap users
# notice first. Pinned to a tag so two people building months apart predict with the same code.
ARG WITH_ADMET=1
ARG ADMELAB_VERSION=v0.3.0

# Its own environment because admet-ai pulls torch and chemprop, whose versions conflict with
# the docking dependencies. That conflict is the reason core/design.py talks to it over a
# subprocess instead of importing it. The CPU wheel is deliberate: these are small models and
# the CUDA build would add roughly 4 GB for no gain.
RUN if [ "$WITH_ADMET" = "1" ]; then \
        git clone --depth 1 --branch "$ADMELAB_VERSION" \
            https://github.com/DiegoAnyG/admelab.git /opt/adme \
     && mamba create -y -n admet python=3.12 \
     && /opt/conda/envs/admet/bin/pip install --no-cache-dir \
            torch --index-url https://download.pytorch.org/whl/cpu \
     && /opt/conda/envs/admet/bin/pip install --no-cache-dir -e /opt/adme \
     && mamba clean -afy; \
    fi

# core/design.py wants a directory that contains admelab/, which the clone above satisfies the
# same way a hand-installed copy does outside the container.
ENV POLISCREEN_ADME_PYTHON=/opt/conda/envs/admet/bin/python \
    POLISCREEN_ADME_ROOT=/opt/adme

# Optional engines, off unless asked for. Both are Linux-only, which is why the container can
# carry them and the native Windows installer never will. They land in $HOME/poliscreen_tools,
# where PoliScreen looks by default, so nothing else has to be configured.
ARG WITH_ADCP=0
ARG WITH_GNINA=0

# Only the installer scripts, so these two layers survive an edit anywhere else in the source.
COPY scripts ${POLISCREEN_HOME}/scripts

# The ADFR installer accepts the Scripps academic licence on the spot, so an image built with
# this flag is for whoever built it and must not be redistributed. parche_adfr.py fixes the
# math domain error that otherwise aborts preparation on phosphoproteins and on structures
# with nearly collinear atoms.
RUN if [ "$WITH_ADCP" = "1" ]; then \
        bash scripts/get_adcp.sh && python scripts/parche_adfr.py; \
    fi

# About 4.5 GB with the CUDA 12 and cuDNN 9 libraries the "static" binary still needs. It does
# nothing without an NVIDIA GPU handed to the container; see docker-compose.gpu.yml.
RUN if [ "$WITH_GNINA" = "1" ]; then \
        bash scripts/get_gnina.sh; \
    fi

# Last, so that editing it rebuilds this layer and the check below, and nothing else.
COPY . ${POLISCREEN_HOME}
RUN pip install --no-deps -e .

# Everything that was asked for has to be usable before the image counts as built. Checking only
# that the imports work let an image ship in which the ADMET engine had been requested and was not
# reachable, so the interface reported it missing at the moment it was needed -- on a machine whose
# owner had asked for it and watched the build succeed. A build that cannot deliver what was
# requested should fail while someone is still watching it.
RUN vina --version && obabel -V && plip --version >/dev/null 2>&1 || true \
 && python -c "import rdkit, py3Dmol, streamlit, poliscreen; print('PoliScreen ready')" \
 && if [ "$WITH_ADMET" = "1" ]; then \
        python -c "from poliscreen.core.design import AdmelabBridge; b = AdmelabBridge(); \
assert b.available(), 'ADMET was requested but the bridge cannot reach it: python=%s root=%s' % (b.python, b.root); \
print('ADMET engine ready:', b.python)"; \
    fi \
 && if [ "$WITH_GNINA" = "1" ]; then \
        python -c "from poliscreen.core.docking import gnina_exe; e = gnina_exe(); \
assert e, 'gnina was requested but no executable was found'; print('gnina ready:', e)"; \
    fi

VOLUME ["/data"]
EXPOSE 8501

CMD ["streamlit", "run", "/opt/poliscreen/src/poliscreen/ui/streamlit_app.py", \
     "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
