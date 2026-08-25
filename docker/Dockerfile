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
 && apt-get install -y --no-install-recommends curl ca-certificates git unzip \
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

# Reading CAVER and CaverDock output. On by default and cheap: it is pure standard library, so
# --no-deps is the whole story and nothing can conflict with the docking stack.
ARG WITH_TUNNELS=1
ARG CAVER_TRANSLATE_VERSION=main
RUN if [ "$WITH_TUNNELS" = "1" ]; then \
        pip install --no-cache-dir --no-deps \
            "git+https://github.com/DiegoAnyG/caver-translate.git@${CAVER_TRANSLATE_VERSION}"; \
    fi

# CAVER itself: finding the tunnels, which is the interactive half -- seconds per structure, not
# the minutes CaverDock takes. It is GPL-3, so unlike CaverDock it may be redistributed, and it is
# a Java program, so it runs the same here as on any other machine. That is the whole reason it can
# be in this image and CaverDock cannot: CaverDock is a native binary built against Ubuntu 20.04's
# libraries, and its MPI daemons cannot be made to run on this base.
#
# There is no stable download URL -- caver.cz serves it behind a form -- so the zip arrives one of
# two ways, and the build says which it used:
#
#   --build-arg CAVER_URL=https://.../caver_3.0.2.zip   (with CAVER_SHA256 to pin it)
#   or drop caver_*.zip into installer/vendor/ before building
#
# Version 3.0.2 on purpose: it reproduces CaverWeb to every digit on 8HTB, where 3.0.3 BETA returns
# eight tunnels for the server's six.
ARG WITH_CAVER=1
ARG CAVER_URL=
ARG CAVER_SHA256=
COPY installer/vendo[r] /tmp/caver_vendor
RUN if [ "$WITH_CAVER" = "1" ]; then \
        zip=""; \
        if [ -n "$CAVER_URL" ]; then \
            curl -fsSL "$CAVER_URL" -o /tmp/caver.zip \
         && if [ -n "$CAVER_SHA256" ]; then echo "${CAVER_SHA256}  /tmp/caver.zip" | sha256sum -c -; fi \
         && zip=/tmp/caver.zip; \
        else \
            zip="$(find /tmp/caver_vendor -maxdepth 1 -name 'caver_*.zip' -print -quit 2>/dev/null || true)"; \
        fi; \
        if [ -n "$zip" ]; then \
            mkdir -p /opt/caver && unzip -q "$zip" -d /opt/caver && rm -f /tmp/caver.zip \
         && mamba install -y -n poliscreen openjdk && mamba clean -afy \
         && echo "CAVER installed from $zip"; \
        else \
            echo "WITH_CAVER=1 but no zip was given: set CAVER_URL or put caver_*.zip in installer/vendor/."; \
            echo "The image will still read tunnel results; it just cannot compute them."; \
        fi; \
    fi

# core/caver.py searches for caver.jar rather than trusting one path, but naming it here means the
# search never runs and a broken unpack is noticed at build time instead of in the interface.
ENV POLISCREEN_CAVER=/opt/caver

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
    fi \
 && if [ "$WITH_TUNNELS" = "1" ]; then \
        python -c "from poliscreen.core import tunnels; \
assert tunnels.available(), 'tunnel reading was requested but caver-translate is not importable'; \
print('tunnel reading ready:', tunnels.version())"; \
    fi \
 && if [ -d /opt/caver ]; then \
        python -c "from poliscreen.core import caver; \
assert caver.caver_available(), 'a CAVER zip was unpacked but no runnable caver.jar came out of it'; \
print('tunnel search ready:', caver.caver_jar())"; \
    fi

# The base image already sets tini as its entrypoint, and `init: true` puts Docker's own tini in
# front of it, so the inner one starts at PID 7 and says so on every start. It is right to say it:
# a process that is not PID 1 and not registered as a subreaper cannot reap what the docking
# leaves behind. Registering it costs one variable, silences the warning and keeps both halves
# doing their job.
ENV TINI_SUBREAPER=1

VOLUME ["/data"]
EXPOSE 8501

CMD ["streamlit", "run", "/opt/poliscreen/src/poliscreen/ui/streamlit_app.py", \
     "--server.address=0.0.0.0", "--server.port=8501", "--server.headless=true"]
