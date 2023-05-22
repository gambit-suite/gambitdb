# Software installation, no database files
FROM mambaorg/micromamba:jammy as app_base

ARG GAMBITDB_SOFTWARE_VERSION="0.0.1"
ARG GAMBITDB_GIT_TAG=v${GAMBITDB_SOFTWARE_VERSION}
ARG GAMBITDB_SRC_URL=https://github.com/theiagen/gambitdb/archive/refs/tags/${GAMBITDB_GIT_TAG}.tar.gz

LABEL base.image="mambaorg/micromamba:0.27.0"
LABEL dockerfile.version="1"
LABEL software="GAMBITDB"
LABEL software.version=${GAMBITDB_SOFTWARE_VERSION}
LABEL description="Create databases for Gambit"
LABEL website="https://github.com/theiagen/gambitdb"
LABEL license="https://github.com/theiagen/gambitdb/blob/master/LICENSE"
LABEL maintainer1="Michelle Scribner"
LABEL maintainer.email1="michelle.scribner@theiagen.com"
LABEL maintainer2="Andrew Page"
LABEL maintainer.email2="andrew.page@theiagen.com"

# Environment
ENV GAMBIT_DB_PATH=./
ENV LC_ALL=C.UTF-8

# Install mamba environment
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # Subsequent RUN commands use environment

RUN mamba install -y gambit

# Install GAMBITDB package
RUN pip install ${GAMBITDB_SRC_URL} && \
  micromamba clean -a -y

USER root
RUN mkdir $GAMBITDB_DB_PATH /data && \
  chown $MAMBA_USER:$MAMBA_USER /data
USER $MAMBA_USER
WORKDIR /data

# Make sure conda, python, and GAMBITDB are in the path
ENV PATH="/opt/conda/bin:${PATH}"

# With database files added
FROM app_base AS app

# Run test
FROM app as test

COPY test.sh .
RUN bash test.sh
