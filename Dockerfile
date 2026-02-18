# syntax=docker/dockerfile:1.4
FROM bolognabiocomp/deepmito:latest

LABEL org.opencontainers.image.title="CoMR"
LABEL org.opencontainers.image.description="Containerized runtime for the CoMR Snakemake pipeline"
LABEL org.opencontainers.image.source="https://github.com/theLabUpstairs/CoMR"

SHELL ["/bin/bash", "-c"]

ARG MAMBA_ROOT_PREFIX=/opt/micromamba
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV MAMBA_ROOT_PREFIX=/opt/micromamba

RUN sed -i 's|deb.debian.org/debian|archive.debian.org/debian|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org/debian-security|archive.debian.org/debian-security|g' /etc/apt/sources.list && \
    printf 'Acquire::Check-Valid-Until "false";\nAcquire::AllowInsecureRepositories "true";\n' > /etc/apt/apt.conf.d/99archive && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        curl ca-certificates bzip2 build-essential make gnat \
        cpanminus && \
    mkdir -p /x86_64-conda-linux-gnu && \
    ln -sfn / /x86_64-conda-linux-gnu/sysroot && \
    printf '#include <locale.h>\n' > /usr/include/xlocale.h && \
    rm -rf /var/lib/apt/lists/*

RUN cpanm --notest Math::Cephes Perl6::Slurp Inline::C Inline && \
    cpanm --notest --local-lib=/opt/software/perl5 Math::Cephes Perl6::Slurp Inline Inline::C

RUN set -euxo pipefail && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp && \
    mv /tmp/bin/micromamba /usr/local/bin/micromamba && \
    rm -rf /tmp/bin

# Create micromamba environment with all runtime dependencies
COPY environment.yml /tmp/environment.yml
RUN --mount=type=cache,target=/opt/micromamba/pkgs \
    micromamba create -y -n comr -f /tmp/environment.yml && \
    micromamba clean -a -y

ENV PATH=/usr/local/bin:/opt/micromamba/envs/comr/bin:$PATH \
    CONDA_DEFAULT_ENV=comr

ENV SOFTWARE_ROOT=/opt/software
RUN mkdir -p ${SOFTWARE_ROOT}

# Install Mitofates and Mitoprot II (fetched ahead of time)
COPY third_party/MitoFates ${SOFTWARE_ROOT}/MitoFates
COPY third_party/mitoprotII ${SOFTWARE_ROOT}/mitoprotII
RUN set -euxo pipefail && \
    chmod 755 ${SOFTWARE_ROOT}/MitoFates/MitoFates.pl && \
    chmod -R 755 ${SOFTWARE_ROOT}/MitoFates && \
    mkdir -p ${SOFTWARE_ROOT}/MitoFates/bin/modules/_Inline && \
    chmod 777 ${SOFTWARE_ROOT}/MitoFates/bin/modules/_Inline && \
    ln -sf /usr/bin/gcc /bin/x86_64-conda-linux-gnu-gcc && \
    cd ${SOFTWARE_ROOT}/mitoprotII && \
    gnatmake mitoprot.adb -o mitoprot && \
    chmod -R 755 ${SOFTWARE_ROOT}/mitoprotII

# Install helper wrappers
COPY containers/bin /usr/local/bin
RUN chmod +x /usr/local/bin/*.sh && \
    ln -sf /usr/local/bin/snakemake.sh /usr/local/bin/snakemake

WORKDIR /opt/CoMR

# Copy pipeline source into the image (users can still bind-mount to override)
COPY . /opt/CoMR

# Default to an interactive shell so snakemake commands can be executed directly
ENTRYPOINT ["/bin/bash"]
