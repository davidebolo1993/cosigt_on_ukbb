FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install base dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    pkg-config \
    cmake \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    pigz \
    parallel \
    python3-distutils \
    python3-dev \
    libjemalloc-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Go 1.19.1
RUN wget -q https://go.dev/dl/go1.19.1.linux-amd64.tar.gz && \
    tar -xf go1.19.1.linux-amd64.tar.gz -C /usr/local && \
    rm go1.19.1.linux-amd64.tar.gz
ENV PATH="/usr/local/go/bin:${PATH}"

# Install Rust and Cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN rustup update
RUN cargo install --locked maturin

# Install samtools 1.22
RUN wget -q https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -jxf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    cd .. && rm -rf samtools-1.22*

# Install htslib 1.22
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -jxf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j $(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.22*

# Install BWA-MEM2 2.2.1
RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    mv bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/ && \
    rm -rf bwa-mem2-2.2.1*

# Install bedtools 2.31.0
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static && \
    chmod +x bedtools.static && \
    mv bedtools.static /usr/local/bin/bedtools

# Install meryl
RUN wget -q https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xf meryl-1.4.1.Linux-amd64.tar.xz && \
    rm meryl-1.4.1.Linux-amd64.tar.xz
ENV PATH="/meryl-1.4.1/bin:${PATH}"

# Install kfilt (Go)
RUN git clone https://github.com/davidebolo1993/kfilt && \
    cd kfilt && \
    go mod init kfilt && \
    go mod tidy && \
    go build && \
    cp kfilt /usr/local/bin/ && \
    cd .. && rm -rf kfilt

# Install ODGI (CMake build) - use correct path
RUN git clone --recursive https://github.com/pangenome/odgi.git && \
    cd odgi && \
    cmake -H. -DCMAKE_BUILD_TYPE=Generic -Bbuild && \
    cmake --build build -- -j $(nproc) && \
    cp bin/odgi /usr/local/bin/odgi && \
    cd .. && rm -rf odgi

# Install gfainject (Rust)
RUN git clone https://github.com/AndreaGuarracino/gfainject && \
    cd gfainject && \
    cargo install --force --path . && \
    cp target/release/gfainject /usr/local/bin/ && \
    cd .. && rm -rf gfainject

# Install gafpack (Rust)
RUN git clone https://github.com/pangenome/gafpack.git && \
    cd gafpack && \
    cargo install --force --path . && \
    cp target/release/gafpack /usr/local/bin/ && \
    cd .. && rm -rf gafpack

# Install panplexity (Rust)
RUN git clone https://github.com/AndreaGuarracino/panplexity && \
    cd panplexity && \
    cargo install --force --path . && \
    cp target/release/panplexity /usr/local/bin/ && \
    cd .. && rm -rf panplexity

# Install cosigt (Go)
RUN git clone https://github.com/davidebolo1993/cosigt.git && \
    cd cosigt && \
    go mod init cosigt && \
    go mod tidy && \
    go build && \
    cp cosigt /usr/local/bin/ && \
    cd .. && rm -rf cosigt

# Install Miniconda and create R environment
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

# Accept Terms of Service for Anaconda channels - otherwise it fails
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create conda R environment with required packages (minimal for our scripts)
RUN conda update -y conda && \
    conda create -y -n renv -c conda-forge \
    r-base=4.3.3 \
    r-data.table=1.16.4 \
    r-rjson=0.2.23 \
    r-reshape2=1.4.4 \
    r-dbscan=1.2.2 \
    r-cluster

# Activate conda environment by default
RUN echo "source activate renv" > ~/.bashrc
ENV PATH=/miniconda/envs/renv/bin:${PATH}

# Copy all scripts and R files into the image
COPY process_region_docker.sh /usr/local/bin/process_region_docker.sh
COPY cluster.r /usr/local/bin/cluster.r
COPY coverage_outliers.r /usr/local/bin/coverage_outliers.r
RUN chmod +x /usr/local/bin/process_region_docker.sh
RUN chmod +x /usr/local/bin/*.r

WORKDIR /work

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

