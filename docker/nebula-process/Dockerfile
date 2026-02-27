FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    ca-certificates \
    git \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/miniconda/bin:${PATH}

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

RUN conda update -y conda && \
    conda create -y -n renv \
        -c conda-forge \
        -c bioconda \
        r-base=4.3.3 \
        r-rcpp \
        r-rcppparallel \
        r-matrix \
        r-nloptr \
        r-rfast \
        r-trust \
        r-parallelly \
        r-future \
        r-dofuture \
        r-foreach \
        r-dorng \
        r-seurat \
        r-spam \
        r-remotes \
        r-cluster \
        r-readr \
        r-optparse \
        r-dplyr \
        r-biocmanager \
        bioconductor-singlecellexperiment

ENV PATH=/miniconda/envs/renv/bin:${PATH}

RUN Rscript -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages('nebula') \
"
COPY nebula_run.r /usr/local/bin/nebula_run.r
RUN chmod +x /usr/local/bin/nebula_run.r

WORKDIR /work

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

