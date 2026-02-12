FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    ca-certificates \
    git \
    pigz \
    && rm -rf /var/lib/apt/lists/*

# Install meryl
RUN wget -q https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz && \
    tar -xf meryl-1.4.1.Linux-amd64.tar.xz && \
    cp meryl-1.4.1/bin/* /usr/local/bin/ && \
    rm -rf meryl-1.4.1*

WORKDIR /work

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
