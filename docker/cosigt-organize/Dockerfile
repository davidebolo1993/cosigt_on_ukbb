FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install minimal dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    coreutils \
    && rm -rf /var/lib/apt/lists/*

# Copy organize script
COPY organize.sh /usr/local/bin/organize.sh
RUN chmod +x /usr/local/bin/organize.sh

# Copy entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

WORKDIR /work

ENTRYPOINT ["/entrypoint.sh"]
