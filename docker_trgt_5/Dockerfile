FROM ubuntu:22.04

ARG TRGT_VERSION=5.0.0

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        tabix && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt/trgt && \
    curl -L \
        "https://github.com/PacificBiosciences/trgt/releases/download/v${TRGT_VERSION}/trgt-v${TRGT_VERSION}-x86_64-unknown-linux-gnu.tar.gz" \
        -o /tmp/trgt.tar.gz && \
    tar -xzf /tmp/trgt.tar.gz -C /opt/trgt && \
    find /opt/trgt -type f -name trgt -exec chmod +x {} \; -exec ln -s {} /usr/local/bin/trgt \; && \
    rm -f /tmp/trgt.tar.gz && \
    trgt --version
