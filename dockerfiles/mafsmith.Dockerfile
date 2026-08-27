FROM rust:1.89-slim AS builder

ENV MAFSMITH_VERSION="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        make \
        gcc \
        libc6-dev \
    && rm -rf /var/lib/apt/lists/*

RUN cargo install --git https://github.com/nf-osi/mafsmith --tag v${MAFSMITH_VERSION} mafsmith

FROM ubuntu:24.04

COPY --from=builder /usr/local/cargo/bin/mafsmith /usr/local/bin/mafsmith
