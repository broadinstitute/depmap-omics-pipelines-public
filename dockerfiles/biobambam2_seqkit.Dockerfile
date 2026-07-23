FROM biobambam2:2.0.87-release-20180301132713

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update \
    && apt-get -y install --no-install-recommends \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# seqkit: prebuilt static binary (the upstream Dockerfile builds from source
# with Go, which is far heavier). Pinned by version and sha256.
ENV SEQKIT_VERSION=2.13.0
ENV SEQKIT_SHA256=7d686de448464fada1b1988e2e07d693bec68768312da62846bc0e2b502bfc46
RUN curl -fsSL "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz" -o seqkit.tar.gz \
    && echo "${SEQKIT_SHA256}  seqkit.tar.gz" | sha256sum -c - \
    && tar -xzf seqkit.tar.gz -C /usr/local/bin seqkit \
    && rm seqkit.tar.gz \
    && seqkit version
