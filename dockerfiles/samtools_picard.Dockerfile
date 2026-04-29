FROM us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools:production

ENV DEBIAN_FRONTEND=noninteractive

RUN curl -fsSL https://packages.adoptium.net/artifactory/api/gpg/key/public \
        | gpg --dearmor -o /usr/share/keyrings/adoptium.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/adoptium.gpg] https://packages.adoptium.net/artifactory/deb bullseye main" \
        > /etc/apt/sources.list.d/adoptium.list \
    && apt-get -y update \
    && apt-get -y install --no-install-recommends --no-install-suggests \
        temurin-17-jre \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PICARD_VERSION="3.4.0"

RUN curl -SL https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar \
    -o /usr/local/bin/picard.jar

ENV BWA_VERSION="0.7.18"

RUN curl -SL https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz \
    -o /tmp/bwa.tar.gz \
    && tar xvf /tmp/bwa.tar.gz -C /usr/local/src --remove-files \
    && mv /usr/local/src/bwa-* /usr/local/src/bwa \
    && cd /usr/local/src/bwa \
    && make \
    && cp bwa /usr/local/bin/bwa
