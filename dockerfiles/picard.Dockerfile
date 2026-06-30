FROM eclipse-temurin:17-jdk AS build
LABEL stage=buildStage

ENV PICARD_VERSION="3.4.0"

# Install git for building
RUN apt-get update && \
    apt-get --no-install-recommends install -y git && \
    rm -rf /var/lib/apt/lists/*

# Clone Picard repository
RUN git clone --branch ${PICARD_VERSION} --depth 1 https://github.com/broadinstitute/picard.git /usr/picard
WORKDIR /usr/picard

# Build Picard
RUN ./gradlew clean printVersion shadowJar -Dhtsjdk.version=4.2.0

FROM eclipse-temurin:17-jdk AS final
MAINTAINER Broad Institute DSDE <dsde-engineering@broadinstitute.org>

# Install R
RUN apt-get update && \
    apt-get --no-install-recommends install -y r-base && \
    apt-get clean autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /usr/picard
COPY --from=build /usr/picard/build/libs/picard.jar /usr/picard/

RUN mkdir /usr/working
WORKDIR /usr/working
