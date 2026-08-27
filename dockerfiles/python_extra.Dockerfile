ARG PYTHON_TAG=3
FROM python:${PYTHON_TAG}

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
    apt-get -y install --no-install-recommends curl && \
    rm -rf /var/lib/apt/lists/*
