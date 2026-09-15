FROM rocker/r-ver:4.4.0

# install system deps for common R packages and Bioconductor
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libuv1-dev \
    libwebp-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    ca-certificates \
    curl \
    git \
    make \
    && rm -rf /var/lib/apt/lists/*

# install bedtools
RUN curl -SL "https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static" \
    -o /usr/bin/bedtools && \
    chmod +x /usr/bin/bedtools

# restore R packages
WORKDIR /app
COPY ./renv.lock ./
RUN Rscript -e "install.packages('remotes')" \
    && Rscript -e "remotes::install_version('renv', version = '1.0.7')" \
    && Rscript -e "renv::restore()"

COPY ./segment.R ./
COPY ./plot_segments.R ./
COPY ./plot_heat_map.R ./
