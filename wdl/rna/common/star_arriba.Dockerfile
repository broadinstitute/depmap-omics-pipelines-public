FROM us-central1-docker.pkg.dev/depmap-omics/terra-images/samtools:production

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
    g++ \
    pigz \
    r-base \
    wget \
    zlib1g-dev

RUN wget --quiet "https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip" && \
    unzip STAR_2.7.11b.zip && \
    rm STAR_2.7.11b.zip && \
    mv STAR_2.7.11b/Linux_x86_64_static/STAR /usr/local/bin/

RUN wget --quiet "https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz" && \
    tar -xzf arriba_v2.4.0.tar.gz && \
    rm arriba_v2.4.0.tar.gz && \
    cd arriba_v2.4.0 && \
    make && \
    cd .. && \
    chmod 755 arriba_v2.4.0/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz && \
    mv arriba_v2.4.0/arriba /usr/local/bin/

RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org', quiet=TRUE)" && \
    R -e "BiocManager::install(c('tximport', 'jsonlite'))"

RUN wget --quiet "http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz" && \
    tar xvzf gffread-0.12.7.Linux_x86_64.tar.gz && \
    rm gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/* /usr/local/bin/
