FROM ubuntu:24.04

# ----- start of hermit-specific requirements ------
RUN apt-get update && \
    apt-get install -y openssh-server sudo wget bzip2 ca-certificates curl software-properties-common
RUN mkdir /run/sshd && \
    (deluser ubuntu || true) && \
    adduser ubuntu --disabled-password --uid 2000 --gecos "" && \
    echo 'ubuntu ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers 

RUN mkdir /home/ubuntu/.ssh && chmod go-rx /home/ubuntu/.ssh && chown ubuntu:ubuntu /home/ubuntu/.ssh
# ----- end of hermit-specific requirements ------


RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|http://mirrors.kernel.org/ubuntu/|g' /etc/apt/sources.list.d/ubuntu.sources

RUN apt-get update && apt-get install -y \
    git \
    npm \
    screen \
    vim \
    nano \
    curl \
    gnupg \
    apt-transport-https \
    ca-certificates \
    software-properties-common

RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list

RUN apt-get update && apt-get install -y google-cloud-cli


ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    chown -R ubuntu:ubuntu $CONDA_DIR

ENV PATH=$CONDA_DIR/bin:$PATH

USER root

RUN git clone -b releases/0.9.0b1 https://github.com/ConesaLab/IsoAnnot.git /opt/IsoAnnot && \
    chown -R ubuntu:ubuntu /opt/IsoAnnot

USER ubuntu
RUN conda init bash

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

WORKDIR /opt/IsoAnnot/IsoAnnot
RUN conda env create -f snakemake.yaml -n isoannot_env

RUN mkdir -p /opt/IsoAnnot/IsoAnnot/software

RUN ln -s /opt/IsoAnnot/IsoAnnot /home/ubuntu/IsoAnnot

RUN echo "conda activate isoannot_env" >> /home/ubuntu/.bashrc
RUN conda run -n isoannot_env pip install jupyter notebook pandas ipykernel

USER root
EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]