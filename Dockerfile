FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq -y update \
    && apt-get -qq -y update \
    && apt-get -qq -y install vim \
    default-jre \
    ncbi-blast+ \ 
    automake \
    autoconf \
    default-jre \ 
    build-essential \
    cd-hit \
    mafft \      
    samtools \           
    curl \
    wget \
    bedtools \
    infernal \
    python3-pip \
    python3-dev \
    git \
    libz-dev

RUN umask 022

WORKDIR /usr/local/bin

RUN ln -s /usr/bin/python3 python

WORKDIR /tmp

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf - \
    && mv ./minimap2-2.24_x64-linux/minimap2 /bin/

RUN git clone https://github.com/RemiAllio/MitoFinder.git \
    && cd MitoFinder \
    && ./install.sh

RUN rm -rf /var/lib/apt/lists/*

WORKDIR /bin

RUN git clone https://github.com/marcelauliano/MitoHiFi.git [github.com]

RUN useradd -m mu

USER mu

RUN pip3 --no-cache-dir install --upgrade pip

RUN pip3 install biopython \
    pandas \
    Pillow \
    matplotlib \
    entrezpy \
    dna-features-viewer \
    bcbio-gff

WORKDIR /tmp

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

RUN printf '\nyes\n\n' | bash Miniconda3-latest-Linux-x86_64.sh

ARG CONDA_DIR=/home/mu/miniconda3

RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.bashrc

RUN $CONDA_DIR/bin/conda install mamba -n base -c conda-forge
