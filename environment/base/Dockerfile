FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq -y update \
    && apt-get -qq -y install \
    autoconf \
    automake \
    bedtools \
    build-essential \
    curl \
    default-jre \ 
    git \
    infernal \
    libopenjp2-7 \
    libtiff5 \
    libz-dev \
    mafft \
    ncbi-blast+ \
    python3-dev \
    python3-pip \
    samtools \
    wget \
    vim \
    && rm -rf /var/lib/apt/lists/*
# this will install python3.6.9 and python2.7.17
# and pip3 9.0.1  under /usr/bin/  

WORKDIR /opt

# /opt/minimap2-2.24_x64-linux
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf  - --no-same-owner

# /opt/hifiasm-0.16.1
RUN curl -L https://github.com/chhylp123/hifiasm/archive/refs/tags/0.16.1.tar.gz \
    | tar -xzvf - \
    && cd hifiasm-0.16.1 \
    && make \
    && wget -P /usr/local/src https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 /usr/local/src/get-pip.py \
    && python2 -m pip --no-cache-dir install biopython==1.70
# This will install pip2 20.3.4 under /usr/local/bin
# and biopython 1.70 under /usr/local/lib/python2.7/dist-packages

# /opt/MitoFinder
RUN git clone https://github.com/RemiAllio/MitoFinder.git \
    && cd MitoFinder \
    && ./install.sh \
    && sed -i 's/\/usr\/bin\/python/\/usr\/bin\/env python/' mitofinder
# This step will change the shebang line in mitofinder to use /usr/bin/env python
# The original shebang line is #!/usr/bin/python2.7 and no effect in this stage
RUN chmod -R 755 /opt/MitoFinder

RUN git clone https://github.com/weizhongli/cdhit.git && \
    cd cdhit && \
    make MAX_SEQ=10000000 

RUN mkdir -p /opt/wrappers

COPY mitos_wrapper.sh /opt/wrappers/runmitos.py
COPY mitofinder_wrapper.sh /opt/wrappers/mitofinder
# This will overwrite the original two scripts
# and change the python installed two different conda environments which are not created yet

RUN chmod -R 755 /opt/wrappers

ARG CONDA_DIR=/opt/conda

RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" \
    && bash Miniforge3-$(uname)-$(uname -m).sh -bfp $CONDA_DIR \
    && rm -rf Miniforge3-$(uname)-$(uname -m).sh


RUN $CONDA_DIR/bin/conda create -n mitos_env -c bioconda -c conda-forge mitos=2.1.0 "r-base>4"
# this will create a conda environment named mitos_env
# with python3.12.7 and r-base 4.1.1, and biopython 1.76
RUN $CONDA_DIR/bin/conda create -n mitofinder_env python=2.7
# this will create a conda environment named mitofinder_env
# with python2.7.15

RUN $CONDA_DIR/bin/conda clean -a

RUN mkdir -p /opt/databases

WORKDIR /opt/databases

RUN curl -L https://zenodo.org/record/4284483/files/refseq89m.tar.bz2?download=1 | tar -jxvf - \
    && curl -L https://zenodo.org/record/4284483/files/refseq89f.tar.bz2?download=1 | tar -jxvf -

WORKDIR /opt

# this will install pip3 under /usr/local/bin
# pip 21.3.1 from /usr/local/lib/python3.6/dist-packages/pip
# which is different from the original /usr/bin/pip3 installed above
RUN pip3 --no-cache-dir install --upgrade pip \
    && pip3 --no-cache-dir install biopython \
    pandas \
    Pillow \
    matplotlib \
    entrezpy \
    dna_features_viewer \
    bcbio-gff
# this will install all packages under /usr/local/lib/python3.6/dist-packages
# include biopython 1.79
