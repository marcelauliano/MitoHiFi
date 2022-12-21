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
    libz-dev \
    libopenjp2-7 \
    libtiff5 \
    && apt-get -y clean \
    && rm -rf /var/lib/apt/lists/* /var/tmp/*

RUN umask 022

WORKDIR /usr/local/bin

RUN ln -s /usr/bin/python3 python

WORKDIR /opt

# /opt/minimap2-2.24_x64-linux
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf -

# /opt/MitoFinder
RUN git clone https://github.com/RemiAllio/MitoFinder.git \
    && cd MitoFinder \
    && ./install.sh

WORKDIR /bin

RUN git clone https://github.com/marcelauliano/MitoHiFi.git
    && pip3 --no-cache-dir install --upgrade pip \
    && pip3 --no-cache-dir install biopython \
    pandas \
    Pillow \
    matplotlib \
    entrezpy \
    dna_features_viewer \
    bcbio-gff

RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.16.1.tar.gz \
    && wget -P ~/.local/lib https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 ~/.local/lib/get-pip.py \
    && python2 -m pip install biopython==1.70 \
    && tar -xzvf 0.16.1.tar.gz \
    && cd hifiasm-0.16.1 && make

RUN echo "#!/usr/bin/env python3" | cat - /bin/MitoHiFi/mitohifi.py | \
                tee /bin/MitoHiFi/mitohifi.py
RUN echo "#!/usr/bin/env python" | cat - /bin/MitoHiFi/findFrameShifts.py | \
                tee /bin/MitoHiFi/findFrameShifts.py
RUN echo "#!/usr/bin/env python" | cat - /bin/MitoHiFi/fixContigHeaders.py | \
                tee /bin/MitoHiFi/fixContigHeaders.py

RUN chmod -R 755 /bin

# /opt/conda
RUN wget -P /usr/local/src https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash /usr/local/src/Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && /opt/conda/bin/conda install -n base conda-libmamba-solver

# /opt/conda/envs/mitos_env
RUN /opt/conda/bin/conda create -n mitos_env --experimental-solver=libmamba -c bioconda -y mitos \
    && mkdir /opt/conda/envs/mitos_env/bin/clean \
    && ln -s ../runmitos.py /opt/conda/envs/mitos_env/bin/clean/runmitos.py \
    && sed -i "/runmitos\.py/s#python2\", \"#/opt/conda/envs/mitos_env//bin#" /opt/MitoHiFi/parallel_annotation_mitos.py

RUN useradd -m mu

USER mu

WORKDIR /tmp

ENV PATH /opt/minimap2-2.24_x64-linux/:${PATH}
ENV PATH /opt/MitoFinder/:${PATH}
ENV PATH /bin/hifiasm-0.16.1/:${PATH}
ENV PATH /bin/MitoHiFi/:${PATH}
ENV PATH /opt/conda/envs/mitos_env/bin/clean/:${PATH}

