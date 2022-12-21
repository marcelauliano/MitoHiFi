FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq -y update \
    && apt-get -qq -y update \
    && apt-get -qq -y install \
    autoconf \
    automake \
    bedtools \
    build-essential \
    cd-hit \
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

RUN umask 022

WORKDIR /usr/local/bin

RUN ln -s /usr/bin/python3 python

WORKDIR /tmp

RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
    | tar -jxvf - \
    && mv ./minimap2-2.24_x64-linux/minimap2 /bin/ \
    && rm -r /tmp/*

WORKDIR /bin

RUN git clone https://github.com/marcelauliano/MitoHiFi.git

RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.16.1.tar.gz \
    && wget -P ~/.local/lib https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 ~/.local/lib/get-pip.py \
    && python2 -m pip install biopython==1.70 \
    && tar -xzvf 0.16.1.tar.gz \
    && cd hifiasm-0.16.1 && make

RUN echo "#!/usr/bin/env python3" | cat - /bin/MitoHiFi/mitohifi.py | \
    tee /bin/MitoHiFi/mitohifi.py
RUN echo "#!/usr/bin/env python3" | cat - /bin/MitoHiFi/findFrameShifts.py | \
    tee /bin/MitoHiFi/findFrameShifts.py
RUN echo "#!/usr/bin/env python3" | cat - /bin/MitoHiFi/fixContigHeaders.py | \
    tee /bin/MitoHiFi/fixContigHeaders.py

RUN mkdir /bin/wrappers

COPY mitos_wrapper.sh /bin/wrappers/runmitos.py

COPY mitofinder_wrapper.sh /bin/wrappers/mitofinder

RUN chmod -R 755 /bin

RUN useradd -m mu

USER mu

RUN python3 -m pip --no-cache-dir install --upgrade pip

RUN python3 -m pip install \
    bcbio-gff \
    biopython \
    dna-features-viewer \
    entrezpy \
    matplotlib \
    pandas \
    Pillow

WORKDIR /tmp

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && printf '\nyes\n\n' | bash Miniconda3-latest-Linux-x86_64.sh \
    && rm -r /tmp/*

ENV CONDA_DIR=/home/mu/miniconda3

RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.bashrc

ENV PATH /bin/wrappers:/bin/MitoHiFi:/bin/hifiasm-0.16.1:${PATH}

RUN $CONDA_DIR/bin/conda install -n base conda-libmamba-solver

RUN $CONDA_DIR/bin/conda create -n mitos_env --experimental-solver=libmamba -c bioconda -y mitos

RUN $CONDA_DIR/bin/conda create -n mitofinder_env --experimental-solver=libmamba -c bioconda -c conda-forge -y mitofinder

RUN $CONDA_DIR/bin/conda clean -a
