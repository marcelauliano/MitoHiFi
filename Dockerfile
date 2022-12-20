# Base Image
FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq -y update \
    && apt -qq -y update \
    && apt -qq -y install vim \
    && apt-get -qq -y install default-jre \
    && apt-get -qq -y install ncbi-blast+ \ 
    && umask 022 \
    && apt-get install -y python3-pip python3-dev \
    && cd /usr/local/bin \
    && ln -s /usr/bin/python3 python \
    && pip3 --no-cache-dir install --upgrade pip \
    && pip3 install biopython \
    && pip3 install pandas \
    && pip3 install matplotlib \
    && pip3 install entrezpy \
    && pip3 install dna-features-viewer \
    && pip3 install bcbio-gff \
    && apt-get -qq -y install automake autoconf \
    && apt -qq -y install default-jre \ 
    && apt-get -qq -y install build-essential \
    && apt-get -qq -y install cd-hit \
    && apt-get -qq -y install mafft \      
    && apt-get -qq -y install samtools \           
    && apt-get -qq -y install curl \
    && apt-get -qq -y install wget \
    && apt-get -qq -y install bedtools \
    && apt-get -qq -y install infernal \
    && curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf - \
    && mv ./minimap2-2.24_x64-linux/minimap2 /bin/ \
    && cd /bin/ \
    && apt-get -qq -y install git \
    && git clone https://github.com/RemiAllio/MitoFinder.git \
    && cd MitoFinder \
    && ./install.sh  \
    && cd /bin/ \
    && apt-get -qq -y install libz-dev \
    && rm -rf /var/lib/apt/lists/* 
    
RUN useradd -m mu

USER mu

WORKDIR /tmp

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

RUN printf '\nyes\n\n' | bash Miniconda3-latest-Linux-x86_64.sh

ARG CONDA_DIR=/home/blobtoolkit/miniconda3

RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.bashrc

RUN $CONDA_DIR/bin/conda install mamba -n base -c conda-forge

RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.16.1.tar.gz \
    && wget -P ~/.local/lib https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 ~/.local/lib/get-pip.py \
    && python2 -m pip install biopython==1.70 \
    && tar -xzvf 0.16.1.tar.gz \
    && cd hifiasm-0.16.1 && make


RUN cd /bin/ \
    && git clone https://github.com/marcelauliano/MitoHiFi.git

ENV PATH /bin/MitoFinder/:${PATH}
ENV PATH /bin/hifiasm-0.16.1/:${PATH}
ENV PATH /bin/MitoHiFi/:${PATH}

RUN echo "#!/usr/bin/env python" | cat - /bin/MitoHiFi/mitohifi.py | \
                tee /bin/MitoHiFi/mitohifi.py
RUN echo "#!/usr/bin/env python" | cat - /bin/MitoHiFi/findFrameShifts.py | \
                tee /bin/MitoHiFi/findFrameShifts.py
RUN echo "#!/usr/bin/env python" | cat - /bin/MitoHiFi/fixContigHeaders.py | \
                tee /bin/MitoHiFi/fixContigHeaders.py

RUN chmod -R 755 /bin

RUN conda activate mitos
