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
    && rm -rf /var/lib/apt/lists/* \
    && wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.16.1.tar.gz \
    && wget -P ~/.local/lib https://bootstrap.pypa.io/pip/2.7/get-pip.py \
    && python2 ~/.local/lib/get-pip.py \
    && python2 -m pip install biopython==1.70 \
    && tar -xzvf 0.16.1.tar.gz \
    && cd hifiasm-0.16.1 && make
    
COPY ./MitoHiFi/mitohifi.py /bin/
RUN echo "#!/usr/bin/env python" | cat - /bin/mitohifi.py | tee /bin/mitohifi.py
COPY ./MitoHiFi/gfa2fa /bin/
COPY ./MitoHiFi/alignContigs.py /bin/
COPY ./MitoHiFi/circularizationCheck.py /bin/
COPY ./MitoHiFi/cleanUpCWD.py /bin/
COPY ./MitoHiFi/filterfasta.py /bin/
COPY ./MitoHiFi/getMitoLength.py /bin/
COPY ./MitoHiFi/getReprContig.backup.py /bin/
COPY ./MitoHiFi/getReprContig.py /bin/
COPY ./MitoHiFi/parse_blast.py /bin/
COPY ./MitoHiFi/rotate.py /bin/
COPY ./MitoHiFi/rotation.py /bin/
COPY ./MitoHiFi/findMitoReference.py /bin/
COPY ./MitoHiFi/createCoveragePlot.py /bin/
COPY ./MitoHiFi/compareGenesLists.py /bin/
COPY ./MitoHiFi/rotation_mitos.py /bin/
COPY ./MitoHiFi/rotate_genbank.py /bin/
COPY ./MitoHiFi/reverse_complement.py /bin/
COPY ./MitoHiFi/plot_coverage.py /bin/
COPY ./MitoHiFi/plot_coverage_final_mito.py /bin/
COPY ./MitoHiFi/plot_annotation.py /bin/
COPY ./MitoHiFi/plot_annotation.GFF.py /bin/
COPY ./MitoHiFi/parse_blast.py /bin/
COPY ./MitoHiFi/parallel_annotation.py /bin/
COPY ./MitoHiFi/make_genome.py /bin/
COPY ./MitoHiFi/gff_to_gbk.py /bin/
COPY ./MitoHiFi/get_mitos_stats.py /bin/
COPY ./MitoHiFi/getGenesList.py /bin/
COPY ./MitoHiFi/get_depth.py /bin/
COPY ./MitoHiFi/fix_MitoFinder_headers.py /bin/
COPY ./MitoHiFi/fix_improper_gff.py /bin/
COPY ./MitoHiFi/fixContigHeaders.py /bin/
COPY ./MitoHiFi/findFrameShifts.py /bin/

COPY ./MitoHiFi/fetch.py /bin/
COPY ./MitoHiFi/fetch_mitos.py /bin/

RUN chmod -R 755 /bin
CMD ["/bin/bash"]
