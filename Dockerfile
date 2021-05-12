FROM ubuntu:18.04
#FROM ubuntu:16.04
#FROM ubuntu:18.10


# Run apt update and install samtools
#RUN apt-get update -qq \
# && apt-get install -y \
# libcurl3 \
# samtools \
# && rm -rf /var/lib/apt/lists/*

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq -y update \ 
    && apt-get -qq -y install ncbi-blast+ \ 
    && umask 022 \
    && apt-get install -y python3-pip python3-dev \
    && cd /usr/local/bin \
    && ln -s /usr/bin/python3 python \
    && pip3 --no-cache-dir install --upgrade pip \
    && pip3 install biopython \
    && pip3 install pandas \
    && pip3 install entrezpy \
    #&& apt-get -qq -y install python-biopython \
    #&& apt-get -qq -y install python-pandas \
    && apt-get -qq -y install automake autoconf \
    && apt -qq -y install default-jre \ 
    && apt-get -qq -y install build-essential \
    && apt-get -qq -y install cd-hit \
    && apt-get -qq -y install mafft \	   
    && apt-get -qq -y install samtools \	   
    && apt-get -qq -y install curl \
    && curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf - \
    && mv ./minimap2-2.17_x64-linux/minimap2 /bin/ \
#    && apt-get -qq -y install minimap \	   
#    && apt-get -qq -y install python \
    && cd /bin/ \
    && apt-get -qq -y install git \
    && git clone https://github.com/RemiAllio/MitoFinder.git \
    && cd MitoFinder \
    && ./install.sh  \
    && cd /bin/ \
    && apt-get -qq -y install wget \
    && apt-get -qq -y install libz-dev \
    && rm -rf /var/lib/apt/lists/* \
    && wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.14.2.tar.gz \
    && tar -xzvf 0.14.2.tar.gz \
    && cd hifiasm-0.14.2 && make

# && curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" \
   # && python get-pip.py \
   # && pip install biopython \
   # && pip install pandas
# Create a dev user
#RUN useradd -ms /bin/bash mitohifi
#WORKDIR /home/mitohifi
#USER mitohifi

ENV PATH /bin/MitoFinder/:${PATH}
ENV PATH /bin/hifiasm-0.14.2/:${PATH}

#RUN mkdir /bin/scripts
#COPY ./MitoHiFi/scripts/gfa2fa /bin/
#COPY ./MitoHiFi/mitohifi-v3-fromCirc_v02.11.3-marcela.py /bin/mitohifi_v2.py
COPY ./MitoHiFi/mitohifi_v2.py /bin/
RUN echo "#!/usr/bin/env python" | cat - /bin/mitohifi_v2.py | tee /bin/mitohifi_v2.py
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
COPY ./MitoHiFi/findFrameShits.py /bin/
COPY ./MitoHiFi/fixContigHeaders.py /bin/
#COPY ./MitoFinder/ /bin/MitoFinder/

RUN chmod -R 755 /bin
#ENV SHELL /bin/bash
##RUN chmod 755 /bin/MitoFinder/*
#RUN echo $PATH
#RUN ls -hl /usr//lib/python2.7
#RUN find /usr/lib | grep SeqIO
#RUN ls -hl /usr/lib/python2.7/dist-packages/
#RUN ls -hl /usr/lib/python2.7/dist-packages/Bio/SeqIO/
#RUN ls -hl /usr/lib/python2.7/dist-packages/Bio/
#RUN python --version
#RUN ls -hl /lib/
CMD ["/bin/bash"]
