FROM ghcr.io/marcelauliano/mitohifi-base:b3.2.3

COPY . /opt/MitoHiFi

WORKDIR /opt/MitoHiFi

RUN cd src \
    && chmod +x *.py \
    && sed -i '1i#!/usr/bin/env python3' mitohifi.py findFrameShifts.py fixContigHeaders.py \
    && sed -i '1 s/python\>/python3/' *.py \
    && sed -i 's/"python2",\s//' parallel_annotation_mitos.py \
    && sed -i 's/MITOS\/data/\/opt\/databases/' parallel_annotation_mitos.py

WORKDIR /tmp

ENV CONDA_DIR=/opt/conda

ENV PATH=/opt/wrappers:/opt/cdhit/:/opt/hifiasm-0.16.1/:/opt/MitoHiFi/src/:/opt/MitoFinder/:/opt/minimap2-2.24_x64-linux/:${PATH}
