#!/bin/bash

(
eval "$($CONDA_DIR/bin/conda shell.bash hook)"
conda activate mitofinder &&
PYTHONPATH= /home/mu/miniconda3/envs/mitofinder/bin/mitofinder "$0"
conda deactivate
)
