#!/bin/bash

(
eval "$(/opt/conda/bin/conda shell.bash hook)"
conda activate mitofinder_env &&
mitofinder "$@"
conda deactivate
)
