#!/bin/bash
# run_pipeline.sh

# activate conda env
source /home/hkaufm49/anaconda3/etc/profile.d/conda.sh
conda activate single_cell_r_env

Rscript scripts/workflow.R
