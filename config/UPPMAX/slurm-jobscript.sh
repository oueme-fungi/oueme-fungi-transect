#!/bin/bash -l
# properties = {properties}
module load texlive &&
conda activate oueme1 &&
{exec_job}
