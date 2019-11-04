#!/bin/bash -l
# SBATCH --mail-type=ALL
# properties = {properties}
module load texlive &&
{exec_job}
