#!/bin/bash -l

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=1:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=muse_qa
#SBATCH --output=../joblog/muse_qa.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=10 # Number of cpu cores on one node
#SBATCH --mem-per-cpu=12GB




export SINGULARITY_LOCALCACHEDIR="/temp_work/$(whoami)/temp_build" \
       SINGULARITY_CACHEDIR="/temp_work/$(whoami)/temp_build" \
       SINGULARITY_TMPDIR="/temp_work/$(whoami)/temp_build"

singularity run -B ..:/project:rw -B /lab-share/Psych-Jalbrzikowski-PNR-e2/Public/code/LunaR_files/Rpackages:/Rpackages /home/$(whoami)/Luna.simg Rscript run_QA.R museEannotDenoised
