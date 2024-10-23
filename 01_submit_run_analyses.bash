#!/bin/bash -l

#SBATCH --partition=bch-compute # queue to be used
#SBATCH --time=2:00:00 # Running time (in hours-minutes-seconds)
#SBATCH --job-name=muse_denoised_eannot
#SBATCH --output=../joblog/muse_denoised_eannot.txt # Name of the output file
#SBATCH --nodes=1 # Number of compute nodes
#SBATCH --ntasks=10 # Number of cpu cores on one node
#SBATCH --mem-per-cpu=12GB



export SINGULARITY_LOCALCACHEDIR="/temp_work/$(whoami)/temp_build" \
       SINGULARITY_CACHEDIR="/temp_work/$(whoami)/temp_build" \
       SINGULARITY_TMPDIR="/temp_work/$(whoami)/temp_build"

##################################
#EDIT THESE FOR DIFFERING ANALYSES
##################################
read_in=0
analysis_name=museAnnot
list_name=muse

if [ $read_in == "1" ]; then
	Rscript read_in_only.R
	
	echo "activating spyder env"
	conda activate spyder-env
	echo "running YASA"
	python3 01a_YASA_autostaging.py
	python3 01b_YASA_autostaging.py
	python3 01c_YASA_autostaging.py
	echo "activating base"
	conda activate base
	echo "making highconf files"
	Rscript make_highest_conf_annots.R
fi

Rscript run_analyses.R ${analysis_name} ${list_name}

chmod -R 775 ../derivatives 2>/dev/null

Rscript extract_alleps.R

chmod -R 775 ../joblog 2>/dev/null

Rscript run_QA.R ${analysis_name}
