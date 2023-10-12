#! /bin/bash
#$ -S /bin/bash -V
#$ -j y
#$ -cwd
#$ -o /mnt/yassamri/iEEG/sandra/analysis_pipeline_final/grid_outputs
#$ -pe openmp 15
#$ -p -2
#$ -l arch=linux-x64
#$ -l h='!Bane'
#$ -q yassa.q,shared.q 

# This is an example job A002_extract_baseline_and_trials A003_condition_wise_spectral_analysis A003a_condition_wise_spectral_analysis_cond_spec_prestim A002b_LFP_from_SU_data
# -pe openmp 3 is the number of CPUs
# -o is the where the standard output will go
echo "$HOSTNAME" 
pwd

/mnt/yassamri/Software/MATLAB/R2018a/bin/matlab -nodisplay -nosplash -nodesktop -r "A003a_condition_wise_spectral_analysis_cond_spec_prestim; quit;"
