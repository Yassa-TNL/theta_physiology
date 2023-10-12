#! /bin/bash
#$ -S /bin/bash -V
#$ -j y
#$ -cwd
#$ -o /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/grid_outputs
#$ -pe openmp 15
#$ -p -2
#$ -l arch=linux-x64
#$ -q yassa.q,shared.q 

# This is an example job
# -pe openmp 3 is the number of CPUs
# -o is the where the standard output will go

echo "$HOSTNAME"
pwd

/tmp/yassamri/Software/MATLAB/R2018a/bin/matlab -nodisplay -nosplash -nodesktop -r "tv_timedomain_GC; quit;"
