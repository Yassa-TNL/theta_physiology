#! /bin/bash
#$ -S /bin/bash -V
#$ -j y
#$ -cwd
#$ -o /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/grid_outputs
#$ -pe openmp 3
#$ -p -2
#$ -l arch=linux-x64
#$ -l h='!Revan&!Bane'
#$ -q yassa.q

# This is an example job
# -pe openmp 3 is the number of CPUs
# -o is the where the standard output will go
# ,shared.q 
echo "$HOSTNAME" 
pwd

/tmp/yassamri/Software/MATLAB/R2018a/bin/matlab -nodisplay -nosplash -nodesktop -r "TransferEntropy; quit;"
