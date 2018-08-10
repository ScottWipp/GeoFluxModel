#!/bin/bash
#SBATCH -t 10
#SBATCH --ntasks=19 
#SBATCH --mem-per-cpu=4000	
#SBATCH --share				
#SBATCH --account=schmerr-lab-hi
#SBATCH --mail-type=ALL
#SBATCH --job-name=geonu_test
#SBATCH -L matlab
. ~/.profile
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "run('/lustre/swipp/code/CODE_LITHO1_CalcAbund_v2_Cluster.m'); exit"

