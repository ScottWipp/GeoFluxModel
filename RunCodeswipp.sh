#!/bin/bash
#SBATCH -t 20
#SBATCH --ntasks=1 
#SBATCH --mem-per-cpu=4000	
#SBATCH --share				
#SBATCH --account=schmerr-lab-hi
#SBATCH --mail-type=ALL
#SBATCH --job-name=geoneutrino_test
. ~/.profile
module load matlab
matlab -nodisplay -nosplash -nodesktop -r "run('/lustre/swipp/code/CODE_LITHO1_CalcAbund_v2_Cluster.m')"

