#!/bin/sh
#$ -S /bin/sh
#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
##SBATCH --qos=long
#SBATCH --partition=amilan
#SBATCH --output=cellbender_jkic_%a.out
#SBATCH --error=cellbender_jkic_%a.err
##SBATCH --mail-type=ALL
##SBATCH --mail-user=jun.inamo@cuanschutz.edu

export TMP="/scratch/alpine/jinamo@xsede.org"
export TMPDIR="/scratch/alpine/jinamo@xsede.org"
export TEMP="/scratch/alpine/jinamo@xsede.org"

# conda create -n cellbender python=3.7
module load anaconda
conda activate cellbender
# pip install cellbender

seq_libs=(List LB183 LB189 LB214 LB215 LB216 LB219 LB220 LB221)
sample_id=${seq_libs[$SLURM_ARRAY_TASK_ID]}

cellbender remove-background --help

DIR=/pl/active/fanzhanglab/jinamo/repertoire/SjS_Takeshita/JKiC/
data_type="GEX"
cd ${DIR}/data/${sample_id}/${data_type}/

inputh5="raw_feature_bc_matrix.h5"
outputh5="cellbender_feature_bc_matrix.h5"

cellbender remove-background \
        --input ${inputh5} \
        --output ${outputh5} \
        --fpr 0.01 \
        --epochs 150 \
        --cpu-threads 4 \
        --debug

# sbatch scripts/repertoire/cellbender_jkic.sh
