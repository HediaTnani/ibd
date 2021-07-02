#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=standard
#SBATCH --mem=8192
#SBATCH --account=tumi

module purge
module load anaconda/2020.11-py3.8

cd /project/tumi/park/ibd/gimme
source ../../env/bin/activate

python3 gimme_9.py