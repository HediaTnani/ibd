#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --partition=standard
#SBATCH --mem=24576
#SBATCH --account=tumi

module purge
module load anaconda/2020.11-py3.8

python3 -m pip install --user virtualenv

cd ..
python3 -m venv env
source env/bin/activate
cd ibd

pip3 install cobra
pip3 install numpy
pip3 install pandas
pip3 install -e git+https://github.com/gregmedlock/driven@devel#egg=driven

python3 gimme_protect.py