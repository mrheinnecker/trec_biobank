#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --job-name=dinov2
#SBATCH --array=40-42   # nb of parallel tasks (both bounds included)
#SBATCH --output=dino_%a.txt
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --gres=gpu:1
#SBATCH --mem=20000

# Load and activate Conda
module purge
# Load the Anaconda module that matches your needs
module load Mamba/4.14.0-0
# Activate conda in the current environment

source <(/g/easybuild/x86_64/Rocky/8/znver2/software/Mamba/4.14.0-0/bin/conda shell.bash hook)


# Activate the conda environment required for your calculations
conda deactivate
conda activate /g/schwab/marco/conda_microsam

# Inside your conda env
python -c "import certifi,sys; print(certifi.where())"
# Suppose it prints /g/schwab/marco/conda_microsam/lib/python3.13/site-packages/certifi/cacert.pem

export SSL_CERT_FILE=/g/schwab/marco/conda_microsam/lib/python3.13/site-packages/certifi/cacert.pem
export REQUESTS_CA_BUNDLE="$SSL_CERT_FILE"
export CURL_CA_BUNDLE="$SSL_CERT_FILE"


# Run your script
python /g/schwab/marco/repos/trec_biobank/MainLine/Embedding/dinov2.py $SLURM_ARRAY_TASK_ID
