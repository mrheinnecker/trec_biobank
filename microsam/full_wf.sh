#!/bin/bash

## run gregoires full workflow

sbatch /g/schwab/marco/repos/trec_biobank/MainLine/Segmenting/microsamsegmenting.sh

sbatch /g/schwab/marco/repos/trec_biobank/MainLine/Segmenting/maskmerge.sh

sbatch /g/schwab/marco/repos/trec_biobank/MainLine/Extraction/getimages.sh

sbatch /g/schwab/marco/repos/trec_biobank/MainLine/Embedding/dinov2.sh


module load Mamba/4.14.0-0
# Activate conda in the current environment

source <(/g/easybuild/x86_64/Rocky/8/znver2/software/Mamba/4.14.0-0/bin/conda shell.bash hook)

conda activate /g/schwab/marco/conda_microsam


## move to directory where scripts are in
cd /g/schwab/marco/repos/trec_biobank/MainLine/Classifying


## from here scripts need to be executed from within python (maybe ill change this)
## if i want to use the pretrained parameters from gregoire, i can simply go for the 
## pretrained script from gregoire and use the pretrained labelled data to predict the cluster of new organelles


## if i want to retrain my clusters (but then i also have to redo the manual labelling)
## there is the knn_new.py script.... which i havent really tested for now

python3


