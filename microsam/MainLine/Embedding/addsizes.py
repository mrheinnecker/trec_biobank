# -*- coding: utf-8 -*-
"""
Created on Thu May 15 10:27:44 2025

@author: TEAM
"""
import numpy as np
from PIL import Image
import sys
import pandas as pd
import zarr
n = sys.argv[1]

            
SCALE=3
def add_sizes(n,destination):
    MODEL_TYPE="vit_b"
    instances =np.asarray(Image.open("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/mask_"+MODEL_TYPE+"_merged.png")).copy()
    masks,nbs = np.unique(instances,return_counts=True)
    df = pd.read_pickle(destination)
    print(len(nbs),len(df))
    masks=[f"c{n}o{i}.png" for i in masks]
    mapping = dict(zip(masks, nbs))
    df['size'] = df['image_name'].map(mapping)
    print(df[['image_name','size']])
    df.to_pickle(destination)
    
def add_avg(n,destination):
    MODEL_TYPE="vit_b"
    path = '/g/schwab/Karel/Mobie_project_dinoflagellate/data/VSM20_A1_AM1/images/ome-zarr/VSM20_A1_AM1_'+"0"*(3-len(n))+n+'.ome.zarr'
    dataset = zarr.open_group(path, mode = 'r')
    image=dataset['s3']
    mask =np.asarray(Image.open("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/mask_"+MODEL_TYPE+"_merged.png")).copy()
    
    H = min(image.shape[0], mask.shape[0])
    W = min(image.shape[1], mask.shape[1])
    image = image[:H, :W]
    mask = mask[:H, :W]
    
    labels = np.unique(mask)
    avg_intensities = {}

    for label in labels:
        if label == 0:
            continue
        region_pixels = image[mask == label]
        avg_intensity = region_pixels.mean() if region_pixels.size > 0 else np.nan
        avg_intensities[f"c{n}o{i}.png"] = avg_intensity
        
    df = pd.read_pickle(destination)
    print(f'Cell nb {n} is done')
    df['average'] = df['image_name'].map(avg_intensities)
    df.to_pickle(destination)

if __name__ == "__main__":
    add_sizes(n)
    add_avg(n)