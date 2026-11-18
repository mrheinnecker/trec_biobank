# -*- coding: utf-8 -*-
"""
Created on Wed May 14 13:56:02 2025

@author: TEAM

Paths correspond to their location in slurm
"""
import numpy as np
from PIL import Image, ImageOps 
import zarr
from micro_sam.automatic_segmentation import get_predictor_and_segmenter, automatic_instance_segmentation
import sys
import os

print(sys.argv)
n = sys.argv[1]


path = "/g/schwab/Karel/Mobie_project_dinoflagellate/data/VSM20_A1_AM1/images/ome-zarr/VSM20_A1_AM1_"+"0"*(3-len(n))+n+".ome.zarr"
print(path)
dataset = zarr.open_group(path, mode = 'r')

MODEL_TYPE = "vit_b"
tilenb = 2
img = dataset["s3"]
#img = Image.open(path)
#img = img.resize((2048,2048)) 
img = np.asarray(img)
h,w=img.shape
print(img.shape)

os.mkdir("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n)

for i in range(tilenb):
    for j in range(tilenb):
        #print(f"Chunk size: {img.chunks}")
        
        image = img#np.asarray(img[max(0,int(i*h/tilenb-20)):min(h,int((i+1)*h/tilenb+20)),max(0,int(j*w/tilenb-20)):min(w,int((j+1)*w/tilenb+20))])
        
        
        # Load the Segment Anything for Microscopy model.
        #print(help(get_predictor_and_segmenter))
        predictor, segmenter = get_predictor_and_segmenter(model_type=MODEL_TYPE,is_tiled=False,device="cuda")
        
        # Run automatic instance segmentation (AIS) on our image.
        instances = automatic_instance_segmentation(predictor=predictor, segmenter=segmenter, input_path=image,ndim=2)#,tile_shape=(int(h/2+40),int(w/2+40)),halo=(40,40))
        mask = Image.fromarray(instances.astype(np.uint16))
        #image = Image.fromarray(image.astype(np.uint16))
        # Save the image
        #image.save("/g/schwab/GregoireMichelDeletie/Ines/ref_"+n+".png")
        mask.save(f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/mask_{MODEL_TYPE}_{i}{j}.png")
        # Visualize the image and corresponding instance segmentation result.
        print(np.max(mask))
        
        

