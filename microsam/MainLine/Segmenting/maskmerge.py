# -*- coding: utf-8 -*-
"""
Created on Fri May 16 09:08:59 2025

@author: TEAM
"""

import numpy as np
from PIL import Image
import sys
#from unmask import masksplit
from maskclean import clean_masks,merge_masks


tilenb = 1
MODEL_TYPE="vit_b"
n = sys.argv[1]
masks=[[0 for j in range(tilenb)]for i in range(tilenb)]

for i in range(tilenb):
    for j in range(tilenb):
        masks[i][j] = np.asarray(Image.open(f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/mask_{MODEL_TYPE}_{i}{j}.png")).copy()

# Load an example image from the 'scikit-image' library.
#image = cells3d()[30, 0]
PROCESS=True
# Load the Segment Anything for Microscopy model.

# Run automatic instance segmentation (AIS) on our image.



masks=clean_masks(masks)
if PROCESS and tilenb>1:
    onemask=merge_masks(masks,20)
    print('one mask : ',onemask.shape)
else:
    onemask=masks[0][0]
'''
print(type(instances))
print(instances.shape)
images = masksplit(image,instances)
for i,img in enumerate(images):
    img = np.clip(img / 256, 0, 255)
    grayscale_image = Image.fromarray(img.astype(np.uint16))
    # Save the image
    grayscale_image.save("C:\\Users\\TEAM\\Desktop\\Gregoire\\maskstore\\organelle_"+str(i)+".png")
'''
if PROCESS:
    tosave = Image.fromarray(onemask.astype(np.uint16))
    tosave.save(f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/mask_{MODEL_TYPE}_merged.png")
#mask = Image.fromarray(instances.astype(np.uint8))
# Save the image
#mask.save("C:\\Users\\TEAM\\Desktop\\Gregoire\\mask.png")
# Visualize the image and corresponding instance segmentation result.
"""
v = napari.Viewer()
v.add_image(image, name="Image")
if not PROCESS:
    for i in range(tilenb):
        for j in range(tilenb):
            print(max(0,int(i*h/tilenb-20)))
            t=(max(0,int(i*h/tilenb-20)),max(0,int(j*w/tilenb-20)))
            v.add_labels(masks[i][j],translate = t, name=MODEL_TYPE+str(i)+str(j)+".png")
if PROCESS:
    v.add_labels(onemask, name=MODEL_TYPE+".png")

napari.run()
"""