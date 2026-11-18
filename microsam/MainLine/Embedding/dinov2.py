# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 13:20:48 2025

@author: TEAM
"""
import sys
print(sys.executable)
import os

import torch
import torch.nn as nn
import torchvision
from torchvision.transforms import Compose, Resize, CenterCrop, ToTensor, Normalize
from PIL import Image
import pandas as pd
from addsizes import add_sizes,add_avg
#from SimCLR import SimCLR
from pl_bolts.models.self_supervised import SimCLR
current_dir = os.path.dirname(__file__)
utils_dir = os.path.join(current_dir, '..','..', 'Other')
sys.path.append(os.path.abspath(utils_dir))
print(sys.path)

# Now you can import the module
from finetuning import DinoFineTuneModel,Dinov2EncoderWrapper

print(sys.argv)
n = sys.argv[1]

print('GPU availability :',torch.cuda.is_available())


    
def load_dinov2_model(pkl_name):
    dinov2 = torch.hub.load('facebookresearch/dinov2', 'dinov2_vits14')
    
    
    
        
    encoder = Dinov2EncoderWrapper(dinov2)
    

    model = DinoFineTuneModel(encoder, nn.Identity())
    model = model.cuda()

    # REPO_DIR = "/g/schwab/marco/models/dinov2"

    # # 2) Make sure your checkpoint is available offline.
    # #    Put the pretrain file into:  $TORCH_HOME/hub/checkpoints/
    # #    Example:
    # #       export TORCH_HOME=/g/schwab/marco/torch_cache
    # #       mkdir -p $TORCH_HOME/hub/checkpoints
    # #       # copy dinov2_vits14 *.pth there (name should match what hubconf expects)
    # os.environ.setdefault("TORCH_HOME", "/g/schwab/marco/torch_cache")

    # # 3) Load from local repo (no internet)
    # dinov2 = torch.hub.load(REPO_DIR, 'dinov2_vits14', source='local')

    # # 4) Your original wrapping
    # encoder = Dinov2EncoderWrapper(dinov2)
    # model = DinoFineTuneModel(encoder, nn.Identity())

    # # 5) Move to device safely
    # device = 'cuda' if torch.cuda.is_available() else 'cpu'
    # model = model.to(device)

    if pkl_name!="originalDino":
    
        checkpoint_path = rf"/g/schwab/marco/projects/trec_biobank/slurm_outputs/{pkl_name}.pth"
        state_dict = torch.load(checkpoint_path, map_location='cuda')
        
        state_dict = state_dict['model_state_dict'] if 'model_state_dict' in state_dict else state_dict
    
        from collections import OrderedDict
        new_state_dict = OrderedDict()
        
        for k, v in state_dict.items():
            if k.startswith('module.'):
                new_key = k[len('module.'):]
            if new_key.startswith('projection_head.') or new_key.startswith('head.')or new_key.startswith('projector.'):
                continue
            new_state_dict[new_key] = v
        
        
        model.load_state_dict(new_state_dict)
    
    
    model.eval()
    return model
def load_SimCLR():
    weight_path = 'https://pl-bolts-weights.s3.us-east-2.amazonaws.com/simclr/bolts_simclr_imagenet/simclr_imagenet.ckpt'
    simclr = SimCLR.load_from_checkpoint(weight_path, strict=False)

    simclr_resnet50 = simclr.encoder
    simclr_resnet50.eval()
    return simclr_resnet50
# Preprocessing pipeline
transform = Compose([
    Resize((252, 252)),
    ToTensor(),
    Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
])

def extract_features(path, model,device):
    # Load image
    img = Image.open(path).convert('RGB')

    # Preprocess
    x = transform(img).unsqueeze(0).to(device)  # add batch dim

    # Extract features
    with torch.no_grad():
        features = model(x)

    return features

def lootfolderpanda(path, model,device):
    images=[]
    image_names=[]
    for x in os.listdir(path):
        if x.endswith(".png"):
            images.append(Image.open(os.path.join(path,x)).convert('RGB'))
            image_names.append(x)
    print(len(images))
    #print(len(os.listdir(path)))
    print(len(image_names))
    batch_tensors = [transform(img) for img in images]  # list of tensors [C, H, W]
    batch = torch.stack(batch_tensors, dim=0).to(device)
    print(len(batch))
    with torch.no_grad():
        features = model(batch)
    print(len(features))
    #print(features[0])
    outputs_np = features.detach().cpu().numpy()#replace with features[0] when using SimCLR
    print(len(outputs_np))
    df = pd.DataFrame(outputs_np)
    df['image_name'] = image_names
    return df

model_name="dino"
if model_name=="dino":
    pkl_name="originalDino"
elif model_name=="finetuned":
    pkl_name="teacherStudent"
elif model_name=="finetuned2":
    pkl_name="teacherStudent2"

if __name__ == "__main__":
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = load_dinov2_model(pkl_name).to(device) # only works if the model has the same architecture as DINOv2
    #model = load_SimCLR().to(device) #change the corresponding line in lootfolderpandas to use this
    path = "/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/maskstoreContext"


    df = lootfolderpanda(path, model,device)
    destination=f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/{pkl_name}.pkl"
    df.to_pickle(destination)
    add_sizes(n,destination)
    add_avg(n,destination)
    "/g/schwab/marco/projects/trec_biobank/maskstore"
    #print(len(df))
    print("Feature vector shape:", df.shape)