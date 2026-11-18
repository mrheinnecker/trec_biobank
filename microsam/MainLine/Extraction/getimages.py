# -*- coding: utf-8 -*-
"""
Created on Thu May 15 10:27:44 2025

@author: TEAM
"""
import numpy as np
import napari
from PIL import Image
import zarr,pprint,os,sys
import scipy
n = sys.argv[1]


def masksplit(dataset,instances,scale=3):
    expansion = 1 #expansion >1 means that the output images will be bigger to include a larger context
    h,w=instances.shape
    images=[]
    contextualized=[]
    nbs=[]
    indivmasks=[]
    for x in range(h):
        if x%(h/10)==0:
            print('row :',x)
        for y in range(w):
            if instances[x,y]!=0:
                while not len(indivmasks)>=instances[x,y]:
                    indivmasks.append([[] for _ in range(h)])
                indivmasks[instances[x,y]-1][x].append(y)
    #print(indivmasks)
    tot= len(indivmasks)
    print("I say :",len(indivmasks))
    print("Sir np says :",len(np.unique(instances))-1)
    #indivmasks=indivmasks[90:]
    pbs=0
    for smth,mask in enumerate(indivmasks):
        #print('mask '+str(smth+1)+'/'+str(tot))
        for i in range(h):
            if mask[i]!=[]:
                #print(instances[i][mask[i][0]])
                #nbs.append(instances[i][mask[i][0]])
                break
        else: #empty mask failcase
            print('skip')
            continue
        for j in range(1,h+1):
            if mask[-j]!=[]:
                break
        else:
            continue
        mini=w;maxi=0
        for k in range(i,h-j+1):
            if mask[k]!=[]:
                mini=min(mini,mask[k][0])
                maxi=max(maxi,mask[k][-1])
        if (h-j+1-i)*(maxi-mini)<500:
            image=dataset['s0']
        elif (h-j+1-i)*(maxi-mini)<2000:
            image=dataset['s1']
            #print(f's1 for mask of dimensions :{h-j+1-i}*{maxi-mini}')
        elif (h-j+1-i)*(maxi-mini)<500000:
            image=dataset['s2']
            #print(f's2 for mask of dimensions :{h-j+1-i}*{maxi-mini}')
        else: 
            image=dataset['s3']
            #print(f's3 for mask of dimensions :{h-j+1-i}*{maxi-mini}')
            
        #print(i,h-j,mini,maxi)
        if (h-j-i+1)*(maxi-mini+1)<100:
            print("Too small ",smth, ", size :",(h-j-i+1)*(maxi-mini+1))
            pbs+=1
            continue
        if  i>0 and mini>0 and maxi+1<w and j>1:
            upscale =image.shape[0]/instances.shape[0]
            #print(f'New dimensions :{round((h-j-i+1)*upscale)}*{round((maxi-mini+1)*upscale)}')
            
            img=np.zeros((round((h-j-i+1)*upscale),round((maxi-mini+1)*upscale)))#both bounds are included\
            
            imgc=np.zeros((round((h-j-i+1)*upscale*expansion),round((maxi-mini+1)*upscale*expansion)))#both bounds are included
            deltax,deltay=(len(imgc[0])-len(img[0]))//2,(len(imgc)-len(img))//2
            small_mask=np.full((h-j-i+1,maxi-mini+1),0.)
            off =round(mini*upscale)
            for k in range(round((h-j-i+1))):
                #print(valid_cols.size)
                if len(mask[round(k+i)]) == 0:
                    continue
                valid_cols = np.arange(round(mask[round(k+i)][0]),round(mask[round(k+i)][-1])+2).astype(int)
                #print(valid_cols)
                row=np.array(mask[round(k+i)])
                valid_cols[np.isin(np.round(valid_cols).astype(int), row,invert=True)]=-1
                valid_cols=np.unique(valid_cols)
                scipy.ndimage.zoom
                #print(len(mask[round(k/upscale+i)])-len(valid_cols),)
                
                if len(valid_cols)==0:
                    continue
                if valid_cols[0]==-1:
                    valid_cols=valid_cols[1:]
                    if len(valid_cols)==0:
                        continue
                if  valid_cols[-1]-off==len(img[0]):
                    valid_cols=valid_cols[:-1]
                    if len(valid_cols)==0:
                        continue
                col=1
                if valid_cols[0]-off>=len(img[0]):
                    continue
                while valid_cols[-col]-off>=len(img[0]):
                    valid_cols[-col]=valid_cols[0]
                    col+=1
                valid_cols= np.unique(valid_cols)

                #assert valid_cols[-1]-off<len(img[0]),f" Modifying {valid_cols[-1]-off}th column in table with size {len(img[0])}"
                small_mask[k,valid_cols-mini]=1.
                #print(start,off,width)
                #print(k, start-off,start+width-off)
                #print(img[k, start-off:start+width-off])
            bigmask = scipy.ndimage.zoom(small_mask,upscale, order=3)
            bigmask = scipy.ndimage.gaussian_filter(bigmask, sigma=5, mode='constant', cval=0.0)
            bigmask=bigmask>0.5
            bigmask=scipy.ndimage.binary_erosion(bigmask,iterations=max(int((upscale-1)/2),1))
            startx,starty=round(mini*upscale),round(i*upscale)+1
            lengthx,lengthy=round((maxi-mini+1)*upscale),round((h-i-j+1)*upscale)
            deltax=min(startx,deltax,len(image[0])-startx-lengthx)
            deltay=min(starty,deltay,len(image)-starty-lengthy)
            img[0:lengthy, 0:lengthx] = image[starty:starty+lengthy, startx:startx+lengthx]
            imgc[0:lengthy+2*deltay, 0:lengthx+2*deltax] = image[starty-deltay:starty+lengthy+deltay, startx-deltax:startx+lengthx+deltax]
            img[~bigmask]=0
            compare=False
            if compare:
                img2=np.zeros((round((h-j-i+1)*upscale),round((maxi-mini+1)*upscale)))#both bounds are included
                img2[0:lengthy, 0:lengthx] = image[starty:starty+lengthy, startx:startx+lengthx]
                naive = scipy.ndimage.zoom(small_mask,upscale, order=0)>0.5
                img2[~naive]=0
                images.append(img2)
            images.append(img)
            contextualized.append(imgc)
            nbs.append(smth+1)
    #print(nbs)
    print("Numer of errors",pbs)
    return contextualized,nbs #return 'images' intead of 'contextualised' for images with a black background
            
SCALE=3
def main():
    path = '/g/schwab/Karel/Mobie_project_dinoflagellate/data/VSM20_A1_AM1/images/ome-zarr/VSM20_A1_AM1_'+"0"*(3-len(n))+n+'.ome.zarr'
    dataset = zarr.open_group(path, mode = 'r')
    MODEL_TYPE="vit_b"
    if not os.path.exists("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/maskstore"):
        os.mkdir("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/maskstore")
    if not os.path.exists("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/maskstoreContext"):
        os.mkdir("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/maskstoreContext")
    #img = dataset["s"+str(SCALE)]
    instances =np.asarray(Image.open("/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_"+n+"/mask_"+MODEL_TYPE+"_merged.png")).copy()
    images,nbs = masksplit(dataset,instances,SCALE)
    for i,img in enumerate(images):
        img = np.clip(img / 256, 0, 255)
        grayscale_image = Image.fromarray(img.astype(np.uint8))
        # Save the image
        grayscale_image.save(f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/maskstoreContext/c{n}o{nbs[i]}.png")
        #print("C:\\Users\\TEAM\\Desktop\\Gregoire\\maskstore\\organelle_"+str(nbs[i])+".png")

if __name__ == "__main__":
    main()