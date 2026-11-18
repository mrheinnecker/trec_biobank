# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 10:33:01 2025

@author: TEAM
"""

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN,HDBSCAN,KMeans, AgglomerativeClustering,OPTICS
from sklearn.neighbors import KNeighborsClassifier,LocalOutlierFactor
from sklearn.model_selection import train_test_split,StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score,confusion_matrix,ConfusionMatrixDisplay
import numpy as np
import umap
from getlabels import get_labeled,matplot_display
import os
import matplotlib.pyplot as plt

labeled_data=get_labeled()

model_name="dino"

if model_name=="dino":
    pkl_name="originalDino"
elif model_name=="finetuned":
    pkl_name="teacherStudent"
elif model_name=="finetuned2":
    pkl_name="teacherStudent2"
elif model_name=="SimCLR":
    pkl_name="SimCLR"
dflist=[]
for n in range(40,42):
    df = pd.read_pickle(f"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cell_nb_{n}/{pkl_name}.pkl")
    df['cellnb'] = n
    #print(df.shape)
    dflist.append(df)
df = pd.concat(dflist)
print(df)
#print(df)
#print(df.iloc[12783])
df.set_index(["cellnb", "image_name"], inplace=True)

## this will not wor if i dont have the labelled data

df["cluster"]=None
for i,key in enumerate(labeled_data.keys()):
    list_of_indices=labeled_data[key]
    list_of_indices=[(int(a),f"organelle_{b}.png") for a,b in list_of_indices]
    #print(list_of_indices)
    #print(df)
    df.loc[list_of_indices,"cluster"]=key
print(df['cluster'][df['cluster'].notna()])

X_sizes = df[["size"]]
Y=df['cluster']
X = df.drop(["size",'cluster'], axis=1)
print(X.columns)
print("Table size :", X.shape)
X=X.values
# Assume X is your data matrix (shape: [n_samples, n_features])

# 1. Standardize features (important for PCA)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 2. Apply PCA then UMAP to reduce dimensions
pca = PCA(n_components=0.9)
X_pca = pca.fit_transform(X_scaled)
print(f"Number of components selected: {pca.n_components_}")
#umap_model = umap.UMAP(n_neighbors=15, n_components=50,random_state=12)
#X_umap=X_pca
#X_umap = umap_model.fit_transform(X_pca)

size_scaled = scaler.fit_transform(np.log1p(X_sizes))
X_pca = scaler.fit_transform(X_pca)



X_final = np.hstack((X_pca, size_scaled))
print(X_final)
print("Table size :", X_final.shape)


# 3. Cluster 
df_final=pd.concat([pd.DataFrame(X_final, index=df.index),Y],axis=1)




labeled=df_final[df_final['cluster'].notna()]
print(labeled)
unlabeled=df_final[~df_final['cluster'].notna()]
method="KNN"


plt.plot([x/len(size_scaled) for x in range(len(size_scaled))], sorted(size_scaled), color='b')
plt.plot([x/len(size_scaled[df_final['cluster'].notna()]) for x in range(len(size_scaled[df_final['cluster'].notna()]))], sorted(size_scaled[df_final['cluster'].notna()]), color='r')
Sthreshold = np.percentile(size_scaled[df_final['cluster'].notna()], 5)
print("Size threshold : ",Sthreshold)
plt.plot([0,1],[Sthreshold,Sthreshold],color='g')

plt.savefig(f'/g/schwab/marco/projects/trec_biobank/slurm_outputs/SizeRep.png')
plt.clf()
print(f"Figure was saved successfully : SizeRep.png")


labeled=labeled[size_scaled[df_final['cluster'].notna()]>Sthreshold]
unlabeled=unlabeled[size_scaled[~df_final['cluster'].notna()]>Sthreshold]

n=5

X=labeled.drop('cluster',axis=1).values
y=labeled['cluster'].values

print("X :",len(X))

evaluation=False
ktuning=False

K=5

if evaluation:
    skf = StratifiedKFold(n_splits=n, shuffle=True, random_state=2)
    ks = []
    training = []
    test = []
    scores = []
    full_training=[]
    if ktuning:
        start,end=2,30
    else:
        start,end=K,K+1
    
    for k in range(start, end):
        clf = KNeighborsClassifier(n_neighbors = k,metric="euclidean")
        clf.fit(X, y)
    
        training_score = clf.score(X, y)
        if k==K:
            conf_matrix = confusion_matrix(y, clf.predict(X))*0

        ks.append(k)
        full_training.append(training_score)
    plt.plot(ks, full_training, color='r')
    i=0
    for train_index, test_index in skf.split(X,y ):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        
        training.append([])
        test.append([])
        scores.append({})
        
        for k in range(start, end):
            clf = KNeighborsClassifier(n_neighbors = k,metric="euclidean")
            clf.fit(X_train, y_train)
        
            training_score = clf.score(X_train, y_train)
            test_score = clf.score(X_test, y_test)
            if k==K:
                conf_matrix +=confusion_matrix(y_test, clf.predict(X_test))

        
            training[i].append(training_score)
            test[i].append(test_score)
            scores[i][k] = [training_score, test_score]
        plt.scatter(ks, training[i], color='k')
        plt.scatter(ks, test[i], color='g')
        i+=1
    if ktuning:
        avg_training=[sum([training[j][i] for j in range(n)])/n for i in range(len(training[0]))]
        avg_test=[sum([test[j][i] for j in range(n)])/n for i in range(len(test[0]))]
        plt.plot(ks, avg_training, color='k')
        plt.plot(ks, avg_test, color='g')
        plt.savefig(f'/g/schwab/marco/projects/trec_biobank/slurm_outputs/ScoreK.png')
        plt.clf()
        print(f"Figure was saved successfully : ScoreK.png")
    conf_matrix_avg = conf_matrix / conf_matrix.sum(axis=1, keepdims=True)* 100
    disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix_avg, display_labels=clf.classes_)
    disp.plot(cmap=plt.cm.Blues, values_format=".0f")
    plt.title(f"Confusion Matrix_{model_name}")
    plt.xticks(rotation=90)
    plt.savefig(f'/g/schwab/marco/projects/trec_biobank/slurm_outputs/ConfMatrix_{model_name}_k{K}.png')
    plt.clf()



#clf = KNeighborsClassifier(n_neighbors = k,weights='uniform')
#clf.fit(X, y)

clf = KNeighborsClassifier(n_neighbors=K,metric="euclidean")
clf.fit(X, y)

labeled["prediction"]=clf.predict(X)
Xu=unlabeled.drop('cluster',axis=1).values
unlabeled["prediction"]=clf.predict(Xu)

# Get distances to the 5 nearest neighbors for each test point
distances, indices = clf.kneighbors(X)
distances2, indices2 = clf.kneighbors(Xu)

# Compute average distance for each test point
average_distances = distances.mean(axis=1)
average_distances2 = distances2.mean(axis=1)



'''
print("X :",len(X))
lof = LocalOutlierFactor(n_neighbors=k, contamination=0.15,novelty=True)
lof.fit(X)


res = lof.predict(X)
print("res :",len(res))
scores = -lof.decision_function(X)
inliers= scores<0
print("Inlier : ",sum(inliers))
#print(scores)
'''
#print("Inlier : ",sum(inliers))
plt.plot([x/len(average_distances) for x in range(len(average_distances))], sorted(average_distances), color='r')
plt.plot([x/len(average_distances2) for x in range(len(average_distances2))], sorted(average_distances2), color='b')
threshold = np.percentile(average_distances, 95)
print("Threshold : ",threshold)
plt.plot([0,1],[threshold,threshold],color='g')
plt.savefig(f'/g/schwab/marco/projects/trec_biobank/slurm_outputs/ScoreLOF.png')
plt.clf()
print(f"Figure was saved successfully : ScoreLOF.png")
# Suppose you define your own threshold based on top 5% scores
#threshold = np.percentile(scores, 95)#1.195
#estimated_outliers = scores > threshold
#estimated_contamination = np.mean(estimated_outliers)
#print("threshold ",threshold)
#print("estimated_contamination ",estimated_contamination)


print(X.shape)
print(X)

map2d=False

if map2d:

    labeled.loc[average_distances>threshold,"prediction"]="None"
    unlabeled.loc[average_distances2>threshold,"prediction"]="None"
    
    print(np.unique(unlabeled["prediction"],return_counts=True))
    
    
    umap_model = umap.UMAP(n_neighbors=10, n_components=2,random_state=74,min_dist=0.1)
    X_final=umap_model.fit_transform(X_final)
    df_final2D=pd.concat([pd.DataFrame(X_final, index=df.index),Y],axis=1)
    labeled2D=df_final2D[df_final2D['cluster'].notna()]
    unlabeled2D=df_final2D[~df_final2D['cluster'].notna()]
    
    labeled2D=labeled2D[size_scaled[df_final['cluster'].notna()]>Sthreshold]
    unlabeled2D=unlabeled2D[size_scaled[~df_final['cluster'].notna()]>Sthreshold]
    
    labeled2D["prediction"]=labeled["prediction"]
    unlabeled2D["prediction"]=unlabeled["prediction"]
    
    cmap=plt.cm.get_cmap('hsv', len(np.unique(labeled2D["prediction"]))+1)
    for i,key in enumerate(np.unique(labeled2D["prediction"])):
        if key == "None":
            print("\n\n\n******\n\n\n")
            new_df=unlabeled2D[unlabeled2D["prediction"]==key]
            plt.scatter(new_df.iloc[:,0],new_df.iloc[:,1],s=0.5,c='black')
            new_df=labeled2D[labeled2D["prediction"]==key]
            plt.scatter(new_df.iloc[:,0],new_df.iloc[:,1],s=0.5,c='black', label=key)
            continue
        new_df=unlabeled2D[unlabeled2D["prediction"]==key]
        plt.scatter(new_df.iloc[:,0],new_df.iloc[:,1],s=0.5,c=cmap(i))
        new_df=labeled2D[labeled2D["prediction"]==key]
        plt.scatter(new_df.iloc[:,0],new_df.iloc[:,1],s=0.5,c=cmap(i), label=key)
    plt.legend()
    
    
    
    plt.savefig(f'/g/schwab/marco/projects/trec_biobank/slurm_outputs/umapmap.png')
    print(f"Figure was saved successfully : umapmap.png")
    final_df=pd.concat([unlabeled2D,labeled2D])
    print(final_df)
    final_df.reset_index()
    
    final_df.to_pickle(rf"/g/schwab/marco/projects/trec_biobank/slurm_outputs/cluster_table_{model_name}_{method}.pkl")
    '''
    absc=Xu[average_distances2<threshold][:,0]
    ordi=Xu[average_distances2<threshold][:,1]
    plt.scatter(absc,ordi,s=0.5,color='g')
    absc=Xu[average_distances2>threshold][:,0]
    ordi=Xu[average_distances2>threshold][:,1]
    plt.scatter(absc,ordi,s=0.5,color='black')
    absc=X[:,0]
    ordi=X[:,1]
    plt.scatter(absc,ordi,s=0.5,color='b')
    absc=X[average_distances<threshold][:,0]
    ordi=X[average_distances<threshold][:,1]
    plt.scatter(absc,ordi,s=0.5,color='r')
    '''
