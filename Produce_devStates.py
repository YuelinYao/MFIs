#test.py
#coding=utf-8
import sys

diffCutoff = sys.argv[1]
devStates_path = sys.argv[2]
trainDat_path = sys.argv[3]
minStateDeviation = sys.argv[4]
minNoCells = sys.argv[5]
stateDevAlpha=sys.argv[6]

#print(type(minStateDeviation))

import pandas as pd
import numpy as np
#import igraph as ig
#import argparse
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cut_tree
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.offsetbox import OffsetImage,AnnotationBbox
#import io
#from PIL import Image
#import seaborn as sns
#sns.set(style='white')
#sns.set_palette('colorblind')
from sklearn.metrics import pairwise_distances

devStates = pd.read_csv(devStates_path, dtype=str, index_col=0)
trainDat=pd.read_csv(trainDat_path)





if len(devStates)==0:
    print('no deviating states, terminating...')
    sys.exit()

if len(devStates)==1:
    print('Only one deviation state, terminating...')
    sys.exit()

devStates = devStates[['genes', 'state', 'enrichment', 'pval_corrected']]
#devStates.columns = ['genes', 'state', 'enrichment', 'pval_corrected']
devStates['enrichment']=pd.to_numeric(devStates['enrichment'])
devStates['pval_corrected']=pd.to_numeric(devStates['pval_corrected'])
devStates = devStates[devStates["enrichment"] > float(minStateDeviation)]
devStates = devStates[devStates["pval_corrected"] < float(stateDevAlpha)]



# Binreps is the binary represenations of the interactions: binReps[i] is 1 if cell i is in the deviating state, 0 otherwise.  
binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1
#n = len(binReps)

devStates["No.Cells"]=np.sum(binReps,axis=1)
devStates = devStates[devStates["No.Cells"] > float(minNoCells)]


# Labels that combine the genes and their states---once as list, once as string with newlines
#labsWithStates = devStates.apply(lambda x: [''.join(g) for g in list(zip(x['genes'].split('_'), ['+' if int(s)==1 else '-' for s in x['state']]))], axis=1)
#labsWithStates_str = labsWithStates.apply(lambda x: '\n'.join(x)).values

# Updated Binreps 
binReps = np.array(devStates.apply(lambda x: (trainDat[x['genes'].rsplit('_')]==[int(g) for g in list(str(x['state']))]).all(axis=1), axis=1))*1


# Save binReps
binReps_pd = pd.DataFrame(binReps)
binReps_pd.index=devStates['genes']+devStates['state']



# linkage defines the distances between the binReps, using the Dice-distance: https://en.wikipedia.org/wiki/Sørensen–Dice_coefficient
linked_full = linkage(binReps, 'average', metric='dice')
print(linked_full.shape)



# modularity score
def modularity_score(adjMat, cluster_labels, verbose=False):
    '''
    Calculate the modularity score for a given clustering of a graph.
    '''
    n = len(adjMat)
    m = np.sum(adjMat) / 2.0
    k = np.sum(adjMat, axis=1)
    q = 0.0
    for i in range(n):
        for j in range(n):
            if (cluster_labels[i] == cluster_labels[j]) & (i!=j):
                q += (adjMat[i][j] - k[i]*k[j]/(2.0*m))       
    if verbose:
        print(f'n:{n}')
        print(f'm:{m}')
        print(f'k:{k}')
        print(f'q:{q}')
    return q/(2.0*m)


# function to calculate the cluster labels for a given cutoff
cutAt = lambda x: fcluster(linked_full, x, criterion = 'distance')
pairwiseDists = pairwise_distances(binReps, metric='dice')

# Modularity calculation for a range of N cutoffs (set to 50 for round numbers, but also forced to be round to 2 decimal places)
N = 50
cutoffs = np.linspace(0.01, 0.99, N)
cutoffs = np.round(cutoffs, 2)

modScores = [modularity_score((1-pairwiseDists), cutAt(d)) for d in cutoffs]

diffCutoff_name=diffCutoff
diffCutoff_max = cutoffs[np.argmax(modScores)]

if diffCutoff == "Optimal":
    diffCutoff = diffCutoff_max
    print(f'Using cutoff of {diffCutoff} to maximise modularity score')

print(f'Max modularity score: {max(modScores)}')
print(f'Optimal dice distance: {diffCutoff_max}')
print(f'Using cutoff of {diffCutoff} to cut')

modularity_scores_path="./"+str(minStateDeviation)+'_'+str(minNoCells)+'_'+str(stateDevAlpha)+"_"+"modularity_scores.csv"
pd.DataFrame(zip(cutoffs, modScores), columns=['Cutoff', 'Modularity score']).to_csv(modularity_scores_path)


devStates['cluster'] = fcluster(linked_full, diffCutoff, criterion = 'distance')
#devStates['cluster'] = fcluster(linked_full, diffCutoff, criterion = 'distance')

path='./'+str(diffCutoff_name)+"_"+str(minStateDeviation)+'_'+str(minNoCells)+'_'+str(stateDevAlpha)+'_devStates.csv'
pd.DataFrame(devStates).to_csv(path)

path2='./'+str(diffCutoff_name)+"_"+str(minStateDeviation)+'_'+str(minNoCells)+'_'+str(stateDevAlpha)+'_binReps.csv'
binReps_pd.to_csv(path2)



#print(devStates['cluster'].shape)
print(len(set(devStates['cluster'])))
print("done")





