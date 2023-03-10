#test.py
#coding=utf-8
import sys

diffCutoff = sys.argv[1]
devStates_path = sys.argv[2]
trainDat_path = sys.argv[3]
minStateDeviation = sys.argv[4]
minNoCells = sys.argv[5]
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


devStates = pd.read_csv(devStates_path, dtype=str, index_col=0)
trainDat=pd.read_csv(trainDat_path)





if len(devStates)==0:
    print('no deviating states, terminating...')
    sys.exit()

if len(devStates)==1:
    print('Only one deviation state, terminating...')
    sys.exit()

devStates.columns = ['genes', 'state', 'dev', 'pval']
devStates['dev']=pd.to_numeric(devStates['dev'])
devStates = devStates[devStates["dev"] > float(minStateDeviation)]



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


# linkage defines the distances between the binReps, using the Dice-distance: https://en.wikipedia.org/wiki/S??rensen???Dice_coefficient
linked_full = linkage(binReps, 'average', metric='dice')


print(linked_full.shape)

devStates['cluster'] = fcluster(linked_full, diffCutoff, criterion = 'distance')
#devStates['cluster'] = fcluster(linked_full, diffCutoff, criterion = 'distance')

path='./'+str(diffCutoff)+"_"+str(minStateDeviation)+'_'+str(minNoCells)+'_devStates.csv'
pd.DataFrame(devStates).to_csv(path)

#print(devStates['cluster'].shape)
print(len(set(devStates['cluster'])))
print("done")





