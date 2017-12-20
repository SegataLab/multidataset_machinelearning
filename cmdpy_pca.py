#!/usr/bin/env python

import argparse as ap
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from sklearn import decomposition
from sklearn import preprocessing as prep
from sklearn import manifold
from sklearn.metrics import pairwise_distances
import pandas as pd
import numpy as np

def read_params(args):
    distances = ['lbraycurtis', 'sbraycurtis', 'canberra', 'chebyshev'\
	       , 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski'\
	       , 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto'\
	       , 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath'\
	       , 'sqeuclidean', 'yule', 'cityblock', 'cosine', 'euclidean'\
	       , 'l1', 'l2', 'manhattan', 'braycurtis', 'precomputed'] 
    p = ap.ArgumentParser()
    arg = p.add_argument
    arg('-if','--stdin',default=sys.stdin,type=str)
    arg('-of','--output_file',default=None,type=str)
    arg('-d','--dist',default='braycurtis',type=str,choices=distances)
    arg('-a','--algorithm',default='mds',type=str,choices=['mds','pca','nmds'])
    arg('-z','--feature_identifier',default='s__',type=str,choices=['k__','s__','PWY','UniRef90'])
    arg('-si','--sample_id',default='sampleID',type=str)
    arg('-nc','--num_classes',default=2,type=int)
    arg('-cn','--class_name',default='study_condition',type=str)
    return vars(p.parse_args())

def pca(f):
    pca = decomposition.PCA(n_components=2)
    return pca.fit_transform(f), pca.explainedA_variance_ratio_

def mds(dist):
    mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-9, dissimilarity='precomputed')
    return mds.fit(dist).embedding_
     
def nmds(dist):
    nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
    return nmds.fit_transform(dist)

def compute_distances(self, data, metric):
    if metric == 'precomputed':
        return data
    elif metric == 'lbraycurtis':
        ldata = np.matrix([[(math.log(1.0+float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
        return pairwise_distances(ldata, metric="braycurtis")
    elif metric == 'sbraycurtis':
        sdata = np.matrix([[(math.sqrt(float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
        return pairwise_distances(sdata, metric="braycurtis")
    else:
        return pairwise_distances(data, metric=metric)

def load_input(stdin): 
    return pd.read_csv(stdin, sep='\t', header=None, index_col=0)

def edited_(f, sid, fi, ci, nc):
    feats = [i for i in f.index.tolist() for fii in fi.split(':') if fii in i]
    return f.loc[[sid,ci]+['color'+str(n+1) for n in range(nc)]+['shape'+str(n+1) for n in range(nc)]] + feats, :]
    
if __name__ == '__main__':
    par = read_params(sys.argv)
    e = edited_(load_input(par['stdin']), par['sample_id'], par['feature_identifier'], par['class_name'])
    ### print e
