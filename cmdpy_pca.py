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
    p = ap.ArgumentParser()
    arg = p.add_argument
    arg('-if','--stdin',default=sys.stdin,type=str)
    arg('-of','--output_file',default=None,type=str)
    arg('-d','--dist',default='braycurtis',type=str)
    arg('-a','--algorithm',default='mds',type=str,choices=['mds','pca','nmds'])
    arg('-z','--feature_identifier',default='s__',type=str,choices=['k__','s__','PWY','UniRef90'])
    arg('-si','--sample_id',default='sampleID',type=str)
    return vars(p.parse_args())
     
def load_input(stdin): 
    return pd.read_csv(stdin, sep='\t', header=None, index_col=0)

def edit_(f, sid, fi, np):
    feats = [i for i in f.index.tolist() if i fi.split(':')]
    for n in range(np):
        nn = str(n+1)
        print f.loc[[sid,'color'+nn,'shape'+nn]+feats]

if __name__ == '__main__':
    par = read_params(sys.argv)
    edit_(load_input[par['stdin'], 'sampleID', par['feature_identifier'], par['num_plot']])
