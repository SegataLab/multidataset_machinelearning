#!/usr/bin/env python

#  USAGE: 
#  python cmdpy_dataset.py CM_periimplantitis.study_condition:peri-implantitis,col:color:crimson,col:shape:o,col:define:periimplantitis CM_periimplantitis.study_condition:control,col:color:blue,col:shape:o,col:define:control --shrink | python cmdpy_betadiversity.py -a mds -sc coords.txt -of pi_plot 


import argparse as ap
import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from sklearn import decomposition
from sklearn import preprocessing as prep
from sklearn import manifold
from sklearn.metrics import pairwise_distances
import seaborn as sns
import pandas as pd
import numpy as np


class beta_diversity(object):

    def __init__(self, params=None):

        if not isinstance(params, dict):  self.args = self.read_params(sys.argv)

        else: self.args = params

        self.stdin = pd.read_csv(self.args['stdin'], sep='\t', header=None, index_col=0)
 
        self.shapes_types = '.,ov^<>1234sp*hH+xDd|_'
        self.coordinates, self.metadata, self.explained_var = self.beta_div()
        self.sample_in_class, self.atts_of_sample, self.sample_to_class = self.guess_classes(self.metadata)
        self.legend = [ll.split() for ll in list(set([' '.join(l) for l in self.atts_of_sample.values()]))] 

        self.sample_and_coordinates = self.stack_sample_coordinates() \
	    if not self.args['load_coordinates'] \
	    else pd.read_csv(self.args['load_coordinates'], sep='\t', header=0, index_col=0)        

        if self.args['save_coordinates']: self.save_coordinates(self.args['save_coordinates'])

        if self.args['mask']: self.mask()

        if self.args['explain']: self.explain()

        if not (not self.args['no_ordination']): self.scatter_plot()

        if self.args['boxplot']: print 'mia nonna'         



    def read_params(self, args):

        p = ap.ArgumentParser()

        distances = ['lbraycurtis','sbraycurtis','canberra','chebyshev','correlation','dice','hamming','jaccard','kulsinski'\
	    ,'mahalanobis','matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener','sokalsneath'\
	    ,'sqeuclidean','yule','cityblock','cosine','euclidean','l1','l2','manhattan','braycurtis','precomputed'] 

        arg = p.add_argument
        arg(	'-if', '--stdin', default=sys.stdin, type=str)
        arg(	'-of', '--stdout', default=None, type=str)
        arg(	'-sc', '--save_coordinates', default=None, type=str)
        arg(	'-lc', '--load_coordinates', default=None, type=str)
        arg(	'-d', '--distance', default='braycurtis', type=str, choices=distances)
        arg(	'-a', '--algorithm', default='mds', type=str, choices=['mds','pca','nmds','boxp'])
        arg(	'-z', '--feature_identifier', default='s__', type=str, choices=['k__','s__','PWY','UniRef90'])
        arg(	'-si', '--sample_id', default='sampleID', type=str)
        arg(	'-ci', '--classes_id', default='define', type=str)
        arg(	'-m', '--mask', default=[], nargs='+', type=str)

        arg(	'-ll', '--legend_loc', default='lower center', type=str)
        arg(	'-fmt', '--format', default='png', type=str, choices=['svg','png'])


        arg(	'--no_ordination', action='store_false')
        arg(	'--boxplot', action='store_true')
        arg(	'-ex', '--explain', action='store_true')

        arg(	'--title', type=str, default='Ordination Plot')
        arg(	'--dot_size', type=float, default=10)


        return vars(p.parse_args())



    def explain(self):
        print ' sample_in_class: a dict with classes as keys and samples for any keys listed as values.' 
        print ' atts_of_sample: a dict with samples as keys and [class, color, shape] as values.'
        print ' sample_to_class: a dict with samples as keys and class for each as value.'



    def pca(self, f):
        pca = decomposition.PCA(n_components=2)
        return pca.fit_transform(f), pca.explained_variance_ratio_



    def mds(self, d):
        try:
            mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-9, dissimilarity='precomputed')
            return mds.fit(d).embedding_
        except ValueError:
            print 'You Have NaNs in the data: here the coordinates'
            exit(1)     



    def nmds(self, d):
        nmds = manifold.MDS(n_components=2, metric=False, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
        return nmds.fit_transform(d)


    def compute_distances(self, data, metric):
        if metric == 'precomputed': return data
        elif metric == 'lbraycurtis':
            ldata = np.matrix([[(math.log(1.0+float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            return pairwise_distances(ldata, metric='braycurtis')
        elif metric == 'sbraycurtis':
            sdata = np.matrix([[(math.sqrt(float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            return pairwise_distances(sdata, metric='braycurtis')
        else:
            return pairwise_distances(data, metric='braycurtis')


    def load_input(self, stdin): return pd.read_csv(stdin, sep='\t', header=None, index_col=0)


    def edit_(self, f, sid, feat_id):
        feats = [i for i in f.index.tolist() for fii in feat_id.split(':') if fii in i]
        return f.loc[[sid,self.args['classes_id'],'color','shape']+feats, :]


    def guess_classes(self, mdf):
        class_to_samples = dict([(c, mdf.loc[self.args['sample_id'], mdf.loc[self.args['classes_id']].isin([c])].tolist()) for c in mdf.loc[self.args['classes_id'], :].tolist()])
        sample_to_attributes = dict([(s, [a[0] for a in mdf.loc[[self.args['classes_id'],'color','shape'], mdf.loc[self.args['sample_id']].isin([s])].values.tolist()]) \
	    for s in mdf.loc[self.args['sample_id'], :].tolist()])
        sample_to_class = dict([(s,c) for s,c in zip(mdf.loc[self.args['sample_id'], :].tolist(), mdf.loc[self.args['classes_id'], :].tolist())])
        return class_to_samples, sample_to_attributes, sample_to_class            


    def transformation(self, way, f, feat_id, distance):
        data = f.loc[[i for i in f.index.tolist() if (feat_id in i)]].T
        metadata = f.loc[[i for i in f.index.tolist() if (not feat_id in i)]]
        T_func = self.pca if way == 'pca' else (self.mds if way == 'mds' else self.nmds)
        if T_func == self.pca: transformed, exp_var = T_func(data)
        else: transformed, exp_var = T_func(self.compute_distances(data,distance)), None
        return transformed, metadata, exp_var


    def beta_div(self):
        edt = self.edit_(self.stdin, self.args['sample_id'], self.args['feature_identifier'])
        transform, metadata, exp_var = self.transformation(self.args['algorithm'], edt, self.args['feature_identifier'], self.args['distance'])
        return transform, metadata, exp_var


    def stack_sample_coordinates(self):
        toreturn = pd.DataFrame({self.args['sample_id']: self.metadata.loc[self.args['sample_id'], :].tolist(), 'x1': self.coordinates[:, 0], 'x2': self.coordinates[:, 1]})
        toreturn.index = toreturn[self.args['sample_id']]
        del toreturn[self.args['sample_id']]
        return toreturn


    def save_coordinates(self, coordinates_f): self.sample_and_coordinates.to_csv(coordinates_f, sep='\t', header=True, index=True)

    
    def mask(self):
        self.sample_and_coordinates = self.sample_and_coordinates.drop([i for i in self.sample_and_coordinates.index.tolist() if self.sample_to_class[i] in self.args['mask']])
        for i in self.sample_to_class.keys(): 
            if self.sample_to_class[i] in self.args['mask']: del self.atts_of_sample[i]        


    def scatter_plot(self, ax=None):

        sns.set_style('darkgrid')
        fig = False

        if not bool(ax):

            fig, ax = plt.subplots(figsize=(8,6))

            if self.args['algorithm'] == 'mds':

                ax.set_xlim(-0.8, 0.8)
                ax.set_ylim(-0.8, 0.8)
                ax.xaxis.set_ticks(np.arange(-0.6, 0.8, 0.2))
                ax.yaxis.set_ticks(np.arange(-0.6, 0.8, 0.2))

        for c in self.sample_in_class.keys():

            present = [s for s in self.sample_in_class[c] if s in set(self.sample_and_coordinates.index.tolist())]
            present_sample_frame = self.sample_and_coordinates.loc[present]

            scatterp = sns.regplot(x='x1', y='x2', data=present_sample_frame, ax=ax, scatter=True, fit_reg=False, scatter_kws={'s': self.args['dot_size']}\
		, label=self.atts_of_sample[present[0]][0], marker=self.atts_of_sample[present[0]][2], color=self.atts_of_sample[present[0]][1])

        if bool(fig):

            plt.legend(bbox_to_anchor=(0., 1.02, 1., 1.102), loc=3, ncol=3, mode="expand", borderaxespad=1., fontsize=8)
            plt.subplots_adjust(top=0.8)
            plt.suptitle(self.args['title'], fontsize=8)
            plt.savefig(self.args['stdout']+'.'+self.args['format'], dpi=400) 

        else:
            return scatterp


if __name__ == '__main__':

    bd = beta_diversity()
     
    ## print bd.sample_and_coordinates
