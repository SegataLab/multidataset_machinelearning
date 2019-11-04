#!/usr/bin/env python



import argparse
import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
from scipy import stats



class plot_of_distance_from_centroid(object):



    def __init__(self):
     
        self.args = self.read_args(sys.argv)
        print self.args

        ### this file is passed through cmdpy_dataset plus the score option
        self.stdin = pd.read_csv(self.args['stdin'], sep='\t', header=None, index_col=0)
        self.avg_distance_from_each_class = pd.read_csv(self.args['avg_dist_from_centroids'], sep='\t', header=0, index_col=0).astype('float')

        self.avg_distance_from_each_class[self.args['index_column']] = list(map(float, [f for _,f in \
					sorted(zip(	\
						   [self.avg_distance_from_each_class.index.tolist().index(x) \
						    for x in self.avg_distance_from_each_class.index.tolist()]	
						  , self.stdin.loc[self.args['index_column'], :]))]))

        
        


        #if self.args['touched_distance']

        #self.touch_distance()

        print self.avg_distance_from_each_class, ' weeeei'    


    
    



    def read_args(self, args):

        p = argparse.ArgumentParser()
        arg = p.add_argument
        arg('-if', '--stdin', default=sys.stdin, type=str)
        arg('-avgd', '--avg_dist_from_centroids', type=str)
        arg('-ic', '--index_column', type=str, default='Bug-Complex-Abundance')
        
        arg('-c', '--centroids', nargs=2, type=str)
        arg('-tc', '--touched_distance', default='similarity_from_groups', type=str)
        arg('-of', '--output_file', type=str)
        arg('-fmt', '--format', default='png')

        #print vars(p.parse_args()), '  da dentro '

        return vars(p.parse_args())




    def touch_distance(self):

        self.avg_distance_from_each_class[self.args['touched_distance']] = \
			[\
			  ((-1 if self.avg_distance_from_each_class.loc[sample, self.args['centroids'][0]]\
			   > self.avg_distance_from_each_class.loc[sample, self.args['centroids'][1]] else 1) \
  			   * max(\
				(self.avg_distance_from_each_class.loc[sample, self.args['centroids'][0]]/\
				 self.avg_distance_from_each_class.loc[sample, self.args['centroids'][1]])\
                            , (self.avg_distance_from_each_class.loc[sample, self.args['centroids'][1]]/\
                              self.avg_distance_from_each_class.loc[sample, self.args['centroids'][0]])))\
			  for sample in self.avg_distance_from_each_class.index.tolist()]
        



    def plot_correlation(self):

        sns.set(style="white", color_codes=True)
        fig,ax = plt.subplots(1,1, figsize=(6,6))
      
        print self.avg_distance_from_each_class


        #single_scatters = [sns.regplot(x=self.args['touched_distance'], y=self.args['index_column'], ax=ax, scatter=True, fit_reg=False\
        #                , data=self.avg_distance_from_each_class\
        #                , columns=['x1','x2'], index=[p]), scatter_kws={'s': self.args['dot_size']} \
        #                , label='', marker=self.atts_of_sample[present[0]][2], color=cmap(self.grads[p])) for p in present]

        ax = sns.jointplot(data = self.avg_distance_from_each_class, x = self.args['touched_distance'], y = self.args['index_column'], kind = 'kde')
        plt.savefig(self.args['output_file'] + '.' + self.args['format'], dpi=400)

        

    
    
if __name__ == '__main__':

    p = plot_of_distance_from_centroid()
    #p.plot_correlation()

    #print p.args, ' dal main di fuori..'


