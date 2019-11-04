#!/usr/bin/env python


import pandas as pd
import argparse as ap
import numpy as np
import sys
import itertools
from cmdpy_betadiversity import beta_diversity
from utils import project_color_code 
from scipy.spatial.distance import cdist as pairwise
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'



class centroids(object):


    def read_params(self):
        par = ap.ArgumentParser()
        arg = par.add_argument

        arg('-if', '--stdin', default=sys.stdin, type=str)
        arg('-of', '--stdout', default='correlation_plot' )

        arg('-d', '--distance', default='braycurtis' )
        arg('--feature_identifier', default='s__' )
        arg('--sample_id', default='sampleID' )
        arg('--classes_id', default='define' )

        arg('--format', default='png' )
        arg('--cmap', default='plasma' )
        arg('--dot_size', default=60 )
        arg('--title', default='correlation between index and distance from centroids' )
        arg('--facecolor', default='white')

        arg('-sc', '--save_centroids', default=None, type=str)        
        arg('-lc', '--load_centroids', default=None, type=str)
        arg('-f', '--from', nargs=1, default=None, type=str, \
           help='The class that have to be consider for the Figure.')

        arg('-bf', '--beta_file', default=None, type=str)
        arg('-go', '--gradient_on', default='Bug-Complex-Abundance', type=str)
        arg('--intra_individual', action='store_true')
        arg('--try_any_metric', action='store_true')

        arg('-cc', '--class_center', choices=['centroids', 'average'], default='centroids', type=str)
        arg('-cb', '--color_by', choices=['gradient', 'classes'], default='classes')

        arg(	'--despine', action='store_true')
        return vars(par.parse_args())

    

    def __init__(self):

        self.args = self.read_params()

        null_params = ['save_coordinates', 'load_coordinates', 'explain', 'annot', 'text_on', 'alpha_diversity'\
                     , 'squared_plot', 'cbar_title', 'gradient_only_for', 'boxplot'\
                     , 'intra_individual', 'save_beta_diversity', 'avg_dist_from_defined_classes']

        copy_params = ['stdin', 'stdout', 'distance', 'feature_identifier', 'sample_id', 'save_centroids'\
                     , 'classes_id', 'format', 'facecolor', 'cmap', 'dot_size', 'title', 'gradient_on', 'class_center']

        params = {'no_ordination': True, 'legend_loc': 'lower_center'\
                , 'p_values': 'above', 'algorithm': 'mds'}

        for f in null_params: params[f] = False
        for f in copy_params: params[f] = self.args[f] 
        for f in ['mask', 'p_values_only_for']: params[f] = []

        self.betadiv_object = beta_diversity(params) 
        #print self.betadiv_object.args

        self.centroids = (self.betadiv_object.centroids \
            if not self.args['load_centroids'] \
            else pd.read_csv(self.args['load_centroids'], \
            sep='\t', header=0, index_col=0)).astype('float')

        self.data = self.betadiv_object.data.astype('float')
        self.metadata = self.betadiv_object.metadata
        self.indexes = self.betadiv_object.indexes

        self.sample_in_class = self.betadiv_object.sample_in_class
        self.atts_of_sample = self.betadiv_object.atts_of_sample
        self.sample_to_class = self.betadiv_object.sample_to_class

        self.external_index = self.betadiv_object.key_column        


        if not self.args['intra_individual']: 
              self.couples_by_class = dict([(k, list(map(list, itertools.combinations(self.sample_in_class[k], 2)))) for k in self.sample_in_class])
        else: self.couples_by_class = dict([(k, [ss for ss in list(map(list, itertools.combinations(self.sample_in_class[k], 2))) \
              if self.subjects[ss[0]]==self.subjects[ss[1]]]) for k in self.sample_in_class])

        if self.args['beta_file']: 
            self.beta_matrix = pd.read_csv(self.args['beta_file'], sep='\t', header=0, index_col=0)
            self.compute_average_distance()

        #if self.args['between']: 
        #    self.linear_distance_from_groups()

        if not self.args['try_any_metric']:
            self.compute_distance_from_centroid(self.args['distance'])
            self.linear_distance_from_groups()

        else:

            for distance in ['lbraycurtis','sbraycurtis','canberra','chebyshev','correlation','dice','hamming','jaccard','kulsinski'\
            ,'mahalanobis','matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener','sokalsneath'\
            ,'sqeuclidean','yule','cityblock','cosine','euclidean','l1','l2','manhattan','braycurtis','precomputed']:

                print '====================================='
                print '=================   ================='
                print 'DISTANCE: ', distance

                try: 
                    self.compute_distance_from_centroid(distance)
                    self.linear_distance_from_groups()
                except:
                    print 'NOT COMPUTABLE go on...'

                print '====================================='
                print '====================================='
                 




    def compute_average_distance(self):

        self.average_betadiv = dict([(clss, np.mean(\
                  [self.beta_matrix.get_value(couple[0],couple[1]) \
                   for couple in self.couples_by_class[clss]])) \
                      for clss in (self.centroids.columns.tolist() \
                          if len(self.args['from'])==0 \
                          else self.args['from'])])   



    def compute_distance_from_centroid(self, metric):

        observations = self.data.values

        try:
            self.distances_from_centroids = dict([(sample\
                , dict([(clss, pairwise(  \
                      [observations[int(self.indexes[sample])-1]]\
                     ,[self.centroids[clss].tolist()]\
                     , metric=metric)[0][0])\
                  for clss in (self.centroids.columns.tolist() \
                      if len(self.args['from'])==0 \
                      else self.args['from'])])) \
                for sample in self.metadata.loc[self.args['sample_id'], :].tolist()])

        except ValueError:
            raise KeyError('You are probably trying to compute the distance'
                           'between samples of classes that wre enot present'
                           'when the centroid have been computed. For safety,'
                           'try re-compute centroids and see wtaht happens.')
    


    #def sign(self, a): return -1 if a[0] > a[1] else 1 



    def linear_distance_from_groups(self):

        #self.sign(self.proximity_to_centroids[sample]) * 
        #

        self.proximity_to_centroids = dict([(sample, [(x) \
            for x in [self.distances_from_centroids[sample][clss] \
                for clss in (self.centroids.columns.tolist()\
                if len(self.args['from'])==0 else self.args['from'])]]) \
            for sample in self.metadata.loc[self.args['sample_id'], :].tolist()])

        #for s in self.proximity_to_centroids: print s, self.proximity_to_centroids[s]

        self.position_wrt_centroids = dict([(sample, self.proximity_to_centroids[sample][0]) \
            for sample in self.metadata.loc[self.args['sample_id'], :].tolist()])

        #for s in self.position_wrt_centroids: print s, self.position_wrt_centroid[s]

        #exit(1)

        #self.position_wrt_centroids = dict([(sample, \
        #    self.proximity_to_centroids[sample]) \
        #        for y in self.proximity_to_centroids[sample]] ) 
        #        for sample in self.metadata.loc[self.args['sample_id'], :].tolist()]) 

        #print 'il nuoo e questo : ', self.position_wrt_centroids

        #print self.metadata

        self.plot_data = pd.DataFrame(dict([(k,v) for k,v in zip([self.args['sample_id'], \
              'proximity_from_' + self.args['class_center'], self.args['gradient_on'],\
              'color', 'shape', 'define'], \
              [[sample for sample in self.metadata.loc[self.args['sample_id'], :].tolist()], \
               [self.position_wrt_centroids[sample] for sample in self.metadata.loc[self.args['sample_id'], :].tolist()], \
               [self.external_index[sample] for sample in self.metadata.loc[self.args['sample_id'], :].tolist()], \
                self.metadata.loc['color',:].tolist(), \
                self.metadata.loc['shape',:].tolist(), \
                self.metadata.loc['define',:].tolist()])]))

        array = self.plot_data['proximity_from_' + self.args['class_center']].tolist()
        self.plot_data['proximity_from_' + self.args['class_center']] =  (np.array(array, dtype=np.float64))##/np.sum(array) ) #* self.average_betadiv.values()[0]
        #- np.min(array)) / (np.max(array) - np.min(array))

        self.pearson = stats.pearsonr(self.plot_data['proximity_from_' + self.args['class_center']].tolist(), self.plot_data[self.args['gradient_on']].tolist())
        self.spearman = stats.spearmanr(self.plot_data['proximity_from_' + self.args['class_center']].tolist(), self.plot_data[self.args['gradient_on']].tolist())
         
        #print np.corrcoef(x=np.array([self.plot_data['proximity_from_' + self.args['class_center']].tolist()\
        #                            , self.plot_data[self.args['gradient_on']].tolist()], dtype=np.float64))
        #print stats.pearsonr(self.plot_data['proximity_from_' + self.args['class_center']].tolist(), self.plot_data[self.args['gradient_on']].tolist())
        #print stats.spearmanr(self.plot_data['proximity_from_' + self.args['class_center']].tolist(), self.plot_data[self.args['gradient_on']].tolist())

        #print self.average_betadiv, ' EHEMMEMMEM'

        #exit(1)



    def plot_correlation(self):
        
        sns.set(style=self.args['facecolor']) 
        cmap = vars(matplotlib.cm)[self.args['cmap']] if self.args['color_by']=='gradient' else None


        fig,ax = plt.subplots(1,1, figsize=(6,6))
        #sns.regplot(data=self.plot_data, \
        #            y='proximity_from_' + self.args['class_center'], \
        #            x=self.args['gradient_on'], \
	#            scatter=False, fit_reg=True, ci=95, \
        #            n_boot=1000, color='dodgerblue', \
        #            ax=ax)
        
        if not cmap:
            for clss in set(self.plot_data['define'].tolist()):

                plot_data = self.plot_data[self.plot_data['define'] == clss]
                #print 'data now selected: ', plot_data

                sns.regplot(x=self.args['gradient_on'], y='proximity_from_' + self.args['class_center']\
           	    	 , data=plot_data\
          		 , ax=ax, scatter=True, fit_reg=False\
          		 , scatter_kws={'s': self.args['dot_size']\
          		 , 'alpha': 1.}  \
          		 , label = plot_data['define'].tolist()[0] \
         		 , marker = plot_data['shape'].tolist()[0] \
         		 , color = plot_data['color'].tolist()[0] )

            plt.legend(bbox_to_anchor=(0., 1.02, 1., 1.102), loc=3, ncol=3, mode="expand", borderaxespad=1., fontsize=8) 

        else:
            for i,s in enumerate(self.metadata.loc[self.args['sample_id'], :].tolist()):
                sns.regplot(x=self.args['gradient_on'], y='proximity_from_' + self.args['class_center'], ax=ax, scatter=True, fit_reg=False\
                      , data=pd.DataFrame(data={'proximity_from_' + self.args['class_center']: \
                        self.plot_data['proximity_from_' + self.args['class_center']].tolist()[i]\
                      , self.args['gradient_on']: self.plot_data[self.args['gradient_on']].tolist()[i]}\
                      , columns=['proximity_from_' + self.args['class_center'], self.args['gradient_on']], index=[i])\
                      , label='', marker='o', color=cmap(self.plot_data[self.args['gradient_on']].tolist()[i]), scatter_kws={'alpha': 1.})


        sns.regplot(data=self.plot_data, \
            y='proximity_from_' + self.args['class_center'], \
            x=self.args['gradient_on'], \
            scatter=False, fit_reg=True, ci=85, \
            n_boot=1000, color='goldenrod', scatter_kws={'alpha': 0.6},\
            ax=ax)

        ax.annotate(('Pearson: ' + str(self.pearson[0]) + '*\n' + 'Spearman: '+str(self.spearman[0])+'*'), (.8, 1.), size=8)

        plt.suptitle(self.args['title'], fontsize=10)
        sps = [spine.set_linewidth(0.5) for spine in ax.spines.itervalues()]
        if self.args['despine']: sns.despine(right=True, top=True)

        plt.subplots_adjust(top=0.8)
        plt.savefig(self.args['stdout'] + '.' + self.args['format'], dpi=400)




if __name__ == '__main__':
   c = centroids()
   c.plot_correlation()
   print 'Tutto a Posto'

   #print c.average_betadiv
   #print c.centroids
   ##print c.distances_from_centroids
   #print c.plot_data

   #for k in c.position_wrt_centroids: print k, c.position_wrt_centroids[k], c.external_index[k]
   #print ' there are %i samples', len(c.distances_from_centroids.keys())
   #exit(0)
