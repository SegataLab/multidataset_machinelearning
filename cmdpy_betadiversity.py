#!/usr/bin/env python

#  USAGE EXAMPLES: 

#  python cmdpy_dataset.py CM_periimplantitis.study_condition:peri-implantitis,col:color:crimson,col:shape:o,col:define:periimplantitis CM_periimplantitis.study_condition:control,col:color:blue,col:shape:o,col:define:control --shrink | python cmdpy_betadiversity.py -a mds -sc coords.txt -of pi_plot 

## b) python cmdpy_dataset.py CM_periimplantitis.study_condition:mucositis,col:color:black,col:shape:o,col:define:mucositis CM_periimplantitis.study_condition:peri-implantitis,col:color:blue,col:shape:o,col:define:periimplantitis CM_periimplantitis.study_condition:healthy,col:color:green,col:shape:o,col:define:healthy --taxon g__ | python cmdpy_betadiversity.py -avgd mettiqua_latua_betadiversity -of general_trial_output -z g__

## meaningof of b): call dataset with three classe belonging to CM_perimolnatitis data: periimplantitis, mucositis, heathy: each has a shape, a color and a definition. The selected taxon is genus (g__). No other specification means: metaphlan feature types and metadata attached. This table is passed to betadiversities. The selected operattion is avgd, that means saving on a file a dataframe with, for any sample, the average beta diversity from any class 



import argparse as ap
import sys, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from sklearn import decomposition
from sklearn import preprocessing as prep
from sklearn import manifold
import matplotlib.gridspec as gridspec
from sklearn.metrics import pairwise_distances
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
from scipy import stats



class beta_diversity(object):


    def __init__(self, params=None):

        self.is_a_comb = lambda c, _l_ : bool(c in _l_) 
        self.they_overlap = lambda pos_a1, pos_a2, pos_b1, pos_b2 : bool( min([pos_a2, pos_b2]) - max([pos_a1, pos_b1]) )

        if not isinstance(params, dict):  self.args = self.read_params(sys.argv)
        else: self.args = params

        self.p_value = lambda group_a, group_b : stats.ttest_ind(group_a, group_b, axis=0, equal_var=self.args['welch'])[1]
        self.stdin = pd.read_csv(self.args['stdin'], sep='\t', header=None, index_col=0)

        if self.args['gradient_on']: 
            self.grads = dict([(s, g) for s,g in zip(\
                self.stdin.loc[self.args['sample_id'],:].tolist()\
              , self.stdin.loc[self.args['gradient_on'],:].astype('float').tolist())])

        else: self.grads = None
        self.subjects = dict([(s, sub) for s,sub in zip(self.stdin.loc[self.args['sample_id'],:].tolist(), self.stdin.loc['subjectID', :].tolist())])

        if self.args['alpha_diversity'] or self.args['class_center'] == 'average':

            f = self.edit_(self.stdin, self.args['sample_id'], self.args['feature_identifier'])
            self.data = f.loc[[i for i in f.index.tolist() if (self.args['feature_identifier'] in i)]].T
            self.metadata = f.loc[[i for i in f.index.tolist() if (not self.args['feature_identifier'] in i)]]

            self.sample_in_class, self.atts_of_sample, self.sample_to_class = self.guess_classes(self.metadata)
            self.legend = [ll.split() for ll in list(set([' '.join(l) for l in self.atts_of_sample.values()]))]
            
            if self.args['alpha_diversity']: self.plot_alpha_diversity(self.data, self.metadata)
            if self.args['class_center'] == 'average': 

                if self.args['gradient_on']: 

                    self.key_column = dict([(sample, \
                            float(self.stdin.loc[self.args['gradient_on']\
                          , self.stdin.loc[self.args['sample_id']].isin([sample]).tolist()].tolist()[0])) \
                            for sample in self.metadata.loc[self.args['sample_id'], :].tolist()])

                self.compute_centroids(self.data, self.metadata)

                #if self.args['passing_centroids']:
                #         for i in self.centroids.index.tolist(): print '\t'.join([i]+list(map(str, self.centroids.loc[i])))


                if self.args['save_centroids']: self.save_centroids()  # and (self.args['class_center'] == 'average'): self.save_centroids() 
            #if self.args['save_centroids'] and (not self.args['class_center']):
                #self.compute_centroids(self.data, self.metadata)
                #self.save_centroids()
            #    raise KeyError('Specify if you want class_center computed as AVERAGE OR as CENTROIDS')                 
    

        else:
            self.shapes_types = '.,ov^<>1234sp*hH+xDd|_'
            if self.args['save_beta_diversity']: self.args['boxplot'] = True

            self.coordinates, self.metadata, self.explained_var = self.beta_div()
            self.sample_in_class, self.atts_of_sample, self.sample_to_class = self.guess_classes(self.metadata)

            if not self.args['intra_individual']:
                self.legend = sorted([ll.split() for ll in list(set([' '.join(l) for l in self.atts_of_sample.values()]))], key = lambda w : w[0])
                self.couples_by_class = dict([(k, list(map(list, itertools.combinations(self.sample_in_class[k], 2)))) for k in self.sample_in_class])  
            else:
                self.legend = sorted([ll.split() for ll in list(set([' '.join(l) for l in self.atts_of_sample.values()]))], key = lambda w : w[0])
                self.couples_by_class = dict([(k, [ss for ss in list(map(list, itertools.combinations(self.sample_in_class[k], 2))) \
                  if self.subjects[ss[0]]==self.subjects[ss[1]]]) for k in self.sample_in_class]) ## self.stdin.loc[self.args['sample_id'],:].tolist()])

            self.sample_and_coordinates = self.stack_sample_coordinates() \
                	if not self.args['load_coordinates'] \
                	else pd.read_csv(self.args['load_coordinates'], sep='\t', header=0, index_col=0)

            if self.args['save_coordinates']: self.save_coordinates(self.args['save_coordinates'])
            if self.args['mask']: self.mask()
            if self.args['explain']: self.explain()

            if (self.args['boxplot']) or \
	           (self.args['avg_dist_from_defined_classes']) or \
		   (self.args['save_beta_diversity']) or \
                   (self.args['class_center'] == 'centroids'):

                self.dist_mat = pd.DataFrame(data=self.coordinates, columns=self.metadata.loc['sampleID', :], index=self.metadata.loc['sampleID', :])

                if self.args['save_beta_diversity']: self.dist_mat.to_csv(self.args['save_beta_diversity'], sep='\t', header=True, index=True)
                if self.args['avg_dist_from_defined_classes']: self.avg_dist_from_defined_classes() 
                if self.args['boxplot']: 
                    if self.args['stdout']: self.box_plot()
               
                if self.args['class_center'] == 'centroids': 

                    f = self.edit_(self.stdin, self.args['sample_id'], self.args['feature_identifier'])
                    self.data = f.loc[[i for i in f.index.tolist() if (self.args['feature_identifier'] in i)]].T

                    #print ' qui e pieno di finocchi'
                    self.compute_centroids(self.data, self.metadata)

                    if self.args['save_centroids']: self.save_centroids() 
                #if self.args['save_centroids']: self.save_centroids()
                    #    and (self.args['class_center'] == 'centroids'): 
                    #    self.compute_centroids(self.data, self.metadata)
                #    self.save_centroids()
                #if self.args['save_centroids'] and (self.args['class_center']):
                #    raise KeyError('Specify if you want class_center computed as AVERAGE OR as CENTROIDS')

                    if self.args['gradient_on']: 
                        self.key_column = dict([(sample, \
                           float(self.stdin.loc[self.args['gradient_on'], self.stdin.loc[self.args['sample_id']].isin([sample]).tolist()].tolist()[0])) \
                           for sample in self.metadata.loc[self.args['sample_id'], :].tolist()])

                    #if self.args['passing_centroids']:
                    #     for i in self.centroids.index.tolist(): print '\t'.join([i]+list(map(str, self.centroids.loc[i])))

            else:

                if not (not self.args['no_ordination']): 
                    self.scatter_plot()



    def compute_centroids(self, data, metadata):
        data = data.astype('float')
        self.indexes = dict([(name, number) for number,name in zip(metadata.columns.tolist(),metadata.loc['sampleID', :].tolist())])

        if self.args['class_center'] == 'average':
            self.centroids = (pd.DataFrame(\
                dict([(clss, \
                data.loc[[self.indexes[name] for name in self.sample_in_class[clss]], :].apply(np.mean, axis=0).tolist()) \
                for clss in self.sample_in_class]), index=data.columns.tolist()) ) # if self.args['class_center'] == 'average' \

        else:  
            classes = [clss for clss in self.sample_in_class.keys()]
            tot_dist_from_the_others = pd.DataFrame(dict([(sample, \
                [np.sum(self.dist_mat.loc[sample, self.sample_in_class[clss]].tolist()) \
                for clss in classes]) \
                for sample in self.stdin.loc[self.args['sample_id'], :].tolist()]), \
                    index=classes).T

            centroids = [(tot_dist_from_the_others.index.tolist()[\
                          tot_dist_from_the_others[clss].tolist().index(\
                          np.min(tot_dist_from_the_others[clss].tolist()))])\
                         for clss in classes]

            centroids_first = pd.DataFrame(dict([(sample, \
                data.loc[self.indexes[sample], :].tolist()) \
                for sample in centroids])\
               ,index=data.columns.tolist())

            self.centroids = pd.DataFrame(dict([(\
                self.atts_of_sample[sample][0], \
                data.loc[self.indexes[sample],:].tolist()) \
                for sample in centroids])\
               ,index=data.columns.tolist())

            #print 'THE BIG KAHUNA'
            #print centroids_first, ' LA PRIMA VERSA PER VDERE SE LE CLASSI RITORNANO'             
            #exit(1)
            #print self.centroids, ' LA SECONDA VERA PER VEDERE  CENTROIDI FIGHI DURI COME POCHI PURI'
            


        
    def save_centroids(self): self.centroids.to_csv(self.args['save_centroids'], sep='\t', header=True, index=True)
        


    def avg_dist_from_defined_classes(self):

        centroids = [cent for cent in self.sample_in_class.keys()]

        self.dist_from_defined_classes = pd.DataFrame(dict([(sample, \
            [np.mean(self.dist_mat.loc[sample, self.sample_in_class[cent]].tolist()) \
            for cent in centroids]) \
            for sample in self.stdin.loc[self.args['sample_id'], :].tolist()]), \
          index = centroids).T

        self.dist_from_defined_classes.to_csv(self.args['avg_dist_from_defined_classes'], sep='\t', header=True, index=True)



    def read_params(self, args):

        p = ap.ArgumentParser()

        distances = ['lbraycurtis','sbraycurtis','canberra','chebyshev','correlation','dice','hamming','jaccard','kulsinski'\
	    ,'mahalanobis','matching','minkowski','rogerstanimoto','russellrao','seuclidean','sokalmichener','sokalsneath'\
	    ,'sqeuclidean','yule','cityblock','cosine','euclidean','l1','l2','manhattan','braycurtis','precomputed'] 

        arg = p.add_argument

        arg(	'-if', '--stdin', default=sys.stdin, type=str)
        arg(	'-of', '--stdout', default=None, type=str)
        arg(	'-adiv', '--alpha_diversity', action='store_true')

        arg(	'-sc', '--save_coordinates', default=None, type=str)
        arg(	'-lc', '--load_coordinates', default=None, type=str)
        arg(	'-d', '--distance', default='braycurtis', type=str, choices=distances)
        arg(	'-a', '--algorithm', default='mds', type=str, choices=['mds','pca','nmds','boxp'])
        arg(	'-z', '--feature_identifier', default='s__', type=str, choices=['k__','s__','PWY','UniRef90','g__'])
        arg(	'-si', '--sample_id', default='sampleID', type=str)
        arg(	'-ci', '--classes_id', default='define', type=str)
        arg(	'-m', '--mask', default=[], nargs='+', type=str)

        arg(	'-ll', '--legend_loc', default='lower center', type=str)
        arg(	'-fmt', '--format', default='png', type=str, choices=['svg','png'])
        arg(	'--no_ordination', action='store_false')
        arg(	'--boxplot', action='store_true')
        arg(	'-ex', '--explain', action='store_true')
        arg(	'--annot', action='store_true')

        arg(	'--title', type=str, default='Ordination Plot')
        arg(	'--dot_size', type=float, default=10)
        arg(	'-txt', '--text_on', type=str, default=None)
   
        arg(	'--alpha', type=float, default=1.0)     
        arg(	'--facecolor', type=str, choices=['white', 'whitegrid', 'darkgrid'])
        arg(	'--squared_plot', action='store_true')

        colors = ['RdYlBu', 'plasma', 'inferno', 'winter', 'copper', 'viridis', 'YlGnBu', 'YlGnBu_r', 'cool', 'cool_r']
        arg(	'-cm', '--cmap', default='viridis', choices=colors+[c+'_r' for c in colors])
        arg(	'-cn', '--cbar_title', default='')

        arg(	'-go', '--gradient_on', default=None, type=str, help='must be a column in the data.')
        arg(	'-gf', '--gradient_only_for', default=None, type=str, help='definition of a class in case you want gradient only for that class.')

        arg(	'--intra_individual', action='store_true')
        arg(	'-sb', '--save_beta_diversity', default=None, type=str)
        arg(	'-avgd', '--avg_dist_from_defined_classes', default=None, type=str)
        
        arg(	'-cc', '--class_center', type=str, default=None, choices=['centroids', 'average']) 
        ## classes are defined by 'define' column in the whole tool
        arg(	       '--save_centroids', default=None, type=str)
        #arg(	       '--print_centroids')
        #arg(	       '--passing_centroids', action='store_true', help='flag to print centroids in stdout')

        arg(	'--p_values', type=str, default='above', choices=['above', 'below'])
        arg(	'--p_values_only_for', type=str, nargs='+', default=[])
        arg(	'--welch', action='store_false')

        return vars(p.parse_args())



    def explain(self):
        print ' sample_in_class: a dict with classes as keys and samples for any keys listed as values.' 
        print ' atts_of_sample: a dict with samples as keys and [class, color, shape] as values.'
        print ' sample_to_class: a dict with samples as keys and class for each as value.'
        print ' legend: ', self.legend



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

        if metric == 'precomputed': 
            return data

        elif metric == 'lbraycurtis':

            #print type(data), data.shape, data[0:6]
            ldata = data.apply(lambda col : [( (1.+np.log(float(x)) if float(x) > 0.0 else 0.0) if str(x)[0].isdigit() else x) for x in col], axis=1) #np.matrix([[(np.log(1.0+float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            #print ldata.shape, ldata[0:10]
            return pairwise_distances(ldata, metric='braycurtis')

        elif metric == 'sbraycurtis':
            sdata = np.matrix([[(math.sqrt(float(l)) if float(l) > 0.0 else 0.0) for l in d] for d in data])
            return pairwise_distances(sdata, metric='braycurtis')

        else:
            #print data
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

        ###for c in class_to_samples: print c, class_to_samples[c], ' UUUUUUU' 
        return class_to_samples, sample_to_attributes, sample_to_class            



    def transformation(self, way, f, feat_id, distance):

        data = f.loc[[i for i in f.index.tolist() if (feat_id in i)]].T
        metadata = f.loc[[i for i in f.index.tolist() if (not feat_id in i)]]
        T_func = (self.pca if way == 'pca' else (self.mds if way == 'mds' else self.nmds))\
            if not any([self.args['boxplot']\
                      , self.args['avg_dist_from_defined_classes']\
                      , self.args['save_beta_diversity']\
                      , self.args['class_center']])\
            else 'just_the_matrix'

        if T_func != 'just_the_matrix':
            if T_func == self.pca: 
                transformed, exp_var = T_func(data)
            else: 
                transformed, exp_var = T_func(self.compute_distances(data,distance)), None
            return transformed, metadata, exp_var
        else:
            return self.compute_distances(data,distance), metadata, None
            #return pd.DataFrame(self.compute_distances(data,distance)\
            #     , columns=metadata.loc[self.args['sample_id'], :].tolist()\
            #     , index=metadata.loc[self.args['sample_id'], :].tolist()), metadata, None



    def beta_div(self):
        edt = self.edit_(self.stdin, self.args['sample_id'], self.args['feature_identifier'])
        transformed, metadata, exp_var = self.transformation(self.args['algorithm'], edt, self.args['feature_identifier'], self.args['distance'])
    
        #if self.args['save_beta_diversity'] and self.args['algorithm'] in ['mds','nmds']:
        #    print transfod
        return transformed, metadata, exp_var

 
   # def save_beta_matrix(self):
   #      
   #      self.transformation('mds')
   #
   #      edt = self.edit_(self.stdin, self.args['sample_id'], self.args['feature_identifier'])
   #      beta_div = self.compute_distances(edt, self.args['distance'])
   #      print beta_div


    def stack_sample_coordinates(self, distmat=None):
        if not distmat:
            toreturn = pd.DataFrame({self.args['sample_id']: self.metadata.loc[self.args['sample_id'], :].tolist(), 'x1': self.coordinates[:, 0], 'x2': self.coordinates[:, 1]})
            toreturn.index = toreturn[self.args['sample_id']]
            del toreturn[self.args['sample_id']]
        else:
            toreturn = pd.DataFrame(distmat, columns=self.metadata.loc[self.args['sample_id'], :].tolist(), index=self.metadata.loc[self.args['sample_id'], :].tolist())
        return toreturn


    def save_coordinates(self, coordinates_f): self.sample_and_coordinates.to_csv(coordinates_f, sep='\t', header=True, index=True)

    

    def mask(self):
        self.sample_and_coordinates = self.sample_and_coordinates.drop([i for i in self.sample_and_coordinates.index.tolist() if self.sample_to_class[i] in self.args['mask']])
        for i in self.sample_to_class.keys(): 
            if self.sample_to_class[i] in self.args['mask']: del self.atts_of_sample[i]        



    def scatter_plot(self, ax=None):

        sns.set_style(self.args['facecolor']) 
        cmap = vars(matplotlib.cm)[self.args['cmap']] if self.grads else None
        fig, scatterp, single_scatters, legend_done = False, [], [], False

        if not bool(ax):
            if not self.grads:
                fig, ax = plt.subplots(figsize=(8,(8 if self.args['squared_plot'] else 6))) #(8,6))
            else:
                fig = plt.figure(figsize=(9,6)) 
                gs = gridspec.GridSpec(1,2, width_ratios=[24,1])
                ax = plt.subplot(gs[0,0])
                ax_cbar = plt.subplot(gs[0,1])

            if self.args['algorithm'] == 'mds':
                ax.set_xlim(-0.8, 0.8)
                ax.set_ylim(-0.8, 0.8)
                ax.xaxis.set_ticks(np.arange(-0.6, 0.8, 0.2))
                ax.yaxis.set_ticks(np.arange(-0.6, 0.8, 0.2))

        for c in self.sample_in_class.keys():

            present = [s for s in self.sample_in_class[c] if s in set(self.sample_and_coordinates.index.tolist())]

            #print c, ' is the class'
            #print present, ' are the present ones'
            #print self.sample_in_class
            #print set(self.sample_and_coordinates.index.tolist()), ' quest e la condizione'
            #print self.atts_of_sample
            #print self.atts_of_sample[present[0]], ' atts of sample con present [con zero]'


            if not self.grads:            
                present_sample_frame = self.sample_and_coordinates.loc[present]

                scatterp = sns.regplot(x='x1', y='x2', data=present_sample_frame, ax=ax, scatter=True, fit_reg=False\
                    , scatter_kws={'s': self.args['dot_size']\
                    , 'alpha': self.args['alpha']}  \
                    , label=self.atts_of_sample[present[0]][0]  \
                    , marker=self.atts_of_sample[present[0]][2] \
                    , color=self.atts_of_sample[present[0]][1])

                sps = [spine.set_linewidth(0.5) for spine in ax.spines.itervalues()]

            else:
    
                if not self.args['gradient_only_for']:   

                    single_scatters = [sns.regplot(x='x1', y='x2', ax=ax, scatter=True, fit_reg=False\
                        , data=pd.DataFrame(data=self.sample_and_coordinates.loc[p].values.reshape([1,2])\
                        , columns=['x1','x2'], index=[p]), scatter_kws={'s': self.args['dot_size']} \
                        , label='', marker=self.atts_of_sample[present[0]][2], color=cmap(self.grads[p])) for p in present]
                else:

                    if self.args['gradient_only_for'] != c:
                    #   pass
                        #scatterp += [sns.regplot(x='x1', y='x2', ax=ax, scatter=True, fit_reg=False\
                        #    , data=pd.DataFrame(data=self.sample_and_coordinates.loc[p].values.reshape([1,2])\
                        #    , columns=['x1','x2'], index=[p]), scatter_kws={'s': self.args['dot_size']} \
                        #    , label='', marker='o', color=cmap(self.grads[p])) for p in present]

                    #else:
                        scatterp += [sns.regplot(x='x1', y='x2'\
                            , data=self.sample_and_coordinates.loc[present]\
                            , ax=ax, scatter=True, fit_reg=False\
                            , scatter_kws={'s': self.args['dot_size']\
                            , 'alpha': self.args['alpha']}  \
                            , label=self.atts_of_sample[present[0]][0]  \
                            , marker=self.atts_of_sample[present[0]][2] \
                            , color=self.atts_of_sample[present[0]][1])]


                    else: ## self.args['gradient_only_for']:
                        #if not legend_done:
                        #    plt.legend(bbox_to_anchor=(0., 1.02, 1., 1.102), loc=3, ncol=3, mode="expand", borderaxespad=1., fontsize=8)
                        #    legend_done = True

                        single_scatters += [sns.regplot(x='x1', y='x2', ax=ax, scatter=True, fit_reg=False\
                            , data=pd.DataFrame(data=self.sample_and_coordinates.loc[p].values.reshape([1,2])\
                            , columns=['x1','x2'], index=[p]), scatter_kws={'s': self.args['dot_size']} \
                            , label='', marker='o', color=cmap(self.grads[p])) for p in present]
                    
                sps = [spine.set_linewidth(1.) for spine in ax.spines.itervalues()]
                norm = matplotlib.colors.Normalize(vmin=min(self.grads.values()), vmax=max(self.grads.values()))   ##vmin=0.,vmax=1.)
                cbar = matplotlib.colorbar.ColorbarBase(\
                    ax_cbar, cmap=self.args['cmap'], norm=norm, extend='neither'\
                  , filled=True, drawedges=False, orientation='vertical')

                #cbar.set_ticks(range(-3, 2, 1)) #[0.0, .25, .5, .75, 1.0])
                #cbar.set_ticklabels(list(map(str, range(-3, 2, 1)))) #np.arange(min(self.grads.values()), max(self.grads.values()), 0.5)) ) )   #['.0', '.25', '.5', '.75', '1.0'])
                ax_cbar.yaxis.set_ticks_position('left')
                cbar.outline.set_linewidth(1.)
                ax_cbar.set_ylabel(self.args['cbar_title'] if self.args['cbar_title'] else self.args['gradient_on'], size=10)
                ax_cbar.tick_params(labelsize=10)

            if self.args['annot']:
                for sample,x,y in zip(present_sample_frame.index.tolist(), present_sample_frame['x1'].tolist(), present_sample_frame['x2'].tolist()):	
                    ax.annotate(sample, (float(x) - 0.08, float(y)), size=3)            

        if bool(fig):

            if not self.grads: 
                plt.legend(bbox_to_anchor=(0., 1.02, 1., 1.102), loc=3, ncol=3, mode="expand", borderaxespad=1., fontsize=8)

            #elif self.args['gradient_only_for']: 
            #    patches = [ plt.scatter([],[], marker='o', s=self.args['dot_size'], color=c, label=lab) \
            #        for c,lab in zip([leg[1] for leg in self.legend], [leg[0] for leg in self.legend]) \
            #        if lab != self.args['gradient_only_for']]
            #    ax.legend(handles=[p[0] for p in patches], bbox_to_anchor=(0., 1.02, 1., 1.102), loc=3, ncol=3, mode="expand", borderaxespad=1., fontsize=8)

            plt.subplots_adjust(top=0.8)
            plt.suptitle(self.args['title'], fontsize=10)
            plt.savefig(self.args['stdout']+('_ANNOT' if self.args['annot'] else '')+'.'+self.args['format'], dpi=400) 

            return 'Got.'
        else:

            return scatterp




    def betadiv_statistical_support(self, text):
        print '\n============================================='
        print 'Statistical Support ' + ('about ' + self.args['title'] if self.args['title'] else '')      
        print '=============================================\n'
        for group_a,group_b in self.combinations_of_classes:
            statistic,p_value = stats.ttest_ind(self.groups[group_a], self.groups[group_b], axis=0, equal_var=False)
            print ' - '.join(list((group_a,group_b))), 
            print ' P-PVALUE = ' + str(p_value), 
            print 'sign. TRUE' if (p_value < 0.05) else ''
        print '\n=============================================\n'




    def box_plot(self):
       
        sns.set_style(self.args['facecolor'])  
        fig, ax = plt.subplots(figsize=(8,(8 if self.args['squared_plot'] else 6)))

        #for c in self.couples_by_class: self.couples_by_class[c], '   copro di mille balene'
        #print self.args['p_values_only_for']

        class box_plot_object(object):
            def __init__(self, class_, color_, cat_var_, couple_of_samples, dist_mat):
                self.class_ = class_
                self.color_ = color_
                self.cat_var_ = cat_var_
                self.beta_diversity = dist_mat.get_value(couple_of_samples[0], couple_of_samples[1])

        data = pd.DataFrame(\
                [[ob.class_, ob.color_, ob.cat_var_, ob.beta_diversity]   \
                  for ob in [\
                  box_plot_object(cl,co,ct,cp,self.dist_mat) for cl,co,ct,cp in zip(\
                    list(itertools.chain.from_iterable(\
                       [[c[0] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[c[1] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[c[2] for i in range(len(self.couples_by_class[c[0]]))] for c in self.legend if len(self.couples_by_class[c[0]]) ]))\
                  , list(itertools.chain.from_iterable(\
                       [[couple for couple in self.couples_by_class[c[0]]] for c in self.legend if len(self.couples_by_class[c[0]])]))\
                  ) ]], columns=['', 'color', 'group_by', 'Beta-Diversity']) 

        ax_w = sns.swarmplot(data=data, x='', y='Beta-Diversity', ax=ax\
          , hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by', dodge=True, s=2, color='black', alpha=1.)
        ax_b = sns.boxplot(data=data, x='', y='Beta-Diversity', linewidth=1., ax=ax\
          , hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by', palette=dict([(c[0], c[1]) for c in self.legend]))

        data.columns = ['class', 'color', 'group_by', 'Beta-Diversity']
        self.groups = dict([(cl, data.loc[data['class'].isin([cl]), 'Beta-Diversity'].tolist()) for cl in list(set(data['class'].tolist()))])
        self.combinations_of_classes = list(map(tuple, itertools.combinations(list(set(data['class'].tolist())), 2)))
               
        self.distributions = []
        for class_ in data['class'].tolist(): 
            if class_ not in self.distributions: self.distributions.append(class_) 

        self.indices = dict([(d,e) for e,d in enumerate(self.distributions)])
        self.couple_levels = dict()
        drop = data['Beta-Diversity'].max()*0.02

        if self.args['p_values'] == 'above': H = data['Beta-Diversity'].max() + drop #* 2
        else: H = data['Beta-Diversity'].min() - drop ##* 2

        level = H
        #cover = dict([(dat,level) for dat in self.distributions])
        distances = dict([(couple, (max([self.indices[couple[0]], self.indices[couple[1]]]) \
                                  - min([self.indices[couple[0]], self.indices[couple[1]]]))) \
                                    for couple in self.combinations_of_classes])
        self.done_ = dict([(couple, False) for couple in self.combinations_of_classes])
        ##prevs = {} #min(distances.values())

        distance = 1
        for i in range(1, len(self.distributions)): ##list(itertools.chain.from_iterable(range(prev, len(self.distributions)), ):
            cup = self.distributions[i], self.distributions[i-1]   
            if self.is_a_comb(cup, self.combinations_of_classes): print 'coppia OK'
            else:
                cup = (cup[1],cup[0])
                if self.is_a_comb(cup, self.combinations_of_classes): print 'reversed coppia OK'
            sign = self.p_value_on_plot(cup, drop, level)
        
        level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2)
        #print distances, '  guarda questo'

        for distance in range(2, max(distances.values())+1):
            for any_comb in self.combinations_of_classes:
                if distances[any_comb] == distance: 
                   sign =  self.p_value_on_plot(any_comb, drop, level)  
                   if sign: 
                       level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2 )   

        plt.setp(ax.get_xticklabels(), rotation=38, ha='right')
        plt.subplots_adjust(bottom=.3) 
        plt.suptitle(self.args['title'], fontsize=10)
        plt.savefig(self.args['stdout']+'.'+self.args['format'], dpi=400)

        self.betadiv_statistical_support(self.args['text_on'])


    def p_value_on_plot(self, cup, drop, level, just_zerotwo=False):
        if not self.done_[cup]:  
            self.done_[cup] = True
            p_ = self.p_value(self.groups[cup[0]], self.groups[cup[1]])    
            if p_ < 0.05:
                plt.plot(\
                   [self.indices[cup[0]] + ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[0]] + ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[1]] - ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))\
                  , self.indices[cup[1]] - ((drop if not just_zerotwo else just_zerotwo) if self.indices[cup[0]] < self.indices[cup[1]] else -(drop if not just_zerotwo else just_zerotwo))]\
                 , [level, level + drop*0.5, level + drop*0.5, level] if self.args['p_values'] == 'above' \
                    else [level, level - drop*0.5, level - drop*0.5, level], lw=1., color='k')
                
                plt.text((self.indices[cup[0]]+self.indices[cup[1]])*0.5, (level if self.args['p_values'] == 'below' else level - drop)\
                        ,'p. '+str(p_), ha='center', va='bottom', color='k', fontsize=5)
                return p_
        return False



    def plot_alpha_diversity(self, f, meta):

        sns.set_style(self.args['facecolor'])  
        #fig, ax = plt.subplots(figsize=(8,(8 if self.args['squared_plot'] else 6)))

        samples = dict([(e+1,n) for e,n in enumerate(meta.loc['sampleID', :].tolist())])

        self.simple_richness = dict([(samples[i], np.count_nonzero(f.loc[i].astype(float).tolist())) for i in f.index.tolist()])
        self.log_richness = dict([(k,np.log(v)) for k,v in self.simple_richness.items()])
        self.shannon_richness = dict([(samples[i], -np.sum(np.array(f.loc[i].astype(float).tolist(), dtype=np.float64) * 0.01\
          * np.nan_to_num(np.log(np.array(f.loc[i].astype(float).tolist(), dtype=np.float64) * 0.01)))) for i in f.index.tolist()])
        self.ginisimpson_richness = dict([(samples[i], 1.-np.sum(np.power(np.array(f.loc[i].astype(float).tolist(), dtype=np.float64) * 0.01\
	  , 2))) for i in f.index.tolist()])
					
        data = pd.DataFrame( [ (self.simple_richness[sam]\
                   , self.log_richness[sam]\
                   , self.shannon_richness[sam]\
                   , self.ginisimpson_richness[sam]\
                   , self.atts_of_sample[sam][2]\
                   , self.atts_of_sample[sam][0])  for sam in samples.values()]\
             , columns = [ 'Richness', 'Log Richness', 'Shannon Richness', 'Gini-Simpson Richness', 'group_by', '' ] )
        
        self.combinations_of_classes = list(map(tuple, itertools.combinations(list(set(data[''].tolist())), 2)))

        #print self.combinations_of_classes

        self.distributions = []
        for class_ in data[''].tolist():
            if class_ not in self.distributions: self.distributions.append(class_)
        #print self.distributions
        
        self.indices = dict([(d,e) for e,d in enumerate(self.distributions)])
        distances = dict([(couple, (max([self.indices[couple[0]], self.indices[couple[1]]]) \
                       - min([self.indices[couple[0]], self.indices[couple[1]]]))) \
                       for couple in self.combinations_of_classes]) 

        #for k in self.indices: print k, self.indices[k]


        self.done_ = dict([(couple, False) for couple in self.combinations_of_classes])

        for type_of_richness in ['Richness', 'Log Richness', 'Shannon Richness', 'Gini-Simpson Richness']:
            data_ = data[[type_of_richness, 'group_by', '']]
            self.groups = dict([(cl, data.loc[data[''].isin([cl]), type_of_richness].tolist()) for cl in self.distributions])         

            print 

            #fig = plt.figure(figsize=(8,6))

            fig, ax = plt.subplots(figsize=(8,(8 if self.args['squared_plot'] else 6)))

            ax_w = sns.swarmplot(data=data_, x='', y=type_of_richness, hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by'\
	     , dodge=True, s=4, color='black', alpha=0.7, ax=ax)
            ax_b = sns.boxplot(data=data_, x='', y=type_of_richness, hue=None if len(list(set(data['group_by'].tolist())))==1 else 'group_by'\
	     , palette=dict([(c[0], c[1]) for c in self.legend]), linewidth=1., ax=ax)

            drop = data[type_of_richness].max()*0.02
            if self.args['p_values'] == 'above': H = data[type_of_richness].max() + drop
            else: H = data[type_of_richness].min() - drop
            level = H   

            distance = 1             
            for i in range(1, len(self.distributions)): ##list(itertools.chain.from_iterable(range(prev, len(self.distributions)), ):
                cup = self.distributions[i], self.distributions[i-1]
                if self.is_a_comb(cup, self.combinations_of_classes): print 'coppia OK'
                else:
                    cup = (cup[1],cup[0])
                if self.is_a_comb(cup, self.combinations_of_classes): print 'reversed coppia OK'
                sign = self.p_value_on_plot(cup, drop, level, 0.01)
                #print sign, ' hko appen trovato ', type_of_richness

            level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2)
            for distance in range(2, max(distances.values())+1):
                for any_comb in self.combinations_of_classes:
                    if distances[any_comb] == distance:
                       sign =  self.p_value_on_plot(any_comb, drop, level, 0.01)
                       ###print sign, ' hko appen trovato ', type_of_richness
                       if sign:
                           level = (level + drop * 2) if self.args['p_values'] == 'above' else (level - drop * 2 )
       
            plt.setp(ax.get_xticklabels(), rotation=38, ha='right')
            plt.subplots_adjust(bottom=.3)
            plt.suptitle(self.args['title'], fontsize=10)
            plt.savefig(self.args['stdout']+type_of_richness.replace(' ','')+'.'+self.args['format'], dpi=400)

        data.columns = ['Richness', 'Log Richness', 'Shannon Richness', 'Gini-Simpson Richness', 'group_by', 'class']
        self.groups = dict([(class_, dict([(type_, np.array([given_dict[s] for s in samples.values() \
          if s in self.sample_in_class[class_]], dtype=np.float64)) for type_,given_dict in zip(\
            ['Richness', 'Log Richness', 'Shannon Richness', 'Gini-Simpson Richness']\
              , [self.simple_richness, self.log_richness, self.shannon_richness, self.ginisimpson_richness])]))\
                for class_ in list(set(data['class'].tolist()))])

        self.alphadiv_statistical_support(self.args['text_on']) 



    def alphadiv_statistical_support(self, text):

        for richness in ['Richness', 'Log Richness', 'Shannon Richness', 'Gini-Simpson Richness']:

            print '\n==============================================================================' 
            print 'Statistical Support for Alpha-Diversity'  + self.args['title'] if self.args['title'] else ''
            print '\tINDEX = %s' %richness 
            print '================================================================================\n'

            for group_a,group_b in self.combinations_of_classes:
                statistic,p_value = stats.ttest_ind(self.groups[group_a][richness], self.groups[group_b][richness], axis=0, equal_var=False)                
                #print self.groups[group_a][richness], p_value
                #print self.groups[group_b][richness], p_value
                print ' - '.join(list((group_a,group_b))),
                print ' P-PVALUE = ' + str(p_value),
                print 'sign. TRUE' if (p_value < 0.05) else ''
                
            print '******************************************'




if __name__ == '__main__':

    bd = beta_diversity()     
    ## print bd.sample_and_coordinates
