#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import utils
import glob
import copy
import os
import subprocess as sp
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import rankdata
from scipy.spatial import distance
from scipy.cluster import hierarchy
import subprocess


## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do heat_map -db metaphlan -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -cm RdYlBu_r -nif 5

## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do heat_map -db pathways -g0 nt:1000 -g1 nsl:5 -g2 c:entropy -nif 5


class feature_heatmap(object):

    dataselect = '/python ../../cmdpy/datasetSelection2.py'
    save_fig = '../Images'
    getfeatnumber = {10:13, 20:24, 30:35, 40:46, 50:57, 60:68} 


    def __init__(self, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python, n_imp_feat, lodo, char_size_hm, cmap, fig_fmt, path='/CM/data/meta/'):
        self.n_imp_feat = n_imp_feat		
        self.char_size_hm = char_size_hm	
        self.lodo = lodo
        self.which_python = which_python	
        self.utils = utils.usefullfuncs(datasets, mixed_taxa=False)
        self.path = path
        self.datasets = datasets
        self.db = db
        self.feat_iden = self.utils.features[db]
        self.algo = algo
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2        
        self.complete_ranking = '../ml/resultofcrossvalidation_ANYonANY_features:metaphlan_experimenttype:standard_rf_grid0:%s_grid1:%s_grid2:%s_FEATURE_RANKING.txt' %(self.grid0, self.grid1, self.grid2)

        self.cmap = cmap
        self.fig_fmt = fig_fmt
        self.segregate_features = True   	
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        prob = self.tests[0].split(':')
        self.bnm = lambda n : n.split('/')[-1]
        self.feat_files = dict([(ds, self.utils.get_transfer_resultnames([ds, ds], self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2)) for ds in self.datasets])
        self.lodo_files = dict([(ds, self.utils.get_lodo_resultnames(ds, self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2)) for ds in self.datasets])
        self.metadata = self.read_metadata().T
        self.phylogenetic_mapper = dict()
        self.most_relevant_features = list()

        self.rename_spp = lambda s_ : s_.replace('unclassified','spp.') if s_.endswith('unclassified') else s_
        self.frames = self.concat_frames(lodo) 	## all the features taken fromt eh random forest result
        self.signs_ = self.signs()
        self.data = self.datamatrix()
        self.data = self.data.loc[self.pos + self.neg]

        self.data_hot = self.data.loc[self.pos].astype('int')
        self.data_cold = self.data.loc[self.neg].astype('int')
        self.data_comparison = self.define_comparison()
        self.data_comparison_hot = self.data_comparison.loc[self.pos].astype('int')
        self.data_comparison_cold = self.data_comparison.loc[self.neg].astype('int')
        self.data_hot.set_index([[k.replace('_',' ') for k in self.data_hot.index.tolist()]], inplace=True)
        self.data_cold.set_index([[k.replace('_',' ') for k in self.data_cold.index.tolist()]], inplace=True)
        self.data_comparison_hot.set_index([[k.replace('_',' ') for k in self.data_comparison_hot.index.tolist()]], inplace=True)
        self.data_comparison_cold.set_index([[k.replace('_',' ') for k in self.data_comparison_cold.index.tolist()]], inplace=True)

        self.linkage = hierarchy.linkage(distance.pdist(self.data.T, metric='braycurtis'), method='average')
        fg,ax_hc = plt.subplots(1,1)
        sns.set_style('white')
        hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
        hc = hierarchy.dendrogram(self.linkage, ax=ax_hc, above_threshold_color='b', orientation='top')
        hierarchy.set_link_color_palette(None)
        new_order = [ds if (not ds in self.utils.data_aliases) else self.utils.data_aliases[ds] for ds in [self.datasets[i] for i in hc['leaves']]]
        self.data_hot = self.data_hot[new_order]
        self.data_cold = self.data_cold[new_order]
        self.annotate = True
        self.heatmap_clustering(new_order, color='warm')



    def select_nlargest(self):
        features = self.data.nlargest(self.n_imp_feat, self.datasets[0]).index.tolist()
        for data in self.datasets[1:]: features += self.data.nlargest(self.n_imp_feat, data).index.tolist()
        return list(set(features))


    def heatmap_clustering(self, sorted_data, color='warm'): 
        vmin_, vmax_ = self.mid_point, 1.
        norm_ = matplotlib.colors.Normalize(self.mid_point, 1)
        revnorm = matplotlib.colors.Normalize(1, self.mid_point)
        fig = plt.figure(figsize=(9,14))
        gs = gridspec.GridSpec(3,3, width_ratios=[10,1,1], height_ratios=[2,10,10])

        self.ax_hclus = plt.subplot(gs[0, 0])
        self.ax_hot = plt.subplot(gs[1, 0])
        self.ax_cold = plt.subplot(gs[2, 0])
        self.ax_comparison_hot = plt.subplot(gs[1, 1])
        self.ax_comparison_cold = plt.subplot(gs[2, 1])
        self.ax_cbar_hot = plt.subplot(gs[1, 2])
        self.ax_cbar_cold = plt.subplot(gs[2, 2])
        self.ax_hclus.axis('off')

        #self.ax_hot.axis('off')
        #self.ax_cold.axis('off')
        #self.ax_cbar_hot.axis('off')
        #self.ax_cbar_cold.axis('off')

        clus = hierarchy.dendrogram(self.linkage, ax=self.ax_hclus, above_threshold_color='b', orientation='top')
        cbar_hot = matplotlib.colorbar.ColorbarBase(self.ax_cbar_hot, cmap='YlOrBr', norm=revnorm, extend='min')
        cbar_cold = matplotlib.colorbar.ColorbarBase(self.ax_cbar_cold, cmap='Blues_r', norm=norm_, extend='max', ticks=range(1,15,3))

        hm_hot = sns.heatmap(self.data_hot, annot=True if self.annotate else False, ax=self.ax_hot, fmt='g'\
                , vmin=1, vmax=self.mid_point, cmap='YlOrBr_r', cbar=False, xticklabels=[], yticklabels=self.data_hot.index.tolist())
        hm_cold = sns.heatmap(self.data_cold, annot=True if self.annotate else False, ax=self.ax_cold, fmt='g'\
                , vmin=1, vmax=self.mid_point, cmap='Blues_r', cbar=False, xticklabels=sorted_data, yticklabels=self.data_cold.index.tolist())
        
        hm_comp_hot = sns.heatmap(self.data_comparison_hot, annot=True if self.annotate else False, ax=self.ax_comparison_hot, fmt='g'\
                , vmin=1, vmax=self.mid_point, cmap='YlOrBr_r', cbar=False, xticklabels=[], yticklabels=[])
        hm_comp_cold = sns.heatmap(self.data_comparison_cold, annot=True if self.annotate else False, ax=self.ax_comparison_cold, fmt='g'\
                , vmin=1, vmax=self.mid_point, cmap='Blues_r', cbar=False, xticklabels=['Global Ranking'], yticklabels=[])


        for lab in self.ax_hot.get_yticklabels(): lab.set_style('italic')
        for lab in self.ax_cold.get_yticklabels(): lab.set_style('italic')
        cbar_hot.set_ticks(range(self.mid_point, 1, -11)) 
        cbar_hot.set_ticklabels(list(map(str, range(1, self.mid_point +1, 11))))
        cbar_cold.set_ticks(range(1, self.mid_point +1, 11)) 

        ##cbar_cold.set_ticklabels()
        ##self.ax_cbar_hot.set_ylabel('Feature Importance Rank', size=16)
        self.ax_cbar_hot.yaxis.set_ticks_position('right')
        self.ax_cbar_cold.yaxis.set_label_position('right')
        plt.subplots_adjust(bottom=0.13, left=0.4)

        pos_hclust = self.ax_hclus.get_position()
        pos_hot = self.ax_hot.get_position()
        pos_cold = self.ax_cold.get_position()
        pos_chot = self.ax_cbar_hot.get_position()
        pos_ccold = self.ax_cbar_cold.get_position()
        pos_comp_hot = self.ax_comparison_hot.get_position()
        pos_comp_cold = self.ax_comparison_cold.get_position()

        pos_hclust2 = [pos_hclust.x0, pos_hclust.y0 - 0.06, pos_hclust.width, pos_hclust.height]
        pos_hot2 = [pos_hot.x0, pos_hot.y0 - 0.03, pos_hot.width, pos_hot.height]
        pos_cold2 = [pos_cold.x0, pos_cold.y0, pos_cold.width, pos_cold.height]
        pos_chot2 = [pos_chot.x0, pos_chot.y0 - 0.03, pos_chot.width, pos_chot.height]
        pos_ccold2 = [pos_ccold.x0, pos_ccold.y0, pos_ccold.width, pos_ccold.height]     
        pos_comp_hot2 = [pos_comp_hot.x0, pos_comp_hot.y0 - 0.03, pos_comp_hot.width, pos_comp_hot.height]
        pos_comp_cold2 = [pos_comp_cold.x0, pos_comp_cold.y0, pos_comp_cold.width, pos_comp_cold.height]

        self.ax_hclus.set_position(pos_hclust2)
        self.ax_hot.set_position(pos_hot2)
        self.ax_cold.set_position(pos_cold2)
        self.ax_cbar_hot.set_position(pos_chot2)
        self.ax_cbar_cold.set_position(pos_ccold2)
        self.ax_comparison_hot.set_position(pos_comp_hot2)
        self.ax_comparison_cold.set_position(pos_comp_cold2)
        plt.setp(self.ax_cold.get_xticklabels(), rotation=70, ha='right')
        plt.setp(self.ax_comparison_cold.get_xticklabels(), rotation=70, ha='right')
        plt.suptitle('Random Forest Feature Ranking')
        if self.db == 'metaphlan': 
            plt.savefig('%s/FeatureHeatmap_%s.%s' %(self.save_fig, color, self.fig_fmt), dpi=600)
        elif self.db == 'pathways':
            plt.savefig('%s/PathwaysHeatmap_%s.%s' %(self.save_fig, color, self.fig_fmt), dpi=600)



    def compute_LODO_vector(self):
        mean, scores, div = [], [], 0
        start = 2 if self.n_imp_feat not in self.getfeatnumber else self.getfeatnumber[self.n_imp_feat]
        for i,ds in enumerate(self.datasets):
            score, coordinate, n = self.utils.lodo_(ds, self.db, self.tests[0], self.algo, self.grid, self.feat, start)
            mean.append(score * n)
            div += n
            scores.append(score)
        scores.append(sum(mean)/float(div))
        return scores



    def compute_CV_vector(self):
        scores = [self.utils.transfer_([ds,ds], self.db, self.tests[0], self.algo, self.grid, self.feat\
		, 2 if self.n_imp_feat not in self.getfeatnumber else self.getfeatnumber[self.n_imp_feat])[0] for i,ds in enumerate(self.datasets)]
        return scores + [self.utils.cross_validation_(self.db, self.tests[0], self.algo, self.grid, self.feat\
		, 2 if self.n_imp_feat not in self.getfeatnumber else self.getfeatnumber[self.n_imp_feat])]



    def define_comparison(self):
        c = self.features_(False, False, line='line', skip=0, filename=True)
        if self.db == 'metaphlan': c = c.set_index([[self.rename_spp(i[3:]) for i in c.index.tolist()]]) 
        c = c.loc[self.pos+self.neg]
        return c



    def datamatrix(self):
        f = self.frames
        for piece in f: self.most_relevant_features += piece.index.tolist()[:self.n_imp_feat]
        self.most_relevant_features = list(set(self.most_relevant_features))
        min_ = np.min([np.min(piece) for piece in f])
        red = [piece.loc[self.most_relevant_features] for piece in f]
        f = red[0]
        for r in red[1:]: f = f.join(r, how='outer')
        f['sum'] = f.sum(axis=1)
        f = f.fillna(min_).astype('int')
        f = f.sort_values('sum', ascending=True, axis=0)
        del f['sum']
        self.neg = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) \
		for ft in f.index.tolist() if self.signs_[ft] < 0.], key=lambda el : el[1], reverse=True)]
        self.pos = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) \
		for ft in f.index.tolist() if self.signs_[ft] > 0.], key=lambda el : el[1], reverse=False)]
        self.mid_point = int(np.ceil(np.median(f.values) / 10.) * 10 )
        f = f.loc[self.pos + self.neg]        
        if self.db == 'metaphlan': 
            self.pos = [self.rename_spp(i[3:]) for i in self.pos]
            self.neg = [self.rename_spp(i[3:]) for i in self.neg]
            f.set_index([[self.rename_spp(i[3:]) for i in f.index.tolist()]], inplace=True)
        f.columns = [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in f.columns.tolist()]
        return f



    def concat_frames(self, lodo): return [self.features_(ds, lodo) for ds in self.datasets]


    def get_annotation_matrix(self, frame): return frame.apply(lambda f : [('+' if n else '') for n in f.tolist()], axis=1)
    def signs(self): return dict([(feat, self.sign_of_a_feature(self.phylogenetic_mapper[feat])) for feat in self.phylogenetic_mapper.keys()])
    def the_greater_is_the_second_in(self, tup): return True if (tup[1]>tup[0]) else False
    def mean_on_dataset(self, ds, s__): return self.subset(ds, self.tests[0].split(':')[1:])[s__].astype('float').mean()
    def mean_on_subset(self, ds, s__, c): return self.subset(ds, c)[s__].astype('float').mean()
    def sign_of_a_feature(self, f): return 1 if (not self.the_greater_is_the_second_in(self.mean_per_class(f))) else -1


    def mean_per_class(self, s__): 
        return tuple(self.metadata[self.metadata[self.tests[0].split(':')[0]].isin([c])][s__].astype('float').mean() for c in self.tests[0].split(':')[1:])

   
    def features_(self, ds, lodo, line='line', skip=0, filename=False):

        if not filename: f = open(self.feat_files[ds]) if not lodo else open(self.lodo_files[ds])
        else: f = open(self.complete_ranking)
        while (not line.startswith('Feature')): line, skip = f.readline(), skip+1
        f.close()
        d = pd.read_csv(((self.feat_files[ds] if not lodo else self.lodo_files[ds]) if not filename else self.complete_ranking), sep='\t', header=None, index_col=0, skiprows=skip)[[1,2]]
        for f in d[1]: self.phylogenetic_mapper[f.split('|')[-1]] = f 

        d[1] = d[1].apply(lambda n : n.split('|')[-1])
        d[2] = d[2].apply(lambda n : float(n))
        if not filename: d.columns = ['feat', str(ds)]
        else: d.columns = ['feat', 'comparison']
        d.set_index('feat', inplace = True)        

        if not filename: d[ds] = [i+1 for i,f in enumerate(d[ds].tolist())]
        else: d['comparison'] = [i+1 for i,f in enumerate(d['comparison'].tolist())]
        return d



    def subset(self, ds, c): 
        sub = self.metadata[self.metadata['dataset_name'].isin([ds]) & self.metadata[self.tests[0].split(':')[0]].isin([c] if isinstance(c, str) else c)]
        sub.reset_index(inplace=True)
        return sub



    def read_metadata(self):
        command = [self.which_python+self.dataselect]+[ds+'.'+self.tests[0] for ds in self.datasets]
        if isinstance(self.db, list):
            for dab in self.db: command += ['--'+self.utils.databases[dab]]
        elif isinstance(self.db, str):
            command += ['--'+self.utils.databases[self.db]]
        command += ['-of ./complete_metadata_'+('_'.join(self.db) if isinstance(self.db, list) else self.db)+'.csv']
        #if os.path.exists('complete_metadata_'+('_'.join(self.db) if isinstance(self.db, list) else self.db)+'.csv'): pass
        #else: sp.call(' '.join(command), shell=True)
        return pd.read_csv('complete_metadata_'+('_'.join(self.db) if isinstance(self.db, list) else self.db)+'.csv', sep='\t', header=None, index_col=0, low_memory=False)        



if __name__=='__main__':
    print 
