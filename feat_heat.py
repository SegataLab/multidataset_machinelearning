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


## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do heat_map -db metaphlan -g0 nt:1000 -g1 nsl:5 -g2 c:entropy -cm RdYlBu_r -nif 5

## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do heat_map -db pathways -g0 nt:1000 -g1 nsl:5 -g2 c:entropy -nif 5


class feature_heatmap(object):

    dataselect = '/python ../../cmdpy/datasetSelection2.py'
    save_fig = '../Images'
    getfeatnumber = {10:13, 20:24, 30:35, 40:46, 50:57, 60:68} 


    def __init__(self, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python, n_imp_feat, lodo, char_size_hm, cmap, fig_fmt, path='/CM/data/meta/'):
        self.n_imp_feat = n_imp_feat		## num of feature to look at for each
        self.char_size_hm = char_size_hm	## char size on the heatmap
        self.lodo = lodo			## a boolean
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
        self.min_rank = 15.
        self.cmap = cmap
        self.fig_fmt = fig_fmt
        self.segregate_features = True   	## no idea
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        prob = self.tests[0].split(':')
        self.bnm = lambda n : n.split('/')[-1]	# basename  ### the following two collecting names of the files for any datasets
        self.feat_files = dict([(ds, self.utils.get_transfer_resultnames([ds, ds], self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2)) for ds in self.datasets])
        self.lodo_files = dict([(ds, self.utils.get_lodo_resultnames(ds, self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2)) for ds in self.datasets])
        self.metadata = self.read_metadata().T
        self.phylogenetic_mapper = dict()
        self.most_relevant_features = list()
        self.frames = self.concat_frames(lodo) 	## all the features taken fromt eh random forest result

        #print self.frames
        #print ' sono al punto di uscita previsto'
        #exit(1)

        #self.signs_ = self.signs()
        #self.frames['total'] = self.frames.sum(axis=1)  ### features sorted from the most releveant descending
        #self.frames = self.frames.sort_values(by='total', axis=0, ascending=False)

        #print '\t' + '\t'.join([c for c in self.frames.columns.tolist()])
        #for i in self.frames.index.tolist()[:20]: print '\t'.join([i]+list(map(str, self.frames.loc[i])))
        #print 
        #print 
        #print self.frames['total']
        #exit(1)

        self.data = self.datamatrix()

        exit(1)

        self.data.columns = [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.data.columns.tolist()]
        self.data = self.data.loc[self.pos + self.neg]

        self.data_hot = self.data.loc[self.pos]
        self.data_cold = self.data.loc[self.neg]



        self.linkage = hierarchy.linkage(distance.pdist(self.data.T, metric='braycurtis'), method='average')
        fg,ax_hc = plt.subplots(1,1)
        sns.set_style('white')
        hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
        hc = hierarchy.dendrogram(self.linkage, ax=ax_hc, above_threshold_color='b', orientation='top')
        hierarchy.set_link_color_palette(None)  # reset to default after use
        new_order = [ds if (not ds in self.utils.data_aliases) else self.utils.data_aliases[ds] for ds in [self.datasets[i] for i in hc['leaves']]]
        self.data_hot = self.data_hot[new_order]
        self.data_cold = self.data_cold[new_order]
        self.data_hot = self.data_hot.loc[self.pos]
        self.data_cold = self.data_cold.loc[self.neg]
        self.annotate = True
        self.heatmap_clustering(new_order, color='warm')



    def select_nlargest(self):
        features = self.data.nlargest(self.n_imp_feat, self.datasets[0]).index.tolist()
        for data in self.datasets[1:]: features += self.data.nlargest(self.n_imp_feat, data).index.tolist()
        return list(set(features))


    def heatmap_clustering(self, sorted_data, color='warm'): 
        vmin_, vmax_ = self.min_rank, 1.
        norm_ = matplotlib.colors.Normalize(self.min_rank, 1)
        revnorm = matplotlib.colors.Normalize(1, self.min_rank)

        fig = plt.figure(figsize=(9,14))
        gs = gridspec.GridSpec(3,2, width_ratios=[10,1], height_ratios=[2,10,10])
        self.ax_hclus = plt.subplot(gs[0, 0])
        self.ax_hot = plt.subplot(gs[1, 0])
        self.ax_cold = plt.subplot(gs[2, 0])
        self.ax_cbar_hot = plt.subplot(gs[1, 1])
        self.ax_cbar_cold = plt.subplot(gs[2, 1])
        self.ax_hclus.axis('off')
        #self.ax_hot.axis('off')
        #self.ax_cold.axis('off')
        #self.ax_cbar_hot.axis('off')
        #self.ax_cbar_cold.axis('off')
        self.data_hot.index = [k.replace('_',' ') for k in self.data_hot.index.tolist()]
        self.data_cold.index = [k.replace('_',' ') for k in self.data_cold.index.tolist()]

        clus = hierarchy.dendrogram(self.linkage, ax=self.ax_hclus, above_threshold_color='b', orientation='top')
        cbar_hot = matplotlib.colorbar.ColorbarBase(self.ax_cbar_hot, cmap='YlOrBr', norm=revnorm, extend='min')
        cbar_cold = matplotlib.colorbar.ColorbarBase(self.ax_cbar_cold, cmap='Blues_r', norm=norm_, extend='max', ticks=range(1,15,3))
        hm_hot = sns.heatmap(self.data_hot, annot=self.data_hot.values if self.annotate else False, ax=self.ax_hot\
		, vmin=1, vmax=15, cmap='YlOrBr_r', cbar=False, xticklabels=[], yticklabels=self.data_hot.index.tolist())
        hm_cold = sns.heatmap(self.data_cold, annot=self.data_cold.values if self.annotate else False, ax=self.ax_cold\
		, vmin=1, vmax=15, cmap='Blues_r', cbar=False, xticklabels=sorted_data, yticklabels=self.data_cold.index.tolist())
        for lab in self.ax_hot.get_yticklabels(): lab.set_style('italic')
        for lab in self.ax_cold.get_yticklabels(): lab.set_style('italic')
        cbar_hot.set_ticks(range(15,1,-3))
        #cbar_cold.set_ticks(range(1,15,3))
        cbar_hot.set_ticklabels(map(str, range(1,15,3)))
        #cbar_cold.set_ticklabels(map(str, range(1,15,3)))
        ##self.ax_cbar_hot.set_ylabel('Feature Importance Rank', size=16)
        self.ax_cbar_hot.yaxis.set_ticks_position('right')
        self.ax_cbar_cold.yaxis.set_label_position('right')

        plt.subplots_adjust(bottom=0.13, left=0.4)
        pos_hclust = self.ax_hclus.get_position()
        pos_hot = self.ax_hot.get_position()
        pos_cold = self.ax_cold.get_position()
        pos_chot = self.ax_cbar_hot.get_position()
        pos_ccold = self.ax_cbar_cold.get_position()
        pos_hclust2 = [pos_hclust.x0, pos_hclust.y0 - 0.06, pos_hclust.width, pos_hclust.height]
        pos_hot2 = [pos_hot.x0, pos_hot.y0 - 0.03, pos_hot.width, pos_hot.height]
        pos_cold2 = [pos_cold.x0, pos_cold.y0, pos_cold.width, pos_cold.height]
        pos_chot2 = [pos_chot.x0, pos_chot.y0 - 0.03, pos_chot.width, pos_chot.height]
        pos_ccold2 = [pos_ccold.x0, pos_ccold.y0, pos_ccold.width, pos_ccold.height]     
        self.ax_hclus.set_position(pos_hclust2)
        self.ax_hot.set_position(pos_hot2)
        self.ax_cold.set_position(pos_cold2)
        self.ax_cbar_hot.set_position(pos_chot2)
        self.ax_cbar_cold.set_position(pos_ccold2)

        plt.setp(self.ax_cold.get_xticklabels(), rotation=70, ha='right')
        #plt.savefig('%s/FeatureHeatmap_%s.png' %(self.save_fig, color), dpi=600)
        #plt.savefig('%s/FeatureHeatmap_%s_noannot.png' %(self.save_fig, color), dpi=600)
        plt.suptitle('Random Forest Feature Ranking')
        if self.db == 'metaphlan': plt.savefig('%s/FeatureHeatmap_%s.%s' %(self.save_fig, color, self.fig_fmt), dpi=600)
        elif self.db == 'pathways': plt.savefig('%s/PathwaysHeatmap_%s.%s' %(self.save_fig, color, self.fig_fmt), dpi=600)


    def heatmap_simple(self):
        fig, (ax_hmap, ax_cbar) = plt.subplots(1,2, figsize=(8,12), gridspec_kw={'width_ratios':[8,1]})
        vmin_, vmax_ = np.min(self.data.values), np.max(self.data.values)
        norm_ = matplotlib.colors.Normalize(vmin=vmin_, vmax=vmax_)
        cbar = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=self.cmap, norm=norm_, extend='both', filled=True, drawedges=True)
        #cbar.outline.set_edgecolor('black')
        #cbar.outline.set_linewidth(3.) 
        #cbar.dividers.set_color('black')
        #cbar.dividers.set_linewidth(.12)
        xticks = [(nm if not nm in self.utils.data_aliases else self.utils.data_aliases[nm]) for nm in self.data.columns.tolist()]
        yticks = [s.replace('_',' ') for s in self.data.index.tolist()] #[s.replace('_',' ') for s in [('$\it{%s}$' %s_) for s_ in self.data.index.tolist()]]
        hm = sns.heatmap(self.data, ax=ax_hmap, cmap=self.cmap, cbar=False, vmin=vmin_, vmax=vmax_, xticklabels=xticks, yticklabels=yticks, linewidths=.0, linecolor='black', center=0.0)
        for label in ax_hmap.get_yticklabels(): label.set_style('italic') 

        ax_cbar.set_ylabel('Feature Importance Rank', size=16) ##, fontweight='bold')
        hm.set_xticklabels(hm.get_xticklabels(), fontsize=9)  ##, fontweight='bold')
        r_ticks = np.arange(vmin_,vmax_,5)
        if r_ticks[-1] != vmax_: r_ticks = list(r_ticks) + [vmax_]
        cbar.outline.set_clip_on(r_ticks)
        cbar.set_ticks(r_ticks)
        cbar.set_ticklabels([str(int(t)) for t in r_ticks]) ##, fontweight='bold')
        ax_cbar.tick_params(labelsize='large')
        ax_hmap.xaxis.set_ticks_position('top')
        plt.setp(ax_hmap.get_xticklabels(), rotation=60, ha='left') 
        plt.subplots_adjust(top=1., left=0.16)  ##, bottom=0.3)
        plt.tight_layout()
        
        pos1 = ax_cbar.get_position()
        pos2 = [pos1.x0*0.95, pos1.y0, pos1.width*1.5, pos1.height]
        ax_cbar.set_position(pos2)

        #plt.suptitle('Feature Importance Ranking')
        if self.db == 'metaphlan':
            plt.savefig('%s/FeatureHeatmap.png' %(self.save_fig), dpi=600) 
        elif self.db == 'pathways':
            plt.savefig('%s/PathwaysHeatmap.png' %(self.save_fig), dpi=600)


    def heatmap_pretty(self):
        fig, (ax0, ax1) = plt.subplots(2,1, figsize=(14,14), gridspec_kw={'width_ratios':None, 'height_ratios':[1, 8]})
        xticks = [(nm if not nm in self.utils.data_aliases else self.utils.data_aliases[nm]) for nm in self.data.columns.tolist()]
        vmin_, vmax_ = np.min(self.data.values), -np.min(self.data.values)
        annotations = self.get_annotation_matrix(self.frames).values
        _vector = self.compute_LODO_vector() if self.lodo else self.compute_CV_vector()
        ccmap = sns.diverging_palette(210, 12, sep=6, as_cmap=True, center='dark')
        ax1.tick_params(labelsize=int(self.char_size_hm))
        hm = sns.heatmap(self.data, ax=ax1, fmt='', cmap=ccmap, cbar=False, vmin=vmin_, vmax=vmax_\
		, xticklabels=False, linewidths=.0, linecolor='black', annot=annotations, annot_kws={'size': int(self.char_size_hm)})
        cbar0 = ax1.figure.colorbar(ax1.collections[0], extend='both', drawedges=True)
        cbar0.set_ticks([vmin_, 0.0, vmax_]) 
        cbar0.set_ticklabels(['%s-leaning' %self.tests[0].split(':')[2], 'non-detected', '%s-leaning' %self.tests[0].split(':')[1]])
        hml = sns.heatmap(pd.DataFrame(data=np.array([_vector[:-1]], dtype=np.float64), columns=self.datasets) \
                , ax=ax0, cbar_kws={'orientation': 'vertical', 'aspect':.8, 'shrink':0.5, 'extend':'both', 'ticks': [0.5, 0.9], 'extendfrac':0.35, 'drawedges':True}\
		, annot=True, fmt='.2f', cmap='hot', cbar=True, vmin=0.5, vmax=0.9, linewidths=.5, linecolor='black', square=True)
        plt.subplots_adjust(left=0.525, bottom=0.3, hspace=0.005) ###fig.autofmt_xdate(rotation=45, ha='right')
        ax0.xaxis.tick_top()
        ax0.set_yticks([0.5])
        ax0.set_yticklabels(['leave-one-dataset-out auc' if self.lodo else 'cross-validation auc'], rotation=0)
        fig.autofmt_xdate(rotation=45, ha='left')
        plt.suptitle(('LODO' if self.lodo else 'single-cross-validation')+'-feature-importance-map, %s, most %i important feats' \
		%(('_'.join(self.db) if isinstance(self.db, list) else self.db), self.n_imp_feat))
        plt.savefig(('%s/PAN_featureHeatmap_' %self.save_fig)+('LODO' if self.lodo else 'CrossValidations')+'_%s_%i_features.png' %(self.db, self.n_imp_feat), dpi=600)


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


    def datamatrix(self):
        f = self.frames

        for d in f: 
            print "===================================="
            print d
            print "===================================="

        #exit(1)

        ## following compiles a list of the most relevant features of all the datsets
        for piece in f: self.most_relevant_features += piece.index.tolist()[:self.n_imp_feat]
        self.most_relevant_features = list(set(self.most_relevant_features))

        ### QUI FAI IL MERGE
            
        print self.most_relevant_features, len(self.most_relevant_features)

        exit(1)     


        self.most_relevant_features = list(set(self.most_relevant_features))
        f = f.loc[self.most_relevant_features]

        ## now features have to ranked, all togheter, without separating the cold from the hot ones
        #for ds in self.datasets: f.loc[:, ds] = rankdata(f.loc[:, ds]) 
        for ds in self.datasets:

            print "==================="
            print "==================="
            print 'parto qui:...'
            print ds

            print f.loc[:, ds]
            print rankdata(f.loc[:, ds])
            print "======="
            print "======="
            print "======="

        exit(1)


        ## following compiles two list, each with the positive or the negative features
        ## then attach them their sum across all dataset, sort them by this sum 
        ## and exclude the sum from the result list.
	## CORRECTED <= (the sum is probably meaninless if the rank is not executed BEFORE, 'cos dataets are different)
        ## NOTABLY: reverse=True: best feature will have a higher value which results in a lower rank, therefore the right order is reversed
        neg = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) for ft in f.index.tolist() if self.signs_[ft] < 0.], key=lambda el : el[1], reverse=True)]
        pos = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) for ft in f.index.tolist() if self.signs_[ft] > 0.], key=lambda el : el[1], reverse=True)]
        f = f.loc[pos + neg]

        ## this can be eventually eliminated 
        self.min_rank = min(np.max(f.loc[pos, :].values), np.max(f.loc[neg, :].values)) 
        self.min_rank = 15 if self.min_rank > 20 else self.min_rank
        
        ## this not
        rename_spp = lambda s_ : s_.replace('unclassified','spp.') if s_.endswith('unclassified') else s_ 
        if self.db == 'metaphlan': 
            pos = [rename_spp(i[3:]) for i in pos]
            neg = [rename_spp(i[3:]) for i in neg]
            f.set_index([[rename_spp(i[3:]) for i in f.index.tolist()]], inplace=True)

        self.pos = pos
        self.neg = neg

        ## THIS WS A MISTAKE
        #f.loc[pos] = f.loc[pos].apply(lambda dist : (self.min_rank - (1. + (((dist - np.min(dist)) / (np.max(dist) - np.min(dist))) * (self.min_rank - (1.))))) , axis=0) + 1.
        #f.loc[neg] = f.loc[neg].apply(lambda dist : (self.min_rank - (1. + (((dist - np.min(dist)) / (np.max(dist) - np.min(dist))) * (self.min_rank - (1.))))) , axis=0) + 1.

        #f = f.loc[self.pos+self.neg]
        
        f['SUM'] = f.sum(axis=1)        
        #del f['SUM']

        print f.loc[pos]
        print '================='
        print f.loc[neg]

        exit(1)        
        return f


    def concat_frames(self, lodo): return [self.features_(ds, lodo) for ds in self.datasets]
    #feat_ = self.features_(self.datasets[0], lodo)
    #for ds in self.datasets[1:]: feat_ = feat_.merge(self.features_(ds, lodo), right_index=True, left_index=True, how='outer').fillna(0.0)
    #return feat_


    def get_annotation_matrix(self, frame): return frame.apply(lambda f : [('+' if n else '') for n in f.tolist()], axis=1)
    def signs(self): return dict([(feat, self.sign_of_a_feature(self.phylogenetic_mapper[feat])) for feat in self.phylogenetic_mapper.keys()])
    def the_greater_is_the_second_in(self, tup): return True if (tup[1]>tup[0]) else False
    def mean_on_dataset(self, ds, s__): return self.subset(ds, self.tests[0].split(':')[1:])[s__].astype('float').mean()
    def mean_on_subset(self, ds, s__, c): return self.subset(ds, c)[s__].astype('float').mean()
    def sign_of_a_feature(self, f): return 1 if (not self.the_greater_is_the_second_in(self.mean_per_class(f))) else -1



    def mean_per_class(self, s__): 
        return tuple(self.metadata[self.metadata[self.tests[0].split(':')[0]].isin([c])][s__].astype('float').mean() for c in self.tests[0].split(':')[1:])


   
    def features_(self, ds, lodo, line='line', skip=0):

        f = open(self.feat_files[ds]) if not lodo else open(self.lodo_files[ds])
        while (not line.startswith('Feature')): line, skip = f.readline(), skip+1

        d = pd.read_csv(self.feat_files[ds] if not lodo else self.lodo_files[ds], sep='\t', header=None, index_col=0, skiprows=skip)[[1,2]]
        for f in d[1]: self.phylogenetic_mapper[f.split('|')[-1]] = f 
        d[1] = d[1].apply(lambda n : n.split('|')[-1])
        d[2] = d[2].apply(lambda n : float(n))

        d.columns = ['feat', str(ds)]
        d.set_index('feat', inplace = True)        

        d[ds] = [i+1 for i,f in enumerate(d[ds].tolist())]
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
