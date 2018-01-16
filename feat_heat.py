#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import utils
import glob
import os
import subprocess as sp
import pandas as pd
import numpy as np
from scipy.stats import rankdata

class feature_heatmap(object):
    dataselect = '/python ../../cmdpy/datasetSelection2.py'
    save_fig = '../Images'
    getfeatnumber = {10:13, 20:24, 30:35, 40:46, 50:57, 60:68} 

    def __init__(self, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python, n_imp_feat, lodo, char_size_hm, path='/CM/data/meta/'):
        self.n_imp_feat = n_imp_feat		# num of feature to look at for each
        self.char_size_hm = char_size_hm	# char size on the heatmap
        self.lodo = lodo			# a boolean
        self.which_python = which_python	
        self.utils = utils.usefullfuncs(datasets)
        self.path = path
        self.datasets = datasets
        self.db = db
        self.feat_iden = self.utils.features[db]
        self.algo = algo
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2
        self.segregate_features = True   	# no idea
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
        self.signs_ = self.signs()
        self.data = self.datamatrix()        	#self.heatmap_pretty()   		## main
        self.heatmap_simple()			## real main


    def select_nlargest(self):
        features = self.data.nlargest(self.n_imp_feat, self.datasets[0]).index.tolist()
        for data in self.datasets[1:]: features += self.data.nlargest(self.n_imp_feat, data).index.tolist()
        return list(set(features))


    def heatmap_simple(self):
        fig,ax = plt.subplots(figsize=(10,10))#, gridspec_kw=grid_kws)
        ##self.data = self.data.loc[self.select_nlargest()]

        xticks = [(nm if not nm in self.utils.data_aliases else self.utils.data_aliases[nm]) for nm in self.data.columns.tolist()]
        vmin_, vmax_ = np.min(self.data.values), np.max(self.data.values)
        #ccmap = sns.diverging_palette(0,255, sep=77, as_cmap=True, center='dark')
        hm = sns.heatmap(self.data, ax=ax, cmap='rainbow', cbar=True, vmin=vmin_, vmax=vmax_, xticklabels=xticks, linewidths=.2, linecolor='white', center=0.0)
        hm.set_xticklabels(hm.get_xticklabels(), fontsize=10) # rotation=45, rotation_mode='anchor')  

        #cbar = ax.figure.colorbar(ax.collections[0])
        #cbar.set_ticks([vmin_, 0.0, vmax_])
        #cbar.set_ticklabels(['%s-leaning' %self.tests[0].split(':')[2], 'non-detected', '%s-leaning' %self.tests[0].split(':')[1]])

        plt.subplots_adjust(left=0.525, bottom=0.3)
        fig.autofmt_xdate(rotation=45, ha='right')
        plt.suptitle('Feature Importance Ranking')
        plt.savefig('%s/FeatureHeatmap.png' %(self.save_fig), dpi=600) 


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
        for ds in f.columns.tolist(): self.most_relevant_features += f.nlargest(self.n_imp_feat, ds).index.tolist()
        self.most_relevant_features = list(set(self.most_relevant_features))

        for ds in f.columns.tolist(): f.loc[:, ds] = rankdata(np.array(f.loc[:, ds].tolist(), dtype=np.float64))

        f = f.apply(lambda row : row*self.signs_[row.name], axis=1)
        f = f.apply(lambda r : [(0.0 if (x==0.0) else x) for x in r], axis=0) 

        f = f.loc[self.most_relevant_features]
        sort_on_sum = lambda dataframe, feat: sum(dataframe.loc[feat, :]) 

        negatives = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) for ft in f.index.tolist() if self.signs_[ft]<0], key=lambda el : el[1], reverse=True)]
        positives = [el[0] for el in sorted([(ft, np.sum(f.loc[ft, :].tolist())) for ft in f.index.tolist() if self.signs_[ft]>0], key=lambda el : el[1], reverse=True)]
        
        for ds in self.datasets:
            f.loc[negatives, ds] = -1*(rankdata(-1*np.array(f.loc[negatives, ds].tolist(), dtype=np.float64)))
            f.loc[positives, ds] = rankdata(np.array(f.loc[positives, ds].tolist(), dtype=np.float64))

        f = f.loc[self.most_relevant_features]
        sorted_indices = [ft for ft in f.index.tolist() if np.sum(f.loc[ft, :])>0.] + [ft for ft in f.index.tolist() if np.sum(f.loc[ft, :])<0.] 
        f = f.loc[sorted_indices, :]

        if self.db == 'metaphlan': f.set_index([[i[3:] for i in f.index.tolist()]], inplace=True)
        return f


    def concat_frames(self, lodo):
        feat_ = self.features_(self.datasets[0], lodo)
        for ds in self.datasets[1:]: 
            feat_ = feat_.merge(self.features_(ds, lodo), right_index=True, left_index=True, how='outer').fillna(0.0)
        return self.normalize_importance(feat_)    


    def get_annotation_matrix(self, frame):
        return frame.apply(lambda f : [('+' if n else '') for n in f.tolist()], axis=1)


    def normalize_importance(self, f):
        f = (f-f.min())/(f.max()-f.min())
        return f


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
        d = pd.read_csv(self.feat_files[ds] if not lodo else self.lodo_files[ds], sep='\t', header=None, index_col=0, skiprows=skip)[[1,2]] ## deleted here nrows = n_imp_rows which was the error
        for f in d[1]: self.phylogenetic_mapper[f.split('|')[-1]] = f
        d[1] = d[1].apply(lambda n : n.split('|')[-1])
        d[2] = d[2].apply(lambda n : float(n))
        d.columns = ['feat', str(ds)]
        d.set_index('feat', inplace=True)        
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
        if os.path.exists('complete_metadata_'+('_'.join(self.db) if isinstance(self.db, list) else self.db)+'.csv'):
            pass
        else:
            sp.call(' '.join(command), shell=True)
        md = pd.read_csv('complete_metadata_'+('_'.join(self.db) if isinstance(self.db, list) else self.db)+'.csv', sep='\t', header=None, index_col=0, low_memory=False)        
        #os.remove('tmp_metadata.csv')
        return md


if __name__=='__main__':
    print 
