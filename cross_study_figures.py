#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import matplotlib.gridspec as gridspec
import itertools
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
import utils


## WITH RANDOM FOREST
## meteaphlan plot:				FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 
## python run.py crc --define study_condition:CRC:control -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -db metaphlan -do cross_figures -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -cm hot

## gene-families plot: 
## python run.py crc --define study_condition:CRC:control -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -db genefamilies -do cross_figures -g0 nsl:5 -g1 nt:1000 -g2 df -cm hot

## new attempt on gene families
## python run.py crc --define study_condition:CRC:control -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -db genefamilies -do cross_figures -g0 nsl:5 -g1 nt:1000 -g2 df -cm hot -ft MachineLearning_crc_genefamilies_last_attempt_before_PUB

## pathways
## python run.py crc --define study_condition:CRC:control -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017  -db pathways -do cross_figures -g0 nt:1000 -g1 nsl:5 -g2 c:entropy -cm hot -w 1.5

## markers
##1)  python run.py crc --define study_condition:CRC:control -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -db markers -do cross_figures -g0 nt:1000 -g1 nsl:5 -g2 df -cm hot 
##2)  python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db markers -do cross_figures -g0 c:gini -g1 nsl:5 -g2 df -cm hot  

## WITH LINEAR SVM
## 3) python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db markers -do cross_figures -g0 null -g1 null -g2 null -al lsvm -ft MachineLearning_crc_markers_linSVM



class score(object):
    def __init__(self, sc, cr):
        self.auc = sc
        self.where = cr


class scores(object):
    def __init__(self):
        self.lodo = list()
        self.cross_study = list()
        self.norm_mean = list()
        self.cmean = list()
        self.rmean = list()
        self.lodo_mean = list()


class cross_figures(object):
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    save_folder = '../Images/'


    def __init__(self, datasets, title, defined_problem, which_python, cmap, width, fig_title, fig_fmt):
	self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.title = title
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.couples = list(itertools.combinations_with_replacement(self.datasets, 2))
        self.dataframe_shape = (len(self.datasets)+1, len(self.datasets)+1)
        self.lodo_shape = len(self.datasets)+1 
        self.scores = scores()
        self.cmap = cmap
        self.width = width 
        self.fig_title = fig_title
        self.resort_names = False
        self.fig_fmt = fig_fmt


    def plot_data(self, db, test, algo, grid0, grid1, grid2, start):

        scores_cross, scores_cmean, scores_rmean, score_cvalidation, scores_lodo = [], [], [], [], []
        coors = []
        test_n = self.utils.test_magns(('_'.join(db) if isinstance(db, list) else db), algo, test, grid0, grid1, grid2, start)
        avg_scores = dict([(d, [0,0]) for d in self.datasets])

        if self.resort_names or db == 'metaphlan':
            new_order = [dataset[0] for dataset in sorted([[cross, self.utils.transfer_([cross,cross]\
			, db, algo, test, grid0, grid1, grid2, start)] for cross in self.datasets], key=lambda a : a[1], reverse=True)]
        else: 
            new_order = self.datasets

        norm_mean = 0.0

        for couple in self.couples:
            if self.utils.isonedata(couple): pools, cross_validation = [couple], True 

            else: pools, cross_validation = [couple, [couple[1],couple[0]]], False

            for pool in pools:
                score_, coordinates, n = self.utils.transfer_(pool, db, algo, test, grid0, grid1, grid2, start) #, new_order)
                self.scores.cross_study.append(score(score_, coordinates))

                if cross_validation: 
                    norm_mean += score_ * test_n[pool[0]]
                else:
                    avg_scores[pool[0]][0] += score_ * test_n[pool[1]]           
                    avg_scores[pool[1]][1] += score_

        for i,d in enumerate(new_order):
            self.scores.rmean.append(score(avg_scores[d][1]/float(len(self.datasets)-1), i))   
            self.scores.cmean.append(score(avg_scores[d][0]/float(sum([test_n[ds] for ds in self.datasets if ds!=d])), i)) 

        mean, div = [], 0
        lodo_row = len(self.datasets) + 1 ## changed 2
        norm_mean /= float(sum([test_n[d] for d in self.datasets]))
        self.scores.norm_mean.append(score(norm_mean, [len(self.datasets),len(self.datasets)]))

        for i,d in enumerate(new_order): ## bug found
            score_, coordinate, n = self.utils.lodo_(d, db, algo, test, grid0, grid1, grid2, start) ##, new_order)
            mean.append(score_ * n)
            div += n
            self.scores.lodo.append(score(score_, i))  ###(lodo_row, i))
        self.scores.lodo_mean.append(score(sum(mean)/float(div), coordinate+1))   
        
        #try:
        #global_cross_validation, n_feat, n = self.utils.cross_validation_(db, test, algo, grid0, grid1, grid2, start)
        return new_order, 'all'


    def plotdata(self, db, algo, grid0, grid1, grid2, start):
        test = self.tests[0]
        new_order, n_feat = self.plot_data(db, test, algo, grid0, grid1, grid2, start)
        new_order = [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in new_order]

        cross_data_matrix = np.empty((len(self.datasets),len(self.datasets)), dtype=np.float64)
        for o in self.scores.cross_study: cross_data_matrix[o.where] = o.auc
        self.cross_study_data = pd.DataFrame(cross_data_matrix, columns=new_order, index=new_order)
        self.lodo_data = pd.DataFrame(dict([(left_out,lodo.auc) for lodo,left_out in zip(self.scores.lodo, new_order)]), index=['LODO'])
        self.cmean_data = pd.DataFrame(np.array([[sc.auc] for sc in self.scores.cmean], dtype=np.float64), index=new_order, columns=['Average'])
        self.rmean_data = pd.DataFrame(dict([(cname,rm.auc) for rm,cname in zip(self.scores.rmean, self.lodo_data.columns.tolist())]), index=['Average'])
        self.norm_mean_data = pd.DataFrame(self.scores.norm_mean[0].auc, index=[''], columns=[''])
        self.lodo_mean_data = pd.DataFrame(self.scores.lodo_mean[0].auc, index=[''], columns=[''])

        self.dataset_sorting = new_order
        self.feature_used = n_feat
        self.lodo_data = self.lodo_data[new_order]
        ##self.rmean_data = self.rmean_data[new_order]


    def define_proportions(self):
        pos1_cross = self.ax_cross.get_position() 
        pos1_rmean = self.ax_rmean.get_position() 
        pos1_cmean = self.ax_cmean.get_position() 
        pos1_lodo = self.ax_lodo.get_position() 
        pos1_norm = self.ax_norm_mean.get_position() 
        pos1_cbar = self.ax_cbar.get_position()
        pos1_lodo_mean = self.ax_lodo_mean.get_position()

        pos2_cbar = [pos1_cbar.x0 + 0.015, pos1_cbar.y0 + 0.08, pos1_cbar.width * 0.7, pos1_cbar.height * 0.80]
        pos2_cross = [pos1_cross.x0, pos1_cross.y0, pos1_cross.width + 0.001, pos1_cross.height]
        pos2_rmean = [pos1_cross.x0, pos1_rmean.y0, pos1_cross.width + 0.001, pos1_rmean.height + 0.01]
        #pos2_rmean = [pos1_cross.x0, pos1_rmean.y0, pos1_cross.width, pos1_rmean.height + 0.01]
        pos2_cmean = [pos1_cmean.x0 - 0.01, pos1_cross.y0, pos1_cmean.width + 0.01, pos1_cross.height] ## same height of th cross
        pos2_lodo = [pos1_cross.x0, pos1_lodo.y0 - 0.01, pos1_cross.width + 0.001, pos1_lodo.height + 0.01]
        ##pos2_lodo = [pos1_cross.x0, pos1_lodo.y0 - 0.01, pos1_cross.width, pos1_lodo.height + 0.01]
        pos2_norm = [pos1_norm.x0 - 0.01, pos1_norm.y0, pos1_norm.width + 0.01, pos1_norm.height + 0.01] 
        pos2_lodo_mean = [pos1_lodo_mean.x0 - 0.01, pos1_lodo_mean.y0 - 0.01, pos1_lodo_mean.width + 0.01, pos1_lodo_mean.height + 0.01]

        #pos2_cbar = [pos1_cbar.x0 + 0.015, pos1_cbar.y0 + 0.05, pos1_cbar.width * 0.7, pos1_cbar.height * 0.9]
        #pos2_cross = [pos1_cross.x0, pos1_cross.y0, pos1_cross.width, pos1_cross.height]
        #pos2_rmean = [pos1_rmean.x0, pos1_rmean.y0, pos1_cross.width , pos1_rmean.height + 0.01]
        ##pos2_rmean = [pos1_cross.x0, pos1_rmean.y0, pos1_cross.width, pos1_rmean.height + 0.01]
        #pos2_cmean = [pos1_cmean.x0 - 0.01, pos1_cross.y0, pos1_cmean.width + 0.01, pos1_cross.height] ## same height of th cross
        #pos2_lodo = [pos1_lodo.x0, pos1_lodo.y0 - 0.01, pos1_cross.width, pos1_lodo.height + 0.01]
        ##pos2_lodo = [pos1_cross.x0, pos1_lodo.y0 - 0.01, pos1_cross.width, pos1_lodo.height + 0.01]
        #pos2_norm = [pos1_norm.x0 - 0.01, pos1_norm.y0, pos1_norm.width + 0.01, pos1_norm.height + 0.01]
        #pos2_lodo_mean = [pos1_lodo_mean.x0 - 0.01, pos1_lodo_mean.y0 - 0.01, pos1_lodo_mean.width + 0.01, pos1_lodo_mean.height + 0.01]

        self.ax_cross.set_position(pos2_cross)
        self.ax_rmean.set_position(pos2_rmean)
        self.ax_cmean.set_position(pos2_cmean)
        self.ax_lodo.set_position(pos2_lodo)
        self.ax_norm_mean.set_position(pos2_norm)
        self.ax_cbar.set_position(pos2_cbar)
        self.ax_lodo_mean.set_position(pos2_lodo_mean)



    def plot(self, db, algo, grid0, grid1, grid2, start):
        sns.set(font_scale=.8)
        self.plotdata(db, algo, grid0, grid1, grid2, start)
        fig = plt.figure(figsize=(11,11))
        gs = gridspec.GridSpec(len(self.datasets)+2, len(self.datasets)+2)  ##  , width_ratios=[3, 1, 1], height_ratios=[3, 1])

        self.ax_cross = plt.subplot(gs[0:7, 0:7])
        self.ax_lodo = plt.subplot(gs[-1, 0:7])
        self.ax_rmean = plt.subplot(gs[-2, 0:7])
        self.ax_cmean = plt.subplot(gs[0:7, -2])
        self.ax_norm_mean = plt.subplot(gs[7, 7])
        self.ax_cbar = plt.subplot(gs[:, -1])
        self.ax_lodo_mean = plt.subplot(gs[-1, -2])

       
        self.annot_lodo = self.lodo_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
        self.annot_cross = self.cross_study_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
        self.annot_rmean = self.rmean_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
        self.annot_cmean = self.cmean_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
        self.annot_norm_mean = self.norm_mean_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
        self.annot_lodo_mean = self.lodo_mean_data.astype('float').apply(lambda a : ['%.2f' %aa for aa in a])
       

        self.fig_cross = sns.heatmap(self.cross_study_data, annot=self.annot_cross.values, ax=self.ax_cross, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=True, xticklabels=True, annot_kws={'fontsize': 27})
        self.fig_lodo = sns.heatmap(self.lodo_data, annot=self.annot_lodo.values, ax=self.ax_lodo, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=False, xticklabels=False, annot_kws={'fontsize': 27})
        self.fig_cmean = sns.heatmap(self.cmean_data, annot=self.annot_cmean.values, ax=self.ax_cmean, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=False, yticklabels=False, xticklabels=True, annot_kws={'fontsize': 27})
        self.fig_rmean = sns.heatmap(self.rmean_data, annot=self.annot_rmean.values, ax=self.ax_rmean, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=False, xticklabels=False, annot_kws={'fontsize': 27})
        self.fig_norm_mean = sns.heatmap(self.norm_mean_data, annot=self.annot_norm_mean.values, ax=self.ax_norm_mean, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=False, xticklabels=False, annot_kws={'fontsize': 27})
        self.fig_lodo_mean = sns.heatmap(self.lodo_mean_data, annot=self.annot_lodo_mean.values, ax=self.ax_lodo_mean, fmt='', cmap=self.cmap, cbar=False, vmin=0.5, vmax=1., linewidth=self.width, square=False, annot_kws={'fontsize': 27})   

        self.norm_ = matplotlib.colors.Normalize(vmin=.5,vmax=1.)
        self.cbar = matplotlib.colorbar.ColorbarBase(self.ax_cbar, cmap=self.cmap, norm=self.norm_, extend='min', filled=True, drawedges=True)
        self.cbar.set_ticks(list(np.arange(0.5,1.0,0.1))+[1])
        self.cbar.set_ticklabels(list(map(str,np.arange(0.5,1.0,.1)))+[1.0])

        ## IF YOU WANT TO SEE BLACK BORDERS (add them to the heatmap it's enogh to linewidth 1.)
        #cbar.outline.set_edgecolor('black')
        #cbar.outline.set_linewidth(5.)
        #cbar.dividers.set_color('black')
        #cbar.dividers.set_linewidth(.2)

        self.ax_cbar.tick_params(labelsize=16)
        self.ax_cbar.set_ylabel('AUC', size=16)
        self.ax_cross.xaxis.set_ticks_position('top')
        self.ax_cmean.xaxis.set_ticks_position('top')
        self.ax_cross.set_xticklabels(self.ax_cross.get_xticklabels())
        self.ax_cross.set_yticklabels(self.ax_cross.get_yticklabels())
        self.ax_lodo.set_yticklabels(self.ax_lodo.get_yticklabels())
        self.ax_cbar.set_yticklabels(self.ax_cbar.get_yticklabels())

        #self.ax_cross.set_ylabel('Training Cohorts', fontsize=12)
        #self.ax_cross.set_xlabel('Testing Cohorts', fontsize=12)
        #self.ax_cross.xaxis.set_label_position('top') 

        plt.setp(self.ax_cross.get_xticklabels(), rotation=38, ha='left')
        plt.setp(self.ax_cmean.get_xticklabels(), rotation=38, ha='left')
        plt.setp(self.ax_lodo.get_yticklabels(), rotation=0)
        plt.setp(self.ax_rmean.get_yticklabels(), rotation=0)

        ###plt.subplots_adjust(left=0.125) ##, top=0.80) #, bottom=0.12)
        self.define_proportions()
        ###plt.subplots_adjust(left=0.14)

        if not self.fig_title:
            plt.savefig('%s/MachineLearning_%s_%s_%s_features.%s' \
	    %(self.save_folder, self.title, db, str(self.feature_used if self.feature_used<400 else 'all'), self.fig_fmt), dpi=600)
        else:
            plt.savefig('%s/%s.%s' %(self.save_folder, self.fig_title, self.fig_fmt), dpi=600)


if __name__ == '__main__': 
    print 'auuuuuuuuuuuuuuuu' 
