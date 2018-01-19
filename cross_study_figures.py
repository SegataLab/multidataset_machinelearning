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
import utils

## meteaphlan plot: 
## python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do cross_figures -g0 nt:1000 -g1 nsl:5 -g2 c:entropy
## gene-families plot:
##

class cross_figures(object):
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    save_folder = '../Images/'


    def __init__(self, datasets, title, defined_problem, which_python):
	self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.title = title
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.couples = list(itertools.combinations_with_replacement(self.datasets, 2))
        self.dataframe_shape = (len(self.datasets)+1, len(self.datasets)+1)
        self.lodo_shape = len(self.datasets)+1 


    def plot_data(self, db, test, algo, grid0, grid1, grid2, start):
        scores = []
        coors = []
        test_n = self.utils.test_magns(('_'.join(db) if isinstance(db, list) else db), test, algo, grid0, grid1, grid2, start)
        avg_scores = dict([(d, [0,0]) for d in self.datasets])
        grid00 = grid0 #.replace(':','_')
        grid11 = grid1 #.replace(':','_')
        grid22 = grid2 #.replace(':','_')
        ### REEFINE COORDINATES BASED ON CROSS VALIDATION:
        new_order = [dataset[0] for dataset in sorted([[cross, self.utils.transfer_([cross,cross], db, test, algo, grid00, grid11, grid22, start)] for cross in self.datasets], key=lambda a : a[1], reverse=True)]
        norm_mean = 0.0

        for couple in self.couples:
            if self.utils.isonedata(couple): pools, cross_validation = [couple], True                
            else: pools, cross_validation = [couple, [couple[1],couple[0]]], False
            for pool in pools:
                score, coordinates, n = self.utils.transfer_(pool, db, test, algo, grid00, grid11, grid22, start, new_order)
                scores.append(score)
                if cross_validation: norm_mean += score*test_n[pool[0]]
                coors.append([int(coordinates[0]), int(coordinates[1])])
                if not cross_validation:
                    #print pool[0], ' SU: ', pool[1], ' ==> ', score, ' * ', test_n[pool[1]] , ' = ', score * test_n[pool[1]] , cross_validation, ' EHEHEHEH'
                    avg_scores[pool[0]][0] += score * test_n[pool[1]]           
                    avg_scores[pool[1]][1] += score

        for i,d in enumerate(self.datasets):
            scores.append(avg_scores[d][1]/float(len(self.datasets)-1))   
            coors.append([len(self.datasets), new_order.index(d)]) ## found bug
            scores.append(avg_scores[d][0]/float(sum([test_n[ds] for ds in self.datasets if ds != d])))
            coors.append([new_order.index(d), len(self.datasets)])

        mean, div = [], 0
        lodo_row = len(self.datasets) + 1 ## changed 2
        norm_mean /= float(sum([test_n[d] for d in self.datasets]))
        scores.append(norm_mean)
        coors.append([len(self.datasets),len(self.datasets)])

        for i,d in enumerate(self.datasets):
            score, coordinate, n = self.utils.lodo_(d, db, test, algo, grid00, grid11, grid22, start, new_order)
            mean.append(score * n)
            div += n
            coors.append(coordinate)  ###(lodo_row, i))
            scores.append(score)

        scores.append(sum(mean)/float(div))
        coors.append(coordinate+1)#(lodo_row, i+1))
        cross_validation, n_feat, n = self.utils.cross_validation_(db, test, algo, grid00, grid11, grid22, start)
        return scores, coors, n_feat, new_order


    def plotdata(self, db, algo, grid0, grid1, grid2, start):
        test = self.tests[0]
        all_data, all_coordinates, n_feat, new_order = self.plot_data(db, test, algo, grid0, grid1, grid2, start)
        cross_data, cross_coordinates = all_data[:-(len(self.datasets)+1)], all_coordinates[:-(len(self.datasets)+1)]
        lodo_data, lodo_coordinates = all_data[-(len(self.datasets)+1):], all_coordinates[-(len(self.datasets)+1):]
        cross_table = np.empty(self.dataframe_shape, dtype=np.float64) 
        lodo_table = np.empty([1, self.lodo_shape], dtype=np.float64)

        for sc,cr in zip(lodo_data, lodo_coordinates): lodo_table[0, int(cr)] = sc
        for sc,cr in zip(cross_data, cross_coordinates): cross_table[tuple(cr)] = sc
        new_order = [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in new_order]
        cross_data = pd.DataFrame(cross_table, columns=new_order + ['Average'], index=new_order + ['Average'])
        lodo_data = pd.DataFrame(lodo_table, columns=new_order + ['Average'], index=['LODO'])        

        mask_ = pd.DataFrame(dict([(c,False) for c in cross_data.columns.tolist()]), index=cross_data.index.tolist(), columns=cross_data.columns.tolist())
        mask_.loc['Average', 'Average'] = True
        return mask_, cross_data, lodo_data, n_feat


    def plot(self, db, algo, grid0, grid1, grid2, start):
        sns.set(font_scale=1.2)
        mask, cross_data, lodo_data, nf  = self.plotdata(db, algo, grid0, grid1, grid2, start)
        fig = plt.figure(figsize=(10,12))
        gs = gridspec.GridSpec(len(self.datasets)+1,len(self.datasets)+1)  ##  , width_ratios=[3, 1, 1], height_ratios=[3, 1])
        ax_lodo = plt.subplot(gs[-1:, 0:-1])
        ax_cross = plt.subplot(gs[0:-1, 0:-1])
        ax_cbar = plt.subplot(gs[2:, -1])

        data = cross_data.append(lodo_data).values
        annot_o, annot_t = cross_data.apply(lambda a : ['%.2f' %aa for aa in a]).values, lodo_data.apply(lambda o : ['%.2f' %oo for oo in o]).values			### 'fontweight': 'bold'
        hm_o = sns.heatmap(cross_data, ax=ax_cross, fmt='', cmap='hot', cbar=False, vmin=0.5, vmax=1., linecolor='black', annot=annot_o, xticklabels=True, square=True, annot_kws={'fontsize': 19}, linewidth=0.)
        hm_t = sns.heatmap(lodo_data, ax=ax_lodo, fmt='', cmap='hot', cbar=False, vmin=0.5, vmax=1., linecolor='black', annot=annot_t, xticklabels=False, annot_kws={'fontsize': 19}, square=True, linewidth=0.) 
        norm_ = matplotlib.colors.Normalize(vmin=.5,vmax=1.)
        cbar = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap='hot', norm=norm_, extend='min', filled=True, drawedges=True)
        #cbar.outline.set_edgecolor('black')
        #cbar.outline.set_linewidth(5.)
        #cbar.dividers.set_color('black')
        #cbar.dividers.set_linewidth(.2)
        cbar.set_ticks(list(np.arange(0.5,1.0,0.1))+[1])
        cbar.set_ticklabels(list(map(str,np.arange(0.5,1.0,.1)))+[1.0])

        ax_cbar.tick_params(labelsize=16)
        ax_cbar.set_ylabel('AUC', size=16, fontweight='bold')
        ax_cross.xaxis.set_ticks_position('top')
        ax_cross.set_xticklabels(ax_cross.get_xticklabels(), fontweight='bold')
        ax_cross.set_yticklabels(ax_cross.get_yticklabels(), fontweight='bold')
        #ax_lodo.set_xticklabels(ax_lodo.get_xticklabels(), fontweight='bold')
        ax_lodo.set_yticklabels(ax_lodo.get_yticklabels(), fontweight='bold')
        ax_cbar.set_yticklabels(ax_cbar.get_yticklabels(), fontweight='bold')
        #plt.setp(ax_lodo.get_xticklabels(), rotation=30, ha='right')

        plt.setp(ax_lodo.get_yticklabels(), rotation=30)
        plt.setp(ax_cross.get_yticklabels(), rotation=30)
        plt.setp(ax_cross.get_xticklabels(), rotation=38, ha='left')

        plt.subplots_adjust(left=0.24) ## 0.18        
        pos1 = ax_lodo.get_position() # get the original position 
        pos2 = [pos1.x0, pos1.y0+0.1, pos1.width, pos1.height]
        ax_lodo.set_position(pos2)
        posc1 = ax_cbar.get_position() 
        posc2 = [posc1.x0+0.01, posc1.y0+0.15, posc1.width*0.6, posc1.height-0.1]
        ax_cbar.set_position(posc2) ##, which='original')
        plt.savefig('%s/MachineLearning_%s_%s_%s_features.png' %(self.save_folder, self.title, db, str(nf if nf<400 else 'all')), dpi=600)
        

if __name__ == '__main__': 
    print 'se io ti rullassi di cartoni' 
