#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
#import matplotlib.patches as patches
import itertools
import utils

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
        self.dataframe_shape = (len(self.datasets)+2, len(self.datasets)+1)

    def plot_data(self, db, test, algo, grid, feat, start):
        scores = []
        coors = []
        test_n = self.utils.test_magns(('_'.join(db) if isinstance(db, list) else db), test, algo, grid, feat, start)
        avg_scores = dict([(d, [0,0]) for d in self.datasets])
        for couple in self.couples:
            if self.utils.isonedata(couple):
                pools, cross_validation = [couple], True                
            else:
                pools, cross_validation = [couple, [couple[1], couple[0]]], False
            for pool in pools:
                score, coordinates, n = self.utils.transfer_(pool, db, test, algo, grid, feat, start)
                scores.append(score)
                coors.append(tuple(coordinates))
                if not cross_validation:
                    #print pool[0], ' SU: ', pool[1], ' ==> ', score, ' * ', test_n[pool[1]] , ' = ', score * test_n[pool[1]] , cross_validation, ' EHEHEHEH'
                    avg_scores[pool[0]][0] += score * test_n[pool[1]]           
                    avg_scores[pool[1]][1] += score
        for i,d in enumerate(self.datasets):
            scores.append(avg_scores[d][1]/float(len(self.datasets)-1))   
            coors.append((len(self.datasets), i))
            scores.append(avg_scores[d][0]/float(sum([test_n[ds] for ds in self.datasets if ds != d])))
            coors.append((i, len(self.datasets)))
        mean, div = [], 0
        lodo_row = len(self.datasets) + 1 ## changed 2
        for i,d in enumerate(self.datasets):
            score, coordinate, n = self.utils.lodo_(d, db, test, algo, grid, feat, start)
            mean.append(score * n)
            div += n
            coors.append((lodo_row, i))
            scores.append(score)
        scores.append(sum(mean)/float(div))
        coors.append((lodo_row, i+1))
        cross_validation, n_feat, n = self.utils.cross_validation_(db, test, algo, grid, feat, start)
        scores.append(cross_validation)
        coors.append((lodo_row-1, coordinate+1))
        return scores, coors, n_feat

    def plotdata(self, db, algo, grid, feat, start):
        test = self.tests[0]
        data, coordinates, n_feat = self.plot_data(db, test, algo, grid, feat, start)
        tb = np.empty(self.dataframe_shape, dtype=np.float64)
        for cr,sc in zip(coordinates[:-1], data[:-1]): tb[cr] = sc
        return coordinates[-1], n_feat, pd.DataFrame(data=tb\
		, columns=[(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets]+['Average']\
		, index=[(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets]+['Average','LODO'])

    def plot(self, db, algo, grid, feat, start):
        mask, nf, data = self.plotdata(db, algo, grid, feat, start)
        fig, ax = plt.subplots(figsize=(10,10))       
        mask_ = pd.DataFrame([[False for t in range(data.shape[1])] for y in range(data.shape[0])])
        mask_.iloc[mask] = True
        hm = sns.heatmap(data, cmap='hot', vmin=0.50, vmax=.95, cbar=False, annot=True, fmt='.2f'\
			     , linewidths=2.0, linecolor='white', annot_kws={'fontsize': 20}, mask=mask_.values\
			     , square=True) #, cbar_kws={'extend': 'both', 'drawedges': True}, square=True)
 
        cbar = ax.figure.colorbar(ax.collections[0], extend='both', drawedges=True)
        cbar.set_ticks([0.50, 0.95])
        cbar.set_ticklabels(['0.50', '0.95'])

        #hm.set_xticklabels(hm.get_xticklabels(), fontsize=10) # rotation=45, rotation_mode='anchor') 
        #plt.subplots_adjust(left=0.525, bottom=0.3)
        #fig.autofmt_xdate(rotation=45, ha='center')

        plt.suptitle('cross-analyses, %s' % (('_'.join(db) if isinstance(db, list) else db)))
        plt.savefig('%s/PAN_AUC_heatmap_%s_%s_features.png' %(self.save_folder, db, str(nf if nf<400 else 'all')), dpi=600)

if __name__ == '__main__': 
    print 'non ti sopporto piu m\'hai rotto i coglioni' 
