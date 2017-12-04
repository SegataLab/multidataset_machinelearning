#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils
import itertools
import glob

class training_support(object):
    save_folder = '../Images/'
    def __init__(self, title, datasets, db, defined_problem, algo, grid, feat_sel, which_python, mean_median):
        self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.db = db
        self.title = title
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.feat = feat_sel
        self.algo = algo
        self.grid = grid
        self.h_col = len(self.dataset)/2 if (len(self.datasets)%2==0) else ((len(self.datasets)+1)/2)
        self.func = np.median if mean_median == 'median' else np.mean 

    def _values(self, start):
        self.test_n = self.utils.test_magns(('_'.join(self.db) if isinstance(self.db, list) else self.db), self.tests[0], self.algo, self.grid, self.feat, start)
        means = dict([(ds,[]) for ds in self.datasets])
        for dataset_on_which in self.datasets:
            transfer_pattern = self.utils.get_all_trans_for_one(self.db, self.algo, self.grid, self.feat, dataset_on_which)    
            points = [self.utils._auc(p, start)[0] for p in glob.glob(transfer_pattern) if p.count(dataset_on_which)==1]
            if self.func == np.mean:
                wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
            else:
                wavg = np.median(points)
            means[dataset_on_which].append(points + wavg)
            for i in range(2, (len(self.datasets)-1)):
                std_support_pattern = self.utils.get_std_support_pattern(self.db, self.algo, self.grid, self.feat, dataset_on_which, i)
                points = [self.utils._auc(p, start)[0] for p in glob.glob(std_support_pattern)]
                if self.func == np.mean:
                    wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
                else:
                    wavg = np.median(points)
                means[dataset_on_which].append(points + wavg)
            last_pattern = self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.tests[0], self.algo, self.grid, self.feat)
            means[dataset_on_which].append([self.utils._auc(self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.tests[0], self.algo, self.grid, self.feat), start)[0]])
        return means

    def plot_axes(self):
        fig, ax = plt.subplots(self.h_col, 2, figsize=(12, 12))
        for h in range(self.h_col):
            for w in range(2):
                ax[h, w].set_xticks([-1]+range(len(self.datasets)-1))
                ax[h, w].set_xticklabels(['']+[str(1+x) for x in range(len(self.datasets)-1)])
                ax[h, w].set_ylim(bottom=0.5, top=.9)
                #ax[h, w].set_yticks(np.arange(0.4, 1.0, 0.1))
                #ax[h, w].set_yticklabels(list(map(str, np.arange(0.4, 1.0, 0.1))))
                ax[h, w].set_ylabel('AUC', fontsize=8)
        ax[h, w].set_ylim(bottom=0.5, top=.9)
        #ax[h, w].set_yticks(np.arange(0.4, 1.0, 0.1))
        #ax[h, w].set_yticklabels(list(map(str, np.arange(0.4, 1.0, 0.1))))
        ax[h, w].set_ylabel('AUC', fontsize=8)
        ax[h, w].set_title('Average' if self.func == np.mean else 'Median', fontsize=8)
        return fig, ax

    def plot_(self, start, n_feat):
        _values = self._values(start)
        fig, ax = self.plot_axes()
        for ds in self.datasets:
            h, w = self.datasets.index(ds) % self.h_col, 0 if (self.datasets.index(ds)<=(len(self.datasets)/2)) else 1
            for i in range(len(self.datasets)-1):
                ax[h, w].scatter([i]*len(_values[ds][i][:-1]), _values[ds][i][:-1], color='dodgerblue')
                ax[h, w].scatter([i], _values[ds][i][-1], color='goldenrod', s=(_values[ds][i][-1]*100))
            ax[h, w].set_title(ds if (not ds in self.utils.crc_paper_converter) else self.utils.crc_paper_converter[ds], fontsize=8)
            ax[h, w].scatter([len(self.datasets)-2], _values[ds][-1][0], color='goldenrod', s=(_values[ds][-1][0]*100))
            ax[h, w].plot(range(len(self.datasets)-1), [_values[ds][j][-1] for j in range(len(_values[ds]))], linestyle='-', color='orange', linewidth=8.0, alpha=0.5)
        average, h = [], h+1
        for i in range(len(_values[self.datasets[-1]])):
            average.append([_values[ds][i][-1] for ds in self.datasets])
        for x,y in zip(range(len(self.datasets)-1), average):
            ax[h, w].scatter([x]*len(y), y, color='dodgerblue')
            ax[h, w].scatter([x], np.mean(y), color='goldenrod', s=(np.mean(y)*100))
        ax[h, w].plot(range(len(self.datasets)-1), [self.func(a) for a in average], color='orange', linewidth=8.0, alpha=0.5)
        plt.suptitle('Progressive Learning; %s, %s feats, %s' %(self.db, str(n_feat), 'average' if self.func==np.mean else 'median'))
        plt.savefig('%s/ProgressiveLearning_%s_%s_%s' %(self.save_folder, self.db, str(n_feat), 'average' if self.func==np.mean else 'median'))

if __name__ == '__main__':
    print 'lascia perdere'
