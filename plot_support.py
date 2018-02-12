#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils
import itertools
import pandas as pd
import glob
import seaborn as sns

## to plot progressive singularly for metaphan:
## python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot -db metaphlan -g0 c:entropy -g1 nt:1000 -g2 nsl:5

## to plot progressive singularly for genefamilies:
## python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot -db genefamilies -g0 nt:1000 -g1 nsl:5 -g2 df

class training_support(object):
    save_folder = '../Images/'
    def __init__(self, title, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python, fig_fmt, mean_median):
        self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.db = db
        self.title = title
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.algo = algo
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid = '_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2
        self.h_col = len(self.dataset)/2 if (len(self.datasets)%2==0) else ((len(self.datasets)+1)/2)
        self.fig_fmt = fig_fmt
        self.func = 'median' if mean_median == 'median' else 'mean'
        

    def _values(self, start):
        problem = self.tests[0]
        self.test_n = self.utils.test_magns(('_'.join(self.db) if isinstance(self.db, list) else self.db), self.tests[0], self.algo, self.grid0, self.grid1, self.grid2, start)
        means = dict([(ds,[]) for ds in self.datasets])
        for dataset_on_which in self.datasets:
            transfer_pattern = self.utils.get_all_trans_for_one(self.db, self.algo, self.grid0, self.grid1, self.grid2, dataset_on_which)    
            points = [self.utils._auc(p, start)[0] for p in glob.glob(transfer_pattern) if p.count(dataset_on_which)==1]
            if self.func == 'mean': wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
            else: wavg = [np.median(points)] ## for p in points]
            means[dataset_on_which].append(points + wavg)
            for i in range(2, (len(self.datasets)-1)):
                std_support_pattern = self.utils.get_std_support_pattern(self.db, problem, self.grid0, self.grid1, self.grid2, dataset_on_which, i)
                points = [self.utils._auc(p, start)[0] for p in glob.glob(std_support_pattern)]
                if self.func == 'mean': wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
                else: wavg = [np.median(points)] ## for p in points]
                means[dataset_on_which].append(points + wavg)
            last_pattern = self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2)
            means[dataset_on_which].append([self.utils._auc(self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.tests[0], self.algo, self.grid0, self.grid1, self.grid2), start)[0]])
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
        ax[h, w].set_title('Average' if self.func == 'mean' else 'Median', fontsize=8)
        return fig, ax

    def plot_(self, start, n_feat):
        _values = self._values(start)
        fig, ax = self.plot_axes()
        for ds in self.datasets:
            h, w = self.datasets.index(ds) % self.h_col, 0 if (self.datasets.index(ds)<=(len(self.datasets)/2)) else 1
            for i in range(len(self.datasets)-1):
                ax[h, w].scatter([i]*len(_values[ds][i][:-1]), _values[ds][i][:-1], color='dodgerblue')
                ax[h, w].scatter([i], _values[ds][i][-1], color='goldenrod', s=(_values[ds][i][-1]*100))
            ax[h, w].set_title(ds if (not ds in self.utils.data_aliases) else self.utils.data_aliases[ds], fontsize=8)
            ax[h, w].scatter([len(self.datasets)-2], _values[ds][-1][0], color='goldenrod', s=(_values[ds][-1][0]*100))
            ax[h, w].plot(range(len(self.datasets)-1), [_values[ds][j][-1] for j in range(len(_values[ds]))], linestyle='-', color='orange', linewidth=8.0, alpha=0.5)
        average, h = [], h+1
        for i in range(len(_values[self.datasets[-1]])): average.append([_values[ds][i][-1] for ds in self.datasets])
        for x,y in zip(range(len(self.datasets)-1), average):
            ax[h, w].scatter([x]*len(y), y, color='dodgerblue')
            ax[h, w].scatter([x], np.mean(y) if self.func=='mean' else np.median(y), color='goldenrod', s=((np.mean(y) if self.func=='mean' else np.median(y))*100))
        ax[h, w].plot(range(len(self.datasets)-1), [(np.mean(a) if self.func=='mean' else np.median(a)) for a in average], color='orange', linewidth=8.0, alpha=0.5)
        plt.suptitle('progressive learning; %s, %s feats, %s' %(self.db, str(n_feat), 'average' if self.func=='mean' else 'median'))
        plt.savefig('%s/PAN_progressiveLearning_%s_%s_%s' %(self.save_folder, self.db, str(n_feat), 'average' if self.func=='mean' else 'median'))

    
    def corr_single(self, dataset, _values):
        fig,ax = plt.subplots(figsize=(4,4))
        indexed = [[(el,i+1) for el in e[:-1]] for i,e in enumerate(_values[:-1])] + [(_values[-1][0], i+1)] 
    
        data = list(itertools.chain.from_iterable(indexed[:-1])) + [indexed[-1]]
        X = [i[1] for i in data]
        Y = [d[0] for d in data]

        

        #for i,d in zip(indexes, data): print i,d
        #for i in data: print i

        exit(1)

        indexes = [[i[1] for i in I] for I in indexed]
        #values = 
        print indexes
        exit()    

        #X = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable([[i[0] for i in I] for I in indexed[:-1]]))))  
        #Y = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable([[i[1] for i in I] for I in indexed[:-1]]))))

        #for i,e in enumerate([i[:-1] for i in _values[:-1]]+[_values[-1]]): print i,e
        for i in indexed: print i
        exit()

        #for i in range(len(_values))
 
        #data = list(itertools.chain.from_iterable(list(itertools.chain([[ii[1] for ii in i[:-1]] for i in indexed[:-1]])) +[indexed[-1][1]]))
        print data, len(data)
        exit(1)

        dataframe = pd.DataFrame(data, columns=['d'])
        print data


        """        

        data = itertools.chain([i[:-1] for i in _values[:-1] + _values[-1]])
        print data, len(data)
       
        #data = pd.DataFrame(np.array(itertools.chain([_values[i][:-1] for i in _values]), dtype=np.float64), columns=['d'])
        #print data
        """


    def plot_singleton(self, dataset, _values):
        fig,ax = plt.subplots(figsize=(2,1))
        ax.set_xticks(range(6))
        ax.set_xticklabels(list(map(str, [1+r for r in range(6)])), fontsize=2)
        ax.set_yticks(np.arange(0.6, 1., 0.1))
        ax.set_yticklabels(list(map(str, np.arange(0.6, 1.0, 0.1))), fontsize=2)
        if not dataset in ['VogtmannE_2016', 'HanniganGD_2017']: ax.set_ylim(0.55, 0.95) 
        else: ax.set_ylim(0.55, 0.85) 
        ax.set_ylabel('AUC', fontsize=2)
        ax.set_xlabel('# Training Datasets', fontsize=2)
        ax.set_xlim(-0.2, 5.2)
        for i in range(len(self.datasets)-1):
            ax.scatter([i]*len(_values[i][:-1]), _values[i][:-1], color='dodgerblue', alpha=0.6, s=2.5)
            ax.scatter([i], _values[i][-1], color='goldenrod', alpha=.6, s=2.5) 
        ax.plot(range(len(self.datasets)-1), [_values[j][-1] for j in range(len(_values))], linestyle='-', color='goldenrod', linewidth=1.2, alpha=0.6)
        #ax.plot([0, 5.], [0.8, 0.8], color='firebrick', linestyle='--', alpha=.5, linewidth=.2)

        #sns.despine(left=False, bottom=False)

        plt.subplots_adjust(left=0.15, bottom=0.24)
        plt.suptitle('%s' % dataset if not dataset in self.utils.data_aliases else self.utils.data_aliases[dataset], fontsize=2) 
        sps = [i.set_linewidth(0.1) for i in ax.spines.itervalues()]
        plt.savefig('%s/%s_progressiveLearning_%s_%s.%s' %(self.save_folder, dataset, self.db, 'average' if self.func=='mean' else 'median', self.fig_fmt), dpi=600, edgecolor='black')
	#, facecolor='white', edgecolor='black', transparent=True)


    def all_singletons(self):
        #sns.set_style('darkgrid')
        sns.set_style("darkgrid", {'grid.linewidth': 0.1})
        _values = self._values(start=2)
        for data in self.datasets:
            self.plot_singleton(data, _values[data]) ## 'CM_rescignocrc', _values['CM_rescignocrc'])


if __name__ == '__main__':

    print 'lascia perdere'
