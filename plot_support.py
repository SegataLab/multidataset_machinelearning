#!/usr/bin/env python


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.gridspec as gridspec
import utils
import itertools
import pandas as pd
import glob
import seaborn as sns


## to plot progressive singularly for metaphan:
## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do support_plot -db metaphlan -g0 c:entropy -g1 nt:1000 -g2 nsl:5

## to plot progressive singularly for genefamilies:
## python run.py crc --define study_condition:CRC:control --datasets FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -al rf -do support_plot -db genefamilies -g0 nsl:5 -g1 nt:1000 -g2 df


class training_support(object):
    #save_folder = '../Images/'

    def __init__(self, title, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python, fig_fmt, mean_median, path):
        self.save_folder = (path + ('/' if (not path.endswith('/')) else '') + 'Images/') if path else '../Images/'
        self.utils = utils.usefullfuncs(datasets, path)
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
        self.test_n = self.utils.test_magns(('_'.join(self.db) if isinstance(self.db, list) else self.db), problem, self.algo, self.grid0, self.grid1, self.grid2, start)
        means = dict([(ds,[]) for ds in self.datasets])
        for dataset_on_which in self.datasets:
            transfer_pattern = self.utils.get_all_trans_for_one(self.db, self.algo, problem, self.grid0, self.grid1, self.grid2, dataset_on_which)    
            points = [self.utils._auc(p, start)[0] for p in glob.glob(transfer_pattern) if p.count(dataset_on_which)==1]

            print points, ' ECCCCCOOOOOO......'

            if self.func == 'mean': wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
            else: wavg = [np.median(points)] ## for p in points]
            means[dataset_on_which].append(points + wavg)
            for i in range(2, (len(self.datasets)-1)):
                std_support_pattern = self.utils.get_std_support_pattern(self.db, problem, self.grid0, self.grid1, self.grid2, dataset_on_which, i)
                points = [self.utils._auc(p, start)[0] for p in glob.glob(std_support_pattern)]
                if self.func == 'mean': wavg = [sum([p*self.test_n[dataset_on_which] for p in points])/float(sum([self.test_n[dataset_on_which] for j in range(len(points))]))]
                else: wavg = [np.median(points)] ## for p in points]
                means[dataset_on_which].append(points + wavg)
            last_pattern = self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.algo, problem, self.grid0, self.grid1, self.grid2)
            means[dataset_on_which].append([self.utils._auc(self.utils.get_lodo_resultnames(dataset_on_which, self.db, self.algo, problem, self.grid0, self.grid1, self.grid2), start)[0]])
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

        fig,ax = plt.subplots(figsize=(5,4))

        ax.set_xticks(range(6))
        ax.set_xticklabels(list(map(str, [1+r for r in range(6)])), fontsize=2)
        ax.set_yticks(np.arange(0.5, 1.1, 0.1))
        ax.set_yticklabels(list(map(str, np.arange(0.5, 1.1, 0.1))), fontsize=2)

        #if not dataset in ['VogtmannE_2016', 'HanniganGD_2017']: ax.set_ylim(0.50, 1.) 
        ax.set_ylim(0.5, 1.0) 
        ax.set_ylabel('AUC', fontsize=2)
        ax.set_xlabel('# Training Datasets', fontsize=2)
        ax.set_xlim(-0.2, 5.2)

        for i in range(len(self.datasets)-1):
            ax.scatter([i]*len(_values[i][:-1]), _values[i][:-1], color='dodgerblue', alpha=0.6, s=50)
            ax.scatter([i], _values[i][-1], color='goldenrod', alpha=.6, s=50) 
            
        ax.plot(range(len(self.datasets)-1), [_values[j][-1] for j in range(len(_values))], linestyle='-', color='goldenrod', linewidth=4.5, alpha=0.6)

        plt.subplots_adjust(left=0.15, bottom=0.24)
        plt.suptitle('%s' % dataset if not dataset in self.utils.data_aliases else self.utils.data_aliases[dataset], fontsize=2) 
        sps = [i.set_linewidth(0.1) for i in ax.spines.itervalues()]
        plt.savefig('%s/%s_progressiveLearning_%s_%s.%s' %(self.save_folder, dataset, self.db, 'average' if self.func=='mean' else 'median', self.fig_fmt), dpi=600, edgecolor='black')

	#, facecolor='white', edgecolor='black', transparent=True)
        #plt.close()
        return [_values[i][-1] for i in range(len(_values))]



    def all_singletons(self):

        sns.set_style("darkgrid", {'grid.linewidth': 0.1})
        frame = []
        _values = self._values(start=2)

        for data in self.datasets:
 
            if len(frame)==0:
                frame = pd.DataFrame(self.plot_singleton(data, _values[data]), columns=[data], index=[str(n+1) for n in range(6)])
            else: 
                frame = frame.join(pd.DataFrame(self.plot_singleton(data, _values[data]), columns=[data], index=[str(n+1) for n in range(6)]), how='outer')

        self.frame = frame.T
        self.plot_heatmap()



    def plot_heatmap(self):

        fig = plt.figure(figsize=(5,4))
        gs = gridspec.GridSpec(1,11)
        ax_hm = plt.subplot(gs[:, 0:10])
        ax_cbar = plt.subplot(gs[:, 10])
        self.frame = self.frame.astype('float')

        #print self.frame.index.tolist()
        #print [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets]

        self.frame.index = [(ds if not ds in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets]
  
        sns.heatmap(data=self.frame, cmap='hot', cbar=False, ax=ax_hm, annot=self.frame.values, fmt='.2f', vmin=0.5, vmax=1., linewidth=1., square=True, annot_kws={'fontsize': 9})
        norm_ = matplotlib.colors.Normalize(vmin=.5,vmax=1.)
        cbar = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap='hot', norm=norm_, extend='min', filled=True, drawedges=True)
        ax_cbar.tick_params(labelsize=6)
        ax_cbar.set_ylabel('AUC', size=6)

        ax_hm.set_xlabel('Number of Training Cohorts', fontsize=8)
        ax_hm.set_ylabel('Testing Cohort', fontsize=8)
        ax_hm.set_xticklabels(ax_hm.get_xticklabels(), fontsize=6)
        ax_hm.set_yticklabels(ax_hm.get_yticklabels(), fontsize=6)
        plt.suptitle('Median Training Increasing Progression', fontsize=6)
        plt.subplots_adjust(left=0.22, bottom=0.24)
        cbar_shape1 = ax_cbar.get_position()
        cbar_shape2 = [cbar_shape1.x0, cbar_shape1.y0, cbar_shape1.width * 0.6, cbar_shape1.height]
        ax_cbar.set_position(cbar_shape2)
        plt.savefig('%sMedianTableProgressiveLearning_%s.%s' %(self.save_folder, self.db.upper(), self.fig_fmt), dpi=600, edgecolor='black')        

  



if __name__ == '__main__':

    print 'lascia perdere'
