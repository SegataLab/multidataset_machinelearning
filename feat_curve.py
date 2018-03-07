#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import utils
import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches
matplotlib.rcParams['svg.fonttype'] = 'none'
import itertools

## python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do feat_curve -g0 c:entropy -g1 nt:1000 -g2 nsl:5
## python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -do feat_curve -g0 nsl:5 -g1 nt:1000 -g2 SEL

class feat_curve(object):
    save_folder = '../Images/'

    def __init__(self, title, datasets, db, defined_problem, algo, grid0, grid1, grid2, fig_fmt):
        sns.set(style='white')
        self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.db = db
        self.title = title
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]

        self.lodo = False
        self.include_cross_validation = False

        self.algo = algo
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid = '_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2

        self.curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
		, sorted([v for v in self.get_cross_values(data)], key = lambda ls : ls[0])[1:]) for data in self.datasets])

        self.lodo_curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
		, sorted([l for l in self.get_lodo_values(data)], key = lambda ls: ls[0])[1:]) for data in self.datasets])
	        
        #for c in self.curves: print c, self.curves[c]
        #exit(1)

        if self.include_cross_validation:
            self.cross_validation = sorted([auc for auc in self.utils.all_auc(self.utils.get_cross_validation(\
		self.db, self.algo, self.tests[0], self.grid0, self.grid1, self.grid2))], key = lambda ls : ls[0])[1:]

            self.curves['Cross-Validation'] = self.cross_validation         
            self.lodo_curves['Cross-Validation'] = self.cross_validation

        self.index = self.curves.keys() + (['Cross-Validation'] if self.include_cross_validation else [])

        self.scores = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data], [t[1] for t in self.curves[data]]) for data in self.index])
        self.lodo_scores = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data], [l[1] for l in self.lodo_curves[data]]) for data in self.index])
        self.std = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data], [t[-1] for t in self.curves[data]]) for data in self.index])
        self.datasets = [(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]) for data in self.datasets]
        self.n_feat = [f[0] for f in self.curves[self.datasets[0]]][:-1] + ['all']

        self.symbols = ['blue', 'crimson', 'limegreen', 'dodgerblue', 'goldenrod', 'lightcoral', 'darkcyan', 'black']
        self.signs = 'Px*s,^+o'

        self.markers = self.get_markers()
        self.signers = self.get_signs()

        if not self.lodo:				 
            self.plottable_data = dict([(d, [\
			np.array(self.scores[d], dtype=np.float64) + np.array(self.std[d], dtype=np.float64)\
		      , np.array(self.scores[d], dtype=np.float64)\
		      , np.array(self.scores[d], dtype=np.float64) - np.array(self.std[d], dtype=np.float64)])\
			      for d in self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])])
        else:
            self.plottable_data = dict([(d, [\
                        np.array(self.lodo_scores[d], dtype=np.float64) + np.array(self.std[d], dtype=np.float64)\
                      , np.array(self.lodo_scores[d], dtype=np.float64)\
                      , np.array(self.lodo_scores[d], dtype=np.float64) - np.array(self.std[d], dtype=np.float64)])\
                              for d in self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])])
							
        #for k in self.plottable_data: print k, self.plottable_data[k], ' MMMMMeRR'
        #exit(1)

        scatter, times = [], []
        fig, ax_ = plt.subplots(figsize=(5,4)) 

        for e,data in enumerate(self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])): # if not self.lodo else [])):
            line_sub, line_, line_sup = self.plottable_data[data]
            #times.append( ax_.plot( range(7), line_sub, alpha=0.6, color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], linewidth=1.))
            times.append( ax_.plot( range(7), line_, alpha=0.6, color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], linewidth=1., linestyle='--'))
            #times.append( ax_.plot( range(7), line_sup, alpha=0.6, color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], linewidth=1.))

            scatter.append(sns.regplot(ax=ax_, x='nf', y='auc', data=pd.DataFrame({'nf': range(7), 'auc': self.scores[data] if not self.lodo else self.lodo_scores[data]}), marker='o'\
	        , color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], fit_reg=False, label=data, scatter_kws={'s': 60})) 

            #times.append(sns.tsplot(np.array(self.plottable_data[data], dtype=np.float64), color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]]\
	    #	, alpha=0.5, legend=False))
        
        ax_.set_ylim(0.5,1.0)
        ax_.set_xlim(-0.2, 6.2)
        ax_.xaxis.set_ticks(range(7))
        ax_.set_xticklabels(['2','4','8','16','32','64','all'], fontsize=4)
        ax_.yaxis.set_ticks(np.arange(0.5, 1., 0.1))
        #yticks = [('%.2f' %i if e%2!=0. else '') for e,i in enumerate(np.arange(0.5, 1.0, 0.05))]
        ax_.set_yticklabels(list(map(str, np.arange(0.5, 1., 0.1))), fontsize=4)
        ax_.set_ylabel('AUC', fontsize=4)
        ax_.set_xlabel('# of features', fontsize=4)

        sps = [i.set_linewidth(0.5) for i in ax_.spines.itervalues()]
        leg = plt.legend(ncol=4, loc=9, fontsize=3, frameon=True, edgecolor='black', mode='expand')
        plt.suptitle(('Cross-Validation' if not self.lodo else 'Leave-One-Data-Out') + ' Progressively Increasing Feature Number', fontsize=3)
        plt.subplots_adjust(bottom=0.2)

        if not self.lodo:
            plt.savefig('%s/Feature_Curve%s_%s.%s' %(self.save_folder, ('_WithCross' if self.include_cross_validation else '')\
			, 'cross', fig_fmt), dpi=600, facecolor='w', frameon=False, edgecolor='w')
        else:
            print ' sto per salvare : ', '%s/Feature_Curve_%s_%s.%s' %(self.save_folder, ('WithCross' if self.include_cross_validation else '')\
                        , 'lodo', fig_fmt)
            plt.savefig('%s/Feature_Curve%s_%s.%s' %(self.save_folder, ('_WithCross' if self.include_cross_validation else '')\
			, 'lodo', fig_fmt), dpi=600, facecolor='w', frameon=False, edgecolor='w')


    def get_markers(self):
        markers, j = {}, 0
        if self.title == 'crc': 
            for data in self.datasets  + (['Cross-Validation'] if self.include_cross_validation else []):
                data_ = data if data not in self.utils.data_aliases else self.utils.data_aliases[data]
                markers[data_] = self.symbols[j]
                j += 1
            return markers


    def get_signs(self):
        markers, j = {}, 0
        if self.title == 'crc':
            for data in self.datasets  + (['Cross-Validation'] if self.include_cross_validation else []):
                data_ = data if data not in self.utils.data_aliases else self.utils.data_aliases[data]
                markers[data_] = self.signs[j]
                j += 1
            return markers


    def get_cross_values(self, dataset): return self.utils.all_auc_from_transfer([dataset, dataset], self.db, self.algo, self.tests[0], self.grid0, self.grid1, self.grid2)
    def get_lodo_values(self, dataset): return self.utils.all_auc_from_lodo(dataset, self.db, self.algo, self.tests[0], self.grid0, self.grid1, self.grid2)


if __name__ == '__main__':
    print 'aaaaaaaaaaaarrrrrrgggghhhh!!'
