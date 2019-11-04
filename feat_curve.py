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

## python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -do feat_curve -g0 c:entropy -g1 nt:1000 -g2 special


### yes



class feat_curve(object):

    #save_folder = '../Images/'


    def __init__(self, title, datasets, db, defined_problem, algo, grid0, grid1, grid2, fig_fmt, path):

        sns.set(style='whitegrid')  #'white')
        self.save_folder = (path + ('/' if (not path.endswith('/')) else '') + 'Images/') if path else '../Images/'
        self.utils = utils.usefullfuncs(datasets, path)
        self.datasets = datasets
        self.db = db
        self.title = title
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]

        self.lodo = False
        self.include_cross_validation = True
        self.algo = algo
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid = '_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2

        self.symbols = ['blue', 'crimson', 'limegreen', 'dodgerblue', 'goldenrod', 'lightcoral', 'darkcyan', 'black']
        self.signs = 'Px*s,^+o'
        self.markers = self.get_markers()
        self.signers = self.get_signs()
        self.fig_fmt = fig_fmt

        if self.db[0] == 'genefamilies':
            self.plot_8192_genefamilies()

        elif self.db[0] == 'metaphlan':
            self.plot_64_species()



    def plot_8192_genefamilies(self):

        self.include_cross_validation = False
        self.lodo = True
        self.curves = {}
        self.lodo_curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
                , sorted([l for l in self.get_lodo_values(data)], key = lambda ls: ls[0])) for data in self.datasets])
        self.index = self.lodo_curves.keys() + (['Cross-Validation'] if self.include_cross_validation else [])

        self.lodo_scores = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data], [l[1] for l in self.lodo_curves[data]]) for data in self.index])
        self.datasets = [(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]) for data in self.datasets]
        self.plottable_data = dict([(d, [np.array(self.lodo_scores[d], dtype=np.float64)]) for d in self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])])
        self.cross_validation = False
        self.lodo = True
        self.num_thresholds = 13
        self.x_ticks = ['2','4','8','16','32','64','128','256','512','1024','2048','4096','8192']
        self.print_even_heatmap_down_genes()



    
    def print_even_heatmap_down_metaphlan(self):

        #figb, axb = plt.subplots(figsize=(4,4))				
        data = pd.DataFrame({'nf': self.x_ticks, str(self.datasets[0]): self.scores[self.datasets[0]] if not self.lodo else self.lodo_scores[self.datasets[0]]})
        for dat in self.datasets[1:] + (['Cross-Validation'] if self.include_cross_validation else []):
            data = data.join(pd.DataFrame({str(dat): self.scores[dat] if not self.lodo else self.lodo_scores[dat]}))
 
        data.index = data['nf']
        del data['nf']
        datat = data.T
        datat.to_csv('Metaphlan_FeatureSelection_'+('CrossValidation' if not self.lodo else 'LODO')+'_Numbers.csv', sep='\t', header=True, index=True)

        #print '\t\t' + '\t'.join([str(c) for c in datat.columns.tolist()]) 
        #for i in datat.index.tolist():
        #    print '\t'.join([i] + list(map(str, datat.loc[i].tolist())))
        #print datat        

        #sns.heatmap(data=datat, ax=axb, annot=True, cmap='hot', cbar=True, annot_kws={'fontsize': 4}, square=True)
        #axb.tick_params(labelsize=4)
 
        #plt.suptitle(' Metaphlan Signature Scores ' + ('in Cross-Validation' if not self.lodo else 'in LODO'))
        #plt.subplots_adjust(bottom=0.2, left=0.4)
        #plt.savefig('%sFeature_CurveNUMBERS_%s_%s.%s' %(self.save_folder, self.db[0]\
        #    , 'inCrossValidation' if not self.lodo else 'inLODO', self.fig_fmt)\
        #    , dpi=400, facecolor='w', frameon=False, edgecolor='w')

        largest = np.array([datat['all'].tolist() for i in range(datat.values.shape[1])], dtype=np.float64).T 
        diff = pd.DataFrame(columns=datat.columns.tolist(), index=datat.index.tolist(), data=(largest - datat.values)*1.)
        diff.to_csv('Metaphlan_FeatureSelection_'+('CrossValidation' if not self.lodo else 'LODO')+'_Differences.csv', sep='\t', header=True, index=True)




    def print_even_heatmap_down_genes(self):
        
        #figb, axb = plt.subplots(figsize=(4,4))
        data = pd.DataFrame({'nf': self.x_ticks, str(self.datasets[0]): self.lodo_scores[self.datasets[0]]})
        for dat in self.datasets[1:] + (['Cross-Validation'] if self.include_cross_validation else []):
            data = data.join(pd.DataFrame({str(dat): self.lodo_scores[dat]}))
        data.index = data['nf']
        del data['nf']
        datat = data.T
        datat.to_csv('GeneFamilies_FeatureSelection_Numbers.csv', sep='\t', header=True, index=True)

        #sns.heatmap(data=datat, ax=axb, annot=True, cmap='hot', cbar=True, annot_kws={'fontsize': 3}, square=True)
        #axb.tick_params(labelsize=4)


        #plt.suptitle(' Gene-Families Signature Scores ' + ('in Cross-Validation' if not self.lodo else 'in LODO'))
        #plt.subplots_adjust(bottom=0.2, left=0.4)
        #plt.savefig('%sFeature_CurveNUMBERS_%s_%s.%s' %(self.save_folder, self.db[0]\
        #          , 'inCrossValidation' if not self.lodo else 'inLODO', self.fig_fmt)\
        #          , dpi=400, facecolor='w', frameon=False, edgecolor='w')
 

        largest = np.array([datat['8192'].tolist() for i in range(datat.values.shape[1])], dtype=np.float64).T
        diff = pd.DataFrame(columns=datat.columns.tolist(), index=datat.index.tolist(), data=(largest - datat.values)*1.)
        #print diff
        diff.to_csv('GeneFamilies_FeatureSelection_Differences.csv', sep='\t', header=True, index=True)

        #diff = np.array([[]])     np.zeros([datat.values.shape], dtype=np.float64)
        





    def plot_64_species(self):

        self.curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
                , sorted([v for v in self.get_cross_values(data)], key = lambda ls : ls[0])[1:]) for data in self.datasets])

        if self.db[0] == 'genefamilies':
            self.lodo_curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
		, sorted([l for l in self.get_lodo_values(data)], key = lambda ls: ls[0])) for data in self.datasets]) ## [1:] dopo ls [0] per metaphlan

        elif self.db[0] == 'metaphlan':
            self.lodo_curves = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
                , sorted([l for l in self.get_lodo_values(data)], key = lambda ls: ls[0])[1:]) for data in self.datasets])

        if self.include_cross_validation:
            self.cross_validation = sorted([auc for auc in self.utils.all_auc(self.utils.get_cross_validation(\
		self.db, self.algo, self.tests[0], self.grid0, self.grid1, self.grid2))], key = lambda ls : ls[0])[1:]

            self.curves['Cross-Validation'] = self.cross_validation         
            self.lodo_curves['Cross-Validation'] = self.cross_validation
            self.index = self.lodo_curves.keys() + (['Cross-Validation'] if self.include_cross_validation else [])

        self.scores = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
            , [l[1] for l in self.curves[data]]) for data in self.index ])
        self.lodo_scores = dict([(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]\
            , [l[1] for l in self.lodo_curves[data]]) for data in self.index ])
        
        self.datasets = [(data if data not in self.utils.data_aliases else self.utils.data_aliases[data]) for data in self.datasets]
        self.plottable_data = dict([(d, [np.array(self.lodo_scores[d] if self.lodo else self.scores[d], \
		dtype=np.float64)]) for d in self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])])
        self.cross_validation = True

        #self.lodo = True
        self.num_thresholds = 7 
        self.x_ticks = ['2','4','8','16','32','64','all']

        self.print_even_heatmap_down_metaphlan()

        #self.n_feat = [f[0] for f in self.curves[self.datasets[0]]][:-1] + ['all']
        #if not self.lodo:				 
        #    self.plottable_data = dict([(d, [\
	#		np.array(self.scores[d], dtype=np.float64) + np.array(self.std[d], dtype=np.float64)\
	#	      , np.array(self.scores[d], dtype=np.float64)\
	#	      , np.array(self.scores[d], dtype=np.float64) - np.array(self.std[d], dtype=np.float64)])\
			  #    for d in self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])])
        #else:
        #for k in self.plottable_data: print k, self.plottable_data[k], ' MMMMMeRR'
        #exit(1)




    def do_the_plot(self):

        sns.set(style='ticks')

        scatter, times = [], []
        fig, ax_ = plt.subplots(figsize=(5,4)) 

        for e,data in enumerate(self.datasets + (['Cross-Validation'] if self.include_cross_validation else [])): # if not self.lodo else [])):
            line_ = self.plottable_data[data][0]
            times.append( ax_.plot( range(self.num_thresholds), line_, alpha=0.6, color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], linewidth=1., linestyle='--'))
            scatter.append(sns.regplot(ax=ax_, x='nf', y='auc', data=pd.DataFrame({'nf': range(self.num_thresholds), 'auc': self.scores[data] if not self.lodo else self.lodo_scores[data]}), marker='o'\
	        , color=self.markers[data if data not in self.utils.data_aliases else self.utils.data_aliases[data]], fit_reg=False, label=data, scatter_kws={'s': 60})) 
        
        ax_.set_ylim(0.5,1.0)
        ax_.set_xlim(-0.2, self.num_thresholds - 0.8)
        ax_.xaxis.set_ticks(range(self.num_thresholds))
        ax_.set_xticklabels(self.x_ticks, fontsize=4)
        ax_.yaxis.set_ticks(np.arange(0.5, 1.0, 0.1))

        #ax_.set_ylim(0.5,1.0)

        ax_.set_yticklabels(list(map(str, np.arange(0.5, 1., 0.1))), fontsize=4)
        ax_.set_ylabel('AUC', fontsize=4)
        ax_.set_xlabel('# of features', fontsize=4)
        sps = [i.set_linewidth(0.5) for i in ax_.spines.itervalues()]
        leg = plt.legend(ncol=4, loc=9, fontsize=3, frameon=True, edgecolor='black', mode='expand')
        plt.suptitle(('Cross-Validation' if not self.lodo else 'Leave-One-Data-Out') + ' Progressively Increasing Feature Number' + ('' if not self.db=='genefamilies' else ' (Selected Genefamilies)'), fontsize=3)
        plt.subplots_adjust(bottom=0.2)

        if not self.lodo:
            plt.savefig('%s/Feature_Curve%s_%s.%s' %(self.save_folder, ('_WithCross' if self.include_cross_validation else ('' if self.db=='metaphlan' else 'WithGeneFamilies'))\
			, 'cross', self.fig_fmt), dpi=600, facecolor='w', frameon=False, edgecolor='w')
        else:
            plt.savefig('%s/Feature_Curve%s_%s.%s' %(self.save_folder, ('_WithCross' if self.include_cross_validation else ('' if self.db=='metaphlan' else 'WithGeneFamilies'))\
			, 'lodo', self.fig_fmt), dpi=600, facecolor='w', frameon=False, edgecolor='w')



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

    print('aaaaaaaaaaaarrrrrrgggghhhh!!')
