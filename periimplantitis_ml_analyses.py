#!/usr/bin/env python


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
import argparse
import utils
import itertools
from feat_heat import feature_heatmap
from feat_curve import feat_curve
import seaborn as sns;


class CM_lilt_paper(object):

    save_folder = '/scratchCM/users/paolo.manghi/peri-implantitis/' #'paologhensi_and_paolomanghi/Figure1/'
    titles = ['Microbial Species', 'Microbial Pathways'] #, 'Species and Pathways'] 
    tests = ['Healthy on Periimplantitis', 'Healthy on Mucositis', 'Periimplantitis on Mucositis']
 
    metaphlan_hpi = '%shealthy_vs_periimplantitis_metaml_metaphlan.txt' % save_folder
    metaphlan_hmu = '%shealthy_vs_mucositis_metaml_metaphlan.txt' % save_folder
    metaphlan_pimu = '%speriimplantitis_vs_mucositis_metaml_metaphlan.txt' % save_folder
    #
    pathways_hpi = '%shealthy_vs_periimplantitis_metaml_genefam.txt' % save_folder
    pathways_hmu = '%shealthy_vs_mucositis_metaml_genefam.txt' % save_folder
    pathways_pimu = '%speriimplantitis_vs_mucositis_metaml_genefam.txt' % save_folder
    #
    metaphlan_and_pathways_hpi = '%shealthy_vs_periimplantitis_metaml_pathways.txt' % save_folder
    metaphlan_and_pathways_hmu = '%shealthy_vs_mucositis_metaml_pathways.txt' % save_folder
    metaphlan_and_pathways_pimu = '%speriimplantitis_vs_mucositis_metaml_pathways.txt' % save_folder

    trip_one = [metaphlan_hpi, metaphlan_hmu, metaphlan_pimu]
    trip_two = [pathways_hpi, pathways_hmu, pathways_pimu]
    trip_three = [metaphlan_and_pathways_hpi, metaphlan_and_pathways_hmu, metaphlan_and_pathways_pimu]


    def __init__(self):

        self.utils = utils.usefullfuncs(['CM_periimplantitis'], self.save_folder)
        self.which_one_to_use = 0

        #self.define()

        sns.set_style('white')
        self.fmt = 'svg'
        self.cmap = 'Wistia_r'

        self.header = lambda t : t.split('_learner')[0][15:]
        self.vertical_limits = [[0.5, 1.0], [.5, 1.], [.5, 1.]]

        self.basis = self.titles
        self.title = 'Periimplantitis ML Analyses'

        self.num_thresholds = 7
        self.x_ticks = [ '2','4','8','16','32','64','all' ]
        self.symbols = [ 'deepskyblue', 'crimson', 'limegreen', 'dodgerblue'\
                       , 'goldenrod', 'lightcoral', 'forestgreen', 'mediumblue', 'brown' ]

        self.define()
        #self.plot_64_species()
                    


    def define(self):
        self.my_data = [ r[self.which_one_to_use] for r in \
               [ self.trip_one, self.trip_two ]]
##, self.trip_three ]]
        self.all_my_data = self.trip_one + self.trip_two #+ self.trip_three 
        
        self.basic_auc = dict([(filename, self.utils._auc(filename, 2)) for filename in self.all_my_data])
        self.most_basic_auc = dict([(filename, self.utils._auc(filename, 2)) for filename in self.my_data])

        self.all_auc = dict([(filename, self.utils.all_auc(filename)) for filename in self.all_my_data])
        self.all_most_basic_auc = dict([(filename, self.utils.all_auc(filename)) for filename in self.my_data])

        self.curves = dict([(data, sorted([v for v in self.all_most_basic_auc[data]]\
                    , key = lambda ls : ls[0])[1:]) for data in self.all_most_basic_auc])

        self.plottable_data = dict([(d, [np.array(self.curves[d], dtype=np.float64)]) \
               for d in self.all_most_basic_auc])
        self.markers = self.get_markers()



    def change_test_type(self, tp):
        self.which_one_to_use = tp
        self.define()



    def do_the_plot(self):

        scatter, times = [], []
        fig, ax_ = plt.subplots(figsize=(5,4))

        for e,data in enumerate(self.my_data):

            #print ' perche tutto plottable data con lui saebber: ', self.plottable_data[data]

            line_ = [v[1] for v in self.plottable_data[data][0]]

            #print data, [v[1] for v in self.plottable_data[data][0]], ' IIIIIIII'
            #print 'perche ora questo e my data', self.my_data
            #print line_, ' this is line'
            #print data, ' MENTRE, L HEADER CORRETTO SAREBBE:  ', self.basis[e]

            times.append( ax_.plot( range(self.num_thresholds), line_, alpha=0.6, color=self.markers[data], linewidth=1., linestyle='--'))

            scatter.append(sns.regplot(ax=ax_, x='nf', y='auc', data=pd.DataFrame({'nf': range(self.num_thresholds), 'auc': [a[1] for a in self.curves[data]]}), marker='o'\
                , color=self.markers[data], fit_reg=False, label=self.basis[e], scatter_kws={'s': 60}))

        ax_.set_ylim(self.vertical_limits[self.which_one_to_use][0], self.vertical_limits[self.which_one_to_use][1])
        ax_.set_xlim(-0.2, self.num_thresholds - 0.8)
        ax_.xaxis.set_ticks(range(self.num_thresholds))
        ax_.set_xticklabels(self.x_ticks, fontsize=4)

        ax_.yaxis.set_ticks(np.arange(self.vertical_limits[self.which_one_to_use][0], self.vertical_limits[self.which_one_to_use][1], 0.1))
        ax_.set_yticklabels(list(map(str, np.arange(self.vertical_limits[self.which_one_to_use][0], self.vertical_limits[self.which_one_to_use][1], 0.1))), fontsize=4)

        ax_.set_ylabel('AUC', fontsize=4)
        ax_.set_xlabel('Number of Features Used', fontsize=4)
        sps = [i.set_linewidth(0.5) for i in ax_.spines.itervalues()]
 
        leg = plt.legend(ncol=4, loc=9, fontsize=3, frameon=True, edgecolor='black', mode='expand')

        plt.suptitle(('minimal microbial or pathway signature in the %s prediction' % self.tests[self.which_one_to_use]).title(), fontsize=8)
        plt.subplots_adjust(bottom=0.2)
        plt.savefig('%s%s_Feature_Curve_%s.%s' %(self.save_folder, self.title.replace(' ','_'), self.tests[self.which_one_to_use], self.fmt)\
            , dpi=600, facecolor='w', frameon=False, edgecolor='w')



    

    def auc_heatmap(self):
        sort_ = lambda t : [t[2], t[1], t[0]]   

        data = pd.DataFrame(data=[[self.basic_auc[element][0] for element in triplette] for triplette in \
             [ self.trip_one, self.trip_two ]]\
             , columns = self.tests, index = self.titles)

        fig = plt.figure(figsize=(5,10)) 
        gs = gridspec.GridSpec(1,4)
        self.ax_plot = plt.subplot(gs[:,:-1])
        self.ax_cbar = plt.subplot(gs[:,-1])

        sns.heatmap(data=data, ax=self.ax_plot, vmin=0.5, vmax=1.0, cmap=self.cmap, cbar=False\
                  , linewidths=.0, linecolor='white', square=True, annot=True, fmt='.2f'\
                  , xticklabels=True, yticklabels=True)
 
        self.norm_ = matplotlib.colors.Normalize(vmin=.5, vmax=1.)
        self.cbar = matplotlib.colorbar.ColorbarBase(self.ax_cbar, cmap=self.cmap\
                  , norm=self.norm_, extend='min', filled=True, drawedges=False)

        self.cbar.set_ticks(list(np.arange(0.5,1.0,0.1))+[1.])
        self.cbar.set_ticklabels(list(map(str,np.arange(0.5,1.0,.1)))+['1.0'])

        self.ax_cbar.set_ylabel('Area-Under-the-Curve', size=6)
        self.ax_cbar.tick_params(labelsize=6)
        self.cbar.outline.set_edgecolor('white')
        self.cbar.outline.set_linewidth(.0)

        #plt.subplots_adjust(left=.6)
        plt.subplots_adjust(left=.6, top=0.9, bottom=0.1)

        ax_pos0 = self.ax_plot.get_position()  ##a  x_pos1 = [ax_pos0, ax_pos0, ax_pos0, ]

        cbarpos0 = self.ax_cbar.get_position() 
        cbarpos1 = [cbarpos0.x0, cbarpos0.y0 + 0.32, cbarpos0.width * 0.5, cbarpos0.height * 0.2] ## cbarpos0.height]
        self.ax_cbar.set_position(cbarpos1)
        ##plt.subplots_adjust(left=.6, top=0.9, bottom=0.1)  ##bottom=0.18, left=.6, )#right=5.)

        plt.setp(self.ax_plot.get_yticklabels(), fontsize=6)
        plt.setp(self.ax_plot.get_xticklabels(), rotation=70, ha='right', fontsize=6)
        plt.suptitle(self.title, fontsize=6)
        plt.savefig(self.save_folder+self.title.replace(' ','_')+'_AUCs' + '.' + self.fmt, dpi=400)   


    def get_markers(self):
        markers, j = {}, 0
        for data in self.my_data:
            markers[data] = self.symbols[j]
            j += 1
        return markers


    def select_features(self, howmany, filename):

        feats = list()

        line, skip = '_', 0
        f = open(filename)
        while not line.startswith('Feature'): 
            line, skip=f.readline(), skip+1
        h = 0
        while h<howmany:
            h += 1
            feats.append( ' - '.join( f.readline().split()[0:2] + [self.header(filename)] ) )

        return feats



if __name__ == '__main__':

    cm = CM_lilt_paper()
    cm.auc_heatmap()

 #   cm.do_the_plot()
 
 #   for i in [1,2]:
 #       cm.change_test_type(i)
 #       cm.do_the_plot()

    #cm.best_model_crc_over_control()    


