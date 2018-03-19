#!/usr/bin/env python


import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from utils import lefse_reader, project_color_code
import seaborn as sns
import argparse
import sys


class effect_size_barplot(object):

    def read_params(self, args):
        pa = argparse.ArgumentParser()
        add = pa.add_argument
        add('lefse_file', type=str)
        add('-wi', '--written_inside', action='store_true')
        add('-fmt', '--fig_fmt', type=str, choices=['png', 'svg'], default='png')
        add('-nim', '--number_imp_feat', type=int, default=5)
        add('-wh', '--where', type=str)
        add('--feature_names', type=str, default='Species')
        add('--color_code', type=str,  choices=['periimplantitis', 'crc'])
        add('--save_features', type=str, default=None)
        add('--basic_fields', type=str, default=['sampleID', 'study_condition'], nargs='+')

        add('--title', default='', type=str)
        add('-txt', '--text_on', default=None, type=str, help='Text file with the statistical support.')
        return vars(pa.parse_args())


    def __init__(self, args):        
        self.args = self.read_params(args)
        self.fmt = self.args['fig_fmt']
        self.color_code = project_color_code(self.args['color_code']).color_code


    def save_feats(self, fname, feat_list):

        with open(fname, 'w') as out:
            [out.write(basic + '\n') for basic in self.args['basic_fields']]
            if self.args['feature_names'] != 'Pathways':
                [out.write((f if self.args['feature_names']!='Species' else 's__'+f)+'\n') for f in feat_list]
            else: 
                [out.write((f if self.args['feature_names']!='Pathways' else f.replace('_','-'))+'\n') for f in feat_list]



    def write_statistical_support(self, texton, data):
        data = data[[self.args['feature_names'], 'C', 'p']]
        data.columns = [[self.args['feature_names'], 'class', 'Mann-Withney-Significance']]
        data.to_csv(texton, sep='\t', header=True, index=False)
            

    def set_plot(self, file_name):
        lefse = lefse_reader(file_name)
        lefse.get()
        lefse.select_highest(self.args['number_imp_feat'])
        classes = lefse.classes_
        data_to_plot = pd.DataFrame(\

                                      {'Effect Size':  [e.esize for e in lefse.data]\
                                     , 'feat':  [e.name for e in lefse.data]\
                                     , 'p':    [e.pv for e in lefse.data]\
                                     , 'C':   [e.cl for e in lefse.data]\
                                     , self.args['feature_names']: [(name[:-13]+'_spp' if name.endswith('_unclassified') else name).replace('_',' ')\
				        for name in [(e.name if e.pv >= 0.05 else e.name + ' *') for e in lefse.data]]\
                                     , 'colors': [self.color_code[e.cl] for e in lefse.data]\

                                      })

        data_to_plot.set_index('feat', inplace=True)
        data_to_plot = data_to_plot.loc[lefse.feats]

        if bool(self.args['save_features']):
            self.save_feats(self.args['save_features'], lefse.feats)
        if bool(self.args['text_on']):
            self.write_statistical_support(self.args['text_on'], data_to_plot)

        self.plot(data_to_plot, classes)

        
    def plot(self, data, classes):
         sns.set_style('darkgrid')
         fig,ax_plot = plt.subplots(1,1, figsize=(6,3))
         sns.barplot(data=data, x='Effect Size', y=self.args['feature_names']\
		   , ax=ax_plot, palette=[c for c in data['colors'].tolist()], saturation=0.5)
         ax_plot.yaxis.set_ticks_position('right')
         for lab in ax_plot.get_yticklabels():   ##  for d in dir(lab): print d
             lab.set_style('italic')
             lab.set_size(8)
         plt.subplots_adjust(right=.6, bottom=0.2)
         plt.savefig(self.args['where'] + ('' if self.args['where'].endswith('/') else '/') + '_'.join(classes) \
             + ('' if not self.args['title'] else '_'+self.args['title']) + '_effect_size.' + self.fmt, dpi=400)
        

if __name__ == '__main__':
    a = effect_size_barplot(sys.argv)
    a.set_plot(a.args['lefse_file']) 
