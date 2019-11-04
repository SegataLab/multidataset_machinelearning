#!/usr/bin/env python



import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from utils import project_color_code
import seaborn as sns
import argparse
import sys
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'



class sample_property_barplot(object):

    def read_params(self, args):

        pa = argparse.ArgumentParser()
        add = pa.add_argument
						
        add('-p', '--property', type=str, default='percentage_of_mapped_by_metaphlan2'\
            , help='You have to specify the column in the table you want to plot.')

        add('-wi', '--written_inside', action='store_true')
        add('-fmt', '--fig_fmt', type=str, choices=['png', 'svg'], default='png')
        add('-wh', '--where', type=str, default='')

        cmaps = ['plasma', 'viridis', 'inferno', 'magma']
        add('-cm','--cmap', type=str, default='plasma', choices=cmaps + [m+'_r' for m in cmaps])
        add('-c','--color', type=str, default='dodgerblue')

        ##add('--color_code', type=str,  choices=['periimplantitis', 'crc'], default='None')
        add('--save_features', type=str, default=None)
        add('--basic_fields', type=str, default=['sampleID', 'study_condition'], nargs='+')

        add('--maxv', type=float, default=100.)
        add('--minv', type=float, default=0.)
        add('--horizontal', action='store_true')
        add('--fontsize', type=int, default=4)

        add('--bar_width', type=float, default=0.5)
        add('--bar_height', type=float, default=0.5)
        add('--italic', action='store_true')
        add('--despine', action='store_true')
 
        add('--sort', type=str, help='Filename of a file with the right order')
        add('--facecolor', type=str, default='white')
        add('--title', default='', type=str)

        add('-l', '--log_scale', action='store_true')

        return vars(pa.parse_args())




    def __init__(self, args):        

        self.args = self.read_params(args)
        self.fmt = self.args['fig_fmt']

        ##self.color_code = project_color_code(self.args['color_code']).color_code

        self.data = pd.read_csv(sys.stdin, sep='\t', index_col=0, header=None)
        self.cmap = vars(matplotlib.cm)[self.args['cmap']] 
        self.title = (self.args['title'] if self.args['title']\
                     else self.data.loc['dataset_name'].tolist()[0] + '_' + self.args['property']).replace('Percentage','%')



    def plot(self):

         float_and_maybe_log = lambda a : float(a) if not self.args['log_scale'] else np.log(float(a))

         sns.set_style(self.args['facecolor'])
         fig,ax_plot = plt.subplots(1,1, figsize=(20,4) if not self.args['horizontal'] else (4,20))

         if self.args['horizontal']: ax_plot.set_xlim(self.args['minv'], float_and_maybe_log(self.args['maxv']))
         else: ax_plot.set_ylim(self.args['minv'], float_and_maybe_log(self.args['maxv']))

         #print ax_plot.get_ylim(), ' questo e il limite della x'

         #ax_plot.set_ylim(-(self.args['bar_height']/2), len(ax_plot.patches)+(self.args['bar_height']/2))
         #ax_plot.set_xlim(0,self.args['x_lim'])

         self.data.columns = self.data.loc['sampleID'].tolist()

         if self.args['sort']: 
             f = open(self.args['sort'])
             prev_samples = self.data.columns.tolist() 
             self.new_data = self.data[f.readline().rstrip().split()]
             f.close()

         #y=list(map(float, self.data.loc[self.args['property']].tolist())) \
         #                if not self.args['horizontal'] \
         #                else self.data.loc['sampleID'].tolist()

         #x=self.data.loc['sampleID'].tolist()\
         #                if not self.args['horizontal'] \
         #                else list(map(float, self.data.loc[self.args['property']].tolist()))

         #print y, ' e ipsilon'
         #print x, ' e ix'
 
         #  float_and_maybe_log = lambda a : float(a) if not self.args['log_scale'] else np.log(float(a))

         y_=list(map(float_and_maybe_log, self.new_data.loc[self.args['property']].tolist())) \
                         if not self.args['horizontal'] \
                         else self.new_data.loc['sampleID'].tolist()

         x_=self.new_data.loc['sampleID'].tolist()\
                         if not self.args['horizontal'] \
                         else list(map(float_and_maybe_log, self.new_data.loc[self.args['property']].tolist()))

         sns.barplot(\
	             y=y_\
                   , x=x_\
		   , ax=ax_plot\
                   , color=self.args['color']\
                    )

                   #, palette=[self.cmap(prp) for prp in list(map(float, self.data.loc[self.args['property']].tolist()))])

#, orient=('v' if not self.args['horizontal'] else 'h')) #\
                   #, palette=[c for c in self.data.loc['color'].tolist()], saturation=0.5)

         for bar in ax_plot.patches: 
             if not self.args['horizontal']: bar.set_width(self.args['bar_width'])
             else: bar.set_height(self.args['bar_height'])

         if not self.args['horizontal']:
             ax_plot.xaxis.set_ticks_position('top')
             plt.setp(ax_plot.get_xticklabels(), rotation=90) #, ha='left')
             plt.setp(ax_plot.get_xticklabels(), fontsize=self.args['fontsize'])
         #ax_plot.yaxis.set_ticks_position('right')

         for label in (ax_plot.get_yticklabels() \
             if not self.args['horizontal'] \
             else ax_plot.get_xticklabels()):
                 if self.args['italic']: label.set_style('italic')
                 label.set_size(self.args['fontsize'])

         sps = [spine.set_linewidth(0.5) for spine in ax_plot.spines.itervalues()]
         if self.args['despine']: sns.despine(left=False, bottom=False)
         if self.args['horizontal']: ax_plot.set_xlabel(self.args['property'])
         else: ax_plot.set_ylabel(self.args['property'])

         plt.subplots_adjust(top=0.6)
         plt.savefig(self.title + '.' + self.fmt, dpi=400)
         #self.args['where'] \
         #    + ('' if self.args['where'].endswith('/') else '/')\
         #    + self.title + '.' + self.fmt, dpi=400)
        


if __name__ == '__main__':
    a = sample_property_barplot(sys.argv)
    a.plot()
