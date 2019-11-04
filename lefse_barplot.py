#!/usr/bin/env python3



import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from utils import \
    lefse_reader, project_color_code
import seaborn as sns
import argparse
import sys
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'



class effect_size_barplot(object):

    def read_params(self, args):

        pa = argparse.ArgumentParser()
        add = pa.add_argument

        add('lefse_file', type=str)

        add('-wi', '--written_inside', action='store_true')
        add('-fmt', '--fig_fmt', type=str, choices=['png', 'svg'], default='svg')
        add('-nim', '--number_imp_feat', type=int, default=5, help='If you specify -1, it s like taking all')
        add('-wh', '--where', type=str)

        add('-of', '--outfile', type=str)
       

        add('--feature_names', type=str, default='Species')
        add('--color_code', type=str,  choices=['periimplantitis', 'crc'])
        add('--save_features', type=str, default=None)
        add('--basic_fields', type=str, default=['sampleID', 'study_condition'], nargs='+')

        add('--bar_height', type=float, default=.4)
        add('--x_lim', type=float, default=5)
        add('--x_start', type=float, default=2)

        add('--facecolor', type=str, default='whitegrid')
        add('--title', default='', type=str)

        add('-es', '--effect_size', action='store_true', help='Write in the text also the effect size.')
        add('-txt', '--text_on', default=None, type=str, help='Text file with the statistical support.')

        add('--feature_dictionary', type=str, default=False)
        add('--right_space', action='store_true')
        add('--read_pwy', action='store_true')
 
        add('--only_sig', default=None, type=float)
        add('--split_based_on', default=None, type=str)

        #add('--feature_dictionary_coded', type=str, default=False)

        return vars(pa.parse_args())



    def __init__(self, args):        
        self.args = self.read_params(args)
        self.fmt = self.args['fig_fmt']
        self.color_code = project_color_code(self.args['color_code']).color_code



    def save_feats(self, fname, feat_list):

        with open(fname, 'w') as out:

            #for f in feat_list: print f

            [out.write(basic + '\n') for basic in self.args['basic_fields']]
            #if self.args['feature_names'] != 'Pathways':
            [out.write(     (f if self.args['feature_names']!='Species' else 's__'+f)  \
                 +'\n') for f in feat_list]
            #else: 
            #    [out.write((f if self.args['feature_names']!='Pathways')+'\n') for f in feat_list]



    def write_statistical_support(self, texton, data):
        data = data[[self.args['feature_names'], 'C', 'p'] + ([] if not self.args['effect_size'] else ['Effect Size'])]
        data.columns = [[self.args['feature_names'], 'class', 'Mann-Withney-Significance'] + ([] if not self.args['effect_size'] else ['Effect Size'])]
        data.to_csv(texton, sep='\t', header=True, index=False)



    def set_plot(self, file_name):

        print(self.args['only_sig'], ' CAREFULL; THIS IS ONLY SIG')


 #       if self.args['feature_dictionary']:
 #            with open(self.args['feature_dictionary']) as fd:
 #                func_dict = dict([( l.rstrip().split()[0].replace('-', '_')\
 #                , ' '.join(l.rstrip().split()[1:])\
 #                ) for l in fd.readlines()])

 #       else: func_dict = False

        lefse = lefse_reader(file_name, self.args['only_sig']\
            , ( False if not self.args['feature_dictionary'] else \
              ( dict( [ tuple( \
                  [ line.rstrip().split('\t')[0], ' '.join(line.rstrip().split('\t')[1:])  ]  ) \
              for line in open(self.args['feature_dictionary']).readlines()   ] ) ) )\
            , True, self.args['read_pwy'])
        lefse.get()
        lefse.select_highest(self.args['number_imp_feat'])
        classes = lefse.classes_

        #
        """
        if self.args['feature_dictionary']:
             with open(self.args['feature_dictionary']) as fd:
                 d = dict([( l.rstrip().split()[0].replace('-', '_')\
                 , ' '.join(l.rstrip().split()[1:])\
                 ) for l in fd.readlines()])
        """
        #for k in d: d[k.replace('-', '_')] = d[k]
        #print d.keys()
        #print [e.name for e in lefse.data]

#        print lefse.feats, '  QUESTI SONO INOSTRI NOMI'
        #lefse.feats = [d[ft] for ft in lefse.feats]
        #print lefse.feats, ' dopo'
        #exit(1)

        #if not self.args['only_sig']:
 
        data_to_plot = pd.DataFrame(\
                                      {'Effect Size': [e.esize for e in lefse.data]\
                                     , 'feat': [e.name for e in lefse.data]\
                                     , 'p': [e.pv for e in lefse.data]\
                                     , 'C': [e.cl for e in lefse.data]\
                                     , self.args['feature_names']: [(name[:-13]+'_spp' \
					if (name.endswith('_unclassified') \
					and (self.args['feature_names']=='Species')) else name).replace('_',' ')\
				        for name in [(e.name if e.pv >= 0.05 else e.name + ' *') for e in lefse.data]]\
                                     , 'colors': [self.color_code[e.cl] for e in lefse.data]\
                                      }\
				   )

        #else:

        #print()
        #print(lefse.feats, ' DEVI GUARDARE QU')
        #print()


        if self.args['split_based_on']:
            data_to_plot['feat'] = [ f.split(self.args['split_based_on'])[-1]  for f in [e.name for e in lefse.data] ]
            lefse.feats = [ f.split(self.args['split_based_on'])[-1]  for f in lefse.feats ]


            #print()
            #print(lefse.feats, ' DEVI GUARDARE QU')
            #print()
 

        #if self.args['only_sig']:         
        #    data_to_plot = data_to_plot.loc[[i for i,v in enumerate([e.pv for e in lefse.data]) if v<self.args['only_sig']]]
            

            #print( [ (e,v)  for e,v in enumerate([e.pv for e in lefse.data])]  ) 
            #print( [i for i,v in enumerate([e.pv for e in lefse.data]) if v<self.args['only_sig']] )
            #print( data_to_plot )


        data_to_plot.set_index('feat', inplace=True)

        #print lefse.feats, ' cazzo ci sei'

        data_to_plot = data_to_plot.loc[lefse.feats]

        print(data_to_plot)


        if bool(self.args['save_features']):
            self.save_feats(self.args['save_features'], lefse.feats)
        if bool(self.args['text_on']):
            self.write_statistical_support(self.args['text_on'], data_to_plot)

        self.plot(data_to_plot, classes)


        
    def plot(self, data, classes):
         sns.set_style(self.args['facecolor'])
         #sns.set_context(rc = {'patch.linewidth': 0.4})

         fig,ax_plot = plt.subplots(1,1, figsize=(6,6))
         ax_plot.set_ylim(-(self.args['bar_height']/2), len(ax_plot.patches)+(self.args['bar_height']/2))
         ax_plot.set_xlim(self.args['x_start'], self.args['x_lim'])

         sns.barplot(data=data, x='Effect Size', y=self.args['feature_names']\
		   , ax=ax_plot, palette=[c for c in data['colors'].tolist()], saturation=0.5, edgecolor='black')

         #ax_plot.xaxis.set_ticks_position('right')
         
         #for bar in ax_plot.patches: 
             
         #    print dir(bar)
         #    bar.outline.set_edgecolor('black')
         #    bar.outline.set_linewidth(5.)
         #    bar.set_height(self.args['bar_height'])

         ax_plot.yaxis.set_ticks_position('right')

         #if self.args['feature_dictionary']:
         #    with open(self.args['feature_dictionary']) as fd:
         #        d = dict([( l.rstrip().split()[0]\
         #        , ' '.join(l.rstrip().split()[1:])\
         #        ) for l in fd.readlines()])

         #fnames = [d[ff] for ff in fnames]


         for lab in ax_plot.get_yticklabels():   ##  for d in dir(lab): print d
             if self.args['feature_names'].lower() == 'species':
                 lab.set_style('italic')

             lab.set_size(12 if self.args['number_imp_feat'] < 10 else 7)
             #if self.args['feature_dictionary']:


             #    print lab.get_name(), ' name'
             #    print lab.get_label(), ' label'
             #    lab.set_name(d[lab.get_name()])
             #    print str(lab)
                 #for dd in  dir(lab): print dd
                 #exit(1)
                 
         sps = [spine.set_linewidth(1.) for spine in ax_plot.spines.values()]
         #sns.despine(left=False, bottom=False)


         plt.subplots_adjust(right=(.3 if not self.args['right_space'] else 0.2), bottom=0.2)
         plt.savefig(self.args['outfile'] + '.' + self.args['fig_fmt'], dpi=400) #self.args['where'] \
             #+ ('' if self.args['where'].endswith('/') else '/') + '_'.join(classes) \
             #+ ('' if not self.args['title'] else '_'+self.args['title']) \
             #+ '_effect_size.' + self.fmt, dpi=400)
        


if __name__ == '__main__':
    a = effect_size_barplot(sys.argv)
    a.set_plot(a.args['lefse_file']) 
