#!/usr/bin/env python


import numpy as np
import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import matplotlib.gridspec as gridspec
import itertools
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
import utils, sys



class double_heatmap_periimplantitis_paper(object):
    
    def read_params(self, args):

        par = argparse.ArgumentParser()
        arg = par.add_argument

        arg('files', type=str, nargs='+')

        arg('-c', '--classes', type=str, nargs='+', required=False)
        arg('-cc', '--contralat_classes', type=str, nargs=2\
            , default=['All_Contralateral', 'Own_Contralateral'])

        arg('-of', '--output', type=str)
        arg('-fmt', '--format', type=str, default='svg')
        arg('-tl', '--title', type=str)

        return vars(par.parse_args())


    def start_from_arguments(self):

        self.files = self.args['files']
        self.classes = self.args['classes']


    def start_from_file(self, file):

        with open(file) as f:
            ops = [l.rstrip().split() for l in f.readlines()]

        self.files = ops[0]
        self.classes = ops[1]


    def __init__(self):

        self.utils = utils.usefullfuncs(['datasets'], '')
        self.args = self.read_params(sys.argv)

        if len(self.args['classes'])==len(self.args['files'])//4: self.start_from_arguments()

        elif (len(self.args['files'])==1) and (len(self.args['classes']) == 0):
            self.start_from_file(self.args['files'][0])

        else: raise IndexError('You have to speficy either a single filename for the params '
              'OR a number of files and then -c [len(files) n of class.]')

        self.n_classes = len(self.classes)
        self.cmap_hot = 'hot'
        self.cmap_cold = 'hot'
        self.width = 0.4
        self.pack_data()



    def pack_data(self):
    
        self.hot_auc = pd.DataFrame(\
            np.array([self.utils._auc(f, 2, True)[0] for f in self.files[:(len(self.files)//2)]]\
            , dtype=np.float64).reshape([2, self.n_classes])\
            , columns=self.classes, index=self.args['contralat_classes'])

        self.cold_auc = pd.DataFrame(\
            np.array([self.utils._auc(f, 2, True)[0] for f in self.files[(len(self.files)//2):]]\
            , dtype=np.float64).reshape([2, self.n_classes])\
            , columns=self.classes, index=self.args['contralat_classes'])

        #print self.hot_auc
        #print
        #print self.cold_auc



    def figure(self):

        sns.set(font_scale=.8)
        fig = plt.figure(figsize=(5,5))
        gs = gridspec.GridSpec(self.n_classes-1\
          , self.n_classes-1, width_ratios=[12,1], height_ratios=[1,1])

        self.ax_hot = plt.subplot(gs[0,0])
        self.ax_cold = plt.subplot(gs[1,0])

        self.cax_hot = plt.subplot(gs[0,1])
        self.cax_cold = plt.subplot(gs[1,1])

        self.norm = matplotlib.colors.Normalize(vmin=.5,vmax=1.)

        self.cbar_hot = matplotlib.colorbar.ColorbarBase(ax=self.cax_hot, cmap=self.cmap_hot\
            , norm=self.norm, extend='min', filled=True, drawedges=False)
        self.cbar_cold = matplotlib.colorbar.ColorbarBase(ax=self.cax_cold, cmap=self.cmap_cold\
            , norm=self.norm, extend='min', filled=True, drawedges=False)

        cold = sns.heatmap(self.cold_auc, ax=self.ax_cold, mask=None, square=True\
            , cbar=False, vmin=0.5, vmax=1.0, cmap=self.cmap_cold, annot=True, xticklabels=True, yticklabels=True)
        hot = sns.heatmap(self.hot_auc, ax=self.ax_hot, mask=None, square=True\
            , cbar=False, vmin=0.5, vmax=1.0, cmap=self.cmap_hot, annot=True, xticklabels=True, yticklabels=True)    

        self.cax_hot.set_ylabel('AUC') #, size=16)
        self.cax_cold.set_ylabel('AUC') #, size=16)
        self.ax_hot.xaxis.set_ticks_position('top')
        self.ax_cold.xaxis.set_ticks_position('top')

        plt.subplots_adjust(left=0.225) ##, top=0.80) #, bottom=0.12)             
        plt.suptitle(self.args['title'])
        plt.savefig('%s.%s' %(self.args['output'], self.args['format']), dpi=400)




class dual_heatmap(object):


    def read_params(self, args):

        par = argparse.ArgumentParser()
        arg = par.add_argument
        arg('files', type=str, nargs='+')
        arg('-c', '--classes', type=str, nargs='+', required=False)
        arg('-of', '--output', type=str)
        arg('-fmt', '--format', type=str, default='png')
        arg('-tl', '--title', type=str)

        return vars(par.parse_args())


    def start_from_arguments(self):

        self.files = self.args['files']
        self.classes = self.args['classes']


    def start_from_file(self, file):

        with open(file) as f:
            ops = [l.rstrip().split() for l in f.readlines()]
        self.files = ops[0]
        self.classes = ops[1] 



    def __init__(self):

        self.utils = utils.usefullfuncs(['datasets'], '')
        self.args = self.read_params(sys.argv)

        if len(self.args['classes'])==len(self.args['files'])//2: self.start_from_arguments()

        elif (len(self.args['files'])==1) and (len(self.args['classes']) == 0):
            self.start_from_file(self.args['files'][0])       

        else: raise IndexError('You have to speficy either a single filename for the params '
              'OR a number of files and then -c [len(files) n of class.]')

        self.n_classes = len(self.classes)
        self.cmap_hot = 'hot'
        self.cmap_cold = 'PuBu_r'
        self.width = 0.4       
        self.pack_data()



    def pack_data(self):

        self.hot_auc = [self.utils._auc(f, 2, True)[0] for f in self.files[:(len(self.files)//2)]]
        self.cold_auc = [self.utils._auc(f, 2, True)[0] for f in self.files[(len(self.files)//2):]]

        self.grid = np.zeros((self.n_classes, self.n_classes))
        self.hot_mask = np.zeros(self.grid.shape, dtype=bool)

        self.cold_mask = np.zeros(self.grid.shape, dtype=bool)
        self.hot_grid = np.zeros((self.n_classes, self.n_classes), dtype=np.float64)

        self.cold_grid = np.zeros((self.n_classes, self.n_classes), dtype=np.float64)

        self.little_debug = dict([(f,s) for f,s in zip(self.files, [self.utils._auc(f, 2, True)[0] for f in self.files])])

        for k in self.little_debug: print (k, self.little_debug[k])

        for i in range(self.grid.shape[0]):
            for j in range(self.grid.shape[0]):
                if i>j:
                    self.cold_mask[i,j] = True
                    self.hot_grid[i,j] += self.hot_auc[(i+j)-1]
                elif i<j:
                    self.hot_mask[i,j] = True
                    self.cold_grid[i,j] += self.cold_auc[(j+i)-1]



    def define_proportions(self):

        pos1_ax_hot_cold = self.ax_hot_cold.get_position()
        pos1_cax_hot = self.cax_hot.get_position()
        pos1_cax_cold = self.cax_cold.get_position()
        pos2_ax = [pos1_ax_hot_cold.x0, pos1_ax_hot_cold.y0, pos1_ax_hot_cold.width * 0.7, pos1_ax_hot_cold.height]
        pos2_cax_h = [pos1_cax_hot.x0 + 0.05, pos1_cax_hot.y0 + 0.2, pos1_cax_hot.width * 0.7, pos1_cax_hot.height * 0.60]
        pos2_cax_c = [pos1_cax_cold.x0 - 0.2, pos1_cax_cold.y0 - 0.05, pos1_cax_cold.width * 0.7, pos1_cax_cold.height * 0.60]
        self.ax_hot_cold.set_position(pos2_ax)
        self.cax_hot.set_position(pos2_cax_h)
        self.cax_cold.set_position(pos2_cax_c)



    def central_array(self):

        shape_ = [['' for i in range(self.n_classes)] for i in range(self.n_classes)]
        for i in range(self.n_classes): shape_[i][i] = self.classes[i]
        return shape_
  


    def figure(self):

        sns.set(font_scale=.8)
        fig = plt.figure(figsize=(5,5))
        gs = gridspec.GridSpec(self.n_classes\
          , self.n_classes, width_ratios=[1,9,1], height_ratios=[1,1,1])  

        self.ax_hot_cold = plt.subplot(gs[:,1:-1])
        self.cax_hot = plt.subplot(gs[1:,0])
        self.cax_cold = plt.subplot(gs[:-1,-1])

        annot_mask = [[True for t in range(self.n_classes)] \
            for t in range(self.n_classes)]

        for t in range(self.n_classes):
            for tt in range(self.n_classes):
                if t==tt: annot_mask[t][tt] = False

        annot_mask = np.array(annot_mask, dtype=bool) 
        central_array = np.array(self.central_array())

        self.norm = matplotlib.colors.Normalize(vmin=.5,vmax=1.)

        self.cbar_hot = matplotlib.colorbar.ColorbarBase(ax=self.cax_hot, cmap=self.cmap_hot\
            , norm=self.norm, extend='min', filled=True, drawedges=False)

        self.cbar_cold = matplotlib.colorbar.ColorbarBase(ax=self.cax_cold, cmap=self.cmap_cold\
            , norm=self.norm, extend='min', filled=True, drawedges=False)

        cold = sns.heatmap(self.cold_grid, ax=self.ax_hot_cold, mask=self.cold_mask, square=True\
            , cbar=False, vmin=0.5, vmax=1.0, cmap=self.cmap_cold, annot=True, xticklabels=False, yticklabels=False)
        hot = sns.heatmap(self.hot_grid, ax=self.ax_hot_cold, mask=self.hot_mask, square=True\
            , cbar=False, vmin=0.5, vmax=1.0, cmap=self.cmap_hot, annot=True, xticklabels=False, yticklabels=False)

        ticks = sns.heatmap(np.ones([self.n_classes, self.n_classes])\
            , ax=self.ax_hot_cold, mask=annot_mask, square=True, annot=central_array\
            , cbar=False, vmin=0.5, vmax=1.0, cmap=self.cmap_hot, fmt=''\
            , xticklabels=False, yticklabels=False)

        self.cax_hot.yaxis.set_ticks_position('left')
        self.cax_hot.yaxis.set_label_position('left')
        self.define_proportions()
        plt.suptitle(self.args['title'])
       
        plt.savefig('%s.%s' %(self.args['output'], self.args['format']), dpi=400)




if __name__ == '__main__':
    

    #print 'Doing Trials'
    s = dual_heatmap()
    s.figure()
    #s.pack_data()

    #p = double_heatmap_periimplantitis_paper()
    #p.figure()

    #*********
