#!/usr/bin/env python


import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import utils
import itertools
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
import os


class score_in_triangle(object):

    def __init__(self, sc, cr, headers):
        self.auc = sc
        self.headers = headers if cr[0]>cr[1] else (headers[1],headers[0])
        self.where = cr if cr[0]>cr[1] else (cr[1],cr[0])


class triangular_batch_plot(object):

    ##save_fig = '../Images'
    
    def __init__(self, datasets, title, defined_class, grid_term0, grid_term1, grid_term2, cmap, width, fig_title, fig_fmt, path):

        self.save_folder = (path + ('/' if (not path.endswith('/')) else '') + 'Images/') if path else '../Images/'
        self.utils = utils.usefullfuncs(datasets, path)
        self.datasets = datasets
        self.title = title
        self.defined_class = defined_class
        self.cmap = cmap
        self.width = width
        self.fig_title = fig_title
        self.fig_fmt = fig_fmt
        self.couples = list(itertools.combinations_with_replacement(self.datasets, 2))
        self.grid_term0 = grid_term0
        self.grid_term1 = grid_term1
        self.grid_term2 = grid_term2
        self.headers = dict()
        self.norm = matplotlib.colors.Normalize(0.8, 1.0)
        
        
    def rescue_results(self, db, algo):
        out = list()
        for couple in self.couples + [[c[1], c[0]] for c in self.couples]:
            #print self.utils.batch_(couple, db, algo, self.grid_term0, self.grid_term1, self.grid_term2, 2)[0],
            if not self.utils.isonedata(couple):
                if os.path.exists(self.utils.get_onlyone_batch_resultname(couple, db, algo, 'control', self.grid_term0, self.grid_term1, self.grid_term2)):
                    out.append(  score_in_triangle( self.utils.batch_(couple, db, algo, 'control', self.grid_term0, self.grid_term1, self.grid_term2, 2)[0] \
                              , (self.datasets.index(couple[0]), self.datasets.index(couple[1]) ) \
			      , (couple[0],couple[1])))
                    for coor,name in zip(out[-1].where, out[-1].headers): self.headers[coor] = name
        return out        


    def triangle_heatmap(self, db, algo):
        
        fig = plt.figure(figsize=(8,6)) 
        gs = gridspec.GridSpec(1,2, width_ratios=[15,1])
        self.ax_hm = plt.subplot(gs[0, 0])
        self.ax_cbar = plt.subplot(gs[0, 1])

        scores = self.rescue_results(db, algo) 
        coors, dims = [s.where for s in scores], len(self.datasets) ## np.max(np.max([s.where[0] for s in scores]), np.max([s.where[1] for s in scores])) +1 
        ticks = [(t[1] if t[1] not in self.utils.data_aliases else self.utils.data_aliases[t[1]]) for t in sorted([(c,n) for c,n in self.headers.items()], key=lambda x : x[0])]
        matrix = np.zeros((dims, dims), dtype=np.float64)
        mask = np.ones((dims, dims), dtype=bool)

        for sc in scores: 
            matrix[sc.where] += sc.auc
            mask[sc.where] = False 
        matrix = matrix
        #matrix.columns = list((ds if ds not in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets)
        #matrix.index = list((ds if ds not in self.utils.data_aliases else self.utils.data_aliases[ds]) for ds in self.datasets)
        sns.heatmap(data=pd.DataFrame(matrix), annot=True, fmt='.2f', ax=self.ax_hm, mask=mask, xticklabels=ticks[:-1], yticklabels=['']+ticks[1:], cbar=False, cmap='hot', vmin=0.8, vmax=1.0)
        cbar_hot = matplotlib.colorbar.ColorbarBase(self.ax_cbar, cmap='hot', norm=self.norm, extend='min')

        plt.setp(self.ax_hm.get_xticklabels(), rotation=38, ha='right')
        plt.suptitle(self.fig_title.replace('_',' '))
        plt.subplots_adjust(left=0.2, bottom=0.24)
        plt.savefig('%sBatchHeatmap_%s.%s' %(self.save_folder, 'control' if not self.fig_title else self.fig_title, self.fig_fmt), dpi=600)        




if __name__=='__main__':

    print 'yabadabadu'
                    
