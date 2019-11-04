#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.patches as patches
import itertools
import utils

class cross_analyses_plot(object):
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    save_folder = '../Images/'

    def __init__(self, datasets, title, defined_problem, which_python):
        self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.title = title
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.couples = list(itertools.combinations_with_replacement(self.datasets, 2))

    def named_axes(self):
        fig, ax = plt.subplots(len(self.datasets)+3, len(self.datasets)+1, figsize=(12, 12))
        for i in range(len(self.datasets)+3):
            for j in range(len(self.datasets)+1): 
                ax[i, j].set_xticks([]) 
                ax[i, j].set_yticks([])
        for i,d in enumerate(self.datasets):
            dt = d if not d in self.utils.crc_paper_converter else self.utils.crc_paper_converter[d]
            ax[0, i].set_xlabel(dt, fontsize=15, rotation=45, ha='left', labelpad=90.0)#60. if len(d)<=10 else 90.0)
            ax[i, 0].set_ylabel(dt, fontsize=15, rotation=45, ha='right', labelpad=20.0)
        i += 1
        ax[0, i].set_xlabel('Average', fontsize=15, rotation=45, ha='left', labelpad=20.0)
        [ax[0, y].xaxis.set_label_position('top') for y in range(i+1)]
        lodo_row = i+2
        for j,d in enumerate(self.datasets):
            ax[lodo_row, j].set_xlabel(d if not d in self.utils.crc_paper_converter else self.utils.crc_paper_converter[d], fontsize=15, rotation=45, ha='left', labelpad=80.0)
        ###[ax[lodo_row, j_].xaxis.set_label_position("top") for j_ in range(j+1)]        
        ax[lodo_row, j+1].set_xlabel('Average', fontsize=15, rotation=45, ha='center', labelpad=20.0)
        ####ax[lodo_row, j+1].yaxis.set_label_position('right')
        ax[lodo_row-1, j+1].set_xlabel('CrossValid', fontsize=15, rotation=45, ha='center', labelpad=50.0)
        ax[lodo_row-1, j+1].xaxis.set_label_position('top')
        ax[lodo_row-1, j+1].xaxis.set_label_position('top')
        ax[lodo_row-2, 0].set_ylabel('Average', fontsize=15, rotation=45, ha='right', labelpad=20.0)
        ax[lodo_row, 0].set_ylabel('LODO', fontsize=15, rotation=45, ha='center', labelpad=50.0)
        for i in range(len(self.datasets)): fig.delaxes(ax[lodo_row-1, i])
        fig.delaxes(ax[len(self.datasets), len(self.datasets)])
        ###fig.autofmt_xdate(bottom=0.2, rotation=45, ha='left')
        return fig, ax

    def plot_data(self, db, test, algo, grid, feat, start):
        scores = []
        coors = []
        test_n = self.utils.test_magns(('_'.join(db) if isinstance(db, list) else db), test, algo, grid, feat, start)
        avg_scores = dict([(d, [0,0]) for d in self.datasets])
        for couple in self.couples:
            if self.utils.isonedata(couple):
                pools, cross_validation = [couple], True                
            else:
                pools, cross_validation = [couple, [couple[1], couple[0]]], False
            for pool in pools:
                score, coordinates, n = self.utils.transfer_(pool, db, test, algo, grid, feat, start)
                scores.append(score)
                coors.append(tuple(coordinates))
                if not cross_validation:
                    #print pool[0], ' SU: ', pool[1], ' ==> ', score, ' * ', test_n[pool[1]] , ' = ', score * test_n[pool[1]] , cross_validation, ' EHEHEHEH'
                    avg_scores[pool[0]][0] += score * test_n[pool[1]]           
                    avg_scores[pool[1]][1] += score
        for i,d in enumerate(self.datasets):
            scores.append(avg_scores[d][1]/float(len(self.datasets)-1))   
            coors.append((len(self.datasets), i))
            scores.append(avg_scores[d][0]/float(sum([test_n[ds] for ds in self.datasets if ds != d])))
            coors.append((i, len(self.datasets)))
        mean, div = [], 0
        lodo_row = len(self.datasets) + 2
        for i,d in enumerate(self.datasets):
            score, coordinate, n = self.utils.lodo_(d, db, test, algo, grid, feat, start)
            mean.append(score * n)        
            div += n
            coors.append((lodo_row, i))
            scores.append(score)
        scores.append(sum(mean)/float(div))
        coors.append((lodo_row, i+1))
        cross_validation, n_feat, n = self.utils.cross_validation_(db, test, algo, grid, feat, start)
        scores.append(cross_validation)
        coors.append((lodo_row-1, coordinate+1))
        return scores, coors, n_feat

    def plot(self, db, algo, grid, feat, start):
        test = self.tests[0]
        fig, ax = self.named_axes()
        plt.subplots_adjust(hspace=.03, wspace=.03)
        scores, coors, n_feat = self.plot_data(db, test, algo, grid, feat, start)
        cmap = plt.get_cmap('hot')
        cmap.set_bad('w', 1)
        cax = ax[len(self.datasets), len(self.datasets)].imshow(np.array([[s] for s in scores], dtype=np.float64), cmap=cmap, interpolation='nearest', vmin=0.5, vmax=1.0)
        cbar = fig.colorbar(cax, ax=ax.ravel().tolist(), ticks=[0.5, 0.75, 1.0], extend='min', drawedges=True)
        cbar.ax.set_yticklabels(['0.5', '0.75', '1.0']) 
        if n_feat >= 400:
            n_feat = 'all-features'
        scores_ = [float('%.2f'%s) for s in scores] + [.99]

        for sc,cr in zip(scores, coors):
            ax[cr].patch.set_facecolor(cmap((sc-np.min(scores_))/float(np.max(scores_)-np.min(scores_))))
        for sc,cr in zip(scores, coors):
            ax[cr].text(\
		0.5*(self.left+self.right), 0.5*(self.bottom+self.top), '%.2f'%sc, horizontalalignment='center'\
		, verticalalignment='center', fontsize=30, color='black' if sc>0.6 else 'white') 

		###if ((sc-np.min(scores_))/float(np.max(scores_)-np.min(scores_)))>. else 'white')    

        #plt.suptitle('cross-study-analyses, %s, %s feats' % (('_'.join(db) if isinstance(db, list) else db), str(n_feat)))
        plt.savefig('%s/PAN_cross_study_analyses_%s_%s_%s_%s_%s_%s_features.png' \
		%(self.save_folder, self.title.upper(), '_'.join(db) if isinstance(db, list) else db\
		, grid, 'random forest' if algo=='rf' else 'svm', ('feature per tree' if algo=='rf' else '')+feat, n_feat), dpi=600)
            
if __name__ == '__main__': 
    print ('non ti sopporto piu m\'hai rotto i') 
