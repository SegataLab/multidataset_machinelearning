#!/usr/bin/env python


import sys, argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'
from sklearn.metrics import roc_auc_score, roc_curve
from scipy import stats

class risk_index(object):

    save_folder = '' #'../ml/'

    def __init__(self, args):

        self.args = self.read_params(args)
        self.data = pd.read_csv(self.args['if'], sep='\t', header=None, index_col=0)
        self.pos, self.negs = self.read_teachers(self.args['teach'])
        self.define_teachers(self.pos, self.negs, self.args['cl'])

        if self.args['auc']:
            self.fpr, self.tpr, self.auc = self.get_scores([self.data.loc[self.args['cl']].tolist()], [list(map(float, self.data.loc[self.args['sc']].tolist()))])[0]
            self.plot_auc(self.fpr, self.tpr, self.auc)

        if bool(self.args['plot']):
            self.fpr, self.tpr, self.auc = self.get_scores([self.data.loc[self.args['cl']].tolist()], [list(map(float, self.data.loc[self.args['sc']].tolist()))])[0]
            self.catch_info(self.args['plot'])
            self.plot_plot(self.auc)
            

    def read_params(self, args):

        par = argparse.ArgumentParser()
        add = par.add_argument
        add('-if', type=str, default=sys.stdin)
        add('-of', type=str, default=sys.stdout)
        add('-cl', type=str, default='define')
        add('-sc', type=str, default='Bug-Complex-Abundance')
        add('-teach', type=str, default='teachers.txt')
        add('-fmt', type=str, default='png', choices=['svg','png'])
        add('-tl', '--title', type=str, default='Generic_Risk_Index')

        add('-cm', '--cmap', type=str, default='cool')
        add('--dot_size', type=float, default=16)

        add('--x_title', type=str)
        add('--y_title', type=str)
        add('--auc', action='store_true')
        add('--plot', type=str, default=None)
        return vars(par.parse_args())


    def catch_info(self, filename):
        with open(filename + '.txt') as i:
            self.info = dict([(line.rstrip().split()[0],float(line.rstrip().split()[1])) for line in i.readlines()])
        self.data.loc[filename] = [self.info[s] for s in self.data.loc['sampleID'].tolist()]


    def plot_plot(self, auc=False):
        sns.set_style('white')
        fig, ax = plt.subplots(1,1,figsize=(6,4))
        cmap = vars(matplotlib.cm)[self.args['cmap']]
        scatterp, single_scatters = [], []

        ax.set_ylim((-.05,1.05))
        ax.set_yticks(np.arange(0,1.05,0.1))

        self.grads = dict([(s,float(g)) for s,g in zip(self.data.loc['sampleID'].tolist(), self.data.loc[self.args['sc']].tolist())] )
        plottable = pd.DataFrame(self.data.loc[['sampleID', self.args['plot'], self.args['sc']]]).T
        plottable.columns=['sampleID',self.args['plot'],self.args['sc']]
        plottable.index = plottable['sampleID'].tolist()        
        del plottable['sampleID']
        plottable = plottable.astype('float')

        #vals = plottable[self.args['sc']].tolist()  #  (grad - np.min(grad)) / (np.max() - np.min(plottable[self.args['sc']].tolist()))
        #plottable[self.args['sc']] = (vals - np.min(vals)) / (np.max(vals) - np.min(vals))

        for sam in plottable.index.tolist():
            dd_ = pd.DataFrame(data=plottable.loc[sam].values.reshape([1,2]), index=[sam], columns=['x1','x2'])
            sns.regplot(x='x1', y='x2', ax=ax, scatter=True, fit_reg=False, data=dd_, marker='o', color=cmap(self.grads[sam]))
   
        sns.regplot(x=self.args['plot'], y=self.args['sc'], ax=ax, scatter=False, data=plottable, color='dodgerblue', fit_reg=True)
        ga = [i for i,e in zip(self.data.loc[self.args['plot']].tolist(), self.data.loc[self.args['cl']].tolist()) if e==1]
        gb = [i for i,e in zip(self.data.loc[self.args['plot']].tolist(), self.data.loc[self.args['cl']].tolist()) if e==0]
        gac = list(map(float, [i for i,e in zip(self.data.loc[self.args['sc']].tolist(), self.data.loc[self.args['cl']].tolist()) if e==1]))
        gbc = list(map(float, [i for i,e in zip(self.data.loc[self.args['sc']].tolist(), self.data.loc[self.args['cl']].tolist()) if e==0]))

        p_value_cat = stats.ttest_ind(ga, gb, axis=0, equal_var=False)[1]
        p_value_con = stats.ttest_ind(gac, gbc, axis=0, equal_var=False)[1]

        if self.args['x_title']: ax.set_xlabel(self.args['x_title'], fontsize=9)
        if self.args['y_title']: ax.set_ylabel(self.args['y_title'], fontsize=9)         
        x_annot = 2
        y_annot = 0.94


        if auc:
           ax.annotate('%s AUC = %.2f' %('Predictor' if not self.args['y_title'] else self.args['y_title'].upper(), auc), (x_annot,y_annot), size=7)
           y_annot -= 0.04
        if p_value_cat < 0.05:
           ax.annotate('%s P-value = %s' %('Index' if not self.args['x_title'] else self.args['x_title'].upper(), p_value_cat), (x_annot,y_annot), size=7)
           y_annot -= 0.04
        if p_value_con < 0.05:
           ax.annotate('%s P-value = %s' %('Index' if not self.args['y_title'] else self.args['y_title'].upper(), p_value_con), (x_annot,y_annot), size=7)
           y_annot -= 0.04

        corr, p = stats.pearsonr( plottable[self.args['plot']].tolist(), plottable[self.args['sc']].tolist())
        ax.annotate('Corr.Coeff = %.3f, (p. %.5f)' %(corr,p), (x_annot,y_annot), size=7)

           #if p_value_cat < 0.05: 
           #    ax.annotate('%s AUC = %.2f\n%s P-Value = %s' %('Predictor' if not self.args['y_title'] else self.args['y_title'].upper(), auc\
                           ##,'Cat. Var' if not self.args['x_title'] else self.args['x_title'].upper(), str(p_value_cat)), (2, 0.94), size=9)
           #else:
           #    ax.annotate('%s AUC = %.2f' %('Predictor' if not self.args['y_title'] else self.args['y_title'].upper(), auc), (2, 0.94), size=9)

        sps = [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
        plt.suptitle('A lower PPD is partially correlated with a lower RIM', fontsize=10) 
        plt.savefig('%s.%s' %(self.args['title'], self.args['fmt']), dpi=400)


    def read_teachers(self, filename):
        with open(filename) as fn:
            line0, line1 = tuple([line.rstrip().split() for line in fn.readlines()])
        return (line0[1:], line1[1:]) if line0[0].startswith('p') else (line1[1:],line0[1:]) 


    def define_teachers(self, plus_sample_list, minus_sample_list, define):
        teach = dict([(s,t) for s,t in zip(plus_sample_list + minus_sample_list\
                   , [1 for i in range(len(plus_sample_list))] + [0 for i in range(len(minus_sample_list))])]) 
        self.data = self.data.loc[:, self.data.loc['sampleID'].isin(plus_sample_list+minus_sample_list)]
        self.data.loc[define] = [teach[c] for c in self.data.loc['sampleID'].tolist()]
       

    def get_score(self, teachers, scores):
        scores = (scores - np.min(scores)) / (np.max(scores) - np.min(scores))

        fpr, tpr, thresholds = roc_curve(teachers, scores, pos_label=None, drop_intermediate=True)
        auc = roc_auc_score(teachers, scores, average='macro', sample_weight=None)
        return fpr, tpr, auc


    def get_scores(self, teachers, scores): return tuple([self.get_score(teacher,score) for teacher,score in zip(teachers,scores)])

    def plot_auc(self, fpr, tpr, auc):
        self.data_in_frame = pd.DataFrame(dict([(ax, points) for ax, points in zip(['x','y'], [fpr,tpr])]))
        sns.set_style('white')
        fig, ax_ = plt.subplots(figsize=(5,4))
        ax_.plot(fpr, tpr, alpha=0.6, color='gold', linewidth=1., linestyle='--')
        sns.regplot(ax=ax_, x='x', y='y', data=self.data_in_frame , marker='o', fit_reg=False, scatter_kws={'s': 60})
        ax_.xaxis.set_ticks(list(np.arange(0.0, 1.0, 0.1)) + [1.0])
        ax_.yaxis.set_ticks(list(np.arange(0.0, 1.0, 0.1)) + [1.0])
        ax_.set_xticklabels(list(map(str, list(np.arange(0.0, 1.0, 0.1)) + [1.0])), fontsize=4)
        ax_.set_yticklabels(list(map(str, list(np.arange(0.0, 1.0, 0.1)) + [1.0])), fontsize=4)
        ax_.set_yticklabels(list(map(str, np.arange(0.5, 1., 0.1))), fontsize=4)
        ax_.set_ylabel('Sensitivity (True Positive Rate)', fontsize=4)
        ax_.set_xlabel('1 - Specificity (False Positive Rate)', fontsize=4)
        sps = [i.set_linewidth(0.5) for i in ax_.spines.itervalues()]
        ax_.annotate('AUC = %.2f' %auc, (0.65, 0.25), size=4)        
        plt.suptitle('%s' % self.args['title'].replace('_',' '), fontsize=6)
        plt.savefig('%s%s.%s' %(self.save_folder, self.args['title'], self.args['fmt']), dpi=400)
        


if __name__ == '__main__':

    ri = risk_index(sys.argv)
    ##ri.plot_auc()

