#!/usr/bin/env python

from cross_analyses import cross_analyses_plot ## plot the analyses
from cross_study_figures import cross_figures
from cross_study import cross_study ## writes the analyses
from plot_analyses import plot_analyses ## plot support training analyses
from plot_support import training_support
from feat_heat import feature_heatmap ## plot a feature importance heatmap
from utils import usefullfuncs ## utilities
import argparse as ap
import sys

def read_params(args):
    par = ap.ArgumentParser()
    arg = par.add_argument
    arg('title', type=str)
    arg('-d', '--define', type=str)
    arg('-ds', '--datasets', nargs='+')
    arg('-db', '--databases', nargs='+', default=['metaphlan'], choices=['metaphlan', 'genefamilies', 'pathways', 'pfam'])
    arg('-al', '--algo', type=str, choices=['rf', 'svm', 'enet', 'lasso'], default='rf')
    arg('-wp', '--which_python', type=str, default='/scratchCM/users/paolo.manghi/anaconda2/bin')
    arg('-do', type=str, choices=['cross_figures', 'cross_study', 'cross_plot', 'support_study', 'support_plot', 'heat_map', 'batch_plot'])
    arg('-nif', '--number_features', type=int, default=10)
    return vars(par.parse_args())

if __name__ == '__main__':
    """
    to plot the usual analyses:
	python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do cross_figures   
    to plot the progressive training:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot
    to plot a feature importance heatmap:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do heat_map

    """
    par = read_params(sys.argv)

    if par['do'] == 'cross_study':
        cs = cross_study(par['datasets'], par['define'], par['databases'], par['which_python'])
        cs.write_cross_study_commands('rf')

    elif par['do'] == 'cross_figures': #== 'cross_plot':
        cap = cross_figures(par['datasets'], par['title'], par['define'], par['which_python']) 
        l = [[cap.plot(db, par['algo'], 'entropy', 'std', start) for start in [2, 13, 24, 35]] for db in par['databases']]

    elif par['do'] == 'support_study':
        pa = plotanalyses(par['datasets'], par['define'], par['databases'], par['which_python'])
        pa.write_plot_data()
        pa.write_any_support_analyses()

    elif par['do'] == 'support_plot':
        sp = [[training_support(par['title'], par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], np.mean).plot_(start, n_feat) \
			for start,n_feat in zip([2,13,24,35],['all','30','40','50','60'])] for db in par['databases']] 
        sp = [[training_support(par['title'], par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], np.median).plot_(start, n_feat) \
			for start,n_feat in zip([2,13,24,35],['all','30','40','50','60'])] for db in par['databases']]

    elif par['do'] == 'heat_map':
        fhm = [feature_heatmap(par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], par['number_features']) for db in par['databases']] ## save_heatmap()
	####   datasets, db, defined_problem, algo, grid, feat_sel, which_python, n_imp_feat, path='/CM/data/meta/'):

    elif par['do'] == 'batchplot':
        print 
    
