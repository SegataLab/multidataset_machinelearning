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
    arg('-db', '--databases', nargs='+', default=['metaphlan'], choices=['metaphlan','genefamilies','pathways','pfam'])
    arg('-al', '--algo', type=str, choices=['rf','svm','enet','lasso'], default='rf')
    arg('-wp', '--which_python', type=str\
	, default='/scratchCM/users/paolo.manghi/anaconda2/bin')
    arg('-do', type=str, default='doall', choices=['doall','cross_figures','cross_study','cross_plot','support_study','support_plot','heat_map','batch_plot'])
    arg('-nif', '--number_features', type=int, default=10, help='This is for the heatmap; the feats selected for ml always 10-60.')
    arg('-cs', '--charsize', type=int, default=5, help='Fontsize for feature heatmap and progressive training.')
    ### following arguments are for handling a good grid search
    ### they are three and can be used together, in combo
    ### they are needed both for writing analyses and get figures
    arg('-r','--runs', type=int, default=20)
    arg('-mf','--max_features', default=None, choices=[None, '0.3', 'auto', 'log2']) ## if don set this it will use grid between 0.3 and auto
    arg('-ncores','--num_cores', default=10, type=int)
    arg('-o','--only',type=str,default=None,choices=['lodo','transfer'])
    grid_search_choices = ['nt:100', 'nt:500', 'nt:1000', 'nt:10000', 'nt:5000', 'nsl:1', 'oob', 'nsl:5', 'nsl:10', 'a:rf', 'a:svm', 'df', 'c:entropy', 'c:gini']
    arg('-g0','--grid0', type=str, default='nt:500', choices=grid_search_choices)
    arg('-g1','--grid1', type=str, default='nsl:1', choices=grid_search_choices)
    arg('-g2','--grid2', type=str, default='c:entropy', choices=grid_search_choices)
    arg('-p','--parallel', action='store_true')
    arg('-hv','--how_verbose', type=int, default=0, choices=[0,1,2])
    arg('-bg','--in_background', action='store_true')
    return vars(par.parse_args())

if __name__ == '__main__':
    """
    USAGES:
    to plot the usual analyses:
        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do cross_figures -g0 nt:500 -g1 nsl:5 -g2 c:entropy
	##-g0 c:gini -g1 nsl:1 -g2 df 
    to write the lines for the progressive plot:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -al rf -do support_study -g0 nt:500 -g1 nsl:5 -g2 c:entropy
    to plot the progressive training:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot
    to plot a feature importance heatmap:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do heat_map
    to write down the standard analyses:
	python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -al rf -do cross_study -g0 c:gini -g1 nsl:1 -g2 df -hv 0 -ncores 1 -mf auto -r 1
    """
    par = read_params(sys.argv)
    if par['do'] == 'cross_study':
        cs = cross_study(par['datasets'],par['define'],par['databases'],par['which_python'],par['runs'],par['max_features']\
        ,par['num_cores'],par['grid0'],par['grid1'],par['grid2'],par['only'],par['how_verbose'],par['in_background'])
        cs.write_cross_command_lines()
        cs.write_cross_study_commands(par['algo'])
    elif par['do'] == 'cross_figures':
        cap = cross_figures(par['datasets'], par['title'], par['define'], par['which_python']) 
        [[cap.plot(db, par['algo'], par['grid0'], par['grid1'], par['grid2'], 2) for db in par['databases']]]  ## start) for start in [2,13,24,35,46]] for db in par['databases']]
    elif par['do'] == 'support_study':
        pa = plot_analyses(par['datasets'], par['define'], par['databases'], par['which_python'], par['grid0'], par['grid1'], par['grid2'])
        pa.write_plot_data()
        pa.write_any_support_analyses()
    elif par['do'] == 'support_plot':
        [training_support(par['title'], par['datasets'], db, par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['which_python'], 'median').all_singletons() for db in par['databases']] 
    elif par['do'] == 'heat_map':
        #fhm = [feature_heatmap(par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], par['number_features'], False, par['charsize']) for db in par['databases']]
        fhml = [feature_heatmap(par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], par['number_features'], True, par['charsize']) for db in par['databases']]
        ####   datasets, db, defined_problem, algo, grid, feat_sel, which_python, n_imp_feat, path='/CM/data/meta/'):
    elif par['do'] == 'batchplot':
        print 
    #####***************************************
