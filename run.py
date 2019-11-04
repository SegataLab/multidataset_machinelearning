#!/usr/bin/env python

from cross_analyses import cross_analyses_plot ## plot the analyses
from cross_study_figures import cross_figures
from cross_study import cross_study ## writes the analyses
from plot_analyses import plot_analyses ## plot support training analyses
from plot_support import training_support
from feat_heat import feature_heatmap ## plot a feature importance heatmap
from feat_curve import feat_curve
from batch_study_figures import triangular_batch_plot
from feat_heat_vector import feature_heatmap_vector
from utils import usefullfuncs ## utilities
import argparse as ap
import sys


def read_params(args):

    par = ap.ArgumentParser()
    arg = par.add_argument
    arg('title', type=str)
    arg('--define', type=str)
    arg('-ds', '--datasets', nargs='+')
    arg('-db', '--databases', nargs='+', default=['metaphlan'], choices=['metaphlan','genefamilies','pathways','pfam','markers'])
    arg('-al', '--algo', type=str, choices=['rf','lsvm','svm','enet','lasso'], default='rf')
    arg('-wp', '--which_python', type=str, default='/scratchCM/users/paolo.manghi/anaconda2/bin')
    arg('-do', type=str, default=None, choices=['doall','cross_figures','cross_study','cross_plot','support_study','support_plot','heat_map','heat_map_vector','feat_curve','cross_batch','cross_batch_figure'])
    arg('-nif', '--number_features', type=int, default=10, help='This is for the heatmap.')
    arg('-cs', '--charsize', type=int, default=5, help='Fontsize for feature heatmap and progressive training.')
    arg('-cm','--cmap',type=str,default='RdYlBu')
    arg('-w','--width',type=float,default=0.0)
    arg('-ph', '--path', default=None, type=str)
    ### following arguments are for handling a good grid search
    ### they are three and can be used together, in combo
    ### they are needed both for writing analyses and get figures
    arg('-r','--runs', type=int, default=1)
    arg('-mf','--max_features', default=None, choices=[None, '0.1', '0.01', '0.3', '0.5', '0.6', 'auto', 'log2']) 
    ###### if don set this it will use grid between 0.3 and auto
    arg('-ncores','--num_cores', default=10, type=int)
    arg('-o','--only',type=str,default=None,choices=['lodo','transfer'])
    grid_search_choices = ['nt:100', 'nt:500', 'nt:1000', 'nt:10000', 'nt:5000', 'nsl:1', 'oob', 'nsl:5', 'nsl:10', 'a:rf', 'a:svm', 'df', 'c:entropy', 'c:gini', 'null', 'special']
    arg('-g0','--grid0', type=str, default='nt:500', choices=grid_search_choices)
    arg('-g1','--grid1', type=str, default='nsl:1', choices=grid_search_choices)
    arg('-g2','--grid2', type=str, default='c:entropy', choices=grid_search_choices)
    arg('-p','--parallel', action='store_true')
    arg('-hv','--how_verbose', type=int, default=0, choices=[0,1,2])
    arg('-bg','--in_background', action='store_true')

    ##arg('-fr', '--force_reading', nargs='+', default=[], type=str)
    arg('-ft', '--figtitle', type=str, default=None)
    arg('-fmt', '--figfmt', type=str, choices=['svg','png'], default='png')
    arg('-mx','--mixed_taxa', action='store_true')
    pars = vars(par.parse_args())
    if pars['algo'] == 'lsvm': pars['grid0'], pars['grid1'], pars['grid2'] = 'null', 'null', 'null'
    return pars


if __name__ == '__main__':

    """
    USAGES:
    to plot the usual analyses:
        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do cross_figures -g0 nt:500 -g1 nsl:5 -g2 c:entropy
	##-g0 c:gini -g1 nsl:1 -g2 df 

        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db markers -al rf -do cross_figures -g0 nt:1000 -g1 nsl:5 -g2 df
        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db pathways -al rf -do cross_figures -g0 nt:1000 -g1 nsl:5 -g2 c:entropy
 
        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do cross_figures --only lodo 



      #       (for adenomas)
      #       python run.py adenomas_vs_CRC --define study_condition:CRC:adenoma -ds ZellerG_2014 FengQ_2015 CM_lilt HanniganGD_2017 -db metaphlan -al rf -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -do cross_figures -ft ML_HM_metaphlan_adenoma_vs_CRC -fmt png
      #       python run.py adenomas_cs_control --define study_condition:adenoma:control -ds ZellerG_2014 FengQ_2015 CM_lilt HanniganGD_2017 -db metaphlan -al rf -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -do cross_figures -ft ML_HM_metaphlan_adenoma_vs_control -fmt png


    to write the lines for the progressive plot:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -al rf -do support_study -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -mf 0.3 -hv 0

    to plot the progressive learning:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot -db metaphlan -g0 nt:1000 -g1 nsl:5 -g2 c:entropy

    to plot a feature importance heatmap:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do heat_map -db metaphlan -g0 nt:1000 -g1 nsl:5 -g2 c:entropy 

	(with a scores vectors):
            python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do heat_map_vector -db metaphlan -g0 nt:1000 -g1 nsl:5 -g2 c:entropy 

    to get the feat curve:
        python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do feat_curve -g0 nt:1000 -g1 nsl:5 -g2 c:entropy 
      
    """
  
# this part is more tricky veacsue you could be confounded

#    to write down the standard analyses:

#       (gene families)
#	python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -al rf -do cross_study -g0 nsl:5 -g1 nt:1000 -g2 df -hv 0 -ncores 2 -mf 0.01 -r 1

#       (markers)
#       python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db markers -al rf -do cross_study -g0 nsl:5 -g1 nt:1000 -g2 df -hv 0 -ncores 2 -mf 0.01 -r 1

#       (standard)
#       python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -al rf -do cross_study -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -ncores 2 -mf 0.3 -r 1

#       (for adenomas)
#       python run.py adenomas_vs_CRC --define study_condition:CRC:adenoma -ds ZellerG_2014 FengQ_2015 CM_lilt HanniganGD_2017 -db metaphlan -al rf -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -ncores 2 -mf 0.3 -r 1 -do cross_study #-ft ML_HM_metaphlan_adenoma_vs_CRC
#       python run.py adenomas_cs_control --define study_condition:adenoma:control -ds ZellerG_2014 FengQ_2015 CM_lilt HanniganGD_2017 -db metaphlan -al rf -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -ncores 2 -mf 0.3 -r 1 -do cross_study #-ft ML_HM_metaphlan_adenoma_vs_control 

    """
    to write down the batch analyses:
        python run.py crc_batch --define study_condition:control:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -al rf -do cross_batch -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -mf 0.3 -hv 0 -ncores 1

    to plot the batch trinagle (the only one currently available):
        python run.py crc -ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 -db metaphlan -al rf -g0 c:entropy -g1 nt:1000 -g2 nsl:5 -hv 0 -ncores 1 --figtitle "Batch on CRC" -do cross_batch_figure

    """
    par = read_params(sys.argv)

    if par['do'] == 'cross_study':
        cs = cross_study(par['datasets'],par['define'],par['databases'],par['which_python'],par['runs'],par['max_features']\
            ,par['num_cores'],par['grid0'],par['grid1'],par['grid2'],par['only'],par['how_verbose'],par['in_background'],par['mixed_taxa'])
        cs.write_cross_command_lines()
        cs.write_cross_study_commands(par['algo'])

    elif par['do'] == 'cross_figures':
        cap = cross_figures(par['datasets'], par['title'], par['define'], par['which_python'], par['cmap'], par['width'], par['figtitle'], par['figfmt'], par['path']) 
        [[cap.plot(db, par['algo'], par['grid0'], par['grid1'], par['grid2'], 2) for db in par['databases']]]  ## start) for start in [2,13,24,35,46]] for db in par['databases']]

    elif par['do'] == 'support_study':
        pa = plot_analyses(par['datasets'], par['define'], par['databases'], par['which_python'], par['grid0'], par['grid1'], par['grid2'], par['max_features'], par['num_cores'], par['runs'], par['how_verbose'])
        pa.write_plot_data()
        pa.write_any_support_analyses()

    elif par['do'] == 'support_plot':
        [training_support(par['title'], par['datasets'], db, par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['which_python'], par['figfmt'], 'median'\
		, par['path']).all_singletons() for db in par['databases']] 

    elif par['do'] == 'feat_curve': ## title, datasets, db, defined_problem, algo, grid0, grid1, grid2, which_python
        fc = feat_curve(par['title'], par['datasets'], par['databases'], par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['figfmt'], par['path'])  
        fc.do_the_plot()

    elif par['do'] == 'heat_map':

        [feature_heatmap(par['datasets'], db, par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['which_python'], par['number_features']\
               , False, par['charsize'], par['figfmt'], par['path']) for db in par['databases']]

        #[feature_heatmap(par['datasets'], db, par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['which_python'], num_features\
	#	, i_have_to_use_the_lodo_question, par['charsize'], par['cmap'], par['figfmt']) for db in par['databases'] \
	#	  for i_have_to_use_the_lodo_question, num_features in zip([False, True],[par['number_features'], par['number_features']*3])]

    elif par['do'] == 'heat_map_vector':
        [feature_heatmap_vector(par['datasets'], db, par['define'], par['algo'], par['grid0'], par['grid1'], par['grid2'], par['which_python'], par['number_features']\
                , False, par['charsize'], par['cmap'], par['figfmt']) for db in par['databases']]

        ##fhml = [feature_heatmap(par['datasets'], db, par['define'], par['algo'], 'entropy', 'std', par['which_python'], par['number_features'], True, par['charsize']) for db in par['databases']]
        ####   datasets, db, defined_problem, algo, grid, feat_sel, which_python, n_imp_feat, path='/CM/data/meta/'):

    elif par['do'] == 'cross_batch':
        cs = cross_study(par['datasets'],par['define'],par['databases'],par['which_python'],par['runs'],par['max_features']\
	    ,par['num_cores'],par['grid0'],par['grid1'],par['grid2'],par['only'],par['how_verbose'],par['in_background'],par['mixed_taxa'])
        cs.write_batchdata_command_lines()
        cs.write_batch_study_commands(par['algo'])
 
    elif par['do'] == 'cross_batch_figure':
        bp = triangular_batch_plot(par['datasets'], par['title'], 'the class is still missing, sorry', par['grid0'], par['grid1'], par['grid2'], par['cmap'], par['width'], par['figtitle'], par['figfmt'], par['path'])
        bp.triangle_heatmap(par['databases'], par['algo'])

    ####**************************
