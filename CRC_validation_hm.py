#!/usr/bin/env python


import argparse as ap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import matplotlib.gridspec as gridspec
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'


class validation_hm(object):
    
    def read_params(self):

        p = ap.ArgumentParser( )
        add = p.add_argument
        add( 'tables', type=str, nargs=2 )
        add( '-fmt', '--format', type=str, default='png', choices=['png', 'svg'] )
        add( '-tl', '--title', default='Validation-Heatmap', type=str )
        add( '-of', '--out_file', type=str, default='Validation_Heatmap' )
        add( '-t', action='store_true' )
        add( '-cm', '--cmap', type=str, default='hot' )
 
        add( '-sv', '--simple_validation', default=None, type=str )
        
        return vars(p.parse_args())    



    def __init__( self ):

        self.args = self.read_params( )
        self.cmap = self.args['cmap']
        self.fr0 = pd.read_csv( self.args['tables'][0], sep='\t', header=0, index_col=0 )
        self.fr1 = pd.read_csv( self.args['tables'][1], sep='\t', header=0, index_col=0 )

        print self.fr0, ' il primo'
        print self.fr1, ' il secondo'

        """
        if self.args['simple_validation']:
            self.fr2 = pd.read_csv( self.args['simple_validation'], sep='\t', header=0, index_col=0 )
            fig2, ax_ = plt.subplots(figsize=(3,3))
            sns.heatmap(self.fr2 if (not self.args['t']) else self.fr2.T\
                , vmin=0.5, vmax=1.0, ax=ax_\
                , fmt='.2f', cmap=self.cmap, cbar=True\
                , linewidth=0.4, square=True, xticklabels=True\
                , yticklabels=True, annot_kws={'fontsize': 10}\
                , annot=True)
            plt.s
        """


    def draw_hm(self, fr, ax):
    
        return sns.heatmap(fr if (not self.args['t']) else fr.T\
		, vmin=0.5, vmax=1.0, ax=ax\
		, fmt='.2f', cmap=self.cmap, cbar=False\
		, linewidth=0.4, square=True, xticklabels=True\
		, yticklabels=True, annot_kws={'fontsize': 10}\
		, annot=True)
	


    def draw_heatmap(self):

        sns.set( font_scale=.8 ) ##if self.args['num_tests'] > 5 else . )
        fig = plt.figure( figsize=(6,6) )
        gs = gridspec.GridSpec( 2, 3, width_ratios=[3,3,1], height_ratios=[1,1] )

        ax_0 = plt.subplot( gs[0, :-1] )
        ax_1 = plt.subplot( gs[1, :-1] )
        ax_cbar = plt.subplot( gs[:, -1] )

        self.norm_ = matplotlib.colors.Normalize(vmin=.5,vmax=1.)
        self.cbar = matplotlib.colorbar.ColorbarBase( ax_cbar\
                	, cmap=self.cmap, norm=self.norm_\
			, extend='min', filled=True, drawedges=True)
        self.cbar.set_ticks(list(np.arange(0.5,1.0,0.1))+[1])
        self.cbar.set_ticklabels(list(map(str,np.arange(0.5,1.0,.1)))+[1.0])

        hm0 = self.draw_hm(self.fr0, ax_0)
        hm1 = self.draw_hm(self.fr1, ax_1)

        ax_cbar.tick_params(labelsize=8)
        ax_cbar.set_ylabel('AUC', size=8)
        ax_0.set_ylabel('GI-tract Dis.', size=8)
        ax_1.set_ylabel('Non-GI-tract Dis.', size=8)
       
        #ax_0.set_ylabel('CM_Cohort1 LODO', size=8)
        #ax_1.set_ylabel('CM_Cohort1 \nCross-Validation', size=8)

        ax_0.xaxis.set_ticks_position('top')
        #plt.setp( ax_0.get_xticklabels(), rotation=70, ha='right' )
        #plt.setp( ax_1.get_xticklabels(), rotation=-70, ha='right' )

        ##plt.suptitle( 'Model Indifference Power\n (GI & non-GI diseases)' )
        ##plt.suptitle( 'LODO/CM_Cohort1 Indifference Power' )  
        ##plt.subplots_adjust( bottom = 0.2 )

        p0_cbar = ax_cbar.get_position()
        p0_ax0 = ax_0.get_position()
        p0_ax1 = ax_1.get_position()
        p1_cbar = [ p0_cbar.x0 , p0_cbar.y0 + 0.33, p0_cbar.width * 0.5, p0_cbar.height * 0.4 ]
        p1_ax0 = [ p0_ax0.x0 , p0_ax0.y0 , p0_ax0.width , p0_ax0.height ]
        p1_ax1 = [ p0_ax1.x0 , p0_ax1.y0 + 0.2 , p0_ax1.width , p0_ax1.height ] 
        ax_cbar.set_position( p1_cbar )
        ax_0.set_position( p1_ax0 )
        ax_1.set_position( p1_ax1 )

        plt.savefig('.'.join([self.args['out_file'],self.args['format']]), dpi=600)



if __name__ == '__main__':

    vhm = validation_hm()
    vhm.draw_heatmap() 
