#!/usr/bin/env python


import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
import matplotlib.gridspec as gridspec
from matplotlib import rc
matplotlib.rcParams['svg.fonttype'] = 'none'


class the_most_simple_heatmap(object):

    def read_params( self ):

        p = argparse.ArgumentParser()
        add = p.add_argument
        add ( 'filename' )
        add ( 'outfile' )
        add ( '-fmt', '--format', choices=['png', 'svg'], default='png' )
        add ( '-t', '--transpose', action='store_true' )
        add ( '-cm', '--cmap', default='hot', type=str )
        add ( '--vmin', type=float, default=0.5 )
        add ( '--vmax', type=float, default=1.0 )
        return vars(p.parse_args())


    def __init__( self ):

        self.args = self.read_params()
        self.cmap = self.args['cmap']
        self.frame = pd.read_csv( self.args['filename'], sep='\t', header=0, index_col=0 )
        fig, ax = plt.subplots(figsize=(3,3))
        sns.heatmap(self.frame if (not self.args['transpose']) else self.frame.T\
            , vmin=self.args['vmin'], vmax=self.args['vmax'], ax=ax, fmt='.2f', cmap=self.cmap, cbar=True\
            , linewidth=0.4, square=True, xticklabels=True, yticklabels=True\
            , annot_kws={'fontsize': 4}, annot=True)
        plt.subplots_adjust(left=0.40)
        plt.savefig(self.args['outfile']+'.'+self.args['format'], dpi=600)


if __name__ == '__main__':
    
    tmsh = the_most_simple_heatmap()
