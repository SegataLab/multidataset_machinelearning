#!/usr/bin/env python

import pandas as pd
import argparse, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3\
                          , venn2_circles\
                          , venn3_circles

class euler(object):
    def read_input(self):

        data, names, colors = [], [], []

        with open(self.args['data'], 'r') as i:

            line = True
            lines = []
            while line != '\n':
                line = i.readline()
                if line != '\n': lines += [line]

            #print lines, '  QUESTO EEEEEEE'

            for line in lines:
                line = line.rstrip().split()
                #print line, ' OOOOOOOOOOOOO'
                names.append(str(line[0]))
                colors.append(str(line[1]))
                data.append(list(map( (int if line[2][0].isdigit() else str), line[2:])))

        return data, names, colors



    def __init__(self):

        self.args = self.read_params(sys.argv)
        #self.args = self.read_params(sys.argv)
        self.func = venn3 if not self.args['two_sets'] else venn2

        self.data, self.names, self.colors = self.read_input()
        self.N = len(self.data)



    def read_params(self, args):
        par = argparse.ArgumentParser()
        add = par.add_argument
        
        add('data', type=str\
           , help='a dataframe with variables as index and sets as columns.')
        add('--two_sets', action='store_true')
        add('-tl', '--title', default='Venn_Diagram')
        add('-fmt', '--format', type=str, default='png')

        return vars(par.parse_args())


        
    def plot_with_3(self):

        #c = venn3_circles(subsets = [len(self.groups[0])], color='orangered')
        

        v = venn3(subsets = (6,6,2, 0,1,0,1) , set_labels = self.names ) # (2, 4, 5) )  #(len(self.groups[0]), len(self.groups[1]), len(self.groups[2
        c = venn3_circles(subsets = (6,6,2, 0,3,0,1) ) ## , set_labels = self.names)

#        c.get_patch_by_id('100').set_alpha(.6)
#        v.get_patch_by_id('100').set_alpha(.6)
 
        c[0].set_color('orangered')
        c[0].set_alpha(0.6)
        c[1].set_color('deepskyblue')
        c[1].set_alpha(0.6)
        c[2].set_color('crimson')
        c[2].set_alpha(0.6) 

 #       for id_,color in zip(self.names, ['orangered', 'deepskyblue', 'crimson']):
 #           v.get_patch_by_id(id_).set_alpha(.8)
 #           v.get_patch_by_id(id_).set_color(color)

        #plt.suptitle('Shotgun-Metagenomics is more powerful than enything else, you see:')
        plt.savefig('%s.%s' %(self.args['title'], self.args['format']), dpi=400)


if __name__ == '__main__':
    
    e = euler()
    #for g in e.groups:
    #    print g
    e.plot_with_3()

        
