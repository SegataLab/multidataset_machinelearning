#!/usr/bin/env python

import sys
import itertools
import argparse as ap
import pandas as pd
import glob
import subprocess as sp

class md(object):
    def __init__(self, fn):
        self.basename = lambda n : n.split('/')[-1]
        self.data_path = '/scratchCM/data/meta/'
        self.f = pd.read_csv(fn,sep='\t',header=0,index_col=None).fillna('NA').T
        self.__d = dict()
        self.__dataset_list = dict()
        for s,dt in zip(self.f.loc['sampleID'].tolist(),self.f.loc['dataset_name'].tolist()):
            if s not in self.__d: self.__d[s] = [dt]
            else: self.__d[s].append(dt)
            if dt not in self.__dataset_list: self.__dataset_list[dt] = [s]
            else: self.__dataset_list[dt].append(s) 
        self.__datasets = self.__dataset_list.keys()
        self.__duplicates = dict()  
        for c in itertools.combinations(self.__datasets, 2):
            for s in set(self.__dataset_list[c[0]]) & set(self.__dataset_list[c[1]]):
                if s not in self.__duplicates: self.__duplicates[s] = [c]
                else: self.__duplicates[s].append(c)
        self.stats = self.get_statistics()
        #print 'Statistics Loaded.'
        ####for s in self.stats: print s, self.stats[s]
        #self.print_cross_section()
        #self.correct_run_over_sra('hmp_mapped.csv', 'HMP_2012')

    def print_cross_section(self): print self.f.loc[:, self.f.loc['dataset_name'].isin(['ZellerG_2014'])]
    def save(self): self.f.T.to_csv('curatedMetadata.tsv',sep='\t',header=True,index=False)
    def get_duplicates(self): return self.__duplicates
    def get_datasets(self): return self.__dataset_list
    def get_data_of_sample(self): return self.__d

    def print_duplicated_sample_table(self):
        with open('duplicated.csv', 'w') as dup:
            for k in md.get_duplicates(): 
                row = [k] + [dataset for dataset in md.get_duplicates()[k][0]]
                dup.write('\t'.join(row)+'\n')

    def load_codes(self, mapp_file, dataset_name, to_mapp_against_sra=True):
        ## the 'mapper' is the output of python ncbi_downloader.py dataset_name #code# -m xml (with verbose=2, currently default) 
        ## and its utility is 
        a, line = open(mapp_file), True
        if to_mapp_against_sra:
            sra_to_run = dict()
            id_to_run = dict()
            while line:
                line = a.readline().rstrip().split()
                if line:
                    sra_to_run[line[2]] = line[5:-2] 
                    id_to_run[line[3]] = line[5:-2]
            a.close()    
            return sra_to_run, id_to_run

    def insert_stats(self, dataset_, verbose=0):
        self.f.loc['number_reads', self.f.loc['dataset_name'].isin([dataset_])] = [self.stats[s][1] for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['number_bases', self.f.loc['dataset_name'].isin([dataset_])] = [self.stats[s][0] for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['minimum_read_length', self.f.loc['dataset_name'].isin([dataset_])] = [self.stats[s][2] for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['median_read_length', self.f.loc['dataset_name'].isin([dataset_])] = [self.stats[s][3] for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        if verbose == 1: print self.f.loc[['sampleID', 'number_bases', 'number_reads', 'minimum_read_length', 'median_read_length'], self.f.loc['dataset_name'].isin([dataset_])]
        elif verbose == 2: self.f.loc[['sampleID', 'number_bases', 'number_reads', 'minimum_read_length', 'median_read_length'], self.f.loc['dataset_name'].isin([dataset_])].T.to_csv(\
			dataset_.lower()+'.tmp', sep='\t', header=True, index=False)       

    def correct_run_over_sra(self, mapp_, dataset_):
        sra_to_run = self.load_codes(mapp_, dataset_, True)[0]
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [';'.join(sra_to_run[s]) for s in self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_ncbi_accession(self, mapp_, dataset_):
        id_to_run = self.load_codes(mapp_, dataset_)[1]
        #print len(id_to_run.keys()), len()
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [';'.join(id_to_run[s]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def get_statistics(self):
        stats = dict()
        for s in glob.glob('all_read_stats/*.stats'):
            with open(s) as st:
                stats[self.basename(s)[:-6]] = list(map(lambda x : int(float(x)), [st.readline() for i in range(2)][-1].rstrip().split()[1:]))
        #for s,dt in zip(self.f.loc['sampleID'].tolist(),self.f.loc['dataset_name'].tolist()):
        #    if s not in stats: stats[s] = self.ext_stats(dt,'scripts/all_read_stats/'+s+'.stats')
        return stats

    def ext_stats(self, dataset_name, sample_name):
        sample_ = self.data_path+dataset_name+'/reads/'+sample_name+'.fastq.bz2'
        sp.call(['python','/scratchCM/repos/pyphlan/fna_len.py',sample_,'scripts/all_read_stats/'+sample_name+'.stats','-q','--stat'])
        with open(self.data_path+'scripts/all_read_stats/'+sample_name+'.stats') as stats:
            line = [stats.readline() for i in range(2)][-1].rstrip().split()
        return line

    def give_report(self):
        print '# samples= %i' %len(list(self.d.keys())) 
        print '# datasets= %i' %len(list(set(self.d.values())))
        print '# body-sites= %i' %len(list(set(self.f.loc['body_site'].tolist()))) 
        for c in self.f.loc['sampleID'].tolist():
            print c, 'number of times: ' , self.f.loc['sampleID'].tolist().count(c)
        return ''

if __name__ == '__main__':
    md = md('cMD_metadata.txt')
    #md.print_duplicated_sample_table()

    #md.correct_run_over_sra('hmp_mapped.csv', 'HMP_2012')
    #md.correct_run_over_sra('karlsson_mapping.txt', 'KarlssonFH_2013') 
    #md.correct_run_over_sra('yu_mapped.csv','YuJ_2015')   

    #md.add_ncbi_accession('tett_mapped.csv','TettAJ_2016')

    #md.insert_stats('LomanNJ_2013')
    md.insert_stats('ZellerG_2014', 2)
    #md.
    #md.
    #md.

    

