#!/usr/bin/env python

import sys
import itertools
import argparse as ap
import pandas as pd
import glob
import os
import subprocess as sp

class md(object):
    def __init__(self, fn, dataname_to_add=None):
        self.basename = lambda n : n.split('/')[-1]
        self.data_path = '/scratchCM/data/meta/'
        self.f = pd.read_csv(fn,sep='\t',header=0,index_col=None).fillna('NA').T
        self.found = ['dataset_name','sampleID','subjectID','body_site','body_subsite','antibiotics_current_use','antibiotics_family','study_condition','disease','disease_subtype'\
		     ,'age','infant_age','age_category','gender','visit_number','BMI','country','location','non_westernized','days_from_first_collection','sequencing_platform'\
		     ,'DNA_extraction_kit','PMID','number_reads','number_bases','minimum_read_length','median_read_length'] ## ,'curator','NCBI_accession']
        self.__d = dict()
        self.__dataset_list = dict()
        if 'dataset_name' in self.f.index.tolist():
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
            if not 'curator' in set(self.f.index.tolist()): self.f.loc['curator',:] = 'Paolo_Manghi' 
        else:
            if dataname_to_add: 
                self.insert_name(dataname_to_add)
                self.stats = self.get_statistics(dataname_to_add)
                self.insert_stats(dataname_to_add)
            for field in self.f.index.tolist():
                self.f.loc[field] = [(STR.strip() if isinstance(STR, str) else STR) for STR in self.f.loc[field].tolist()]
            if 'DNA_extraction_kit' in self.f.index.tolist():
                self.f.loc['DNA_extraction_kit'] = [(x if not (x.lower()=='tiangen') else 'Tiangen') for x in self.f.loc['DNA_extraction_kit'].tolist()]
            else:
                self.f.loc['DNA_extraction_kit'] = 'NA'
            self.sams = self.f.loc['sampleID']
        if 'age' in self.f.index.tolist(): self.f.loc['age'] = [(int(a) if str(a)[0].isdigit() else a) for a in self.f.loc['age'].tolist()]
 
    def print_cross_section(self, data): print self.f.loc[:, self.f.loc['dataset_name'].isin([data])]
    def save(self): self.f.T.to_csv('curatedMetadata.tsv',sep='\t',header=True,index=False)
    def get_duplicates(self): return self.__duplicates
    def get_datasets(self): return self.__dataset_list
    def get_data_of_sample(self): return self.__d
    def data_list(self): return self.f.loc['dataset_name', :].tolist()
    def add_curator_first_time(self, curator, dataset_name): self.f.loc['curator', self.loc['dataset_name'].isin([dataset_name])] = curator
    def equal_cols(self, cola, colb): self.f.loc[colb] = self.f.loc[cola].tolist() 
    def sort_vector_column(self,column): self.f.loc[column] = [';'.join(sorted(values.split(';'))) for values in self.f.loc[column].tolist()]

    def capitalize_column(self, dataname, column): self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])] = [c.upper() for c in self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])]]
    def lower_column(self, dataname, column): self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])] = [c.lower() for c in self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])]]
    def cast_column(self, dataname, column, datatype): self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])] = list(map(datatype, self.f.loc[column, self.f.loc['dataset_name'].isin([dataname])].tolist()))


    def append_something_to_all(self, dataname, column, what, after=True):
         self.datum = self.f.loc[:, self.f.loc['dataset_name'].isin([dataname])].copy(deep=True).T         
         self.datum[column] = [(what+'_'+r if not after else r+'_'+what) for r in self.datum[column].tolist()]
         return self.datum.T
        
    def sort_fields(self, fields):
        found = [c for c in self.found if c in fields]
        non_found = [c for c in fields if not c in found + ['curator','NCBI_accession']]
        return found + non_found + ['curator','NCBI_accession']

    def save_curated_table(self, dataset_name):
        if not os.path.isdir('./%s' %dataset_name): os.mkdir(dataset_name)         
        self.f.loc['disease'] = [(d if not (d is 'none') else 'healthy') for d in self.f.loc['disease'].tolist()]
        if 'disease_subtype' in self.f.index.tolist():
            self.f.loc['disease_subtype'] = [(ds if (ds.strip() != 'healthy') else 'NA') for ds in self.f.loc['disease_subtype'].tolist()]
        fc = self.f.loc[:, self.f.loc['dataset_name'].isin([dataset_name])].copy(deep=True).T
        fc_columns = fc.columns.tolist()
        mask = fc.apply(lambda row : row.tolist().count('NA')==len(row.tolist()), axis=0)
        fc = fc.loc[:, ~mask]
        fc_columns = [c for c in fc_columns if ((c in fc.columns.tolist()) & (not c in ['NCBI_accession','dataset_name']))] + ['NCBI_accession']
        if not 'NCBI_accession' in fc.columns.tolist(): fc['NCBI_accession'] = 'NA'
        if not 'DNA_extraction_kit' in fc.columns.tolist(): fc['DNA_extraction_kit'] = 'NA'
        fc = fc[self.sort_fields(fc_columns)]
        if not 'antibiotics_current_use' in fc_columns:
            if 'body_subsite' in fc_columns: fc.insert(4,'antibiotics_current_use','NA')
            else: fc.insert(3,'antibiotics_current_use','NA')
        fc.to_csv(dataset_name+'/'+dataset_name+'_metadata.tsv', sep='\t', header=True, index=False)

    def model_column_on_another(self, dataname, column_to_change, column_to_look_at, value_not_to_copy='NA'):
        f = lambda cz, co, passing_by : [(cz[i] if co[i]==passing_by else co[i]) for i in range(len(cz))]
        self.datum = self.f.loc[:, self.f.loc['dataset_name'].isin([dataname])].copy(deep=True).T
        self.datum[column_to_change] = f(self.datum[column_to_change].tolist(), self.datum[column_to_look_at].tolist(), value_not_to_copy)
        return self.datum.T

    def data_set(self, prettyprint=False): 
        if not prettyprint:
            l= list(set(self.f.loc['dataset_name', :].tolist()))
            return l, len(l)
        else:
            return ' '.join(list(set(self.f.loc['dataset_name', :].tolist()))) 

    def print_duplicated_sample_table(self):
        with open('duplicated.csv', 'w') as dup:
            for k in md.get_duplicates(): 
                row = [k] + [dataset for dataset in md.get_duplicates()[k][0]]
                dup.write('\t'.join(row)+'\n')

    def load_codes(self, mapp_file, dataset_name, to_mapp_against_sra=True):
        ## the 'mapper' is the output of python ncbi_downloader.py dataset_name #code# -m xml (with verbose=2, currently default) 
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
        else:
            id_to_run = dict()
            while line:
                line = a.readline().rstrip().split()
                if line:
                    id_to_run[line[2]] = line[4:-2]
            a.close()
            return id_to_run

    def drop_name(self): self.f.drop('dataset_name',axis=0,inplace=True)
    def insert_name(self, dataset_): self.f.loc['dataset_name'] = dataset_
    def add_curator(self, curator_): self.f.loc['curator'] = curator_

    def insert_stats(self, dataset_, verbose=0):
        self.f.loc['number_reads', self.f.loc['dataset_name'].isin([dataset_])] = [int(self.stats[s][1]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['number_bases', self.f.loc['dataset_name'].isin([dataset_])] = [int(self.stats[s][0]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['minimum_read_length', self.f.loc['dataset_name'].isin([dataset_])] = [int(self.stats[s][2]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        self.f.loc['median_read_length', self.f.loc['dataset_name'].isin([dataset_])] = [int(self.stats[s][3]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]
        if verbose == 1: print self.f.loc[['sampleID', 'number_bases', 'number_reads', 'minimum_read_length', 'median_read_length'], self.f.loc['dataset_name'].isin([dataset_])]
        elif verbose == 2: self.f.loc[['sampleID', 'number_bases', 'number_reads', 'minimum_read_length', 'median_read_length'], self.f.loc['dataset_name'].isin([dataset_])].T.to_csv(\
			dataset_.lower()+'.tmp', sep='\t', header=True, index=False)       

    def correct_run_over_sra(self, mapp_, dataset_):
        sra_to_run = self.load_codes(mapp_, dataset_, True)[0]
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [';'.join(sra_to_run[s]) for s in self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_ncbi_accession(self, mapp_, dataset_):
        id_to_run = self.load_codes(mapp_, dataset_)[1]
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [';'.join(id_to_run[s]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_ncbi_accession_no_sra(self, mapp_, dataset_):
        id_to_run = self.load_codes(mapp_, dataset_, False)
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [';'.join(id_to_run[s]) for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_ncbi_accession_with_excep(self, mapp_, dataset_):
        id_to_run = self.load_codes(mapp_, dataset_)[1]			
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [(';'.join(id_to_run[s]) if s in id_to_run else 'na') for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_accession_easyway(self, mapp_, dataset_):					
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [dict([tuple(line.rstrip().split()) for line in open(mapp_)])[s] \
			for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]

    def add_accession_easyway_with_excep(self, mapp_, dataset_):
        a = open(mapp_)
        available_map = dict([tuple(line.rstrip().split()) for line in a.readlines()])
        a.close()
        self.f.loc['NCBI_accession', self.f.loc['dataset_name'].isin([dataset_])] = [ ('NA' if s not in available_map else available_map[s]) \
		for s in self.f.loc['sampleID', self.f.loc['dataset_name'].isin([dataset_])].tolist()]


    def get_statistics(self, name=None):
        stats = dict()
        if not name: 
            for s in glob.glob('all_read_stats/*.stats'):
                with open(s) as st:
                    stats[self.basename(s)[:-6]] = list(map(lambda x : int(float(x)), [st.readline() for i in range(2)][-1].rstrip().split()[1:]))
            #for s,dt in zip(self.f.loc['sampleID'].tolist(),self.f.loc['dataset_name'].tolist()):
            #    if s not in stats: stats[s] = self.ext_stats(dt,'scripts/all_read_stats/'+s+'.stats')
            return stats
        else:
            for s in self.f.loc['sampleID'].tolist():
                with open('all_read_stats/'+s+'.stats') as st:
                    stats[s] = list(map(lambda x : int(float(x)), [st.readline() for i in range(2)][-1].rstrip().split()[1:]))
            return stats
	
    def correct_column(self, dataname, title, prev, new): 
        self.f.loc[title, self.f.loc['dataset_name'].isin([dataname])] = [(n if n!=prev else new) for n in self.f.loc[title, self.f.loc['dataset_name'].isin([dataname])].tolist()]

    def correct_column_conditional(self, dataname, tp, tc, nw, cp):
        new_col = [(nw if p==cp else cur) for p,cur in zip(self.f.loc[tp, self.f.loc['dataset_name'].isin([dataname])].tolist(), self.f.loc[tc, self.f.loc['dataset_name'].isin([dataname])].tolist())]
        self.f.loc[tc, self.f.loc['dataset_name'].isin([dataname])] = new_col

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
        for c in self.f.loc['sampleID'].tolist(): print c, 'number of times: ' , self.f.loc['sampleID'].tolist().count(c)
        return ''

    def get_report_strain_evolution(self):			
        self.fc = self.f.copy(deep=True).T
        adults = self.fc[(self.fc['days_from_first_collection']!='NA') & (self.fc['age_category']!='child') & (self.fc['age_category']!='newborn') & (self.fc['age_category']!='schoolage')]
        childs = self.fc[(self.fc['days_from_first_collection']!='NA') & ((self.fc['age_category']=='newborn') | (self.fc['age_category']=='child'))]
        print childs.shape, ' childs'
        print adults.shape, ' adults'
        sample_per_tp_child = dict([( sample, len(childs[subj]==childs[childs[sample] & childs[subj]])) for sample,subj in zip(childs['sampleID'].tolist(), childs['subjectID'].tolist())])
        sample_per_tp_adult = dict([( sample, len(adults[subj]==adults[adults[sample] & adults[subj]])) for sample,subj in zip(adults['sampleID'].tolist(), adults['subjectID'].tolist())])
        print sample_per_tp_child
        print sample_per_tp_adult

if __name__ == '__main__':
  ## EXAMPLES: curation of the dataset BackhedF_2015
  ##	       starting from a metadata table 'initial_backhed.csv'

  md = md('initial_backhed.csv','BackhedF_2015')
  md.insert_stats('BackhedF_2015')
  md.add_ncbi_accession_no_sra('backhedf_2015_mapped.csv','BackhedF_2015')
  #md.f = md.append_something_to_all('BackhedF_2015','sampleID',)
  #md.f = md.append_something_to_all('BackhedF_2015','subjectID',)

  md.lower_column('BackhedF_2015','disease')
  md.lower_column('BackhedF_2015','born_method')
  md.cast_to_int('BackhedF_2015', 'infant_age')
  md.cast_to_int('BackhedF_2015', 'days_from_first_collection')
  md.add_curator('Valentina_Giunchiglia')##, 'BackhedF_2015')
  md.save_curated_table('BackhedF_2015')

