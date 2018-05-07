#!/usr/bin/env python

import argparse as ap
import subprocess as sp
import scipy.stats as sts
import sys,os
import numpy as np
import pandas as pd

"""

- This tool is a command-line provider of
- meta-datasets; definition-level is higher 
- wrt CMD R-package. As in the 
- curatedMetagenomicDataset, starting bases 
- are metadata and functional and taxonomical 
- profiles of 9 types from metaphlan2 and humann2.
- The script assumes that for any dataset
- a folder is present containing metadata 
- and profiles in subfolders.

- 1st step: Export main dataset folder in bash_profile.

- NOTE: the tool also contain the converter class 
- to pass from any human profile to a profile
- shaped like a metaphlan profile.
- This is necessary for the script merge_metaphlan_profiles
- to be used and to keep things sorted.

- NOTE: for very large databases (genefamilies, markers)
- there's an option to perform Wilcoxon - rank 
- in order to srhink features.
"""


class metadataset:

    def __init__(self):

        self.args = self.read_params(sys.argv)
        self.imp_fields = \
            ['dataset_name'\
            ,'sampleID','subjectID','body_site'\
            ,'study_condition','disease','age'\
            ,'age_category','gender','country'\
            ,'sequencing_platform'\
            ,'DNA_extraction_kit'\
            ,'PMID']

        self.numerator = 0

        self.profile_species = set()
        self.profile_markab = set()
        self.profile_markpres = set()
        self.profile_cov = set()
        self.profile_pwys = set()
        self.profile_genes = set()
        self.profile_pfam = set()
        self.profile_mirna = set()
        self.metadata = ''
        self.features = list()

    def read_params(self, args):
        p = ap.ArgumentParser()
        add = p.add_argument
        # - DATA TO TAKE
        add( 'datasets', nargs='+', help='E.g.: LiJ_2014.study_condition:control,body_site:stool .'\
           'Another possibility is setting: LiJ_2014.study_condition:control,col:country:FRA'\
           ',col:body_subsite:subgingival_plaque to add a new column to the loaded metadata table.')

	# - FOR SPECIFY DIFFERENT PATHS		
        add( '--base_path', type=str, default='/CM/data/meta/')
        add( '--metadata_path', type=str, default='/CM/data/meta/')
        add( '--metadata_name', type=str, default=None)

        add( '--metaphlan_path', type=str, default='/CM/data/meta/')
        #add( '--pwyrelab_path', type=str, default='')
        #add( '--genefam_path', type=str, default='')
        #add( '--marker_path', type=str, default='')

        add('--change_profile_func', action='store_true') ##handle_profile)

        # - ACTIVATE TYPES OF PROFILES
        add( '--metaphlan', action='store_true')		
        add( '--pwyrelab', action='store_true')
        add( '--genefam', action='store_true')
        add( '--markab', action='store_true')
        add( '--markpres', action='store_true')
        add( '--pwycov', action='store_true')
        add( '--pfam', action='store_true')
        add( '--mirna', action='store_true')

        # - NAMES FOR THE PROFILES FOLDERS
        add( '--metaphlan_folder', type=str, default='metaphlan2/')
        add( '--pwyrelab_folder', type=str, default='pathways/')
        add( '--genefam_folder', type=str, default='gene_families/')
        add( '--pwycoverage_folder', type=str, default='humann2_PWYcov/')
        add( '--markab_folder', type=str, default='marker_abundance/')
        add( '--markpres_folder', type=str, default='marker_presence/')
        add( '--pfam_folder', type=str, default='pfam_profiles/')
        add( '--mirna_folder', type=str, default='miRNA/')

        # - WORD AFTER '_' TO DESIGN THE PROFILE FILE
        add( '--metaphlan_title', type=str, default='_profile')
        add( '--pwy_title', type=str, default='_profile')
        add( '--cov_title', type=str, default='_profile') #, choices=['_profile','_complete_profile'])
        add( '--genefam_title', type=str, default='_profile') #, choices=['_profile','_complete_profile'])
        add( '--pfam_title', type=str, default='_profile')
        add( '--mirna_title', type=str, default='_profile')

        # - SAMPLE A RANDOMIZE POPULATION FOMRM A BIGGER DATA LIST
        add( '-rc', '--randomize_controls', action='store_true')
        add( '-pr', '--proportions', type=str, default='1,3', help='1,2 (1:moiety of cnts wrt cases, 2:moiety of training cnts wrt to test cnts.)')
        add( '-cf', '--control_filename', type=str, default='controls.csv')
        add( '-mulo', '--random_controls_multiple_of', type=int, default=None\
          , help='If not specified, uses as bases for the mutiplicity the allready collected samples.')
        add( '-cc', '--control_conditions', default=None, help='e.g: study_condition:IBD:control, gender:female,study_condition:control')
        add( '-cn', '--control_names', default=['training_controls','testing_controls'], type=str, nargs=2\
          , help='the names to the train and test of the randomized controls;must be two but can be equal: the first is the training.')
        add( '-cw', '--control_word', type=str, default=None, help='add a common field to all the externally derived controls\
             , the added field will be in third position (after dataset_name sampleID) under header "ext_characterization"')
        
        add( '-sc', '--select_columns', type=str, default=[], nargs='+', help='"Selected Columns"')
        add( '-scff', '--select_columns_from_file', type=str, default=None, help='"Select Columns from a File"')
        # - FURTHER METAPHLAN SPECIFICS
        add( '-mx', '--mixed_taxa', action='store_true', help='"With Metaphlan, Uses All the Taxa Together"')
        add( '-ys', '--yes_strain', action='store_true', help='"With Metaphlan, Leave the Strain Level There"')
        add( '-tx', '--taxon', type=str, choices=['k__','p__','c__','o__','f__','g__','s__'], default='s__')
        add( '-sr', '--shrink', action='store_true', help='"With Metaphlan, Uses as Names only The Currenty Selected Taxon"')

        # - only metadata is lasofor prepare the bigger group to randomize after. 	
        add( '-x', '--exclude_samples', nargs='+', type=str, default=[]) ## options to avoid as much as possible 
        add( '-om', '--only_metadata', action='store_true', help='This Option is Usefull Also for Control To Be Ramdomised After.')
        add( '-mgp', '--merge_profile_exec', default='/scratchCM/users/paolo.manghi/multidataset_machinelearning/merge_metaphlan_tables02.py')

        ##add('--top_features', default=None, type=int)
        # - TRANSFORMATIONS
        add( '-pn', '--percentile_normalization', type=str, default=None\
            ,help='percentile normalization: eg -pn study_condtion:control:CRC whatever is in the 1st field is the column to look at. Whats is int eh 2nd \
                  is the reference percentiles distribution (should be some kind of control) whatever is in the last field will be uniformed to the reference: \
                  make usage of add_column option of this script to help yourself.')
        add( '-lg', '--log_transform', action='store_true', help='performs (log(1+anyfeature))')	
        add( '-dm', '--dense_matrix', action='store_true'\
            , help='cut out features being zeros in <zero_threshold> percentage of the samples [ default = 10 per ]')
        add( '-zt', '--zero_threshold', nargs='+', default=[0.1], type=float)

        # - ON FEATURES
        add( '-mf', '--meta_features', action='store_true')
        add( '-fs', '--feature_selection', default=None, type=str \
            , help='"Differs From -scff in That This Only Has Control on Features, -scff Also On Metadata"')

        # - WILKOXON RANK FOR TOO LARGE FEATURE SETS
        add( '-wk', '--wilcoxon', type=str, default=None, help='--wilcoxon study_condition:CRC:control')
        add( '--bonferroni', action='store_true')

        add('--mapped_number', action='store_true')
        add('--mapped_percentage', action='store_true')


        # - OUTPUT SPECIFICS
        add( '-t', '--transpose', action='store_true', help='Output Has Field Names On Columns')
        add( '-b', '--both', action='store_true', help='Outputs Two Results, One Has Fields On Indexes The Other On Columns.')
        add( '-of', '--output_file', type=str, default=None, help='"If Don\'t Specify the -of, There You\'ll Be Able to Sedn the Result As A STDOUT."')
        add( '-gs', '--give_statistics', type=str, default=None\
          , help='Statistics on metadata -gs study_condition:sample_type:disease, or -gs random/randomization/rand/randrandrand/randanything')

        add( '-fo', '--feat_only', action='store_true')
        add( '-fac', '--feat_and_condition', default=[], nargs='+', help='"Alternative to -scff, Completely On Command Line"')

        add( '--grad', type=str, default=[], nargs='+')
        add( '--grad_col', type=str, default='Bug-Complex-Abundance')
        add( '--log_gradient', action='store_true')


        pp = vars(p.parse_args())
        if (not pp['metaphlan']) and (not pp['pwyrelab']) and (not pp['genefam']) and (not pp['mirna']) and\
           (not pp['markab']) and (not pp['markpres']) and (not pp['pwycov']) and (not pp['only_metadata']) and (not pp['pfam']):
            pp['metaphlan'] = True

        pp['proportions'] = list(map(int, pp['proportions'].split(',')))
        if len(pp['proportions'])!=2:
            raise IOError('proportions must contain 2 numbers separated by a comma.')
        num_sel, queries, cols_to_add, datasets = len(pp['datasets']), list(), [[] for j in range(len(pp['datasets']))], list()

        for i,data in enumerate(pp['datasets']): 

            col_added, dataset, params = False, data.split('.')[0], data.split('.')[1].split(',')
            tot_reqs, cols, catch = len(params), list(), list()

            for q in params:

                if (not q.startswith('col:')): catch.append(q)

                else: cols_to_add[i] += [q[4:].split(':')]  

            #col_added = True
            #if not col_added: cols_to_add.append([])
            queries.append(','.join(catch))
            datasets.append(dataset)

        pp['select'], pp['columns'], pp['dataset_input'] = queries, cols_to_add, datasets	

        if (all([x.split(':')[1:]==['control'] for x in queries]) and (not pp['only_metadata'])):

            if pp['output_file']: print 'WARNING: <only_metadata> option is not specififed, even if all the queries are \'controls\'.'
            if pp['output_file']: print 'Ensure this dataset will not be mixed with other via the <randomize_controls> option,',

        if pp['transpose'] and pp['both']: pp['traspose'],pp['both'] = False,True

        return pp



    #def get_percentage_of_reads_mapped(self):
    #    for s in glob.glob('metaphlan2/*/*norm*.tsv'):
    #        with open(s, 'r') as tor:
    #            nn.append(float(tor.readlines()[-1].split()[-1]))
	


    def handle_profile(self, sample, dataset_name):

        if self.args['metaphlan'] and (not sample in self.args['exclude_samples']): 
            self.profile_species.add(self.args['base_path']+dataset_name+'/'+self.args['metaphlan_folder']+sample+'/'+sample+self.args['metaphlan_title']+'.tsv')

        if self.args['pfam'] and (not sample in self.args['exclude_samples']):
            self.profile_pfam.add(self.args['base_path']+dataset_name+'/'+self.args['pfam_folder']+sample+'/'+sample+self.args['pfam_title']+'.txt')

        if self.args['markab'] and (not sample in self.args['exclude_samples']):
            self.profile_markab.add(self.args['base_path']+dataset_name+'/'+self.args['markab_folder']+sample+'/'+sample+'_profile'+'.txt')

        if self.args['markpres'] and (not sample in self.args['exclude_samples']):
            self.profile_markpres.add(self.args['base_path']+dataset_name+'/'+self.args['markpres_folder']+sample+'/'+sample+'_profile'+'.txt')

        if self.args['pwycov'] and (not sample in self.args['exclude_samples']):
            self.profile_cov.add(self.args['base_path']+dataset_name+'/'+self.args['pwycoverage_folder']+sample+'/'+sample+self.args['cov_title']+'.txt')

        if self.args['pwyrelab'] and (not sample in self.args['exclude_samples']): 
            self.profile_pwys.add(self.args['base_path']+dataset_name+'/'+self.args['pwyrelab_folder']+sample+'/'+sample+self.args['pwy_title']+'.tsv')

        if self.args['genefam'] and (not sample in self.args['exclude_samples']): 
            self.profile_genes.add(self.args['base_path']+dataset_name+'/'+self.args['genefam_folder']+sample+'/'+ sample+self.args['genefam_title']+'.tsv')

        if self.args['mirna'] and (not sample in self.args['exclude_samples']):
            self.profile_mirna.add(self.args['base_path']+dataset_name+'/'+self.args['mirna_folder']+sample+'/'+sample+self.args['mirna_title']+'.tsv')


        #if self.args['reads_mapped'] and (not sample in self.args['exclude_samples']):
        #    self.reads_mapped.append(self.args['base_path']+dataset_name+'/'+self.args['metaphlan_folder']+sample+'/'+sample+'_norm_reads.tsv')

            #self.numerator += 1 
            #print os.path.exists(self.args['base_path']+dataset_name+'/'+self.args['mirna_folder']+sample+'/'+sample+self.args['mirna_title']+'.tsv'), 'that ', self.numerator, '  exist'



    def handle_profile_direct_folder(self, sample, dataset_name):

        if self.args['metaphlan'] and (not sample in self.args['exclude_samples']):
            self.profile_species.add( self.args['base_path'] + dataset_name + '/' + self.args['metaphlan_folder'] + sample + self.args['metaphlan_title'] + '.tsv')

        if self.args['pwyrelab'] and (not sample in self.args['exclude_samples']):
            self.profile_pwys.add( self.args['base_path'] + dataset_name + '/' + self.args['pwyrelab_folder'] + sample + self.args['pwy_title'] + '.tsv')



    def create_dataset(self):


        def datasets(path,data_name,selection_criteria,log,random_control=False):
            pf = pd.DataFrame()

            if not random_control:

                f = pd.read_table(path+data_name+'/'+data_name+'_metadata.tsv'\
                    if not self.args['metadata_name']\
                    else self.args['metadata_name'],sep='\t',header=0, index_col=False)
  
                f.insert(0,'dataset_name',data_name)
                ##print data_name, ' questo e f appena inserito i nome del dataset'

            else:

                total_controls_added = int(selection_criteria.split(':')[0][6:])
                bigger_class = (total_controls_added//(self.args['proportions'][1]+1))*self.args['proportions'][1]

                f = pd.read_table(data_name,sep='\t',header=None,index_col=0).T
                f = f.sample(total_controls_added,random_state=1975)

                for data_name,sample in zip(f.dataset_name.tolist(), f.sampleID.tolist()): self.handle_profile(sample,data_name)
                names = [(self.args['control_names'][0] if i<bigger_class else self.args['control_names'][1]) for i in range(total_controls_added)]

                #if self.args['give_statistics'] and self.args['give_statistics'].startswith('rand'):
                    #print 'Randomized dataset :'
                    #print '%i samples added under name: training_controls' %len([n for n in names if n==self.args['control_names'][0]])
                    #print '%i samples added under name: testing_controls' %len([n for n in names if n==self.args['control_names'][1]])

                f['dataset_name'] = names

                if bool(self.args['control_word']):
                    f.insert(2,'ext_characterization',self.args['control_word'])
                selection_criteria = ':'.join(selection_criteria.split(':')[1:])					

            if ((selection_criteria) and (selection_criteria != 'all')):
                pf = pf.append(pd.DataFrame([s.split(':') for s in selection_criteria.split(',')], index=['select']*(selection_criteria.count(',')+1)))

            app = False

            if selection_criteria == 'all': 
                app = True 

            else:

                for i in range(len(pf)):

                    if pf.index[i] == 'select': app, f = True, f[f[pf.iloc[i,0]].isin(pf.iloc[i,1:])]

            if self.args['output_file']:

                if log: print 'i am the log: ', log + '.'


            return f if app else None

        #*******************************************************************

        selection = dict([(str(n+1),self.args['select'][n]) for n in range(len(self.args['select']))])
        f = pd.DataFrame()

        for i in range(len(self.args['select'])):
            n = str(i+1)

            if selection[n]:
                tab = datasets(self.args['metadata_path'],self.args['dataset_input'][i],selection[n],'calling function dataset number '+str(n), False)
                if self.args['columns'][i]: 
                    for cta in self.args['columns'][i]: tab[cta[0]] = cta[1]

            for sample in tab.sampleID.tolist(): 
                if self.args['change_profile_func']:
                    self.handle_profile_direct_folder(sample, self.args['dataset_input'][i])
                else:
                    self.handle_profile(sample, self.args['dataset_input'][i])

            if len(f) == 0: f = tab
            else:  f = f.append(tab) if (not tab is None) else f


        if self.args['randomize_controls']:  
            ### this block append a randomized controls
            ### from a big large dataset of controls
            act_samples = len(f.sampleID.tolist()) if not self.args['random_controls_multiple_of'] else self.args['random_controls_multiple_of']

            if len(f)==0:
                f = datasets(self.args['metadata_path'],self.args['control_filename'],'select'+str(act_samples*int(self.args['proportions'][0]))+\
                (':study_condition:control' if not self.args['control_conditions'] else self.args['control_conditions']),'calling function for controls',True)
            else:
                f = f.append(datasets(self.args['metadata_path'],self.args['control_filename'],'select'+str(act_samples*int(self.args['proportions'][0]))+\
                (':study_condition:control' if not self.args['control_conditions'] else self.args['control_conditions']),'calling function for controls',True))


        f = f.reset_index(drop=True)

        meta = [s for s in f.columns]
        meta1 = [s for s in self.imp_fields if s in meta]
        meta2 = [s for s in meta if s not in self.imp_fields]

        f = f[meta1+meta2].fillna('NA')
        self.metadata = meta1+meta2
        return f


    def merge_profiles(self,profile_list,profile_type):
        chain = ['python',self.args['merge_profile_exec']]+list(profile_list)+['-of','merged_profiles_'+profile_type+'.csv']
        sp.call(chain) 				


    def get_dataset(self):

        def statistics(dataset, columns):
            title='\n'+self.args['output_file'] if self.args['output_file'] else 'dataset'+'\t'.join(list(set(dataset[columns[0]].tolist())))
            stats='\n#samples\t'+str(len(dataset['sampleID'].tolist())) + title
            s=dict([(c, dict([(f, dataset[c].tolist().count(f)) for f in set(dataset[c].tolist())])) for c in columns]) 
            for c in columns: 
                stats += ''.join(['\n%s:\t'%str(c),'\t'.join(['%s: %i\t'%(str(f),s[c][f]) for f in set(dataset[c].tolist())]),'\n'])
            if self.args['output_file']: 
                with open(self.args['output_file']+'.stats','w') as of: of.write(stats)
            else: print stats 	

        #***********************#

        def shrink_taxa(frame):
            columns_n = [(c.split('|')[-1]) for c in frame.columns.tolist()]
            frame.columns = columns_n
            return frame

        #***********************#

        def profile_(fname): 

            #print fname, ' eil filename: esiste? risp => ', os.path.exists(fname)
	
            df = pd.read_csv(fname, header=0, sep='\t', low_memory=(\
						False if \
						    (not self.args['genefam'] \
				                     or (self.args['change_profile_func'] and self.args['pwyrelab'])) \
                                           else True))



            df = df[1:]
            df = df.reset_index(drop=True)
            df = df.T
            df.columns = df.iloc[0,:]
            df = df.iloc[1:,:]
            
            if self.args['change_profile_func']: 
                df.index = [i+'_profile' for i in df.index.tolist()]

            df.index = [prev.split('_'+'complete' if 'complete' in prev else '_prof')[0] for prev in df.index.tolist()]
            return df

        #***********************#

        self.feat = []

        if bool(self.args['only_metadata']): 
            data = self.create_dataset()

        else:
            data = self.create_dataset()
            pro = pd.DataFrame()

        if bool(self.args['metaphlan']):
            if (self.args['metaphlan'] and (self.args['pwyrelab'] or self.args['genefam'] or self.args['markab'] or self.args['markpres'] or self.args['pwycov']\
	        or self.args['mirna'])) and (not self.args['shrink']): 
                self.args['shrink'] = True
            self.merge_profiles(self.profile_species, 'species')
            pro = profile_('merged_profiles_species.csv')

          
            if (self.args['yes_strain'] and self.args['taxon'] == 's__'): 
                pro = pro[[t for t in pro.columns.tolist() if (self.args['taxon'] in t and (not 't__' in t))]] 

            else: 
                #print [t for t in pro.columns.tolist() if (t.split('|')[-1].startswith(self.args['taxon']))]
                if not self.args['mixed_taxa']:
                    pro = pro[[t for t in pro.columns.tolist() if (t.split('|')[-1].startswith(self.args['taxon']))]]
                else:
                    pro = pro[[t for t in pro.columns.tolist() if (t.split('|')[-1][:3] in \
                   (['k__','p__','c__','o__','f__','g__','s__'] \
                    if not self.args['yes_strain'] else \
                    ['k__','p__','c__','o__','f__','g__','s__','t__']))]]

            data = data.merge( shrink_taxa(pro) if self.args['shrink'] else pro, left_on='sampleID', right_index=True, how='left')
            self.feat += pro.columns.tolist()

        if bool(self.args['markab']):
            self.merge_profiles(self.profile_markab, 'markab')
            pro = profile_('merged_profiles_markab.csv')
	    data = data.merge(pro, left_on='sampleID', right_index=True, how='left')

        if bool(self.args['markpres']):
            self.merge_profiles(self.profile_markpres, 'markpres')
            pro = profile_('merged_profiles_markpres.csv')
            data = data.merge(pro, left_on='sampleID', right_index=True, how='left')
            self.feat += pro.columns.tolist()

        """
        if bool(self.args['pwycov']):
            self.merge_profiles(self.profile_cov, 'cover')
            pro = profile_('merged_profiles_cover.csv') 
            pro.columns = [(i+'-PWYC' if 'PWYC' not in i else i) for i in pro.columns.tolist()]
            data = data.merge(pro, left_on='sampleID', right_index=True, how='left')
        """ 

        if bool(self.args['pwyrelab']):
            self.merge_profiles(self.profile_pwys, 'pwys')
            pro = profile_('merged_profiles_pwys.csv')
            pro.columns = [(i+'-PWY' if 'PWY' not in i else i) for i in pro.columns.tolist()]
            data = data.merge(pro, left_on='sampleID', right_index=True, how='left')
            self.feat += pro.columns.tolist()

        """
        if bool(self.args['pfam']):
            self.merge_profiles(self.profile_pfam, 'pfam')
            pro = profile_('merged_profiles_pfam.csv')
            data = data.merge(pro, left_on='sampleID', right_index=True, how='left')
            self.feat += pro.columns.tolist()
        """

        if bool(self.args['genefam']): 
            self.merge_profiles(self.profile_genes, 'genes')
            pro_ = profile_('merged_profiles_genes.csv')
            data = data.merge(pro_, left_on='sampleID', right_index=True, how='left')
            self.feat += pro_.columns.tolist()
        
        if bool(self.args['mirna']):
            self.merge_profiles(self.profile_mirna, 'mirna')
            pro_ = profile_('merged_profiles_mirna.csv')
            data = data.merge(pro_, left_on='sampleID', right_index=True, how='left')
            self.feat += pro_.columns.tolist()

        if self.args['give_statistics'] and (not self.args['give_statistics'].startswith('rand')):
            statistics(data, self.args['give_statistics'].split(':') )

        if self.args['output_file'] and (os.path.exists(self.args['output_file'])): 
            os.remove(self.args['output_file'])
 
        for exte in ['species','pwys','markab','markpres','cover','pfam','mirna']:
            if os.path.exists('merged_profiles_'+exte+'.csv'): os.remove('merged_profiles_'+exte+'.csv')


        if bool(self.args['wilcoxon']):
            column, dist1, dist2 = self.args['wilcoxon'].split(':')
            ind_1 = data[data[column].isin([dist1])].index.tolist()
            ind_2 = data[data[column].isin([dist2])].index.tolist()

            bonf = float(len(self.feat))
            self.feat = [ff for ff in self.feat if \
	         (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<(0.01 \
                                                                    if not self.args['bonferroni'] else (0.01/bonf)))]
            data = data.reset_index(drop=True)


        if bool(self.args['bonferroni']):
            if bool(self.args['output_file']):
                print 'A usable alpha for this dataset (BOFERRONI): ', 0.01/float(len(self.feat)) 
            else:
                bonferroni = open('bonferroni_alpha.txt', 'w')
                bonferroni.write(str(0.01/float(len(self.feat))))   
                bonferroni.close()


        if self.args['meta_features']:
            pwy_ab = np.count_nonzero(data.loc[:, pwy_feat].values.astype(np.float64), axis=1)
            phlan_ab = np.count_nonzero(data.loc[:, phlan_feat].values.astype(np.float64), axis=1)
            gene_ab = np.count_nonzero(data.loc[:, gene_feat].values.astype(np.float64), axis=1)             
            if np.any(pwy_ab): data['PWY_richness'] = pwy_ab
            if np.any(phlan_ab): data[self.args['taxon']+'richness'] = phlan_ab
            if np.any(gene_ab): data['GnFm_richness'] = gene_ab
            data = data.reset_index(drop=True)

		
        def percentile_normalization(distribution_to_use, datas, feats):
            if feats == []: return ''
            ori_ind = datas.index.tolist()
            cn_ind = datas[datas[distribution_to_use.split(':')[0]].isin([distribution_to_use.split(':')[1]])].index.tolist()
            others_ind = datas[datas[distribution_to_use.split(':')[0]].isin([distribution_to_use.split(':')[2]])].index.tolist()
            f = datas[feats].values
            xT = pd.DataFrame(np.array([[sts.percentileofscore(f[cn_ind, i], f[j, i], kind='mean')\
                for j in ori_ind]for i in range(len(feats))], dtype=np.float64).T, columns=feats) 
            return xT



        if bool(self.args['mapped_percentage']) or bool(self.args['mapped_number']):
            mr, mrp = list(), list()
            for ds,s,nr in zip(data['dataset_name'].tolist(), data['sampleID'].tolist(), data['number_reads'].tolist()):
                with open(self.args['base_path']+ds+'/'+self.args['metaphlan_folder']+s+'/'+s+'_norm_reads.tsv') as tor:
                    reads_mapped = float(tor.readlines()[-1].split()[-1])
                mr.append(reads_mapped)                     
                mrp.append((reads_mapped*100.)/float(nr))

            try:            
                data.insert(10, 'number_of_mapped_by_metaphlan2', mr)
                data.insert(10, 'percentage_of_mapped_by_metaphlan2', mrp)
            except IndexError:
                data.insert(4, 'number_of_mapped_by_metaphlan2', mr)
                data.insert(4, 'percentage_of_mapped_by_metaphlan2', mrp)



        if bool(self.args['feature_selection']):
            with open(self.args['feature_selection']) as fs:
                data = data[self.metadata + [feat.rstrip() for feat in fs.readlines()]]

        if self.args['output_file']: print 'Current dataset: %i features.' % len(self.feat)
        if self.args['output_file']: print 'Current dataset: %i samples.' % len(data.sampleID.tolist())

        if bool(self.args['log_transform']): 
            data.loc[:,self.feat]=data.loc[:,self.feat].apply(lambda d: np.log(np.array(map(float,d), dtype=np.float64)+1.),axis=0)

        if bool(self.args['exclude_samples']):
            data = data[~data['sampleID'].isin(set(self.args['exclude_samples']))] 
            data = data.reset_index(drop=True)  ##  print data
 
        return data



    def send_output(self, data):

        if bool(self.args['dense_matrix']):

            is_negative_above_th = lambda feat_row : ((len(feat_row)-np.count_nonzero(feat_row))/float(len(feat_row)))
            stats = open('gene_family_rarefactions/statistics.txt','w')
            for th in self.args['zero_threshold']:
                dt = data.T
                feats = [ft for ft in self.feat if not is_negative_above_th(dt.loc[ft].astype('float'))>th]
                dt.drop(feats, axis=0, inplace=True)
                print 'Saving threshold: %.2f' %th
                stats.write('cutoff %.2f. #feats %i\n' %(th, len(feats)))
                dt.to_csv(self.args['output_file'][:-4]+'_cutat_%.2f'%th+'.csv',sep='\t',header=False,index=True)
            stats.close()

        else:

            if bool(self.args['grad']):

                #print data
                
                if not self.args['log_gradient']:
                    grad = np.sum([data[gr].astype('float').tolist() for gr in self.args['grad']], axis=0)
                else:

                    log_and_inf = lambda xx : np.array([(-np.log(x) if bool(x) else 0.0) for x in xx], dtype=np.float64)

                    #print (np.mean([data[gr].astype('float').tolist() for gr in self.args['grad']], axis=0))
                    ###print np.nan_to_num(np.log(np.mean([data[gr].astype('float').tolist() for gr in self.args['grad']], axis=0)))

                    #print [n for n in (-np.nan_to_num(np.log(np.mean(\
                    #        [data[gr].astype('float').tolist() for gr in self.args['grad']], axis=0)))*10)], ' WOOOW'
                    #exit(1)
                     ####    np.nan_to_num
                    ######################grad = [n for n in (-np.nan_to_num(log_inf(

                    #print self.args['grad']

                    #exit(1)
                    grad = -log_and_inf(np.mean([np.array(data[gr].astype('float').tolist(), dtype=np.float64)*0.01 for gr in self.args['grad']], axis=0))

                    #grad = np.array(map(int, [n for n in (-np.nan_to_num(np.log(np.mean(\
                    #       [data[gr].astype('float').tolist() for gr in self.args['grad']], axis=0)))*10)]), dtype=np.int64)
                    #for g in grad: 
                    #    print type(g), np.isfinite(g), g
                    #print grad, ' appena definito'
                    ##### exit(1)

                #print grad, '  qui ci siamo....'
                #print grad.shape

                #print [data[gr].astype('float').tolist() for gr in self.args['grad']]
                #exit(1) 
   
                ########if not self.args['log_gradient']:
                ########    data[self.args['grad_col']] = (grad - np.min(grad)) / (np.max(grad) - np.min(grad)) #if not self.args['log_gradient'] else grad
                ########else: 
                ########    data[self.args['grad_col']] = grad

                #print data[self.args['grad_col']].tolist()
                

            if self.args['select_columns']:
                data = data[self.args['select_columns']]

            if self.args['select_columns_from_file']:
                fts = open(self.args['select_columns_from_file'], 'r')

                #with open(self.args['select_columns_from_file'], 'r') as fts:
                     #print [feat.rstrip() for feat in fts.readlines()], ' ma che avra di strano...'
               # for f in fts:
                #    print f.rstrip() in data.columns.tolist()
                
 
                data = data[[f.rstrip() for f in fts.readlines()]]
                fts.close()

            #print ' sono uscio da suo if'

            if self.args['feat_only']:
                data = data[['sampleID'] + self.feat]
            
            elif len(self.args['feat_and_condition']) > 0:
                data = data[self.args['feat_and_condition'] + self.feat] 

            if not self.args['output_file']: 



                datat = data.T
                for i in datat.index.tolist(): 

                    print '\t'.join([i]+list(map(str, datat.loc[i])))

            else:

                if self.args['transpose']:

                    data.to_csv(self.args['output_file'].split('.')[0]+'_Headers.csv',sep='\t',header=True,index=False)

                elif self.args['both']:

                    datat = data.T
                    data.to_csv(self.args['output_file'].split('.')[0]+'_Headers.csv',sep='\t',header=True,index=False)
                    datat.to_csv(self.args['output_file'],sep='\t',header=False,index=True)

                else:

                    datat = data.T
                    datat.to_csv(self.args['output_file'],sep='\t',header=False,index=True)


if __name__=='__main__':
	st = metadataset()
	st.send_output(st.get_dataset())	
