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
        self.profile_species = set()
        self.profile_marka = set()
        self.profile_markp = set()
        self.profile_cov = set()
        self.profile_pwys = set()
        self.profile_genes = set()
        self.profile_pfam = set()
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
        add( '--metadata_name', type=str, default=None)
        #add( '--metaphlan_path', type=str, default='')
        #add( '--pwyrelab_path', type=str, default='')
        #add( '--genefam_path', type=str, default='')
        #add( '--marker_path', type=str, default='')

        # - ACTIVATE TYPES OF PROFILES
        add( '--metaphlan', action='store_true')		
        add( '--pwyrelab', action='store_true')
        add( '--genefam', action='store_true')
        add( '--marka', action='store_true')
        add( '--markp', action='store_true')
        add( '--pwycov', action='store_true')
        add( '--pfam', action='store_true')

        # - NAMES FOR THE PROFILES FOLDERS
        add( '--metaphlan_folder', type=str, default='metaphlan2')
        add( '--pwyrelab_folder', type=str, default='pathways')
        add( '--genefam_folder', type=str, default='gene_families')
        add( '--pwycoverage_folder', type=str, default='humann2_PWYcov')
        add( '--marka_folder', type=str, default='marker_abundance')
        add( '--markp_folder', type=str, default='marker_presence')
        add( '--pfam_folder', type=str, default='pfam_profiles')

        # - WORD AFTER '_' TO DESIGN THE PROFILE FILE
        add( '--metaphlan_title', type=str, default='profile')
        add( '--pwy_title', type=str, default='profile', choices=['profile','complete_profile'])
        add( '--cov_title', type=str, default='profile', choices=['profile','complete_profile'])
        add( '--genefam_title', type=str, default='profile', choices=['profile','complete_profile'])
        add( '--pfam_title', type=str, default='profile')

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
        # - FURTHER METAPHLAN SPECIFICS
        add( '-mx', '--mixed_taxa', action='store_true')
        add( '-tx', '--taxon', type=str, choices=['k__','p__','c__','o__','f__','g__','s__'], default='s__')
        add( '-sr', '--shrink', action='store_true')
        # - only metadata is for prepare the bigger group to
        # - randomize after. 	
        add( '-x', '--exclude_samples', nargs='+', type=str, default=[]) ## options to avoid as much as possible 
        add( '-om', '--only_metadata', action='store_true', help='This option is for gen control to ramdomize after.')
        add( '-mgp', '--merge_profile_exec', default='/scratchCM/users/paolo.manghi/cmdpy/merge_metaphlan_tables02.py')

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
            , help='either a series of features to eval, or a file with a list of features')		

        # - WILKOXON RANK FOR TOO LARGE FEATURE SETS
        add( '--not_gene_wilcoxon', action='store_false')
        add( '--pwy_wilcoxon', action='store_true')
        add( '--cov_wilcoxon', action='store_true')
        add( '--taxa_wilcoxon', action='store_true')
        add( '--not_mark_wilcoxon', action='store_false')
        add( '-wk', '--wilcoxon', type=str, default=None, help='--wilcoxon study_condition:CRC:control')

        # - OUTPUT SPECIFICS
        add( '-t', '--transpose', action='store_true', help='Output with fields on headers instead of first column.')
        add( '-b', '--both', action='store_true', help='Outputs the straight output plus the one with fields as headers.')
        add( '-of', '--output_file', type=str, default=None)
        add( '-gs', '--give_statistics', type=str, default=None\
          , help='Statistics on metadata -gs study_condition:sample_type:disease, or -gs random/randomization/rand/randrandrand/randanything')

        pp = vars(p.parse_args())
        if ( not pp['metaphlan'] ) and ( not pp['pwyrelab'] ) and ( not pp['genefam']) and \
           ( not pp['marka']) and ( not pp['markp']) and ( not pp['pwycov']) and ( not pp['only_metadata'] and (not pp['pfam']) ):
           #print 'Automatic settings on taxa relative abundance.'
           pp['metaphlan'] = True
        pp['proportions'] = list(map(int, pp['proportions'].split(',')))
        if len(pp['proportions'])!=2: raise IOError('proportions must contain 2 numbers separated by a comma.')
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
	
    def handle_profile(self, sample, dataset_name):
        if self.args['metaphlan'] and (not sample in self.args['exclude_samples']): 
            self.profile_species.add(self.args['base_path']+dataset_name+'/'+self.args['metaphlan_folder']+'/'+sample+'/'+sample+'_'+self.args['metaphlan_title']+'.txt')
        if self.args['pfam'] and (not sample in self.args['exclude_samples']):
            self.profile_pfam.add(self.args['base_path']+dataset_name+'/'+self.args['pfam_folder']+'/'+sample+'/'+sample+'_'+self.args['pfam_title']+'.txt')
        if self.args['marka'] and (not sample in self.args['exclude_samples']):
            self.profile_marka.add(self.args['base_path']+dataset_name+'/'+self.args['marka_folder']+'/'+sample+'/'+sample+'_profile'+'.txt')
        if self.args['markp'] and (not sample in self.args['exclude_samples']):
            self.profile_markp.add(self.args['base_path']+dataset_name+'/'+self.args['markp_folder']+'/'+sample+'/'+sample+'_profile'+'.txt')
        if self.args['pwycov'] and (not sample in self.args['exclude_samples']):
            self.profile_cov.add(self.args['base_path']+dataset_name+'/'+self.args['pwycoverage_folder']+'/'+sample+'/'+sample+'_'+self.args['cov_title']+'.txt')
        if self.args['pwyrelab'] and (not sample in self.args['exclude_samples']): 
            self.profile_pwys.add(self.args['base_path']+dataset_name+'/'+self.args['pwyrelab_folder']+'/'+sample+'/'+sample+'_'+self.args['pwy_title']+'.txt')
        if self.args['genefam'] and (not sample in self.args['exclude_samples']): 
            self.profile_genes.add(self.args['base_path']+dataset_name+'/'+self.args['genefam_folder']+'/'+ sample+'/'+ sample+'_'+self.args['genefam_title']+'.txt')
 
    def create_dataset(self):
        #*********************************************
        def datasets(path,data_name,selection_criteria,log,random_control=False):
            pf = pd.DataFrame()
            if not random_control:
                f = pd.read_table(path+data_name+'/'+data_name+'_metadata.txt'\
                if not self.args['metadata_name']\
                else self.args['metadata_name'],sep='\t',header=0, index_col=False)
                f.insert(0,'dataset_name',data_name)
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
            if selection_criteria:
                pf = pf.append(pd.DataFrame([s.split(':') for s in selection_criteria.split(',')], index=['select']*(selection_criteria.count(',')+1)))
            app = False
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
                tab = datasets(self.args['base_path'],self.args['dataset_input'][i],selection[n],'calling function dataset number '+str(n), False)
                if self.args['columns'][i]: 
                    for cta in self.args['columns'][i]: tab[cta[0]] = cta[1]
            for sample in tab.sampleID.tolist(): self.handle_profile(sample,self.args['dataset_input'][i])
            if len(f) == 0: f = tab
            else:  f = f.append(tab) if (not tab is None) else f
        if self.args['randomize_controls']:  
            ### this block append a randomized controls
            ### from a big large dataset of controls
            act_samples = len(f.sampleID.tolist()) if not self.args['random_controls_multiple_of'] else self.args['random_controls_multiple_of']
            if len(f)==0:
                f = datasets(self.args['base_path'],self.args['control_filename'],'select'+str(act_samples*int(self.args['proportions'][0]))+\
                (':study_condition:control' if not self.args['control_conditions'] else self.args['control_conditions']),'calling function for controls',True)
            else:
                f = f.append(datasets(self.args['base_path'],self.args['control_filename'],'select'+str(act_samples*int(self.args['proportions'][0]))+\
                (':study_condition:control' if not self.args['control_conditions'] else self.args['control_conditions']),'calling function for controls',True))
        f = f.reset_index(drop=True)
        meta = [s for s in f.columns]
        meta1 = [s for s in self.imp_fields if s in meta]
        meta2 = [s for s in meta if s not in self.imp_fields]
        f = f[meta1+meta2].fillna('na')
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
            df = pd.read_csv(fname, header=0, sep='\t', low_memory=False if not self.args['genefam'] else True)
            df = df[1:]
            df = df.reset_index(drop=True)
            df = df.T
            df.columns = df.iloc[0,:]
            df = df.iloc[1:,:]
            df.index = [prev.split('_'+'complete' if 'complete' in prev else '_prof')[0] for prev in df.index.tolist()]
            return df
        #***********************#
        self.feat = []
        if bool(self.args['only_metadata']): data = self.create_dataset()
        else:
            data = self.create_dataset()
            pro = pd.DataFrame()

        if bool(self.args['metaphlan']):
            if (self.args['metaphlan'] and (self.args['pwyrelab'] or self.args['genefam'] or self.args['marka'] or self.args['markp'] or self.args['pwycov']\
	        )) and (not self.args['shrink']): self.args['shrink'] = True
            self.merge_profiles(self.profile_species, 'species')
            pro = profile_('merged_profiles_species.csv')
            pro = pro[[t for t in pro.columns.tolist() if (self.args['taxon'] in t and (not 't__' in t))]] 
            data = data.merge( shrink_taxa(pro) if self.args['shrink'] else pro, left_on='sampleID', right_index=True, how='left')
            self.feat += pro.columns.tolist()

        """
        if bool(self.args['marka']):
            self.merge_profiles(self.profile_marka, 'marka')
            pro = profile_('merged_profiles_marka.csv')
	    data = data.merge(profile_('merged_profiles_marka.csv'), left_on='sampleID', right_index=True, how='left')
        if bool(self.args['markp']):
            self.merge_profiles(self.profile_markp, 'markp')
            pro = profile_('merged_profiles_markp.csv')
            data = data.merge( profile_('merged_profiles_markp.csv'), left_on='sampleID', right_index=True, how='left')
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
            #print 'merging profiles...'
            self.merge_profiles(self.profile_genes, 'genes')
            #print 'merged finished...'
            pro_ = profile_('merged_profiles_genes.csv')
            data = data.merge(pro_, left_on='sampleID', right_index=True, how='left')
            self.feat += pro_.columns.tolist()

        if self.args['give_statistics'] and (not self.args['give_statistics'].startswith('rand')): statistics(data, self.args['give_statistics'].split(':') )
        if self.args['output_file'] and (os.path.exists(self.args['output_file'])): 
            os.remove(self.args['output_file'])
        for exte in ['species','pwys','marka','markp','cover','pfam']:
            if os.path.exists('merged_profiles_'+exte+'.csv'): os.remove('merged_profiles_'+exte+'.csv')
        if self.args['pwyrelab']: 
            self.feat = [(f+'-PWY' if ( (not 'PWY' in f) and (not self.args['taxon'] in f) and (not 'UniRef90' in f) & (not 'marka' in f) \
                and (not 'markp' in f) and (not 'PWYC' in f) & (not 'PF' in f)) else f) for f in self.feat]	
        data = data.reset_index(drop=True)

        if bool(self.args['wilcoxon']):
            column, dist1, dist2 = self.args['wilcoxon'].split(':')
            ind_1 = data[data[column].isin([dist1])].index.tolist()
            ind_2 = data[data[column].isin([dist2])].index.tolist()
            if not self.args['not_gene_wilcoxon']:				 
                gene_feat = [ff for ff in gene_feat if (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<=0.00001)]
            if not self.args['not_mark_wilcoxon']:
                marka_feat = [ff for ff in marka_feat if (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<=0.00001)]
                markp_feat = [ff for ff in marka_feat if (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<=0.00001)]
            if self.args['taxa_wilcoxon']:
                phlan_feat = [ff for ff in phlan_feat if (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<=0.00001)]
            if self.args['pwy_wilcoxon']:
                pwy_feat = [ff for ff in pwy_feat if (sts.mannwhitneyu(data[ff].iloc[ind_1], data[ff].iloc[ind_2], alternative='two-sided')<=0.00001)]
            self.feat = phlan_feat + pwy_feat + gene_feat + marka_feat + markp_feat + cov_feat
            data = data.reset_index(drop=True)

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

        if bool(self.args['feature_selection']):
            with open(self.args['feature_selection']) as fs:
                f_imp = True if len(fs.readline().rstrip().split())==2 else False 
                line = True
                features = [line.rstrip() for line in fs.readlines()] if not f_imp else dict([tuple(line.rstrip().split()[:2]) for line in fs.readlines()])	
                if isinstance(features, list):
                    features,data,self.feat = [o for o in (set(features) & set(self.feat))], data[self.metadata+features], features
                elif isinstance(features, dict): 
                    for k in features: 
                        if k not in (set(features) & set(feat)): del features[k]
                    data = data[self.metadata + [k for k in features.keys()]]
                    for k in features.keys(): data.loc[:,k].apply(lambda a : a*features[k]*10 )##, axis=1)
                    feat = [k for k in features.keys()]

        data.drop(data.loc[:,self.feat].columns[data.loc[:,self.feat].max().astype(np.float64)==data.loc[:,self.feat]\
	        .min().astype(np.float64)],axis=1,inplace=True)

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
