#!/usr/bin/env python

# grid0, grid1, grid2 are terms present in the various file names
# for parameters to be tuned.

class usefullfuncs(object):
    #result_folder = '../ml/'

    data_aliases = {'CM_rescignocrc':'CM_Cohort2', 'CM_lilt':'CM_Cohort1'}
    #features = dict([(data,feat) for data,feat in zip(\
    #		['metaphlan','pfam','metaphlan_pfam','genefamilies','pathways','markers']\
    #	      , ['s__','PF','s__:PF','UniRef90','PWY','MK'])])

    databases = dict([(db, dbname) for db,dbname in zip(\
		['metaphlan','pathways','genefamilies','markers','pfam']\
	      , ['metaphlan','pwyrelab','genefam','markpres','pfam'])])

    def __init__(self, datasets, mixed_taxa=False):
        self.datasets = datasets        

        self.features = dict([(data,feat) for data,feat in zip(\
                ['metaphlan','pfam','genefamilies','pathways','markers']\
              , [('s__' if not mixed_taxa else 'k__'),'PF','UniRef90','PWY','MK'])])


        #self.problem = defined_problem.split(':')
        #self.tests = [':'.join(self.problem), ':'.join([self.problem[0], self.problem[1]]), ':'.join([self.problem[0], self.problem[2]])]

    def _auc(self, result, start, crossv=False):
        with open(result) as r:
            for s in range(start):
                line = r.readline()
                if s==0:
                    n_samples = int(line.split()[1])
            if crossv:
                n_features = int(line.rstrip().split()[1])
            else:
                n_features = 'x'
            while line[0]!='auc':
                line = r.readline().rstrip().split()
            if not crossv:
                return float(line[1]), n_samples
            else:
                return float(line[1]), n_features, n_samples

    def all_auc(self, result):
        print 'PARSING RESULT NAMED: ', result

        n_f, au, stds = [], [], []
        with open(result) as r:
            line = [' ']
            while (not line[0].startswith('Feature')):
                line = r.readline().rstrip().split()
                if line[0] == '#features': n_f.append(int(line[1]))
                elif line[0] == 'auc': 
                    au.append(float(line[1]))
                    stds.append(float(line[2]))
        return [[n,u,s] for n,u,s in zip(n_f, au, stds)]

    def all_auc_from_lodo(self, leftout, db, algo, test, grid0, grid1, grid2): 
        return self.all_auc(self.get_lodo_resultnames(leftout, db, algo, test, grid0, grid1, grid2))

    def all_auc_from_transfer(self, pool, db, algo, test, grid0, grid1, grid2):
        return self.all_auc(self.get_transfer_resultnames(pool, db, algo, test, grid0, grid1, grid2))

    def isonedata(self, s): return True if s[0]==s[1] else False

    def get_lodo_resultnames(self, leftdataset, database, algo, test, grid0, grid1, grid2):
        return '../ml/resultoflodo_ANYon_'+leftdataset+'_features:'+('_'.join(database) if isinstance(database, list) else database)\
	      +'_experimenttype:'+test+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

        ### ml/resultoflodo_ANYon_ZellerG_2014_features:genefamilies_experimenttype:standard_rf_grid0:c:gini_grid1:nsl:5_grid2:df

    def get_transfer_resultnames(self, pool_of_datasets, db, algo, test, grid0, grid1, grid2):
        return '../ml/resultoftransfer_'+'_ON_'.join(pool_of_datasets)+'_features:'+('_'.join(db) if isinstance(db, list) else db)\
	      +'_experimenttype:'+test+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

    def get_cross_validation(self, database, algo, test, grid0, grid1, grid2):
        return '../ml/resultofcrossvalidation_ANYonANY_features:'+('_'.join(database) if isinstance(database, list) else database)\
              +'_experimenttype:'+test+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

        ### resultofcrossvalidation_ANYonANY_features:metaphlan_experimenttype:standard_rf_grid0:nt:500_grid1:nsl:5_grid2:c:entropy

    def get_onlyone_batch_resultname(self, pool_of_datasets, db, algo, test_word, grid0, grid1, grid2):
        return '../ml/resultofbatch_'+'_ON_'.join(pool_of_datasets)+'_features:'+('_'.join(db) if isinstance(db, list) else db)\
	      +'_experimenttype:batch_on_'+test_word+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

        ###  ../ml/resultofbatch_ZellerG_2014_ON_FengQ_2015_features:metaphlan_experimenttype:batch_rf_grid0:c:entropy_grid1:nt:1000_grid2:nsl:5.txt

    def get_std_support_test(self, database, algo, grid0, grid1, grid2, dataset_on_which, training_level):
        return '../ml/plotresult_ANYon_'+dataset_on_which+'_features:'+(database if isinstance(database, str) else '_'.join(database))\
	      +'_experimenttype:aucimprove_traininglevel:'+str(training_level)+'_unId:*'+'_'+ algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

    def get_all_trans_for_one(self, database, algo, test, grid0, grid1, grid2, dataset_on_which):
        return '../ml/resultoftransfer_*'+'_ON_'+dataset_on_which+'_features:'+('_'.join(database) if isinstance(database, list) else database)\
	      +'_experimenttype:'+test+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

        ###  ml/cross_validation_ANYonANY_features:pathways_experimenttype:standard_rf_grid:gini_featsel:auto.txt

    def transfer_(self, two_datasets, db, algo, test, grid0, grid1, grid2, start, changed_coordinates=[]):
        score, n_samples = self._auc(self.get_transfer_resultnames(two_datasets, db, algo, test, grid0, grid1, grid2), start)
        if not changed_coordinates:
            coors = self.datasets.index(two_datasets[0]), self.datasets.index(two_datasets[1])
        else:
            coors = changed_coordinates.index(two_datasets[0]), changed_coordinates.index(two_datasets[1])
        return score, coors, n_samples

    def lodo_(self, ds, db, algo, test, grid0, grid1, grid2, start, changed_coordinates=[]):
        score, n_samples = self._auc(self.get_lodo_resultnames(ds, db, algo, test, grid0, grid1, grid2), start)
        if not changed_coordinates:
            coors = self.datasets.index(ds)
        else:
            coors = changed_coordinates.index(ds)
        return score, coors, n_samples

    def batch_(self, two_datasets, db, algo, test, grid0, grid1, grid2, start, changed_coordinates=[]):
        score, n_samples = self._auc(self.get_onlyone_batch_resultname(two_datasets, db, algo, test, grid0, grid1, grid2), start)
        if not changed_coordinates:
            coors = self.datasets.index(two_datasets[0]), self.datasets.index(two_datasets[1])
        else:
            coors = changed_coordinates.index(two_datasets[0]), changed_coordinates.index(two_datasets[1])
        return score, coors, n_samples

    def cross_validation_(self, db, algo, test, grid0, grid1, grid2, start):
        return self._auc(self.get_cross_validation(db, algo, test, grid0, grid1, grid2), start, True)

    def get_all_trans_for_one(self, db, algo, test, grid0, grid1, grid2, dataset_on_which):
        return '../ml/resultoftransfer_*'+'_ON_'+dataset_on_which+'_features:'+('_'.join(db) if isinstance(db, list) else db)\
		+'_experimenttype:'+test+'_'+algo+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

    def get_std_support_pattern(self, db, problem, grid0, grid1, grid2, dataset_on_which, training_level):
        return '../ml/plotresult_ANYon_'+dataset_on_which+'_features:'+(db if isinstance(db, str) else '_'.join(db)) \
                + '_experimenttype:aucimprove_traininglevel:'+str(training_level)+'_unId:*'+'_'+problem+'_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2+'.txt'

    def test_magns(self, db, algo, test, grid0, grid1, grid2, start):
        return dict([(source, self._auc(self.get_transfer_resultnames([source, source], db, algo, test, grid0, grid1, grid2), start)[1]) for source in self.datasets])



class effect_size(object):

    def __init__(self, e, n, c, p):
        self.esize = e
        self.name = n
        self.cl = c
        self.pv = p



class lefse_reader(object):

    def __init__(self, lefse_output):
         self.to_read = lefse_output
         self.data = list()
         self.classes_ = set()
         ##self.howmany = howmany


    def select_highest(self, howmany=5):
        feats = list()
        for c in self.classes_:

            if howmany > 0:
                feats += [feat[1] for feat in sorted([(d.esize, d.name) for d in self.data if d.cl==c]\
                    , key = lambda obj : obj[0], reverse=True)][:howmany]
            else:
                feats += [feat[1] for feat in sorted([(d.esize, d.name) for d in self.data if d.cl==c]\
                    , key = lambda obj : obj[0], reverse=True)]

            #print ' inside: ', feats
            #exit(1)

        self.feats = feats
        self.data = [d for d in self.data if d.name in feats]


    def get(self):
        line = True
        with open(self.to_read) as ip:
            while line:
                line = ip.readline().rstrip().split()
                if len(line) == 5:
                    feat = str(line[0]) if not line[0].startswith('s__') else str(line[0][3:])
                    class_ = str(line[2])
                    pvalue = float(line[-1])
                    e_size = float(line[1])
                    self.data.append(effect_size(e_size, feat, class_, pvalue))
                    self.classes_.add(class_)




class project_color_code(object):

    def __init__(self, project_name):

        if project_name == 'periimplantitis':
            self.color_code = dict([(cat,cl) for cat,cl in zip(\
                ['healthy', 'peri-implantitis', 'mucositis']\
              , ['green', 'red', 'deepskyblue']\
               )])

        



if __name__ == '__main__':
    print 'morditi il culo!, for parameters to be tuned. passa parola'
