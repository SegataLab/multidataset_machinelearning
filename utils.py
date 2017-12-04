#!/usr/bin/env python

class usefullfuncs(object):
    data_aliases = {'CM_rescignocrc':'CM_Cohort1', 'CM_lilt':'CM_Cohort2'}
    features = dict([(data,feat) for data,feat in zip(['metaphlan','pfam','metaphlan_pfam'], ['s__','PF','s__:PF'])])

    def __init__(self, datasets):
        self.datasets = datasets        
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

    def isonedata(self, s): return True if s[0]==s[1] else False

    def get_lodo_resultnames(self, leftdataset, database, test, algo, grid, feat):
        return '../ml/resultoflodo_ANYon_'+leftdataset+'_features:'+('_'.join(database) if isinstance(database, list) else database)\
								+'_experimenttype:standard_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt'
        ### ml/resultoflodo_ANYon_ZellerG_2014_features:pathways_metaphlan_experimenttype:standard_rf_grid:gini_featsel:auto.txt

    def get_transfer_resultnames(self, pool_of_datasets, db, test, algo, grid, feat):
        return '../ml/resultoftransfer_'+'_ON_'.join(pool_of_datasets)+'_features:'+('_'.join(db) if isinstance(db, list) else db)\
	      +'_experimenttype:standard_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt'

    def get_cross_validation(self, database, algo, grid, feat):
        return '../ml/cross_validation_ANYonANY_features:'+('_'.join(database) if isinstance(database, list) else database)\
              +'_experimenttype:standard_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt' ## # + rf_grid:gini_featsel:auto.txt

    def get_std_support_test(self, database, algo, grid, feat, dataset_on_which, training_level):
        return '../ml/plotresult_ANYon_'+dataset_on_which+'_features:'+(database if isinstance(database, str) else '_'.join(database))\
	    +'_experimenttype:aucimprove_traininglevel:'+str(training_level)+'_unId:*'+'_'+ algo+'_grid:'+grid+'_featsel:'+feat+'.txt'

    def get_all_trans_for_one(self, database, algo, grid, feat, dataset_on_which):
        return '../ml/resultoftransfer_*'+'_ON_'+dataset_on_which+'_features:'+('_'.join(database) if isinstance(database, list) else database)\
									+'_experimenttype:standard_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt'
        ## ml/cross_validation_ANYonANY_features:pathways_experimenttype:standard_rf_grid:gini_featsel:auto.txt

    def transfer_(self, two_datasets, db, test, algo, grid, feat, start):
        score, n_samples = self._auc(self.get_transfer_resultnames(two_datasets, db, test, algo, grid, feat), start)
        coors = self.datasets.index(two_datasets[0]), self.datasets.index(two_datasets[1])
        return score, coors, n_samples

    def lodo_(self, ds, db, test, algo, grid, feat, start):
        score, n_samples = self._auc(self.get_lodo_resultnames(ds, db, test, algo, grid, feat), start)
        coors = self.datasets.index(ds)
        return score, coors, n_samples

    def cross_validation_(self, db, test, algo, grid, feat, start):
        return self._auc(self.get_cross_validation(db, algo, grid, feat), start, True)

    def get_all_trans_for_one(self, db, algo, grid, feat, dataset_on_which):
        return '../ml/resultoftransfer_*'+'_ON_'+dataset_on_which+'_features:'+('_'.join(db) if isinstance(db, list) else db)\
								+'_experimenttype:standard_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt'

    def get_std_support_pattern(self, db, algo, grid, feat, dataset_on_which, training_level):
        return '../ml/plotresult_ANYon_'+dataset_on_which+'_features:'+(db if isinstance(db, str) else '_'.join(db)) \
                + '_experimenttype:aucimprove_traininglevel:'+str(training_level)+'_unId:*'+'_'+algo+'_grid:'+grid+'_featsel:'+feat+'.txt'

    def test_magns(self, db, test, algo, grid, feat, start):
        return dict([(source, self._auc(self.get_transfer_resultnames([source, source], db, test, algo, grid, feat), start)[1]) for source in self.datasets])

if __name__ == '__main__':
    print 'morditi il culo!, passa parola'
