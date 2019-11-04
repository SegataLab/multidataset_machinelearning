#!/usr/bin/env python

import utils
import itertools

class plot_analyses(object):
    def __init__(self, datasets, defined_problem, databases, which_python, grid0, grid1, grid2, max_feat, ncores, n_iter, how_verbose):
        self.utils = utils.usefullfuncs(datasets)
        self.datasets = datasets
        self.databases = databases
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])]
        self.progressive_counter = 0

        ###self.db_name = dict([(db,rdb) for db,ndb in zip(['metaphlan','pfam','genefamilies'],['metaphlan','',''])])
        self.features = dict([(data,feat) for data,feat in zip(['metaphlan','pfam','genefamilies'], ['s__','PF','UniRef90'])])
        self.grid0 = grid0
        self.grid1 = grid1
        self.grid2 = grid2
        self.grid = '_grid0:'+grid0+'_grid1:'+grid1+'_grid2:'+grid2 
        self.max_feat = max_feat
        self.ncores = ncores
        self.n_iter = n_iter
        self.how_verbose = how_verbose

        self.combinations = [list(map(list, itertools.combinations(self.datasets, i))) for i in range(3, len(self.datasets))] 

    def dataset_name(self, db, pool, test):
        name = ['data/datasetforpredictionplot:', '_AND_'.join(pool), 'experiment:%s' %':'.join(test.split(':')[1:]), 'datasetnumber:'+str(len(pool))]
        name.append('_databases:' + ('_'.join(db) if isinstance(db, list) else db))
        name.append('.csv')
        return ''.join(name) 

    def dataset_line(self, db, pool, test):
        if db == 'genefamilies': db = 'genefam'        
        if isinstance(db, list): #db = [('genefam' if db[i]=='genefamilies' else db[i]) for i in range(len(db))]
            raise NotImplementedError
        line = ['python', '../multidataset_machinelearning/cmdpy_dataset.py', '--'+db ] #('_'.join(db) if isinstance(db, list) else db)]
        for ds in pool: line.append(ds+'.'+test)
        line.append('-of %s' %self.dataset_name(db, pool, test))
        return ' '.join(line)

    def write_plot_data(self):
        test = self.tests[0]
        for db in self.databases:
            datasetcaller = open('../plot_datasetcaller'+('_'.join(db) if isinstance(db, list) else db)+'.sh', 'w')
            datasetcaller.write('#!/bin/bash\n') ##export PATH=%s:$PATH\n' %self.which_python)
            for pool_of_datasets in itertools.chain.from_iterable(self.combinations):
                datasetcaller.write(self.dataset_line(db, pool_of_datasets, test)+'\n')
            datasetcaller.close()

    def output_name(self, db, pool, test_i, algo='rf'):
        out = ['ml/plotresult', 'ANYon_'+pool[test_i]]
        out.append('features:%s' %('_'.join(db) if isinstance(db, list) else db))
        out.append('experimenttype:aucimprove_traininglevel:%i_unId:%i' %(len(pool)-1, self.progressive_counter))
        self.progressive_counter += 1
        out.append(algo+self.grid)
        return '_'.join(out)

    def random_forest_support(self, db, pool, test_i, test, algo='rf'):
        line = [self.which_python+'/python', '../metaml/classification_dev.py', self.dataset_name((db if db!= 'genefamilies' else 'genefam'), pool, test), self.output_name(db, pool, test_i, test)]
        line.append('-r 1 -s roc_auc -l %s -t dataset_name:%s -d 1:%s' %(algo, pool[test_i], ':'.join(test.split(':')[:2])))
        g = lambda term : ('-'+term.split(':')[0]+' '+term.split(':')[1]) if len(term.split(':'))==2 else '-'+term
        line.append(g(self.grid0))
        line.append(g(self.grid1))
        line.append(g(self.grid2))
        line += ['-hv %i -nc 1' % (self.how_verbose)]
        if self.max_feat: line += ['-mf '+self.max_feat]
        if isinstance(db, list): line.append('-z '+ ':'.join([self.features[dbs] for dbs in db]))
        elif isinstance(db, str): line.append('-z ' + self.features[db])
        return ' '.join(line)

    def write_any_support_analyses(self, algo='rf'):
        test = self.tests[0]
        for db in self.databases:
            exp = open('../plotexperiment_'+(db if isinstance(db, str) else '_'.join(db))+'_'+test+'_'+'rf'+self.grid+'.sh','w')
            exp.write('#!/bin/bash\n') ##export PATH=%s:$PATH\n' %self.which_python)
            for pool_of_datasets in itertools.chain.from_iterable(self.combinations):
                for test_i in range(len(pool_of_datasets)):
                 exp.write(self.random_forest_support(db, pool_of_datasets, test_i, test)+'\n')
            exp.close()
      
if __name__ == '__main__':
    print 'bubble, bubble, bubble. always he makes bubbles'
