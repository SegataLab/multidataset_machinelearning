#!/usr/bin/env python

import itertools
import utils_mod

## last attempt on gene families
## python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -al rf -do cross_study -g0 nsl:5 -g1 nt:1000 -g2 df -hv 0 -ncores 2 -mf 0.01 -r 1

## pfam ??
##  for adenomas


class cross_study(object):

    def __init__(self, datasets, defined_problem, databases, which_python, n_iter, max_feat, ncores, grid_term0, grid_term1, grid_term2, restrict, how_verbose, in_background, mixed_taxa):

        self.utils = utils_mod.usefullfuncs(datasets, mixed_taxa)
        self.datasets = datasets
        self.couples = list(itertools.combinations_with_replacement(self.datasets, 2))
        self.databases = databases
        self.mixed_taxa = mixed_taxa
        self.which_python = which_python
        self.problem = defined_problem.split(':')
        self.max_feat = max_feat
        self.ncores = ncores
        self.n_iter = n_iter
        self.grid_term0 = grid_term0
        self.grid_term1 = grid_term1
        self.grid_term2 = grid_term2
        self.restrict = restrict
        self.tests = [':'.join(self.problem), ':'.join([self.problem[0],self.problem[1]]),':'.join([self.problem[0],self.problem[2]])] 

        self.how_verbose = how_verbose
        self.in_back = in_background

        #if all([d in ['metaphlan', 'pfam'] for d in databases]): 
        #    self.databases.append(['metaphlan', 'pfam'])
	#if all([d in ['metaphlan', 'pathways'] for d in databases]): 
        #    self.databases.append(['metaphlan', 'pathways'])


    def dataset_name(self, data, db, test):
        name = ['data/datasetfortransfer:' if len(data)==2 else 'data/datasetforlodo:']
        if (len(data)==len(self.datasets)): name.append('_ALL_')
        else: name.append('_AND_'.join(data)+'_' if not self.utils.isonedata(data) else data[0])
        name.append('experiment:'+test)
        name.append('_datasetnumber:'+str(len(data) if not self.utils.isonedata(data) else 1))
        name.append('_databases:'+ (db if isinstance(db, str) else '_'.join(db)))
        if self.mixed_taxa: name.append('_mixed_taxa')
        name.append('.csv')
        return ''.join(name)


    def dataset_name_batch(self, data, db, test):
        name = ['data/datasetforbatch:', '_AND_'.join(data)+'_', 'experiment:batch:'+test.split(':')[1], '_datasetnumber:2']
        name.append('_databases:'+ (db if isinstance(db, str) else '_'.join(db)))
        name.append('.csv')
        return ''.join(name)


    def dataset_cross_study(self, pool, db, test):
        commandline = [self.which_python+'/python','../cmdpy/cmdpy_dataset.py']
        if self.utils.isonedata(pool): commandline.append(pool[0]+'.'+test)

        else:
            for data in pool: commandline.append(data+'.'+test)

        if isinstance(db, str): 
            commandline.append('--'+self.utils.databases[db])
        elif isinstance(db, list):
            for dab in db: commandline.append('--'+self.utils.databases[dab])

        if self.mixed_taxa: commandline.append('--taxon k__')

        commandline.append('-of '+self.dataset_name(pool, db, test))

        return ' '.join(commandline)+'\n'


    def dataset_batch_study(self, pool, db, test):
        commandline = [self.which_python+'/python','../cmdpy/cmdpy_dataset.py']

        for data in pool: commandline.append(data+'.'+test)

        if isinstance(db, str):
            commandline.append('--'+self.utils.databases[db])
        elif isinstance(db, list):
            for dab in db: commandline.append('--'+self.utils.databases[dab])
       
        commandline.append('-of '+self.dataset_name_batch(pool, db, test))
        return ' '.join(commandline)+'\n'
            

    def write_cross_command_lines(self):
        test = self.tests[0]
        for db in self.databases:
            datasetcaller = open('datasetcaller_'+(db if isinstance(db, str) else '_'.join(db))+('.sh' if not self.mixed_taxa else '_mixed_taxa.sh') ,'w')
            datasetcaller.write('#!/bin/bash\n')
            if not self.restrict == 'lodo':
                for pool in self.couples: datasetcaller.write(self.dataset_cross_study(pool, db, test))
            if not self.restrict == 'transfer':
                datasetcaller.write(self.dataset_cross_study(self.datasets, db, test))
            datasetcaller.close()


    def write_batchdata_command_lines(self):
        test = ':'.join(self.tests[0].split(':')[:2])
        print test, ' test per noi...'

        for db in self.databases:

            datasetcaller = open('../datasetcaller_batch_' + (db if isinstance(db, str) else '_'.join(db))+('.sh'), 'w')
            datasetcaller.write('#!/bin/bash\n')

            for pool in self.couples: 
                if not self.utils.isonedata(pool): datasetcaller.write(self.dataset_batch_study(pool, db, test))

            datasetcaller.close()

            
    def output_name(self, db, algo, data, test_i, test):
        if test_i<0: out = ['ml/resultofcrossvalidation']
        else: out = ['ml/resultoftransfer' if len(data)==2 else 'ml/resultoflodo'] 
        if len(data)==2: out.append(('_ON_'.join(data)) if test_i==1 else ('_ON_'.join([data[u] for u in reversed(range(len(data)))])))
        else:
            if test_i<0: out.append('ANYonANY')
            else: out.append('ANYon_'+data[test_i])
        out.append('features:%s' % (db if isinstance(db, str) else '_'.join(db)))
        out.append('experimenttype:%s' %test)
        out.append(algo)
        out.append('grid0:'+self.grid_term0)
        out.append('grid1:'+self.grid_term1)
        out.append('grid2:'+self.grid_term2)
        if self.mixed_taxa: out.append('mixedtaxa')
        return '_'.join(out)


    def output_name_batch(self, db, algo, data, test_i):
        out = ['ml/resultofbatch', '_ON_'.join(data), ]
        out.append('features:%s' % (db if isinstance(db, str) else '_'.join(db)))
        out.append('experimenttype:batch')
        out.append(algo)
        out.append('grid0:'+self.grid_term0)
        out.append('grid1:'+self.grid_term1)
        out.append('grid2:'+self.grid_term2)
        return '_'.join(out)


    def random_forest_experiment(self, db, data, test_i, test):

        line = [self.which_python+'/python','../metaml/classification_dev.py', self.dataset_name(data,db,test), self.output_name(db,'rf',data,test_i,test)]
        line += [('-hv %i -r %i -s roc_auc -l rf -d 1:' % (self.how_verbose, 20 if self.utils.isonedata(data) else 1)) + ':'.join(test.split(':')[:2])]

        if self.max_feat: line += ['-mf '+self.max_feat]
        if (not self.utils.isonedata(data)) and (not test_i<0): line += ['-t dataset_name:%s' %data[test_i]]
        if self.ncores!=10: line += ['-nc %i' % self.ncores] ## if not self.utils.isonedata(data) else ['-nc %i' % (self.n_iter/5)]

        if isinstance(db, list): line.append('-z '+ ':'.join([self.utils.features[dab] for dab in db]))
        elif isinstance(db, str): line.append('-z ' + self.utils.features[db])

        line.append('-'+' '.join(self.grid_term0.split(':')))
        line.append('-'+' '.join(self.grid_term1.split(':')))
        line.append('-'+' '.join(self.grid_term2.split(':')))       
        return ' '.join(line)


    def random_forest_experiment_batch(self, db, data, test_i, test):
        line = [ self.which_python+'/python','../metaml/classification_dev.py', self.dataset_name_batch(data,db,test), self.output_name_batch(db,'rf',data,test_i)]
        line += [ ('-hv %i -r 1 -s roc_auc -l rf -d 1:' % self.how_verbose ) + 'dataset_name:' + data[1]]

        if self.max_feat: line += ['-mf '+self.max_feat]
        if self.ncores!=10: line += ['-nc %i' % self.ncores] 

        if isinstance(db, list): line.append('-z '+ ':'.join([self.utils.features[dab] for dab in db]))
        elif isinstance(db, str): line.append('-z ' + self.utils.features[db])

        line.append('-'+' '.join(self.grid_term0.split(':')))
        line.append('-'+' '.join(self.grid_term1.split(':')))
        line.append('-'+' '.join(self.grid_term2.split(':')))       
        return ' '.join(line) 


    def lsvm_experiment(self, db, data, test_i, test):
        dataset = self.dataset_name(data, db, test)
        line = [self.which_python +'/python','../metaml/classification_dev.py', self.dataset_name(data, db, test), self.output_name(db,'svm',data,test_i,test)]

        line += [('-hv %i -r %i -s roc_auc -l lsvm -d 1:' % (self.how_verbose, 20 if self.utils.isonedata(data) else 1)) + ':'.join(test.split(':')[:2])]
        if (not self.utils.isonedata(data)) and (not test_i<0): line += ['-t dataset_name:%s' %data[test_i]]

        if isinstance(db, list): line.append('-z '+':'.join([self.utils.features[dab] for dab in db]))
        elif isinstance(db, str): line.append('-z '+self.utils.features[db])

        return ' '.join(line)


    def write_cross_study_commands(self, algo):
        test = self.tests[0]
        proc = 0
        nohup_s = ''#'nohup ' if self.in_back else ''
        nohup_e = ''#' &' if self.in_back else ''
        ##' > nohup_%i.out &' if self.in_back else ''

        if (algo == 'rf'): 
            explanation = '_gridterm0:%s_gridterm1:%s_gridterm2:%s' %('_'.join(self.grid_term0.split(':')), '_'.join(self.grid_term1.split(':')), '_'.join(self.grid_term2.split(':')))
            algorithm = self.random_forest_experiment if algo=='rf' else self.lsvm_experiment
        else: 
            raise NotImplementedError('Only rand-oh forest & support crazy horse-vectors are ready so far, and BTW rand-oh is the only one working, so WTF')

        for db in self.databases:
            if not self.restrict == 'lodo':
                ext = open('../transfer_'+(db if isinstance(db, str) else '_'.join(db))+'_'+test+'_'+algo+explanation+('.sh' if not self.mixed_taxa else '_mixed_taxa.sh'),'w')
                if not self.in_back: ext.write('#!/bin/bash\n')  

                for pool in self.couples:
                    if self.utils.isonedata(pool): 
                        ext.write(nohup_s + algorithm(db, pool, 0, test) + nohup_e)
                        proc += 1
                        #nohup_e = (' > nohup_%i.out &' %proc) if self.in_back else ''
                    else: 
                        ext.write(nohup_s + algorithm(db, pool, 1, test) + nohup_e)
                        proc += 1
                        #nohup_e = (' > nohup_%i.out &' %proc) if self.in_back else ''
                        ext.write('\n')
                        ext.write(nohup_s + algorithm(db, pool, 0, test) + nohup_e) # *!!!*
                        proc += 1
                        #nohup_e = (' > nohup_%i.out &' %proc) if self.in_back else ''
                    ext.write('\n')
                ext.close()

            if not self.restrict == 'transfer':
                exl = open('../lodo_'+(db if isinstance(db, str) else '_'.join(db))+'_'+test+'_'+algo+explanation+('.sh' if not self.mixed_taxa else '_mixed_taxa.sh'),'w')
                if not self.in_back: exl.write('#!/bin/bash\n') 
                for i in range(len(self.datasets)): 
                    exl.write(nohup_s + algorithm(db, self.datasets, i, test) + nohup_e + '\n')
                    proc += 1
                    #nohup_e = (' > nohup_%i.out &' %proc) if self.in_back else ''
                exl.write(nohup_s + algorithm(db, self.datasets, -1, test) + nohup_e + '\n')
                proc += 1
                #nohup_e = (' > nohup_%i.out &' %proc) if self.in_back else ''
                exl.close()


    def write_batch_study_commands(self, algo):
        test = self.tests[0]
        class_ = self.tests[0].split(':')[1]

        print ' dentro la funzione che scrive lgi esperimenti:'
        print test
        print class_

        if (algo == 'rf') | (algo.endswith('svm')):
            explanation = '_gridterm0:%s_gridterm1:%s_gridterm2:%s' %('_'.join(self.grid_term0.split(':')), '_'.join(self.grid_term1.split(':')), '_'.join(self.grid_term2.split(':')))
            algorithm = self.random_forest_experiment_batch
        else:
            raise NotImplementedError('Only rand-oh forest & support crazy horse-vectors are ready so far, and BTW rand-oh is the only one working, so WTF')
        
        for db in self.databases:
            ext = open('../batch_'+class_+'_'+(db if isinstance(db, str) else '_'.join(db))+algo+explanation+('.sh'), 'w')
            ext.write('#!/bin/bash\n')
            for pool in [p for p in self.couples if not self.utils.isonedata(p)]:
                ext.write(algorithm(db, pool, 0, test))
                ext.write('\n')
            ext.close()


if __name__=='__main__':
    print 'pappappararararappappa'
