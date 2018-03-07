import sys
import glob
import os
import argparse as ap

def marker_abundance(sample,name):
    header = '#SampleID\tMarker_Abundance_Analysis\n'
    s = open(name+'_profile.txt','w')
    s.write(header)
    with open(sample) as ab:
        line = ab.readline()
        while bool(line):
            line = ab.readline().rstrip()
            if line:
                s.write('MK-'+line.rstrip()+'\n')   
    s.close()

def marker_presence(sample,name):
    s = open(name+'_profile.txt','w')
    header = '#SampleID\tMarker_Presence_Analysis\n'
    s.write(header)
    with open(sample) as pr:
        line = pr.readline()
        while bool(line):
            line = pr.readline().rstrip()
            if line:
                s.write('MK-'+line.rstrip()+'\n')
    s.close()


#def coverage(sample,name):
#    s = open(name+'_profile.txt', 'w')
#    header = '#SampleId\tHumanN2_Analysis\n'
#    s.write(header)
#    with open(sample) as ab:
#        line = ab.readline().rstrip()
#        while bool(line):
#            line = ab.readline().rstrip()
#            if '|' in line: pass
#            elif not line: pass
#            elif 'UNMAPPED' in line: pass
#            elif 'UNINTEGRATED' in line: pass
#            else:
#                pwy = line.split(':')[0].replace('PWY','PWYC')
#		if (not 'PWYC' in pwy):
#                    pwy = 'PWYC-' + pwy
#                abundance = line.split()[-1]
#                s.write("\t".join([pwy,abundance])+'\n')   
#    s.close()

#def coverage_and_bugs(sample,name):
#    s = open(name+'_complete_profile.txt', 'w')
#    header = '#SampleId\tHumanN2_Analysis\n'
#    s.write(header)
#    line = True
#    with open(sample) as ab:
#        while bool(line):
#            line = ab.readline().rstrip()
#            if not line.count('|'):
#                pass
#            else:
#                if len(line.split(':')) > 1:
#                    pathway = line.split(':')[0].replace('PWY','PWYC')
#                    if (not 'PWYC' in pathway): 
#                        pathway = 'PWYC-' + pathway 
#                    taxon,abundance = tuple(line.split('|')[-1].split())
#                    taxon = taxon.split('.')[-1]
#                    s.write("".join([pathway,'|',taxon,'\t',abundance])+'\n')
#    s.close()

#def pathways_and_bugs(sample,name):
#    s = open(name+'_complete_profile.txt', 'w')
#    header = '#SampleId\tHumanN2_Analysis\n'
#    s.write(header)
#    line = True
#    with open(sample) as ab:
#        while bool(line):
#            line = ab.readline().rstrip()
#            if not line.count('|'):
#                pass
#            else:
#                if len(line.split(':')) > 1:
#                    pathway = line.split(':')[0]
#                    taxon,abundance = tuple(line.split('|')[-1].split())
#                    taxon = taxon.split('.')[-1]
#                    s.write("".join([pathway,'|',taxon,'\t',abundance])+'\n')
#                else:
#                    if line.startswith('UNINTEGRATED|'):
#                        taxon,abundance = tuple(line.split('|')[-1].split())
#                        taxon = taxon.split('.')[-1]
#                        s.write("".join(['UNINTEGRATED_PWY|',taxon,'\t',abundance])+'\n')
#    s.close()


def pathways(sample,name):
    s = open(name+'_profile.txt', 'w')
    basename = lambda x : x.split('/')[-1] 
    header = '#sampleID\t%s\n' %basename(sample)[:-4]
    ss = []
    SUM = 0.
    ss.append(header)
    with open(sample) as ab:
        line = ab.readline().rstrip()
        while bool(line):
            line = ab.readline().rstrip()
            if '|' in line: pass
            elif not line: pass
            elif 'UNMAPPED' in line: pass
            elif 'UNINTEGRATED' in line: pass
            else:
                pwy = line.split(':')[0]
                abundance = float(line.split()[-1])
                ##ss.append("\t".join([pwy,abundance])+'\n')
		ss.append([pwy, abundance])
		SUM += abundance
    s.write(ss[0])
    for a in range(1, len(ss)):
        s.write('\t'.join([ss[a][0], str(ss[a][1]/float(SUM))])+'\n')
    s.close()

def genefamilies(sample, name):
    s = open(name+'_profile.txt', 'w')
    basename = lambda x : x.split('/')[-1]
    header = '#sampleID\t%s\n' %basename(sample)[:-4]
    ss = []
    SUM = 0.
    ss.append(header)
    #s.write(header)
    with open(sample) as ab:
        line = ab.readline()
        while bool(line):
            if not line.startswith("#"):
                if not 'UNMAPPED' in line:#line.startswith("UNMAPPED"):
                  #if not '_unknown\t' in line:
                  if '|' not in line:
                      line = line.split()
                      gene = line[0]
		      abundance = float(line[1])
                      SUM += abundance
                      ss.append([gene, abundance])
                      #line = "GnFm_" + "_".join(line[0].split('_')[1:]) + '\t' + line[1]
                      #s.write(line+'\n')
            line = ab.readline().rstrip()
    s.write(ss[0])
    for a in range(1, len(ss)):
        s.write('\t'.join([ss[a][0], str(ss[a][1]/float(SUM))])+'\n') 
    s.close()

#def genefamilies_and_bugs(sample, name):
#    s = open(name+'_complete_profile.txt', 'w')
#    header = '#SampleId\tGeneFamilies_Analysis\n'
#    s.write(header)
#    with open(sample) as ab:
#        line = ab.readline()
#        while bool(line):
#            if not line.startswith("#"):
#                if not line.startswith("UNMAPPED"):
#                    if '|' in line:
#                        line = line.split()
#                        line = "GnFm_" + "_".join(line[0].split('_')[1:]) + '\t' + line[1]
#                        s.write(line+'\n')
#            line = ab.readline().rstrip()
#    s.close()	
				
def read_params():					
   parser = ap.ArgumentParser()
   add = parser.add_argument
   add('dataset', type=str)

   ###add('--gene_families_norm_folder', type=str, default='gene_families')
   add('--pathways_f', type=str, default='pathways') #'humann2_PWYrelab')
   add('--markerpres_f', type=str, default='marker_presence')
   ##add('--markerab_f', type=str, default='marker_abundance')
   add('--genefamilies_f', type=str, default='gene_families')#'humann2_GeneFam')

   ##add('--other_datasets', nargs='+', default=[])

   add('--origin_path', type=str, default='/scratchCM/tmp_projects/pmanghi-epasolli-humann2profiling/')
   add('--destination_path', type=str, default='/CM/data/meta/')  #'/scratchCM/data/meta/')
   add('--gene_fish', type=str, default='genefamilies_relab' )
   add('--pwy_fish', type=str, default='pathabundance_relab')
   add('--markerpres_fish', type=str, default='marker_presence')
   ##add('--markerab_fish', type=str, default='marker_abundance')

   return vars(parser.parse_args())  


if __name__=='__main__':

    par = read_params()  

    def normal():
        #par = read_params()       
        #os.mkdir(par['destination_path']+par['dataset'])

        try:
            os.mkdir(par['destination_path'] + par['dataset'] + '/' + par['genefamilies_f'])
        except:
            print 'dir ' + par['destination_path'] + par['dataset'] + '/' + par['genefamilies_f'] + '/ present.'

        try:
            os.mkdir(par['destination_path'] + par['dataset'] + '/' + par['pathways_f'])
        except:
            print 'dir ' + par['destination_path']+par['dataset']+'/'+par['pathways_f'] + '/ present.'

        #try:
        #    os.mkdir(par['destination_path']+par['dataset']+'/'+par['markerpres_f'])
        #except:
        #    print 'dir ' + par['destination_path']+par['dataset']+'/'+par['markerpres_f'] + '/ present.'
        
        #os.mkdir(par['destination_path']+par['dataset']+'/'+par['pathways_f'])
        #os.mkdir(par['destination_path']+par['dataset']+'/'+par['coverage_f'])
        ### os.mkdir(par['destination_path']+par['dataset']+'/'+par['markerpres_f'])
        ###os.mkdir(par['destination_path']+par['dataset']+'/'+par['markerab_f'])
       

        gn_samples = par['origin_path'] + par['dataset'] + '/' + par['gene_fish'] + '/*.tsv'
        pwy_samples = par['origin_path'] +par['dataset']+'/'+ par['pwy_fish'] + '/*.tsv'
        #mar_samples = par['origin_path'] +par['dataset']+'/'+ par['markerpres_fish'] + '/*.tsv'

        basename = lambda s : s.split('/')[-1]

        """
        #markab_samples = par['origin_path'] +par['dataset']+'/'+ par['markerab_fish'] + '/*.tsv'
        #markpres_samples = par['origin_path'] +par['dataset']+'/'+ par['markerpres_fish'] + '/*.tsv'
        basename = lambda s : s.split('/')[-1]

        
        for t in glob.glob(cov_samples):
            os.mkdir(par['destination_path']+par['dataset']+'/'+par['coverage_f']+'/'+basename(t)[:-4])
            coverage_and_bugs(t,par['destination_path']+par['dataset']+'/'+par['coverage_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
            coverage(t,par['destination_path']+par['dataset']+'/'+par['coverage_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
        print 'done coverage'
        """

        ##### NORMAL
        for t in glob.glob(pwy_samples):
            print 'sample ', t
            try:
                os.mkdir(par['destination_path']+par['dataset']+'/'+par['pathways_f']+'/'+basename(t)[:-4])        
            except: print 'dir ' + par['destination_path']+par['dataset']+'/'+par['pathways_f']+'/'+basename(t)[:-4] + '/ is present'
            pathways(t, par['destination_path']+par['dataset']+'/'+ par['pathways_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])

        for t in glob.glob(gn_samples):
            try:
                os.mkdir(par['destination_path']+par['dataset']+'/'+par['genefamilies_f']+'/'+basename(t)[:-4])    
            except: print 'dir ' + par['destination_path']+par['dataset']+'/'+par['genefamilies_f']+'/'+basename(t)[:-4] + '/ is present'
            genefamilies(t, par['destination_path']+par['dataset']+'/'+ par['genefamilies_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])


        #for t in glob.glob(mar_samples):
        #    try:
        #        os.mkdir(par['destination_path']+par['dataset']+'/'+par['markerpres_f']+'/'+basename(t)[:-4])    
        #   except: print 'dir ' + par['destination_path']+par['dataset']+'/'+par['markerpres_f']+'/'+basename(t)[:-4] + '/ is present'

        
        #    os.mkdir(par['destination_path']+par['dataset']+'/'+par['pathways_f']+'/'+basename(t)[:-4])
        #####    #pathways_and_bugs(t,par['destination_path']+par['dataset']+'/'+par['pathways_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
        #    pathways(t,par['destination_path']+par['dataset']+'/'+par['pathways_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
        #print 'done pathways'
        #for t in glob.glob(gn_samples):
        #    os.mkdir(par['destination_path']+par['dataset']+'/'+par['genefamilies_f']+'/'+basename(t)[:-4])
        #####    #genefamilies_and_bugs(t,par['destination_path']+par['dataset']+'/'+par['genefamilies_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
        #    genefamilies(t,par['destination_path']+par['dataset']+'/'+par['genefamilies_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
        #print 'done genefamilies'

        #for dab,fo,func in zip([markab_samples, markpres_samples],[par['markerab_f'],par['markerpres_f']],[marker_abundance, marker_presence]):

                
        #fo = par['gene_families_norm_folder']
        #os.mkdir(par['destination_path']+par['dataset']+'/'+fo)
        #func = 
        #for t in glob.glob(gn_samples):
        #    os.mkdir(par['destination_path']+par['dataset']+'/'+fo+'/'+basename(t)[:-4])
        #    func(t, par['destination_path']+par['dataset']+'/'+fo+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])


    def subsampled():
        #par = read_params()
        for dada in par['other_datasets']:
            os.mkdir(par['destination_path'] + par['dataset'] + '/' + dada + '/' + par['genefamilies_f'])
            os.mkdir(par['destination_path'] + par['dataset'] + '/' + dada + '/' + par['pathways_f'])
            gn_samples = par['origin_path'] + dada + '/' + par['gene_fish'] + '/*.tsv'
            pwy_samples = par['origin_path'] + dada + '/' + par['pwy_fish'] + '/*.tsv'
            basename = lambda s : s.split('/')[-1]

            for t in glob.glob(pwy_samples):
                os.mkdir(par['destination_path']+par['dataset'] + '/' + dada + '/' + par['pathways_f']+'/'+basename(t)[:-4])
                pathways(t,par['destination_path']+par['dataset'] + '/' + dada + '/' + par['pathways_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
            print 'done pathways'

            for t in glob.glob(gn_samples):
                os.mkdir(par['destination_path']+par['dataset']+'/'+ dada + '/'+par['genefamilies_f']+'/'+basename(t)[:-4])
                genefamilies(t,par['destination_path']+par['dataset']+'/'+ dada + '/'+par['genefamilies_f']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
            print 'done genefamilies'


    #if not par['other_datasets']:
    normal()
    #else:
    #    subsampled()

    """
    for t in glob.glob(marka_samples):
        os.mkdir(par['destination_path']+par['dataset']+'/'+par['marker_af']+'/'+basename(t)[:-4])
        marker_abundance(t,par['destination_path']+par['dataset']+'/'+par['marker_af']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
    print 'done marker abundance'
    for t in glob.glob(markp_samples):   
        os.mkdir(par['destination_path']+par['dataset']+'/'+par['marker_pf']+'/'+basename(t)[:-4])
        marker_presence(t,par['destination_path']+par['dataset']+'/'+par['marker_pf']+'/'+basename(t)[:-4]+'/'+basename(t)[:-4])
    print 'done marker presence'    
    """

#mkdir -p ${destination_path}${d}/humann2_GeneFam/${sample}
#for t in ${origin_path}${d}/${or_profiles}/*.tsv;
