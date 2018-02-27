#!/user/bin/env/ python

# ******************************************************************************************************
# *                     ***************************************************                            *
# *                       *************    NCBI-DOWNLOADER    ***********                              *                               
# *                 USAGE: python ncbi_downloader.py datasetname bprj_code                             *                               
# *                 other usefull modes: only_download; allready_downloaded (specular, it prevents the *                                                        
# *                 creation of new directories if all the necessary have been selcetd from a previous *
# *                 download.                                                                          *
# ******************************************************************************************************

__author__ = "Francesco Beghini, Paolo Manghi (francesco.beghini@studenti.unitn.it, paolo.manghi@unitn.it)"
__version__ = "0.4"
__date__ = "4 Mar 2017"

import pathos.multiprocessing as mp
import lxml.etree as xmlet
import pandas as pd
import subprocess
import logging
import shutil
import time
import bz2
import os, re
import argparse
import numpy as np
from Bio import Entrez

#**********************************************
Entrez.email = "paolo.manghi@unitn.it" #"paolomanghi1974@gmail.com" ##"paolo.manghi@studio.unibo.it"
#**********************************************

class XML_Parser(object): 
    def __init__(self, bioproject_code, datasetname, strategy, verbose=2):
        self.precomputeds = {"PRJEB8249": 288548, "PRJEB18780": 360524}
        self.metadata = pd.DataFrame()
        try:
            if not str(bioproject_code)[0].isdigit(): 
                wholecode = True
                bioproject_handle = Entrez.efetch(db="bioproject", id=str(bioproject_code), retmode="xml")
                entry = xmlet.parse(bioproject_handle)
                bprj = int(entry.xpath('//RecordSet/DocumentSummary/Project/ProjectID/ArchiveID/@id')[0])
            else:
                wholecode, bprj = False, int(bioproject_code)  
        except:
            if wholecode:
                try: 
                    bprj = int(self.precomputeds[bioproject_code])
                except KeyError:
                    raise IOError("This code %s is not public; "%bioproject_code  +\
                                  "go on the ncbi bioproject database and search "+
                                  " the numerical code in the top left corner.")
                    exit(0)
            else:  
                raise SyntaxError("The code you provided ( %s )failed."%bioproject_code+\
                                  " Are you sure this is a real ncbi"+\
                                  "bioproject code?")
                exit(0)
        self.bprj = bprj
        bioproject_sra_link_handle = Entrez.elink(dbfrom="bioproject", db="sra", id=bprj)
        bioproject_sra_link_result = Entrez.read(bioproject_sra_link_handle)

        try:
            self.SRAs = [link['Id'] for link in bioproject_sra_link_result[0]["LinkSetDb"][0]['Link'] ]
        except IndexError:
            raise IOError("This PRJNA doesn't correspond to any data"+\
                          " apparently, as if is not deposited yet.")
            exit(0)

        self.mappruns, self.mapplibs = self.mapp_sra_runs(strategy,verbose)
        try:
            self.metadata.to_csv(datasetname+"_init_metadata.csv", sep='\t')
        except UnicodeEncodeError: 
            print "Unicode error, passing by"
            self.metadata.reset_index(drop=True) 
            print self.metadata.loc[:15]
	finally: 
            return None

    def mapp_sra_runs(self,Strategy,verbose):
        mapping_run, mapping_layout = {}, {}
        link = Entrez.elink(dbfrom="bioproject", db="sra", id=int(self.bprj), retmode="xml")
        tree = Entrez.read(link)
        if(verbose): print "Processing linked sra"
        wgs, others = 0,0

        for access in self.SRAs: ##[241:]:
            sra_handle = Entrez.efetch(db="sra", id=str(access), retmode="xml") 
            entry = xmlet.parse(sra_handle)
            strategy = entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY')[0].text

            if verbose==3:
                experiment_primary_sra = entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/EXPERIMENT/IDENTIFIERS/PRIMARY_ID')[0].text
            if verbose==2:
                sample_primary_sra = entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/SAMPLE/IDENTIFIERS/PRIMARY_ID')[0].text
            
            if (strategy==Strategy):  ##  "WGS"):
                wgs += 1
                runs = entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN/@accession') #[0]         
                if entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED'): 
                    layout = "Paired"
                else: layout = "Single"
                sample_metadata = self.catchMetadata(entry)  
                Id = entry.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/SAMPLE/@alias')[0]
                ID = Id.replace(" ", "_")
                if ID[0].isdigit(): ID = "SID"+ID
                sample_metadata["Studysampleid"] = ID  
                self.metadata = self.metadata.append(sample_metadata)
                for run in runs: 
                    if not ID in mapping_run: mapping_run[ID] = []
                    mapping_run[ID].append(run)
                    mapping_layout[run] = layout
                if verbose==1: 
                    print "MAPPED\tID:\t%s\t"%ID + "\tTO\t" + " ".join([r for r in mapping_run[ID]]) + " %i " %wgs +  " (%s)" %strategy
		elif verbose==2: 
                    print "MAPPED\tID-SRA:\t%s\t%s\t"%(sample_primary_sra, ID) + "\tTO\t" + " ".join([r for r in mapping_run[ID]]) + " %i " %wgs +  " (%s)" %strategy
                elif verbose==3: 
                    print "MAPPED\tID-SRA:\t%s\t%s\t"%(experiment_primary_sra, ID) + "\tTO\t" + " ".join([r for r in mapping_run[ID]]) + " %i " %wgs +  " (%s)" %strategy
            else:
                others += 1
                if verbose: print "!= FROM DESIRED STRATEGY n. %i (%s)" %(others,strategy)    
        return mapping_run, mapping_layout

    def catchMetadata(self, entrez):
        df = pd.DataFrame()
        df['size'] = entrez.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/RUN_SET/RUN/@size')
        try:
            sample_attributes = entrez.xpath('//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/SAMPLE/SAMPLE_ATTRIBUTES')[0]
            for attribute in sample_attributes: 
                tag = attribute[0].text
                value = attribute[1].text
                df[tag] = value
        #except IndexError: print "missing attributes"
	#except: print "General Error missing Metadata"
        except UnicodeEncodeError: print "Unicode error, passing by"
        finally: return df

#***************************
#***************************

class download(object):
    def __init__(self,asp,vdb,_path,samplename,run,down):
        if down: self.status = self.download_run(_path,samplename,run,asp,vdb)

    def download_run(self,pp,sample,run,asp,vdb):
        if not os.path.exists(pp+sample+"/downloads/"+run+".sra"):
            if not os.path.isdir(pp+sample+"/downloads/"): os.makedirs(pp+sample+"/downloads/")
            result = self.iterate_download_attempt(asp,vdb,run,pp,sample,0,False) 
        else:  #### if the run is present and controll is True, no download should be tented
            result = self.iterate_download_attempt(asp,vdb,run,pp,sample,0,self.check_posterior_md5(vdb,pp+sample+"/downloads/"+run+".sra"))
        if result=="Not Done!": return "Failed"
        else: return "Ok" 

    def call_asp(self,asp,run,path_):
        cmd = "%s/bin/ascp -T -i %s/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra %s"\
              %(asp,asp,run[0:3],run[0:6],run,run,path_)
        try:
            Cmd = subprocess.check_call(cmd.split())
            return True
        except subprocess.CalledProcessError: 
            return False

    def iterate_download_attempt(self,asp,vdb,run,pp,samplename,att=0,controll=False):
        if controll: 
            print "Controll OK"
            return "Done"
        while ((not controll) and (att<5)):
            proc = self.call_asp(asp,run,pp+samplename+"/downloads/")
            if proc: controll = self.check_posterior_md5(vdb,pp+samplename+"/downloads/"+run+".sra") 
            if not controll: time.sleep(20)
            else:
                print "Controll OK"
                return "Done"
            att += 1
        return "Not Done!"    
    
    def check_posterior_md5(self,vdb,path):
        validation = subprocess.Popen([vdb,path],shell=False,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        out = validation.stdout.read()
        md5s = re.findall(r"md5", out)
        oks = re.findall(r"md5\sok", out)
        if (len(oks)==len(md5s)): 
            return True
        else:
            print out,"\nfrom the func check_posterior, the wrong check\n"
            return False
        #return True if(len(oks)==len(md5s))else False

#**** This class is used to execute a download (single)
#******************************************************

class LogFile(object):
    def __init__(self, dictionaryidruns, layouts, valid_p, asp_p, fqd_p, totpaired, totsingle, readspath):
        self.idtoruns = dictionaryidruns ## from Download
        self.layouts = layouts ## from Download
        self.valid_p = valid_p ## vdb-validate
        self.asp_p = asp_p ## aspera
        self.fqd_p = fqd_p ## fastq dump
        self.totsingle = totsingle ## from Download
        self.totpaired = totpaired ## from Download
        self.readspath = readspath ## ''
        
    def write_download(self,path_,run):  ### contains check md5 and 5 attempts in bash
        cmd = "%s/bin/ascp -T -i %s/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra %s\n"\
              %(self.asp_p,self.asp_p,run[0:3],run[0:6],run,run,path_)
        whileblock0 = "for i in 1 2 3 4 5\ndo\n"+cmd+"\n"
        whileblock1 = "%s %s%s.sra 2> %s.validation \nif grep -q \"is consistent\" %s.validation; " %(self.valid_p,path_,run,run,run)
        whileblock2 = "then\n\tbreak\nfi\ndone\n" 
        return whileblock0+whileblock1+whileblock2

    def write_dump(self,path_,run,layout):
        if layout=="Single":
            return "%s %s%s\".sra\"\nmv %s\".fastq\" %s" %(self.fqd_p,path_,run,run,path_)    
        return "\n%s -I --split-files %s%s.sra\nmv%s\"_1.fastq\" %s%s\"_R1.fastq\"\nmv%s\"_2.fastq\" %s%s\"_R2.fastq\"\n" \
        %(self.fqd_p,path_,run,run,path_,run,run,path_,run)    
    
    def write_bzip(self,ID,_path):
        res = ""
        if self.totpaired: 
            for R in ["_R1", "_R2"]: 
                for run in self.idtoruns[ID]: res += "\ncat %s >> %s\n" %(_path+run+R+".fastq",_path+ID+R+".fastq\n")     
        if self.totsingle: 
            for run in self.idtoruns[ID]: res += "\ncat %s >> %s\n" %(_path+run+"_UN.fastq", _path+ID+"_UN.fastq\n")
        res += "\nif [ -f %s%s\"_R1.fastq\" ] ; then \n bzip2 %s%s\"_R1.fastq\"\n mv %s%s\"_R1.fastq.bz2\" %s../\nfi\n" \
            %(_path,ID,_path,ID,_path,ID,_path) 
        res += "\nif [ -f %s%s\"_R2.fastq\" ] ; then \n bzip2 %s%s\"_R2.fastq\"\n mv %s%s\"_R2.fastq.bz2\" %s../\nfi\n" \
            %(_path,ID,_path,ID,_path,ID,_path) 
        res += "\nif [ -f %s%s\"_UN.fastq\" ] ; then \n bzip2 %s%s\"_UN.fastq\"\n mv %s%s\"_UN.fastq.bz2\" %s../\nfi\n" \
            %(_path,ID,_path,ID,_path,ID,_path) 
        return res

#******************************************************

class Download(object):
    def __init__(self,bprj,datasetname,aspera,fastqdump,where,vdb\
                ,ncores,dont_keep_downloads,strategy,not_to_bzip\
                ,dofastA,xml_verbosity,separate,mode="all",samples_to_download=None,split_3=True):

        #************* decisions for the logging file ****************************
        self.mode = mode
	self.xml_verbosity = xml_verbosity # 0 1 2
	self.split_3 = split_3
        self.dofasta = dofastA
        self.separated = separate

        if self.mode=="all" and (samples_to_download!=None): self.mode = "onebyone"
        if samples_to_download!=None: to_download = samples_to_download

        if (self.mode not in ["xml","only_download","controll","onebyone"] and os.path.exists("%s_DOWNLOAD.sh"%datasetname)):
            os.remove("%s_DOWNLOAD.sh"%datasetname) 
        if self.mode in ["all","ghost","allready_downloaded"]:
            logging.basicConfig(filename="%s_DOWNLOAD.sh"%datasetname,level=logging.INFO) 
        self.failed=[]

        #**********  paths ***
        self.strategy = strategy
        self.aspera_p = aspera
        self.vdb_p = vdb
        self.fastqdump_p = fastqdump
        self.downloadfolder = where
        self.datasetname = datasetname
        self.readspath = self.downloadfolder+self.datasetname+"/reads/" 

        if self.mode in ["ghost","controll"]:
            try:
                with open(self.datasetname+"samplelist.txt") as sams:
                    samlist = [row.rstrip() for row in sams.readlines()]
            except IOError: 
                print "A list of the samples which download has to be written\
                is missing, therefore the list of all the samples will be used."
                samlist = False

        #******************* functions for dumping: **************************
        self.dumppaired = lambda fqd,p,r : "%s --split-3 %s%s.sra" %(fqd,p,r)

	#if not self.split_3: ## changed
	#    self.dumppaired = lambda fqd,p,r : "%s --split-files %s%s.sra" %(fqd,p,r)

	self.dumpto_fasta = lambda fqd,p,r : "%s --fasta 0 %s%s.sra" %(fqd,p,r)
        self.dumpsingle = lambda fqd,p,r : "%s %s%s.sra" %(fqd,p,r) ### p is path, r run

        self.mvpaired1 = lambda p,r : ["%s_1.fastq"%r, "%s%s_R1.fastq" %(p,r)]
        self.mvpaired2 = lambda p,r : ["%s_2.fastq"%r, "%s%s_R2.fastq" %(p,r)]
        self.mvsingle = lambda p,r : ["%s.fastq"%r, "%s%s_UN.fastq" %(p,r)]

        #*************************************************************************
        if not os.path.isdir(self.readspath): os.makedirs(os.path.dirname(self.readspath))
        if not self.mode in ["xml","only_download","controll"]: 
            logging.info("\nmkdir %s\n"%self.datasetname+"\nmkdir %s"%self.datasetname+"/reads/\n")

        #****************  xml parsing information (sample id & runs) **************
        xml = XML_Parser(bprj, datasetname, self.strategy, self.xml_verbosity)
        self.idtoruns = xml.mappruns
        self.librarylayouts = xml.mapplibs  ### grasp things
        self.ncores = ncores
        self.dont_keep_downloads = dont_keep_downloads  ### NEW
        self.not_to_bzip = not_to_bzip

        if self.mode=="xml": exit(0) 
        #*************** several specifications of the logging of a downlaoder in .sh
        self.Totpaired_reads, self.Totsingle_reads, self.mean = self.stupid_sts()

        if self.mode in ["all", "allready_downloaded", "ghost", "controll"]: 
            logging.info("\n#Tot. samples : %i\n"%(len(self.idtoruns.keys())))
            logging.info("\n#Tot. paired end reads : %i\n#Tot. single end reads : %i\n"%(self.Totpaired_reads,self.Totsingle_reads)+\
                             "\n#Avg. num of runs per samples : %.1f\n"%(float(np.sum(self.mean)/self.mean.shape[0])))

        if self.mode in ["all", "ghost", "allready_downloaded", "controll"]: 
            self.LOG = LogFile(self.idtoruns, self.librarylayouts, self.vdb_p,\
            self.aspera_p, self.fastqdump_p, self.Totpaired_reads, self.Totsingle_reads, self.readspath)

        else: self.LOG = False  ## xml, only_download

        if self.mode=="ghost":
            samples = samlist if samlist else list(self.idtoruns.keys())
            self.log_all(samples,LOG)
            exit(0)

        elif self.mode=="controll":
            samples = samlist if samlist else list(self.idtoruns.keys()) 
            self.exec_only_controll(samples)
            exit(0)

        elif self.mode in ["onebyone", "onlydownload_onebyone", "onebyone_onlydownload"]: ### to specfically download a single sample
             for sample in to_download:
                 try:
                     for run in self.idtoruns[sample]:
                         run_d = download(self.aspera_p, self.vdb_p, self.readspath, sample, run, True).status
                         if (run_d=="Failed"): print "run " + str(run_d) + " on sample " + str(sample) + " failed "
                 except KeyError:
                         print "There was a Key Error in the dictionary: Key you could use are : " + " ".join(list(self.idtoruns.keys()))

                 if (not self.mode in ["onlydownload_onebyone", "onebyone_onlydownload"]):
                     self.Dump_Sample(sample,self.readspath+sample+"/downloads/")

                 #if (not self.not_to_bzip): self.Bzp_Sample(sample,self.readspath+sample+"/downloads/")   
             exit(0)

        ## main    (only_download, all, allready_download modes)
        else: 
            self.go_with_download()
        
    def stupid_sts(self):
        p,s,c = (0,0,0)
        mean = np.zeros(len(self.idtoruns.keys()), np.float64)
        for ID in self.idtoruns:
            mean[c] += len(self.idtoruns[ID])
            c += 1
            for run in self.idtoruns[ID]:
                if self.librarylayouts[run]=="Paired": p += 1
                else: s += 1
        return p,s,mean

    def log_all(self,samplelist,LOG):
        for ID in samplelist:
            _path = self.readspath+ID+"/downloads/" 
            logging.info("\nmkdir %s\n" %self.readspath+ID)  ###  alldownloads and all dumps is easier to parallelize 
            logging.info("\nmkdir %s\n" %_path)
            for run in self.idtoruns[ID]: logging.info(LOG.write_download(_path,run)) ###  in the future
            for run in self.idtoruns[ID]: logging.info(LOG.write_dump(_path,run,LOG.layouts[run]))
            logging.info(LOG.write_bzip(ID,_path))

    def listFAILURES(self):
        if len(self.failed):
            fails = open(self.datasetname+"failed.txt", "w")
            for E in self.failed:fails.write(E[0]+"\t"+E[1]+"\t"+E[2]+"\n")
            fails.close()

    def exec_only_controll(self,samples):
        print "****** CONTROLL ******"
        atleastone = False
        for sample in samples:
            if not os.path.isdir(self.readspath+sample+"/"):
                print "Sample %s DIR not present." %self.readspath+sample+"/downloads/"
            else: 
                atleastone = True
                for run in self.idtoruns[sample]:
                    run_d = download(self.aspera_p,self.vdb_p,self.readspath,sample,run,False).\
                        check_posterior_md5(self.vdb_p,self.readspath+sample+"/downloads/"+run+".sra") 
                    if run_d: print "Sample %s CONTROLLED IS OK (check=%s)" %(sample,str(run_d))
                    else: print "Sample %s CHECK IS %s, CONTROLL GOT AN ERROR" %(sample,str(run_d)) 
        if not atleastone: 
            print "None of the downloads directories was present, still. \
                They have probably been eliminated during download\
                . If you want to keep the /downloads/ for another download, \
                in the future use the flag '-kd' (--keep_downloads).."

    def exec_download_of_any_run(self,samplename):  ### create a download object with the whole bullshit
        if os.path.exists(self.datasetname+"failed.txt"): os.remove(self.datasetname+"failed.txt")
        for run in self.idtoruns[samplename]:
            run_d = download(self.aspera_p,self.vdb_p,self.readspath,samplename,run,True).status 
            if (run_d=="Failed"): 
                print "run_d results being : %s indeed we append it." %str(run_d)
                self.failed.append((samplename,run,"download")) ### doesn't happen ever

    def handle_physical_download(self):
        n = 0
        all_the_samples = list(self.idtoruns.keys())
        if not self.mode in ["allready_downloaded"]:  ## onlydownload, all
            for ID in list(all_the_samples):  
                if not os.path.isdir(self.readspath+ID+"/"):
                    os.makedirs(self.readspath+ID+"/") 
                self.exec_download_of_any_run(ID)
                n += 1
        else:   ### allready downloaded:  delete from dict useless samples
            for anymissingsample in all_the_samples: ## (otherwise they mess around after)
                if not os.path.isdir(self.readspath+anymissingsample+"/"): 
                    del self.idtoruns[str(anymissingsample)]
        return n

    def go_with_download(self):
        count=self.handle_physical_download()
        if self.mode in ["only_download"]:
            print "Download completed with success (%i samples.)" %count 
            exit(0)

        all_the_samples = list(self.idtoruns.keys())
        ###print 'mi aspetto che ci siano 74 campioni: ', len(all_the_samples) 

        def getfree(ID):
            try:
                _path = self.readspath+ID+"/downloads/"
                if (not self.dofasta): 
                    self.Dump_Sample(ID,_path)

                else:
                    self.Dump_Fasta(ID,_path)
                    subprocess.call(['ls'])           

                if (not self.not_to_bzip): 
                    if not self.separated:
                        self.Bzp_Sample(ID, _path)
                    else:
                        self.Bzp_Sample_Separated(ID, _path)
                return 
            except NameError:
                print 'Wish to known which sample made such mess?? Here: %s' %ID
                return 
            except IOError:
                print 'A sample made a different mess: the RUN CODE was not detected.'
                print 'Don ask me why, not this one, but the other 6000.'
                return

        
        with mp.ProcessingPool(self.ncores) as pp:
            pp.map(getfree, all_the_samples)

        self.listFAILURES()
        self.log_all(all_the_samples, self.LOG)


    def Dump_Fasta(self,ID,_path):
        for r in self.idtoruns[ID]:
            try: subprocess.check_call(self.dumpto_fasta(self.fastqdump_p,_path,r).split())
            except subprocess.CalledProcessError: print 'DIDN\'T manage to download', r 
            shutil.move("%s.fasta" %r, "%s%s.fasta" %(_path,r))    


    def Dump_Sample(self, ID, _path):
        for run in self.idtoruns[ID]:
            if (self.librarylayouts[run] == "Paired"):
                try: subprocess.check_call(self.dumppaired(self.fastqdump_p, _path, run).split())
                except subprocess.CalledProcessError: self.failed.append((samplename, run, "dump"))

                shutil.move(self.mvpaired1(_path, run)[0], self.mvpaired1(_path, run)[1])
                shutil.move(self.mvpaired2(_path, run)[0], self.mvpaired2(_path, run)[1])

            else: ### CHANGED again to "dumpsingle"
                try: subprocess.check_call(self.dumpsingle(self.fastqdump_p, _path, run).split())    
                except subprocess.CalledProcessError: self.failed.append((samplename, run, "dump")) 

                shutil.move(self.mvsingle(_path, run)[0], self.mvsingle(_path, run)[1])


    def Bzp_Sample(self, ID, _path):
        if (self.Totpaired_reads): 
            for R in ["_R1", "_R2"]:
                with bz2.BZ2File(_path+ID+R+".fastq.bz2","w") as outfile:    
                    for run in self.idtoruns[ID]:
                        if (self.librarylayouts[run]=="Paired"):
                            with open(_path+run+R+".fastq") as infile:
                                for line in infile: outfile.write(line)
        if (self.Totsingle_reads):
            with bz2.BZ2File(_path+ID+"_UN.fastq.bz2","w") as outfile:
                for run in self.idtoruns[ID]:
                    if (self.librarylayouts[run]=="Single"):
                        with open(_path+run+"_UN.fastq") as infile:
                            for line in infile: outfile.write(line)
        self.moving_and_removing(ID,_path)


    def moving_and_removing(self, ID, _p):
        for ext in ["_R1.fastq.bz2","_R2.fastq.bz2","_UN.fastq.bz2"]:
            if ((os.path.exists(_p+ID+ext)) and (os.path.getsize(_p+ID+ext)>0)):    
                shutil.move(_p+ID+ext,_p+"../"+ID+ext)
        if self.dont_keep_downloads: # subprocess.call(["tar","-zcvf",_p+"../"+"download_archive.tar.gz",_p])
            shutil.rmtree(_p)


    def Bzp_Sample_Separated(self, ID, _path):
        if (self.Totpaired_reads):
            for R in ["_R1", "_R2"]:
                for enum,run in enumerate(self.idtoruns[ID]):
                    with bz2.BZ2File(_path+ID+"_"+str(enum)+"_"+R+".fastq.bz2", "w") as outfile:
                        if (self.librarylayouts[run] == "Paired"): ## sure?
                            with open(_path+run+R+".fastq") as infile:
                                for line in infile: outfile.write(line)
        if (self.Totsingle_reads):
            for enum,run in enumerate(self.idtoruns[ID]):
                with bz2.BZ2File(_path+ID+"_"+str(enum)+"_"+"_UN.fastq.bz2", "w") as outfile:
                    if (self.librarylayouts[run] == "Single"): ##
                        with open(_path+run+"_UN.fastq") as infile:
                            for line in infile: outfile.write(line)

    #def moving_and_removing_separated(self, ID, _path)
    #    import glob
    #    for ext in ["_R1.fastq.bz2","_R2.fastq.bz2","_UN.fastq.bz2"]:
    #        for enum in glob.glob(_p+ID+"*"+ext)

#**********************************************************************************
#************************ following functions are for handling the general settings

def optionsfromfile(Op):
    try: f = open(Op)
    except IOError:  
        raise SyntaxError(\
        "Ehi man, it's ok, but you want to specify a file but that file is either missing or unspecified.\n\
         File is %s; does it look like a dataset name? This function is called when no PRJNA is present \n\
         and so args.brpj=NULL. In case you want to pass a dataset you have also to specifify a code." %Op)
        exit(1)
    list_of_datasets_to_download = f.readlines()
    list_of = [tuple(line.rstrip().split()) for line in list_of_datasets_to_download]
    f.close()
    return list_of     
 
def main(bprj,dataset,asp,fastqdump,downdir,vdb,ncores,dontkeepdownloads,strategy,dontbzip,fasta,xml_verb,separate,modality,stodownload,split_3):
    if not bprj=="NULL": 
        D=Download(bprj,dataset,asp,fastqdump,downdir,vdb,ncores,dontkeepdownloads,strategy,dontbzip,fasta,xml_verb,separate,modality,stodownload,split_3)
    else: 
        OPs = optionsfromfile(dataset) ### here, dataset is the name of the listfile
        #if stodownload!=None: stodownload=None
        #print 'A list of smples to download was present but even a list of datasets to download was, these options are not compatible.'
        for nm,cd in OPs: D=Download(cd,nm,asp,fastqdump,downdir,vdb,ncores,dontkeepdownloads,strategy,dontbzip,fasta,xml_verb,separate,modality,stodownload,split_3) 

def read_arguments():         
    parser = argparse.ArgumentParser(prog="NCBI WGS dataset downloader", description= " Provide a bprj code and download a dataset."
        "\nUSAGE EXAMPLE: python ncbi_downloader.py ZellerG_2014 266076 \n or provide an option with multiple lines of names and codes")
    parser.add_argument("dataset", type=str, metavar="download singularly or acc to a listfile.",help=" The name of the dataset or of the listfile. A listfile line must be: dataset PRBJNA"
        " FirstAuthSurnameN_YEAR: Castro-NallarE_2015 Either you can pass the name of a file containing a list of dataset and codes ")
    parser.add_argument("bprj", metavar="bioprojectcode",type=str, default="NULL", help="sdt project/dataset code on ncbi: PRJNA255439\n  bioproject code is ass a 6-digits numerical code \n"
        " (eq to the last 6 or not).\n This numerical code (e.g. 266076) can be passed too. Default is NULL for dataset listfile.")
    parser.add_argument("-asp", "--aspera_p",default="/CM/tools/aspera-3.5.6/connect", help="set the place in which to search for aspera")
    parser.add_argument("-fqp", "--fastqdump_p",default="/CM/tools/sratoolkit-2.3.4/bin/fastq-dump.2.3.4", help="set the place in which to search for fastq-dump")
    parser.add_argument("-vdb", "--vdb_validate",default="/CM/tools/sratoolkit-2.3.4/bin/vdb-validate.2.3.4", help="set the path of vdb-validate")
    parser.add_argument("-where", "--download_dir",default="/scratchCM/repos/datasets_metagenomics/",help="the place where the dataset is downloaded")
    parser.add_argument("-dkd", "--dont_keep_downloads",action="store_true",help="when used as flag, it TARs the /downloads/ dir of any sample. By default, /downloads/ is removed.")
    parser.add_argument("-ncores", "--core_number",type=int, default=4, help="the number of cores used in the parallel processing; default is 4")
    parser.add_argument("-st", "--strategy", default="WGS",help="This script is though to download WGS-strategy data\
        but changing the strategy can work even for other types of sequencing (e.g.: AMPLICON)")
    parser.add_argument("-fasta",action="store_true",help="Usable for dumping to get fastA, still trial." )
    parser.add_argument("-ntbz", "--not_to_bzip", action="store_true",help="If specified, samples are not bzipped.")
    parser.add_argument("-nsp3", "--non_split_3", action="store_false")
    parser.add_argument("-m", "--mode",choices=["all","only_download","jd","xml","bash","ghost","controll","allready_downloaded","onebyone","one","sample","onlydownload_onebyone","onebyone_onlydownload"]\
        ,default="all",help=" You can use these argument in case you need to handle a complicate download.") 
    parser.add_argument("-xv","--xml_verbosity",default=2,type=int,choices=[0,1,2,3])
    parser.add_argument("-std","--samples_to_download",default=None,nargs="+",help="Speficied, you download only the samples argumented and no others, independently from modality.")
    parser.add_argument("-s","--separate",action='store_true')

    args = parser.parse_args()    
    return vars(args)

params = read_arguments()
mode = params['mode']
aliases = {"jd":"only_download","bash":"ghost","ghost":"ghost"\
          ,"sample":"onebyone","one":"onebyone"}

if mode in aliases: mode=aliases[mode] 

if (__name__=='__main__'): 
    #print args.bprj, "  lo stai vedendo?  "
    #exit(0)
    main( params['bprj']\
         ,params['dataset']\
         ,params['aspera_p']\
         ,params['fastqdump_p']\
         ,params['download_dir']\
         ,params['vdb_validate']\
         ,params['core_number']\
         ,params['dont_keep_downloads']\
         ,params['strategy']\
         ,params['not_to_bzip']\
         ,params['fasta']\
	 ,params['xml_verbosity']\
         ,params['separate']\
         ,mode\
         ,params['samples_to_download']\
	 ,params['non_split_3'])
