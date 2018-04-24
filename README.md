# Scope and Requirements #

## This repo provides:

* Programs for an extensive machine learning meta-analyses (based on metagenomic profile data) including figures. 
* An easy frame-work for alpha & beta diversities and ordination plots.

### Machine Learning Meta-Analysis:

* Consider having 3 datasets named data_a, data_b, data_c. For any, you have provided a metadata table.
* There are a couple of database possibile. For now consider the most basic, Metaphlan2.
* The ** step 1 ** is to compile a list of command lines

```  
python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -do cross_figures -g0 nt:500 -g1 nsl:5 -g2 c:entropy
```

```
to write the lines for the progressive plot:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db metaphlan -al rf -do support_study -g0 nt:500 -g1 nsl:5 -g2 c:entropy
    to plot the progressive training:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do support_plot
    to plot a feature importance heatmap:
        python run.py crc --define study_condition:CRC:control --datasets ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -al rf -do heat_map
    to write down the standard analyses:
	python run.py crc --define study_condition:CRC:control -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 -db genefamilies -al rf -do cross_study -g0 c:gini -g1 nsl:1 -g2 df -hv 0 -ncores 1 -mf auto -r 1


```








### A comprehensive machine learning meta-analyses based on a meta-dataset like the one mentioned above ###

*
* 
* 
* 

* Downloading and curating public 
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact