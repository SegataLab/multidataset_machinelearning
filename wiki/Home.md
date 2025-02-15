## Who to contact? 
since this tutorial is still in progress, please for any needing feel free to get in touch with me at paolomanghi1974@gmail.com.
I really appreciate it, thanks for the patience.

## Where you are?
This repo provides a series of script for Data-Driven Whole-Genome-Shotgun (WGS) Metagenomics, in particular, it allows the generation of figures relative to the Machine-Learning analysis in the paper: "Metagenomic analysis of colorectal cancer datasets identifies cross-cohort microbial diagnostic signatures and a link with choline degradation."

## What does this repo actually contain?
* 1 program allowing the download of a dataset from NCBI.
* 2 programs for speed up the boring procedure of writing in Bash a whole Machine-Learning (ML) meta-analyses based on several WGS datasets (all or part of them, maybe, just downloaded with the above named program).  
* 3 programs producing figures relative to the different aspect of the ML analyses conducted with the above program (the indeed become quite detailed).
* 4 programs devoted to the basic part of a standard WGS analyses, including the generation of figures such as: any type of ordination plot (literally), beta-diversity box-plots, alpha diversity boxplots.

## Getting Started step n.1
To clone the present repo:
```
hg clone https://<USER>@bitbucket.org/CibioCM/multidataset_machinelearning
```
The analysis listed here are mainly dovoted to generate the figures of the paper. Can be applied in a similar fashions to other series of
datasets/other classification problems.

## NOTE: they are THEY ARE BASED ON TWO OTHER PUBLIC REPOS: cmdpy and metaml. You need to clone them in order to run the analysis!
To clone cmdpy and metaml repos:
```
hg clone https://<USER>@bitbucket.org/CibioCM/metaml
hg clone https://<USER>@bitbucket.org/CibioCM/cmdpy

```
- cmdpy repo contains the tool cmdpy_dataset.py which is used to merge quantitative profiles with metadata
- metaml repo contains the tool classification_thomas-manghi.py which is needed to run the main analysis

## Getting Started step n.2

#### NOTE: these figures have been create from a server in which I'm also maintaining the curatedMetagenomicDataset package!
#### unfortunately, that means that is conceived as an architecture with the following tree-shape:

- /base_folder
- /base_folder/dataset_name(e.g. FengQ_2015)
- /base_folder/dataset_name/dataset_name_metadata.tsv
- /base_folder/dataset_name/metaphlan2
- /base_folder/dataset_name/metaphlan2/sample_#1_name
- /base_folder/dataset_name/metaphlan2/sample_#2_name
- /base_folder/dataset_name/metaphlan2/sample_#1_name/sample_#1_name_profile.tsv
- /base_folder/dataset_name/metaphlan2/sample_#1_name/sample_#2_name_profile.tsv

That means that the analysis you would carry out is based a on place in which you keep your datasets,
each as a separate folder. In each dataset folder you must have the computational profiles and the metadata
for the the datasets subdivided as listed above.
Luckily, it is enough for the base-folder to be a folder in your home for which you know the location and 
and have permission, nothing more!

## Getting Started step n.3
##### once the architecture is set (that means also having dataset-based folders, metadata table and metaphlan2 profiles), please note that the essential script refers to such architecture, indeed do:

```
vi ../cmdpy/cmdpy_dataset.py
```

you'll see the following lines:

```
- NOTE: for very large databases (genefamilies, markers)
- there's an option to perform Wilcoxon - rank 
- in order to srhink features.
"""

BASE_PATH='/CM/data/meta/'

class metadataset:

    def __init__(self):
```

manually change BASE_PATH from "/CM/data/meta/" to "/YOUR/SERVER/BASE-ARCHITECTURE/" considering that this must be the place in which you keep your different datasets (SEE above)
then exit 

```
press ESC, then type:
:wq
```

Now do:

```
vi run.py

```

you should be able to see the following lines:

```
from utils import usefullfuncs ## utilities
import argparse as ap
import sys

BASE_PYTHON='/shares/CIBIO-Storage/CM/scratch/users/paolo.manghi/anaconda3/bin/'

def read_params(args):
```

please, change BASE_PYTHON from "/shares/CIBIO-Storage/CM/scratch/users/paolo.manghi/anaconda3/bin/" to "/YOUR/FOVOURITE/PYTHON.3/bin/"
and exit

```
press ESC, then type:
:wq
```

Now you should be able to run the analysis. If not, please refer to paolomanghi1974@gmail.com.

## Prepare the analysis step n.1
##### the first step is to proceed to the dataset generation (ML == starting datasets) THEREFORE, in order to generate the datasets: do

```
python run.py crc --define study_condition:CRC:control \
     -ds ZellerG_2014 YuJ_2015 FengQ_2015 VogtmannE_2016 CM_rescignocrc CM_lilt HanniganGD_2017 \
	-db metaphlan \
	-do cross_study \
	-g0 nt:500 -g1 nsl:5 -g2 c:entropy
```
	
if run without errors, you should be able to see in the folder you run it to files.sh:

```
ls	

datasetcaller_metaphlan.sh	transfer_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
lodo_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
```

## Running the analysis:
In order to run the analysis make the three wrapper executable:

```
chmod +x datasetcaller_metaphlan.sh
chmod +x transfer_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
chmod +x lodo_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
```

and execute them:

```
./datasetcaller_metaphlan.sh
./transfer_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
./lodo_metaphlan_study_condition:CRC:control_rf_gridterm0:nt_500_gridterm1:nsl_5_gridterm2:c_entropy.sh
```

Please note that these analysis may take more than some minute, depending on the computational sources adopted.
For errors, please refer to paolomanghi1974@gmail.com

##### these .sh are wrappers for all the analysis needed
##### IF YOU INSTALLED THE CMDPY AND METAML REPOS YOU SHOULD BE ABLE TO RUN THESE COMMANDS IN THE ORDER LISTED

### Once the analysis have run, in order to generate the figures you might run:

* figure 2 panel a:
** This analysis consist of a cross prediction matrix consisting of 
** n. datasets * n. datasets tests + n. datasets Leave-One-Dataset-Out predictions

```
python run.py crc --define study_condition:CRC:control \
	-ds FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 \
		-db metaphlan \
		-do cross_figures \
		-g0 c:entropy -g1 nt:500 -g2 nsl:5 
		-cm hot 
		--path Fig_Cross_Prediction/	
```

If you manage to run this command without errors, you should now have generated a figure named

```
MachineLearning_crc_metaphlan_all_features.png
```
which, in the present case should look like:

![MachineLearning_crc_metaphlan_all_features2%20copy.jpg]
(https://bitbucket.org/repo/p4MR7K5/images/3967973586-MachineLearning_crc_metaphlan_all_features2%20copy.jpg)

Now, before you can proceed, you need a complete ranking file in order to plot Figure3.
I order to obtain a complete ranking file you need you need to perform 3 steps:

* duplicate one of the dataset of yoru analysis and change its name (or use another one): it is necessary to have a differently named dataset with a least one sample (metaphlan2 + metadata). Suppose to call it FAKE_DATASET
** run :

```
python ../cmdpy_dataset.py FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 FAKE_DATASET -of dataset_for_global_ranking.csv
python ../metaml/classification_thomas-manghi.py dataset_for_global_ranking.csv dataset_for_global_ranking -d 1:study_condition:CRC -t dataset_name:FAKE_DATASET -l rf -nt 500 -nsl 5 -c entropy -nc 10
```

The above command substantially add a fake-dataset to the meta-cohort, and then perform a feature selection ML algorithm on the meta-cohort adopting the fake-dataset as a test. In this way, the featre-ranking of the overall meta-cohort is obtained.
Now, before proceeding to the analysis, do:

```
vi feat_heat.py 
```
You should see the following lines:
```
BASE_PATH='/CM/data/meta/'
COMPLETE_RANKING_PATH='../ml/resultofcrossvalidation_ANYonANY_features:metaphlan_experimenttype:study_condition:CRC:control_rf_grid0:%s_grid1:%s_grid2:%s_FEATURE_RANKING.txt'


class feature_heatmap(object):

    dataselect = '/python ../../cmdpy/cmdpy_dataset.py'

```
please, change BASE_PYTHON from "/CM/data/meta/" to "/YOUR/SERVER/BASE-ARCHITECTURE/", then
change "COMPLETE_RANKING_PATH" to "YOUR_COMPLETE_RANKING_FILE_NAME" (your ranking file name is the one just obtained, in the previous case you be named "global_ranking.txt").
then
```
press ESC, then type:
:wq
```

* Figure 3 panel a:
** This analysis consist in a feature-ranking of 7 cross-validation
*** Therefore the step n. 1 is to perform these cross-validation
	
```
python run.py crc --define study_condition:CRC:control \
	--datasets \
		FengQ_2015 ZellerG_2014 CM_rescignocrc YuJ_2015 CM_lilt VogtmannE_2016 HanniganGD_2017 \
		-al rf \
		-do heat_map 
		-db metaphlan 
		-g0 c:entropy -g1 nt:1000 -g2 nsl:5 -nif 5 
		--path Feat_Rank_Heatmap/
```

If you manage to run this analysis without errors, you should be able to see file named:
```
FeatureHeatmap_warm.png
```

which, in the present case should look like:

![FeatureHeatmap_warm2.png]
(https://bitbucket.org/repo/p4MR7K5/images/1280546409-FeatureHeatmap_warm2.png)

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