# Scope and Requirements #
 
* multidataset_machinlearning is a collection of Python programs aimed to speed up the most essential steps in a WGS metagenomic analyses.
* it relies upone several commonly used Python libraries: among these numpy, scipy, pandas, scikit-learn, matplotlib, seaborn & biopython are not part of the original native libraries.
* it is not meant to be used as a standalone: the idea is that of having your own dataset profiled through some kind of either taxonomical functional tool devoted to profiling microbial communities. The present package has been designed to ease the work on datasest profiled though Metaphlan2 and Humann2 tools, but other tools can be employed. (see paragraphs below).

## There are essentially 3 goals pursued by the script in the collection: ##
### Downloading and Curating public Whole-Genome-Shotgun cohorts ###

* ncbi_downloader.py & ncbi_downloader_dev.py are the programs for executing the physical download of public datasets from NCBI: these script can be applid to different data types, such as AMPLICON sequencing and Viromes.

```
usage: NCBI WGS dataset downloader [-h] [-asp ASPERA_P] [-fqp FASTQDUMP_P]
                                   [-vdb VDB_VALIDATE] [-where DOWNLOAD_DIR]
                                   [-dkd] [-ncores CORE_NUMBER] [-st STRATEGY]
                                   [-fasta] [-ntbz] [-nsp3]
                                   [-m {all,only_download,jd,xml,bash,ghost,controll,allready_downloaded,onebyone,one,sample,onlydownload_onebyone,onebyone_onlydownload}]
                                   [-xv {0,1,2,3}]
                                   [-std SAMPLES_TO_DOWNLOAD [SAMPLES_TO_DOWNLOAD ...]]
                                   [-s]
                                   download singularly or acc to a listfile.
                                   bioprojectcode

Provide a bprj code and download a dataset. USAGE EXAMPLE: python
ncbi_downloader.py ZellerG_2014 266076 or provide an option with multiple
lines of names and codes

positional arguments:
  download singularly or acc to a listfile.
                        The name of the dataset or of the listfile. A listfile
                        line must be: dataset PRBJNA FirstAuthSurnameN_YEAR:
                        Castro-NallarE_2015 Either you can pass the name of a
                        file containing a list of dataset and codes
  bioprojectcode        sdt project/dataset code on ncbi: PRJNA255439
                        bioproject code is ass a 6-digits numerical code (eq
                        to the last 6 or not). This numerical code (e.g.
                        266076) can be passed too. Default is NULL for dataset
                        listfile.

optional arguments:
  -h, --help            show this help message and exit
  -asp ASPERA_P, --aspera_p ASPERA_P
                        set the place in which to search for aspera
  -fqp FASTQDUMP_P, --fastqdump_p FASTQDUMP_P
                        set the place in which to search for fastq-dump
  -vdb VDB_VALIDATE, --vdb_validate VDB_VALIDATE
                        set the path of vdb-validate
  -where DOWNLOAD_DIR, --download_dir DOWNLOAD_DIR
                        the place where the dataset is downloaded
  -dkd, --dont_keep_downloads
                        when used as flag, it TARs the /downloads/ dir of any
                        sample. By default, /downloads/ is removed.
  -ncores CORE_NUMBER, --core_number CORE_NUMBER
                        the number of cores used in the parallel processing;
                        default is 4
  -st STRATEGY, --strategy STRATEGY
                        This script is though to download WGS-strategy data
                        but changing the strategy can work even for other
                        types of sequencing (e.g.: AMPLICON)
  -fasta                Usable for dumping to get fastA, still trial.
  -ntbz, --not_to_bzip  If specified, samples are not bzipped.
  -nsp3, --non_split_3
  -m {all,only_download,jd,xml,bash,ghost,controll,allready_downloaded,onebyone,one,sample,onlydownload_onebyone,onebyone_onlydownload}, --mode {all,only_download,jd,xml,bash,ghost,controll,allready_downloaded,onebyone,one,sample,onlydownload_onebyone,onebyone_onlydownload}
                        You can use these argument in case you need to handle
                        a complicate download.
  -xv {0,1,2,3}, --xml_verbosity {0,1,2,3}
  -std SAMPLES_TO_DOWNLOAD [SAMPLES_TO_DOWNLOAD ...], --samples_to_download SAMPLES_TO_DOWNLOAD [SAMPLES_TO_DOWNLOAD ...]
                        Speficied, you download only the samples argumented
                        and no others, independently from modality.
  -s, --separate

```

* cure.py is a script containing many utility functionalities aimed to manually curate clinical metadata. Examples of the utility of this curation are the metadata table of the CuratedMetagenomicDataset R package.

##### Flag -m of ncbi_downloader_dev.py can be used to reduce the all pipe to part of it, including the sole exploration of the current state of the dataset on NCBI (-m xml) #####

### Performing explorative analyses on a meta-cohort or a meta-dataset ###

* cmdpy_dataset is a tool for collecting a theoretically unlimited number of datasets: it assumes that for any dataset involved in the meta-analyses, community profiles are available, so as corresponding metadata tables. The profiles for which the tool has been conceived were Metaphlan2 and Humann2 profiles, along the lines of the curatedMetagenomicDataset R package. It has been then applied to other types of profiling. 

```
usage: cmdpy_dataset.py [-h] [--base_path BASE_PATH]
                        [--metadata_path METADATA_PATH]
                        [--metadata_name METADATA_NAME]
                        [--metaphlan_path METAPHLAN_PATH]
                        [--change_profile_func] [--metaphlan] [--pwyrelab]
                        [--genefam] [--markab] [--markpres] [--pwycov]
                        [--pfam] [--mirna]
                        [--metaphlan_folder METAPHLAN_FOLDER]
                        [--pwyrelab_folder PWYRELAB_FOLDER]
                        [--genefam_folder GENEFAM_FOLDER]
                        [--pwycoverage_folder PWYCOVERAGE_FOLDER]
                        [--markab_folder MARKAB_FOLDER]
                        [--markpres_folder MARKPRES_FOLDER]
                        [--pfam_folder PFAM_FOLDER]
                        [--mirna_folder MIRNA_FOLDER]
                        [--metaphlan_title METAPHLAN_TITLE]
                        [--pwy_title PWY_TITLE] [--cov_title COV_TITLE]
                        [--genefam_title GENEFAM_TITLE]
                        [--pfam_title PFAM_TITLE] [--mirna_title MIRNA_TITLE]
                        [-rc] [-pr PROPORTIONS] [-cf CONTROL_FILENAME]
                        [-mulo RANDOM_CONTROLS_MULTIPLE_OF]
                        [-cc CONTROL_CONDITIONS]
                        [-cn CONTROL_NAMES CONTROL_NAMES] [-cw CONTROL_WORD]
                        [-sc SELECT_COLUMNS [SELECT_COLUMNS ...]]
                        [-scff SELECT_COLUMNS_FROM_FILE] [-mx] [-ys]
                        [-tx {k__,p__,c__,o__,f__,g__,s__}] [-sr]
                        [-x EXCLUDE_SAMPLES [EXCLUDE_SAMPLES ...]] [-om]
                        [-mgp MERGE_PROFILE_EXEC]
                        [-pn PERCENTILE_NORMALIZATION] [-lg] [-dm]
                        [-zt ZERO_THRESHOLD [ZERO_THRESHOLD ...]] [-mf]
                        [-fs FEATURE_SELECTION] [-wk WILCOXON] [--bonferroni]
                        [--mapped_number] [--mapped_percentage] [-t] [-b]
                        [-of OUTPUT_FILE] [-gs GIVE_STATISTICS] [-fo]
                        [-fac FEAT_AND_CONDITION [FEAT_AND_CONDITION ...]]
                        [--grad GRAD [GRAD ...]] [--grad_col GRAD_COL]
                        [--log_gradient]
                        datasets [datasets ...]

positional arguments:
  datasets              E.g.: LiJ_2014.study_condition:control,body_site:stool
                        .Another possibility is setting: LiJ_2014.study_condit
                        ion:control,col:country:FRA,col:body_subsite:subgingiv
                        al_plaque to add a new column to the loaded metadata
                        table.

optional arguments:
  -h, --help            show this help message and exit
  --base_path BASE_PATH
  --metadata_path METADATA_PATH
  --metadata_name METADATA_NAME
  --metaphlan_path METAPHLAN_PATH
  --change_profile_func
  --metaphlan
  --pwyrelab
  --genefam
  --markab
  --markpres
  --pwycov
  --pfam
  --mirna
  --metaphlan_folder METAPHLAN_FOLDER
  --pwyrelab_folder PWYRELAB_FOLDER
  --genefam_folder GENEFAM_FOLDER
  --pwycoverage_folder PWYCOVERAGE_FOLDER
  --markab_folder MARKAB_FOLDER
  --markpres_folder MARKPRES_FOLDER
  --pfam_folder PFAM_FOLDER
  --mirna_folder MIRNA_FOLDER
  --metaphlan_title METAPHLAN_TITLE
  --pwy_title PWY_TITLE
  --cov_title COV_TITLE
  --genefam_title GENEFAM_TITLE
  --pfam_title PFAM_TITLE
  --mirna_title MIRNA_TITLE
  -rc, --randomize_controls
  -pr PROPORTIONS, --proportions PROPORTIONS
                        1,2 (1:moiety of cnts wrt cases, 2:moiety of training
                        cnts wrt to test cnts.)
  -cf CONTROL_FILENAME, --control_filename CONTROL_FILENAME
  -mulo RANDOM_CONTROLS_MULTIPLE_OF, --random_controls_multiple_of RANDOM_CONTROLS_MULTIPLE_OF
                        If not specified, uses as bases for the mutiplicity
                        the allready collected samples.
  -cc CONTROL_CONDITIONS, --control_conditions CONTROL_CONDITIONS
                        e.g: study_condition:IBD:control,
                        gender:female,study_condition:control
  -cn CONTROL_NAMES CONTROL_NAMES, --control_names CONTROL_NAMES CONTROL_NAMES
                        the names to the train and test of the randomized
                        controls;must be two but can be equal: the first is
                        the training.
  -cw CONTROL_WORD, --control_word CONTROL_WORD
                        add a common field to all the externally derived
                        controls , the added field will be in third position
                        (after dataset_name sampleID) under header
                        "ext_characterization"
  -sc SELECT_COLUMNS [SELECT_COLUMNS ...], --select_columns SELECT_COLUMNS [SELECT_COLUMNS ...]
                        "Selected Columns"
  -scff SELECT_COLUMNS_FROM_FILE, --select_columns_from_file SELECT_COLUMNS_FROM_FILE
                        "Select Columns from a File"
  -mx, --mixed_taxa     "With Metaphlan, Uses All the Taxa Together"
  -ys, --yes_strain     "With Metaphlan, Leave the Strain Level There"
  -tx {k__,p__,c__,o__,f__,g__,s__}, --taxon {k__,p__,c__,o__,f__,g__,s__}
  -sr, --shrink         "With Metaphlan, Uses as Names only The Currenty
                        Selected Taxon"
  -x EXCLUDE_SAMPLES [EXCLUDE_SAMPLES ...], --exclude_samples EXCLUDE_SAMPLES [EXCLUDE_SAMPLES ...]
  -om, --only_metadata  This Option is Usefull Also for Control To Be
                        Ramdomised After.
  -mgp MERGE_PROFILE_EXEC, --merge_profile_exec MERGE_PROFILE_EXEC
  -pn PERCENTILE_NORMALIZATION, --percentile_normalization PERCENTILE_NORMALIZATION
                        percentile normalization: eg -pn
                        study_condtion:control:CRC whatever is in the 1st
                        field is the column to look at. Whats is int eh 2nd is
                        the reference percentiles distribution (should be some
                        kind of control) whatever is in the last field will be
                        uniformed to the reference: make usage of add_column
                        option of this script to help yourself.
  -lg, --log_transform  performs (log(1+anyfeature))
  -dm, --dense_matrix   cut out features being zeros in <zero_threshold>
                        percentage of the samples [ default = 10 per ]
  -zt ZERO_THRESHOLD [ZERO_THRESHOLD ...], --zero_threshold ZERO_THRESHOLD [ZERO_THRESHOLD ...]
  -mf, --meta_features
  -fs FEATURE_SELECTION, --feature_selection FEATURE_SELECTION
                        "Differs From -scff in That This Only Has Control on
                        Features, -scff Also On Metadata"
  -wk WILCOXON, --wilcoxon WILCOXON
                        --wilcoxon study_condition:CRC:control
  --bonferroni
  --mapped_number
  --mapped_percentage
  -t, --transpose       Output Has Field Names On Columns
  -b, --both            Outputs Two Results, One Has Fields On Indexes The
                        Other On Columns.
  -of OUTPUT_FILE, --output_file OUTPUT_FILE
                        "If Don't Specify the -of, There You'll Be Able to
                        Sedn the Result As A STDOUT."
  -gs GIVE_STATISTICS, --give_statistics GIVE_STATISTICS
                        Statistics on metadata -gs
                        study_condition:sample_type:disease, or -gs
                        random/randomization/rand/randrandrand/randanything
  -fo, --feat_only
  -fac FEAT_AND_CONDITION [FEAT_AND_CONDITION ...], --feat_and_condition FEAT_AND_CONDITION [FEAT_AND_CONDITION ...]
                        "Alternative to -scff, Completely On Command Line"
  --grad GRAD [GRAD ...]
  --grad_col GRAD_COL
  --log_gradient

```

* cmdpy_betadiversity.py: the data collected by cmdpy_dataset.py can be send in STDIN to cmdpy_betadiversity.py (see below).

```
usage: cmdpy_betadiversity.py [-h] [-if STDIN] [-of STDOUT] [-adiv]
                              [-sc SAVE_COORDINATES] [-lc LOAD_COORDINATES]
                              [-d {lbraycurtis,sbraycurtis,canberra,chebyshev,correlation,dice,hamming,jaccard,kulsinski,mahalanobis,matching,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,yule,cityblock,cosine,euclidean,l1,l2,manhattan,braycurtis,precomputed}]
                              [-a {mds,pca,nmds,boxp}]
                              [-z {k__,s__,PWY,UniRef90,g__}] [-si SAMPLE_ID]
                              [-ci CLASSES_ID] [-m MASK [MASK ...]]
                              [-ll LEGEND_LOC] [-fmt {svg,png}]
                              [--no_ordination] [--boxplot] [-ex] [--annot]
                              [--title TITLE] [--dot_size DOT_SIZE]
                              [-txt TEXT_ON] [--alpha ALPHA]
                              [--facecolor {white,whitegrid,darkgrid}]
                              [--squared_plot]
                              [-cm {RdYlBu,plasma,inferno,winter,copper,viridis,YlGnBu,YlGnBu_r,cool,cool_r,RdYlBu_r,plasma_r,inferno_r,winter_r,copper_r,viridis_r,YlGnBu_r,YlGnBu_r_r,cool_r,cool_r_r}]
                              [-cn CBAR_TITLE] [-go GRADIENT_ON]
                              [-gf GRADIENT_ONLY_FOR] [--intra_individual]
                              [--p_values {above,below}]
                              [--p_values_only_for P_VALUES_ONLY_FOR [P_VALUES_ONLY_FOR ...]]
                              [--welch]

optional arguments:
  -h, --help            show this help message and exit
  -if STDIN, --stdin STDIN
  -of STDOUT, --stdout STDOUT
-adiv, --alpha_diversity
  -sc SAVE_COORDINATES, --save_coordinates SAVE_COORDINATES
  -lc LOAD_COORDINATES, --load_coordinates LOAD_COORDINATES
  -d {lbraycurtis,sbraycurtis,canberra,chebyshev,correlation,dice,hamming,jaccard,kulsinski,mahalanobis,matching,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,yule,cityblock,cosine,euclidean,l1,l2,manhattan,braycurtis,precomputed}, --distance {lbraycurtis,sbraycurtis,canberra,chebyshev,correlation,dice,hamming,jaccard,kulsinski,mahalanobis,matching,minkowski,rogerstanimoto,russellrao,seuclidean,sokalmichener,sokalsneath,sqeuclidean,yule,cityblock,cosine,euclidean,l1,l2,manhattan,braycurtis,precomputed}
  -a {mds,pca,nmds,boxp}, --algorithm {mds,pca,nmds,boxp}
  -z {k__,s__,PWY,UniRef90,g__}, --feature_identifier {k__,s__,PWY,UniRef90,g__}
  -si SAMPLE_ID, --sample_id SAMPLE_ID
  -ci CLASSES_ID, --classes_id CLASSES_ID
  -m MASK [MASK ...], --mask MASK [MASK ...]
  -ll LEGEND_LOC, --legend_loc LEGEND_LOC
  -fmt {svg,png}, --format {svg,png}
  --no_ordination
  --boxplot
  -ex, --explain
  --annot
  --title TITLE
  --dot_size DOT_SIZE
  -txt TEXT_ON, --text_on TEXT_ON
  --alpha ALPHA
  --facecolor {white,whitegrid,darkgrid}
  --squared_plot
  -cm {RdYlBu,plasma,inferno,winter,copper,viridis,YlGnBu,YlGnBu_r,cool,cool_r,RdYlBu_r,plasma_r,inferno_r,winter_r,copper_r,viridis_r,YlGnBu_r,YlGnBu_r_r,cool_r,cool_r_r}, --cmap {RdYlBu,plasma,inferno,winter,copper,viridis,YlGnBu,YlGnBu_r,cool,cool_r,RdYlBu_r,plasma_r,inferno_r,winter_r,copper_r,viridis_r,YlGnBu_r,YlGnBu_r_r,cool_r,cool_r_r}
  -cn CBAR_TITLE, --cbar_title CBAR_TITLE
  -go GRADIENT_ON, --gradient_on GRADIENT_ON
                        must be a column in the data.
  -gf GRADIENT_ONLY_FOR, --gradient_only_for GRADIENT_ONLY_FOR
                        definition of a class in case you want gradient only
                        for that class.
  --intra_individual
  --p_values {above,below}
  --p_values_only_for P_VALUES_ONLY_FOR [P_VALUES_ONLY_FOR ...]
  --welch

```

#### cmdpy_betadiversity.py computes alpha and beta diversity, and provides graphical support to the analyses. ####

* this tool computes beta and alpha diversity, and provides graphical and statistical representation of the required result, including full-optional ordination plots (MDS, PCA, NMDS).
* when alpha diversity is required 4 box-plots are returned, according to different richness measueres: Shannon entropy, Gini-Simpson information content, taxon richness and logarithmic taxon richness.
* when beta-diversity is required, the box-plot graphical representation is available so as the ordination-plot representation.

##### cmdpy_betadiversity contains many options to speed up ordination-plot represenation needings #####

* the follwing are the parameters to saving and loading point coordinates among different ordination plots.
* relying o the same sample IDs, different graphical versions of the same points and possible

```
-sc SAVE_COORDINATES, --save_coordinates SAVE_COORDINATES
-lc LOAD_COORDINATES, --load_coordinates LOAD_COORDINATES
```

* according to the same principle

```
-go GRADIENT_ON, --gradient_on GRADIENT_ON
                        must be a column in the data.
-gf GRADIENT_ONLY_FOR, --gradient_only_for GRADIENT_ONLY_FOR
                        definition of a class in case you want gradient only
                        for that class.
```

* are the parameters allowing a uniform color gradient among all or part of the points in the graph


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