#!/bin/bash

bp="/scratchCM/users/paolo.manghi/multidataset_machinelearning/"
bd="CM_periimplantitis"
 
#python "${bp}"cmdpy_dataset.py "${bd}".study_condition:healthy,col:color:green,col:shape:o,col:define:healthy "${bd}".study_condition:peri-implantitis,col:color:red,col:shape:o,col:define:periimplantitis "${bd}".study_condition:mucositis,col:color:deepskyblue,col:shape:o,col:define:mucositis --shrink --grad s__Porphyromonas_gingivalis s__Porphyromonas_endodontalis s__Tannerella_forsythia s__Fretibacterium_fastidiosum s__Fusobacterium_nucleatum  s__Prevotella_intermedia s__Treponema_denticola | python cmdpy_betadiversity.py --save_beta_diversity trial_beta_diversity.csv #--save_centroids trial_centroids.csv

python "${bp}"cmdpy_dataset.py "${bd}".study_condition:healthy,col:color:green,col:shape:o,col:define:healthy "${bd}".study_condition:peri-implantitis,col:color:red,col:shape:o,col:define:periimplantitis "${bd}".study_condition:mucositis,col:color:deepskyblue,col:shape:o,col:define:mucositis "${bd}".study_group:peri-implantitis,study_condition:control,col:color:maroon,col:shape:o,col:define:periimplantitis-contr "${bd}".study_group:mucositis,study_condition:control,col:color:blue,col:shape:o,col:define:mucositis-contr "${bd}".study_group:healthy,study_condition:control,col:color:forestgreen,col:shape:o,col:define:healthy-contr --shrink --grad s__Porphyromonas_gingivalis s__Porphyromonas_endodontalis s__Tannerella_forsythia s__Fretibacterium_fastidiosum s__Fusobacterium_nucleatum  s__Prevotella_intermedia s__Treponema_denticola | python cmdpy_betadiversity.py --save_beta_diversity trial_beta_diversity_2.csv --save_centroids trial_centroids_2.csv
 

python "${bp}"cmdpy_dataset.py "${bd}".study_condition:healthy,col:color:green,col:shape:o,col:define:healthy "${bd}".study_condition:peri-implantitis,col:color:red,col:shape:o,col:define:periimplantitis "${bd}".study_condition:mucositis,col:color:deepskyblue,col:shape:o,col:define:mucositis "${bd}".study_group:peri-implantitis,study_condition:control,col:color:maroon,col:shape:o,col:define:periimplantitis-contr "${bd}".study_group:mucositis,study_condition:control,col:color:blue,col:shape:o,col:define:mucositis-contr "${bd}".study_group:healthy,study_condition:control,col:color:forestgreen,col:shape:o,col:define:healthy-contr --shrink --grad s__Porphyromonas_gingivalis s__Porphyromonas_endodontalis s__Tannerella_forsythia s__Fretibacterium_fastidiosum s__Fusobacterium_nucleatum  s__Prevotella_intermedia s__Treponema_denticola | python cmdpy_centroids.py --load_centroids trial_centroids_2.csv --from periimplantitis --distance cosine --color_by classes


#python "${bp}"cmdpy_dataset.py "${bd}".study_condition:healthy,col:color:green,col:shape:o,col:define:healthy "${bd}".study_condition:peri-implantitis,col:color:red,col:shape:o,col:define:periimplantitis "${bd}".study_condition:mucositis,col:color:deepskyblue,col:shape:o,col:define:mucositis --shrink --grad s__Porphyromonas_gingivalis s__Porphyromonas_endodontalis s__Tannerella_forsythia s__Fretibacterium_fastidiosum s__Fusobacterium_nucleatum s__Prevotella_intermedia s__Treponema_denticola | python cmdpy_centroids.py --load_centroids trial_centroids.csv --from periimplantitis --beta_file trial_beta_diversity.csv --distance correlation --color_by gradient ##--try_any_metric



### s__Fusobacterium_nucleatum


#python "${bp}"cmdpy_dataset.py "${bd}".study_condition:healthy,col:color:green,col:shape:o,col:define:healthy "${bd}".study_condition:peri-implantitis,col:color:red,col:shape:o,col:define:periimplantitis --shrink --grad s__Porphyromonas_gingivalis s__Porphyromonas_endodontalis s__Tannerella_forsythia s__Fusobacterium_nucleatum s__Fretibacterium_fastidiosum s__Prevotella_intermedia s__Treponema_denticola | python cmdpy_betadiversity.py --save_centroids trial_centroids.csv #| python cmdpy_centroids.py



#python cmdpy_betadiversity.py --compute_centroids
#average_distance_from_centroids.py --avg_dist_from_centroids mettiqua_latua_betadiversity --centroids periimplantitis healthy -of plot_dentistry
