# LOTUS

We propose a new computational method called LOTUS to predict new driver genes. LOTUS is a machine-learning based approach which allows to integrate various types of data in a versatile manner, including information about gene mutations and protein-protein interactions. In addition, LOTUS can predict cancer driver genes in a pan-cancer setting as well as for specific cancer types, using a multitask learning strategy to share information across cancer types.

## Submitted paper

Our paper can be downloaded in the folder "paper", along with the supplementary files mentioned in it. If appropriate, please cite us as

```
O. Collier, V. Stoven and J.-P. Vert, LOTUS: a Single- and Multitask Machine Learning
Algorithm for the Prediction of Cancer Driver Genes, Preprint. doi: https://doi.org/10.1101/398537
```

## Algorithm

The codes in R and the data that we used to obtain the results presented in our papers are fully available in the folders "code" and "data", except some heavy matrices that need to be downloaded online. To reproduce our results, you need to copy all files in your choice folder. Then, you can directly run one of the files "run_lotus.R", "evaluate_lotus.R" or "evaluate_lotus_multitask.R".

* To run the script "run_lotus.R", you need at a minimum the source file "lotus.R", a feature matrix (for example "features_2020.txt") and a corresponding PPI kernel (for example "PPIKernel_2020.txt"). Then, choose the arguments in the script:

```R
driver_type <- 'tsg'     # 'tsg' or 'og'
dataset <- 'cosmicv86'        # '2020', 'cosmicv86', 'mutsig' or 'tuson'
ppi <- 'yes'           # 'yes', 'no' or 'only'
```

For example, when driver_type is "tsg", dataset is "tuson" and ppi is "yes", run_lotus.R computes a ranking of all genes listed in features_tuson.txt, except for the known tumor suppressors in tsg_tuson.txt, according to their potential as tumor suppressors, using a training set of known tumor suppressors from tsg_tuson.txt and information from the PPI network. 

To use LOTUS with different datasets, you only have to change the lines where the features, the set of known driver genes and possibly the PPI kernel so that the code fits the name of your datasets. Pay attention to the fact that the feature matrix should contain columns named "Frameshift", "LOF", "Splice" for tumor suppressor prediction, and "Entropy.Score", "Missense.Damaging", "Missense.total" for oncogene prediction. Alternatively, for example if you wish to use other features, you will have to change the corresponding lines in the script.

* To run the script "evaluate_lotus.R", you need at a minimum the source files "lotus.R" and "evaluation.R", a feature matrix (for example "features_2020.txt") and a corresponding PPI kernel (for example "PPIKernel_2020.txt"). Then, choose the arguments in the script as for "run_lotus.R". This script returns the consistency error of LOTUS trained with the chosen data.

* To run the script "run_lotus_multitask.R", you need at a minimum the source file "lotus.R", a feature matrix (for example "features_2020.txt"), a corresponding PPI kernel (for example "PPIKernel_2020.txt"), a similarity matrix between diseases (for example "Kernel_similarity.RData") and lists of known cancer driver genes for every considered disease (for example "tsg_per_Diseases.RData"). Then, choose the arguments in the script:

```R
task <- 1
driver_type <- 'tsg' # 'tsg' or 'og'
version <- 'lotus' # 'lotus' for LOTUS, 'lotus2' for LOTUS2, 'aggregation' for aggregation LOTUS and 'onetask' for onetask LOTUS
```

For example, when task is 1, driver_type is "tsg" and version is "lotus", run_lotus_multitask.R computes a ranking of all genes listed in features_2020.txt (which is the default dataset), except for the known type-specific tumor suppressors in tsg_per_Diseases.RData, according to their potential as tumor suppressors, using LOTUS with a training set of known tumor suppressors from tsg_per_Diseases.RData. "Kernel_similarity.RData" is only necessary for LOTUS2.

To use LOTUS or LOTUS2 with different datasets, similar changes have to be made as for the one-task version of LOTUS (see above).

* To run the script "evaluate_lotus_multitask.R", you need at a minimum the source files "lotus.R" and "evaluation.R", a feature matrix (for example "features_2020.txt"), a corresponding PPI kernel (for example "PPIKernel_2020.txt"), a similarity matrix between diseases (for example "Kernel_similarity.RData") and lists of known cancer driver genes for every considered disease (for example "tsg_per_Diseases.RData"). Then, choose the arguments in the script as for "run_lotus_multitask.R".



