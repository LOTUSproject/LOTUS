
# LOTUS: a Single- and Multitask Machine Learning Algorithm for the Prediction of Cancer Driver Genes

## Abstract 

Cancer driver genes, i.e., oncogenes and tumor suppressor genes, are involved in the acquisition of important functions in tumors, providing a selective growth advantage,
allowing uncontrolled proliferation and avoiding apoptosis. It is therefore important to identify these driver genes, both for the fundamental understanding of cancer and to
help finding new therapeutic targets or biomarkers. Although the most frequently mutated driver genes have been identified, it is believed that many more remain to be discovered, particularly for driver genes specific to some cancer types.

In this paper, we propose a new computational method called LOTUS to predict new driver genes. LOTUS is a machine-learning based approach which allows to
integrate various types of data in a versatile manner, including information about gene mutations and protein-protein interactions. In addition, LOTUS can predict cancer
driver genes in a pan-cancer setting as well as for specific cancer types, using a multitask learning strategy to share information across cancer types.

We empirically show that LOTUS outperforms four other state-of-the-art driver gene prediction methods, both in terms of intrinsic consistency and prediction accuracy, and provide predictions of new cancer genes across many cancer types.

## Supplementary files

* S1 Table - List of cancer types. Cancer types derived from annotations in the 20/20 mutation dataset along with their numbers of associated OG and TSG. 
* S2 Table - Description of cancer types. Descriptors of all cancer types according to their localizations and types that are used to compute the disease kernel used by LOTUS2. 
* S3 Table - TSG and OG rankings for LOTUS with the 20/20, the TUSON and the CGC v86 datasets. Note that the training sets were removed every time. 

