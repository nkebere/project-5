# Project Description

This project is a replication of a Single Cell RNA-seq Anlaysis done on pancreatic cells from human and mice samples by Baron et. al (2016)

N.B. The study this project is referencing used sc-RNA-seq data from 4 human donors and 2 mice strains for single cell anlysis. However, for this project the sc-RNA-seq data associated with a 51-year-old female donor was only used.

Reference:
Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., Ryu, J. H., Wagner, B. K., Shen-Orr, S. S., Klein, A. M., Melton, D. A., & Yanai, I. (2016). A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure. Cell systems, 3(4):346–360.e4. https://doi.org/10.1016/j.cels.2016.08.011


# Contributors

Nitsueh Kebere - nkebere@bu.edu



# Repository Contents

This repository contains all the scripts and packages used in the single cell RNA-seq analysis.


# Installation 

To successfully perfom this analysis install the requried packages provided in the packages-used folder 

# Scripts/codes description

The counts_barcodes folder has shell scripts that calculate the number of reads per distinct barcode for a FASTQ file containing CB+UMI raw sequence. This folder also contains the outputs generated from the count shell scripts 

The whitelist folder contains whitelist.R script that generated the whitelsit barodes for three of the libraries.
It also has the generated barcode whitelists for three of the libraries 

The salmon-count-matrix folder contains all scripts needed to run salmon alevin. It contains shell scripts to build salmon index for reference transcriptome, script to create transcript to gene mapping file, and shell script to submit a salmon alevin job. Additionally, there is a merge_matrix.py file that merges the UMI count matrix generated by salmon for the three libraries into one UMI count matrix 

The sc_analysis.R  is a script for the single-cell analysis pipeline using the sc-RNA-seq data associated with a 51-year-old female donor from Baron et. al.

