# Food production resistome

This repository contains the code used for the analysis of antimicrobial resistance (AMR, resistome) screening in 113 European facilities from 5 different countries (i.e. Austria, Iceland, Ireland, Italy and Spain) producing cheese, meat, vegetables and fish products. Overall, a total of 1,780 samples of raw materials, final products and industrial surfaces collected and subjected to shotgun sequencing (Illumina NovaSeq). The results have been submitted for publication together with all the additional metadata.  

These analyses are part of the activities performed within the [EU-MASTER project](https://www.master-h2020.eu/readmore.html), where a repository with code for the different bioinformatic analyses exists and is hosted here: [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines).  

More precisely, the analysis of the resistome in food industries was one of the main goals of WP5, that ended up with the release of the deliverable [FoodGenVir database](https://zenodo.org/records/8344969). The foodGenVir database contains those MAGs obtained from food and food-related environments samples that were collected in the MASTER project and that harbored antimicrobial resistance (AMR) genes. In total, 1828 samples from food and food-producing environments that were successfully sequenced in MASTER and included in the analyses, which lead to 14,815 high-quality MAGs. The data present in the current release contains the assembled-genome FASTA files from the AMR-carrying MAGs and a MASTER_MAGs_FGV.tsv table, that contains the taxonomic information of each MAG as well as the different AMR genes found in its genome and the antimicrobial compounds that those genes confer resistance to according to the ResFinder database.  

<br>

The following sections contain a thorough description of the different steps conducted for the analysis of the resistome in food-producing companies. Worth to mention that the **workflow and code used here can be adaptaed for its use in any microbiome investigation**


