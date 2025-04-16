# Food production resistome

This repository contains the code used for the analysis of antimicrobial resistance (AMR, resistome) screening in 113 European facilities from 5 different countries (i.e. Austria, Iceland, Ireland, Italy and Spain) producing cheese, meat, vegetables and fish products. Overall, a total of 1,780 samples of raw materials, final products and industrial surfaces collected and subjected to shotgun sequencing (Illumina NovaSeq). The results have been submitted for publication together with all the additional metadata.  

These analyses are part of the activities performed within the [EU-MASTER project](https://www.master-h2020.eu/readmore.html), where a repository with code for the different bioinformatic analyses exists and is hosted here: [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines).  

More precisely, the analysis of the resistome in food industries was one of the main goals of WP5, that ended up with the release of the deliverable [FoodGenVir database](https://zenodo.org/records/8344969). The foodGenVir database contains those MAGs obtained from food and food-related environments samples that were collected in the MASTER project and that harbored antimicrobial resistance (AMR) genes. In total, 1828 samples from food and food-producing environments that were successfully sequenced in MASTER and included in the analyses, which lead to 14,815 high-quality MAGs. The data present in the current release contains the assembled-genome FASTA files from the AMR-carrying MAGs and a MASTER_MAGs_FGV.tsv table, that contains the taxonomic information of each MAG as well as the different AMR genes found in its genome and the antimicrobial compounds that those genes confer resistance to according to the ResFinder database.  

<br>

The following sections contain a thorough description of the different steps conducted for the analysis of the resistome in food-producing companies. Worth to mention that the **workflow and code used here can be adaptaed for its use in any microbiome investigation**

<br>

**Index**  
1. [Installation](#id1)
2. [Sequencing data and QC](#id2)
3. [Assembly-free analysis](#id2)
4. [Assembly-based analysis](#id3)
4.1. [Contig-level analysis](#id31)
4.2. [Metagenome Associated Genomes (MAGs)](#id32)
5. [Visualization](#id4)

<br>

## 1. Installation<a name="id1"></a>

In order to reproduce the analyses conducted in the manuscript, we uploaded a yaml file containing the different software and versions used.  
> Please note that the recreation of complex environments by using conda can fail when multiple software with multiple versions are required. If that happens, please try to "relax" the yaml file by enabling conda to look for different software versions (instead of pushing the one in the yaml, if there is not any incompatibilities) or try to generate different environments with less software in each yaml file.

```bash
conda env create -n master -f master.yml
conda activate master
```

<br>

[Back to index](#idx)

<br>

## 2. Sequencing data and QC<a name="id2"></a>

The sequencing data per sample, their corresponding BioProject and BioSample IDs for public download and their corresponding metadata can be observed in the Supplementary Table 3 of the manuscript (not publicly available yet) and also [here].  
The different accession numbers contain the already QCed data, while the process is described in the corresponding section of [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines/tree/master/02-Preprocessing).  

Briefly:

```bash
conda install preprocessing -c fasnicar
parallel -j NCPU 'preprocess.sh -i {} [other params]' ::: `ls input_folder`
```

where:

- `preprocess.sh` takes one parameter which is the input folder containing the raw reads
- input folder should contains the raw reads

<br>

[Back to index](#idx)

<br>

## 3. Assembly-free analysis<a name="id3"></a>

<br>

[Back to index](#idx)

<br>

## 4. Assembly-free analysis<a name="id4"></a>

<br>

[Back to index](#idx)

<br>

