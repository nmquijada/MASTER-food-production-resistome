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
3. [Assembly-free analysis](#id3)
4. [Assembly-based analysis](#id4)  
4.1. [Contig-level analysis](#id41)  
4.2. [MAG-level analysis](#id42)
5. [Mobile Genetic Elements](#id5)
6. [Visualization](#id6)

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

The sequencing data per sample, their corresponding BioProject and BioSample IDs for public download and their corresponding metadata can be observed in the Supplementary Table 3 of the manuscript (not publicly available yet) and also [here](https://raw.githubusercontent.com/nmquijada/food-production-resistome/refs/heads/main/files/Supplementary_Table_3_MASTER_metadata.tsv).  
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

## 4. Assembly-based analysis<a name="id4"></a>

From the QC reads, an assembly and was performed as thorughlly described in the corresponding section of [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines/tree/master/05-Assembly_pipeline).  
The resulting draft metagenomes (contigs) in FASTA file were either subjected to the [TORMES pipeline](https://github.com/nmquijada/tormes) for further downstream analysis at contig-level or used for metagenome assembled genomes (MAGs) extraction.  

<br>

[Back to index](#idx)

<br>

## 4.1 Contig-level analysis<a name="id41"></a>

1. Define your variables

```bash
WORKDIR=# define your working directory
# Then define the "input directory" containing the contigs FASTA files in different directories regarding the sample. e.g
INDIR=/path/to/02.metagenome_assembly
```

The `metagenome_assembly` directory is expected to have the assemblies sorted by samples, as:

```
|-- /path/to/metagenome_assembly/00.raw_reads
|-- /path/to/metagenome_assembly/01.QC_reads
|-- /path/to/metagenome_assembly/02.metagenome_assembly
```

```
|-- /path/to/metagenome_assembly/02.metagenome_assembly/Sample_01/contigs.fasta
|-- /path/to/metagenome_assembly/02.metagenome_assembly/Sample_02/contigs.fasta
|-- ...
|-- /path/to/metagenome_assembly/02.metagenome_assembly/Sample_n/contigs.fasta
```

2. Create TORMES metadata. There's an useful tutorial on how to do it [here](https://github.com/nmquijada/tormes/wiki/Shortcut-to-generate-the-metadata-file-for-TORMES)

```bash
for SAMPLE in $(ls ${INDIR}); do
  mkdir ${WORKDIR}/${SAMPLE}
  echo -e "${SAMPLE}\tGENOME\t${INDIR}/${SAMPLE}/contigs.fasta\tcontigs from ${SAMPLE}" >> ${WORKDIR}/tormes-metadata.txt
done
sed -i "1iSamples\tRead1\tRead2\tDescription" > ${WORKDIR}/tormes-metadata.txt
# further description and metadata fields can be used, if needed
```

3. Run TORMES
```bash
tormes -m ${WORKDIR}/tormes-metadata.txt -t 64 --no_pangenome --gene_min_id 80 --gene_min_cov 80 --only_gene_prediction --prodigal_options "-p meta" --custom_genes_db "plasmidfinder" --min_contig_len 1000 -o ${WORKDIR}/MASTER-tormes 
```

This command will generate the AMR screening results based on sequence alignment (BLASTN) of the contigs against 3 databases (ResFinder, CARD, Argannot).  
For this study we used RedFinder database, so the results will be stored in:  

`${WORKDIR}/MASTER-tormes/resfinder/`


<br>

[Back to index](#idx)

<br>

## 4.2 MAG-level analysis<a name="id42"></a>

<br>

[Back to index](#idx)

<br>

## 5. Mobile Genetic Elements<a name="id5"></a>

Different software were used for the investigation of MGE and the main steps are described below and with more detail [here](https://github.com/nmquijada/MASTER-food-production-resistome/blob/main/assembly-based_workflow/mobilome_workflow.md)  

### 5.1 Plasmids

<br>

#### PlasmidFinder

The results from this database were performed with TORMES and stored in: 

`${WORKDIR}/MASTER-tormes/custom_genes_db/plasmidfinder`

The results indicate all those contigs where a plasmidic replicon was found.  
Due to the short length of the plasmidic replicons contained in PlasmidFinder, only those hits > 95% identity and coverage where considered

<br>

#### Platon

```bash
# Define database
PLATON_DB=/PATH/TO/Platon-DB/db

# Generate output directory
mkdir ${WORKDIR}/MGE/platon
platon --db ${PLATON_DB} --output ${WORKDIR}/MGE/platon -p ${SAMPLE} --threads 8 --meta ${WORKDIR}/MASTER-tormes/${SAMPLE}.fasta
```

The output will contain all those contigs that are potentially plasmidic in a file `${SAMPLE}.plasmid.fasta` located in `${WORKDIR}/MGE/platon`

<br>

#### geNomad

```bash
# Define database
GENOMAD_DB=/PATH/TO/genomad/genomad_db

# Run geNomad
genomad end-to-end --cleanup ${WORKDIR}/MASTER-tormes/${SAMPLE}.fasta ${WORKDIR}/genomad ${GENOMAD_DB}
```

The results in file `${SAMPLE}_plasmid_summary.tsv` will contain the ID of all those contigs that are potentially plasmidic.  
This file is stored in `${WORKDIR}/genomad/${SAMPLE}_summary/`

<br>

[Back to index](#idx)

<br>

## 6. Visualization<a name="id6"></a>

<br>

[Back to index](#idx)

<br>


