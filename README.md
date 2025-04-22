# Food production resistome

This repository contains the code used for the analysis of antimicrobial resistance (AMR, resistome) screening in 113 European facilities from 5 different countries (i.e. Austria, Iceland, Ireland, Italy and Spain) producing cheese, meat, vegetables and fish products. Overall, a total of 1,780 samples of raw materials, final products and industrial surfaces collected and subjected to shotgun sequencing (Illumina NovaSeq). The results have been submitted for publication together with all the additional metadata.  

These analyses are part of the activities performed within the [EU-MASTER project](https://www.master-h2020.eu/readmore.html), where a repository with code for the different bioinformatic analyses exists and is hosted here: [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines).  

More precisely, the analysis of the resistome in food industries was one of the main goals of WP5, that ended up with the release of the deliverable [FoodGenVir database](https://zenodo.org/records/8344969). The foodGenVir database contains those MAGs obtained from food and food-related environments samples that were collected in the MASTER project and that harbored antimicrobial resistance (AMR) genes. In total, 1828 samples from food and food-producing environments that were successfully sequenced in MASTER and included in the analyses, which lead to 14,815 high-quality MAGs. The data present in the current release contains the assembled-genome FASTA files from the AMR-carrying MAGs and a MASTER_MAGs_FGV.tsv table, that contains the taxonomic information of each MAG as well as the different AMR genes found in its genome and the antimicrobial compounds that those genes confer resistance to according to the ResFinder database.  

<br>

The following sections contain a thorough description of the different steps conducted for the analysis of the resistome in food-producing companies. Worth to mention that the **workflow and code used here can be adaptaed for its use in any microbiome investigation**

<br>

**Index**  
1. [Installation and set up](#id1)
2. [Sequencing data and QC](#id2)
3. [Assembly-free analysis](#id3)
4. [Assembly-based analysis](#id4)  
5. [Mobile Genetic Elements](#id5)
6. [Visualization](#id6)

<br>

## 1. Installation and set up<a name="id1"></a>

In order to reproduce the analyses conducted in the manuscript, we uploaded [different yaml files](https://github.com/nmquijada/MASTER-food-production-resistome/tree/main/files) containing the different software and versions used.  
> Please note that the recreation of complex environments by using conda can fail when multiple software with multiple versions are required. If that happens, please try to "relax" the yaml file by enabling conda to look for different software versions (instead of pushing the one in the yaml, if there is not any incompatibilities) or try to generate different environments with less software in each yaml file.

```bash
conda env create -n <environment> -f <environment>.yml
conda activate <environment>
```

<br>

The reproduction of the code requires that the data is stored hierarchically as follows.   

The first level is `${DATASET}` represent the different company partners. Three main directories per `${DATASET}` are available:  

```
|-- /path/to/${DATASET}/reads/
|-- /path/to/${DATASET}/contigs/
|-- /path/to/${DATASET}/mags/
```

Second, each sample from the dataset is a subdirectory stored in the different directories, depending on the type of data you will encounter (i.e. reads, contigs or MAGs). Here, the total of samples are represented as `${SAMPLE}`:  

```
|-- /path/to/${DATASET}/reads/${SAMPLE}/
|-- /path/to/${DATASET}/contigs/${SAMPLE}/
|-- /path/to/${DATASET}/mags/${SAMPLE}/
```

Within each `${SAMPLE}` sudirectory, you will find the corresponding data:

```
|-- /path/to/${DATASET}/reads/${SAMPLE}/${SAMPLE}_R1.fastq.gz
|-- /path/to/${DATASET}/reads/${SAMPLE}/${SAMPLE}_R2.fastq.gz
|-- /path/to/${DATASET}/contigs/${SAMPLE}/${SAMPLE}.fasta
|-- /path/to/${DATASET}/mags/${SAMPLE}/${SAMPLE}__bin1.fasta
|-- /path/to/${DATASET}/mags${SAMPLE}//${SAMPLE}__bin2.fasta
|-- ...
|-- /path/to/${DATASET}/mags/${SAMPLE}__binN.fasta
```


All these information can be be observed in the Supplementary Table 3 of the manuscript (not publicly available yet) and also [here](https://raw.githubusercontent.com/nmquijada/food-production-resistome/refs/heads/main/files/Supplementary_Table_3_MASTER_metadata.tsv).  


Due to the complexity of this organization, **the main steps and general code to run each process is explained below**, while specific code and examples are given in their dedicated sections.

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

For resistome analysis, filtered reads were aligned versus the ResFinder database (<https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/all.fsa>), using bowtie2 (<https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html>) with *--very-sensitive* *--end-to-end* parameters by:

```         
ruby 01.bowtie2_resfinder folder_filtered_reads
```

The folder_filtered_reads must contain a folder per sample (originated from <https://github.com/SegataLab/MASTER-WP5-pipelines/tree/master/02-Preprocessing>) with 3 fastq.bz2 files inside (R1,R2,unpaired).

The obtained sam files were filtered by an in-house ruby script read_counts_v4.rb, which removes the gene over-estimation occurring when forward and reverse reads are aligned with the same gene:

```         
ruby read_counts_v4.rb sam_files_folder
```

To run correctly this script, you need to build previously the gene_list.txt file (in our case, a modified version of Supplementary Table 2, removing first column and reordering the other columns), with contains 4 columns:

-   Family: antibiotic family to which the ARG confer resistance to

-   Genefam: gene family build by CD-HIT alignment at 90% identity of
    ARG from ResFinder (and manually inspected to give the name
    according to the most abundant gene within each cluster or by the
    first number of the cluster within bla and addA genes)

-   Gene: gene name, removing everything after "\_" from the hit

-   Hit: hit or gene name according to ResFinder database

To obtain the *Family*, *Gene* and *Hit* columns information, we recommend to download the *phenotypes.txt* and *all.fsa* files from ResFinder (<https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/>), to run CD-HIT and following script:

```         
cdhit -i all.fsa -o cdhit.txt -c 0.9
ruby extract_cluster_table.rb phenotype_file.txt cdhit_clstr_file.clstr
```
Where in that case *phenotype_file.txt* is the *phenotype.txt* file from ResFinder datababase repository (*phenotype_20250421.txt* in the attached example); and *cdhit_clstr_file* is the *.clstr* output from cd-hit.

You must curate (and add) the genefam column manually according to cluster and gene columns resulted.



The obtained counts matrix was processed by R-scripts to calculate the counts per million reads (CPM) adding a bacterial marker modification according to the formula:

CPM = (target_genes_reads \* 10\^6 \*Bacterial_Markers_alignment)/total_reads

where CPM is the total counts per million reads value for each gene; 
target_genes_reads are the number of reads that match with the target genes; 
total_reads is the total number of reads obtained on each sample;
and Bacterial_Markers_alignment is the value obtained from viromeQC, which was employed using the script 01.viromeqc.rb by:

```         
ruby 01.viromeqc.rb folder_filtered_reads 
```

The *Bacterial Markers alignment* value indicates the proportion of bacterial DNA on a metagenomic fastq file, thus applying this parameter on the equation we remove those reads not assigned to bacterial taxa.
<br>

[Back to index](#idx)

<br>

## 4. Assembly-based analysis<a name="id4"></a>

From the QC reads, an assembly and was performed as thorughlly described in the corresponding section of [MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines/tree/master/05-Assembly_pipeline).  

Some manual edition of the [main assembly pipeline](https://github.com/SegataLab/MASTER-WP5-pipelines/blob/master/05-Assembly_pipeline/pipeline_assembly.sh) are required to redirect the variables to the location of the data and your installed software in your system, and then simply run:

```
/path/to/pipeline_assembly.sh path/to/${DATASET}/reads/
```

This pipeline will:
1) Run the `run_single_assembly.sh` script to perform individual assemblies of the different QC reads from samples ion `path/to/${DATASET}/reads/`
2) Filter contigs to a mininimum length and store them as `${SAMPLE}_filtered_contigs.fa` in `/path/to/${DATASET}/contigs/${SAMPLE}/`
3) Align the QC reads against the metagenome assemblies using `bowtie2` and `samtools`.
4) Find contig depths by using `jgi_summarize_bam_contig_depths` from `metabat2`. 
5) Binning: compact contigs into bins by using `metabat2`.
6) Verify completeness and contamination of the reconstructed bins and choose which are adequate MAGs by using `checkm2`. MAGs that overcome the process will be stored in: `/path/to/${DATASET}/mags/${SAMPLE}/`

<br>

The first step of donwstream analysis involve including the resulting metagenomes (contigs, `/path/to/${DATASET}/contigs/${SAMPLE}/`) and MAGs (`/path/to/${DATASET}/mags/${SAMPLE}/`) in FASTA file to the [TORMES pipeline](https://github.com/nmquijada/tormes) for CDS prediction and annotation, and AMR and virulence screening.  

Even though [TORMES](https://github.com/nmquijada/tormes) was devised for WGS of single bacterial isolates, it also accepts draft assemblies and/or MAGs with some minor rearrangements. The steps followed are described below and the main improvements are being used to develop TORMES v2.  

### Running TORMES for downstream analysis

1. Define your variables

```bash
WORKDIR=# define your working directory
# Then define the "input directory" containing the contigs FASTA files in different directories regarding the sample. e.g
INDIR=/path/to/${DATASET}/
```

<br> 

2. Create TORMES metadata. There's an useful tutorial on how to do it [here](https://github.com/nmquijada/tormes/wiki/Shortcut-to-generate-the-metadata-file-for-TORMES).  

<br>

For **contigs**:

```bash
for SAMPLE in $(ls ${INDIR}/contigs/); do
  echo -e "${SAMPLE}\tGENOME\t${INDIR}/contigs/${SAMPLE}/${SAMPLE}\tcontigs from ${SAMPLE} in ${DATASET}" >> ${WORKDIR}/tormes-metadata-contigs.txt
done
sed -i "1iSamples\tRead1\tRead2\tDescription" > ${WORKDIR}/tormes-metadata-contigs.txt
# further description and metadata fields can be used, if needed
```

For **MAGs**:

```bash
for SAMPLE in $(ls ${INDIR}/mags/); do
  for BINS in $(ls ${INDIR}/mags/${SAMPLE}/); do
    echo -e "${SAMPLE}\tGENOME\t${INDIR}/contigs/${SAMPLE}/${BIN}\MAGs from ${SAMPLE} in ${DATASET}" >> ${WORKDIR}/tormes-metadata-mags.txt
done
sed -i "1iSamples\tRead1\tRead2\tDescription" > ${WORKDIR}/tormes-metadata-mags.txt
```

<br>

3. Run TORMES
```bash
conda activate tormes-1.3.0
for DATATYPE in contigs mags; do
  tormes -m ${WORKDIR}/tormes-metadata-${DATATYPE}.txt -t ${NCPUS} --no_pangenome --gene_min_id 80 --gene_min_cov 80 --only_gene_prediction --prodigal_options "-p meta" --custom_genes_db "plasmidfinder" --min_contig_len 1000 -o ${WORKDIR}/MASTER-tormes-${DATATYPE}
done
```

Where `${NCPUS}` is the number of threads to use (default = 1)

This command will, among other functions:
1) Assign taxonomy to the different contigs by using `kraken2`. It is recommended to install the latest database from [Kraken2 main repository](https://benlangmead.github.io/aws-indexes/k2)
> Note: the flag `--use-names` from `kraken2` is not implemented in `tormes v1.3.0` but in the current beta version (private). Such flag can be added to the main code without disturbance.
2) Predict CDS by using `prodigal`.
3) Generate the AMR screening based on sequence alignment (`blastn`) of the contigs against 3 databases (ResFinder, CARD, Argannot). For this study we used RedFinder database, so the results will be stored in: `${WORKDIR}/MASTER-tormes-${DATATYPE}/antimicrobial_resistance_genes/resfinder/`
4) Generate virulence creening based on sequence alignment (`blastn`) of the contigs against the VFDB database. Results located in: `${WORKDIR}/MASTER-tormes-${DATATYPE}/virulence_genes/`
5) Search for potential plasmid replicons by sequence alignment (`blastn`) of the contigs against the PlasmidFinder database. Results located in: `${WORKDIR}/MASTER-tormes-${DATATYPE}/custom_genes_db/plasmidfinder/`

> The thresholds used for sequence similarity alignments are set to 80% identity and coverage, but the output can be simply filtered by using awk, for instance:

```
threshold=90
awk -v OFS="\t" -F "\t" -v threshold="$threshold" 'NR == 1 || ($10 > threshold)' ${WORKDIR}/MASTER-tormes-${DATATYPE}/antimicrobial_resistance_genes/resfinder/${SAMPLE}_resfinder.tab | awk -v OFS="\t" -F "\t" -v threshold="$threshold" 'NR == 1 || ($11 > threshold)'
```


<br>

[Back to index](#idx)

<br>


## 5. Mobile Genetic Elements<a name="id5"></a>

Different software were used for the investigation of MGE and the main steps are described below and with more detail [here](https://github.com/nmquijada/MASTER-food-production-resistome/blob/main/assembly-based_workflow/mobilome_workflow.md)  

<br>

### 5.1 Plasmids



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
platon --db ${PLATON_DB} --output ${WORKDIR}/MGE/platon -p ${SAMPLE} --threads 8 --meta /path/to/${DATASET}/contigs/${SAMPLE}/${SAMPLE}.fasta
```

The output will contain all those contigs that are potentially plasmidic in a file `${SAMPLE}.plasmid.fasta` located in `${WORKDIR}/MGE/platon`

<br>

#### geNomad

```bash
# Define database
GENOMAD_DB=/PATH/TO/genomad/genomad_db

# Run geNomad
genomad end-to-end --cleanup /path/to/${DATASET}/contigs/${SAMPLE}/${SAMPLE}.fasta ${WORKDIR}/MGE/genomad ${GENOMAD_DB}
```

The results in file `${SAMPLE}_plasmid_summary.tsv` will contain the ID of all those contigs that are potentially plasmidic.  
This file is stored in `${WORKDIR}/MGE/genomad/${SAMPLE}_summary/`

<br>

### 5.2 Integrons

```
mkdir ${WORKDIR}/MGE/integrons
integron_finder --local-max /path/to/${DATASET}/contigs/${SAMPLE}/${SAMPLE}.fasta --cpu 1 --outdir ${WORKDIR}/MGE/integrons
```

The results file `{WORKDIR}/MGE/integrons/${SAMPLE}/${SAMPLE}.summary` contains the screening for each contig, and different values for consider if they inlcude integrons or not (CALIN, In0, etc.).  
The results file `{WORKDIR}/MGE/integrons/${SAMPLE}/${SAMPLE}.integrons` contains a thorough description if replicons are found in certain positions within each contig.


<br>

[Back to index](#idx)

<br>

## 6. Visualization<a name="id6"></a>

<br>

[Back to index](#idx)

<br>


