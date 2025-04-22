# Assembly-free workflow

For resistome and virulome analysis, filtered reads were aligned versus
the ResFinder database
(<https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/all.fsa>),
using bowtie2 (<https://rnnh.github.io/bioinfo-notebook/docs/bowtie2.html>) with *--very-sensitive* *--end-to-end* parameters by

```         
ruby 01.bowtie2_resfinder folder_filtered_reads
```

The folder_filtered_reads must contain a folder per sample (originated
from
<https://github.com/SegataLab/MASTER-WP5-pipelines/tree/master/02-Preprocessing>)
with 3 fastq.bz2 files inside (R1,R2,unpaired).

The obtained sam files were filtered by an in-house ruby script
read_counts_v4.rb, which removes the gene over-estimation occurring when
forward and reverse reads are aligned with the same gene:

```         
ruby read_counts_v4.rb sam_files_folder
```

To run correctly this script, you need to build previously the
gene_list.txt file (in our case, a modified version of Supplementary
Table 2, removing first column and reordering the other columns), with
contains 4 columns:

-   Family: antibiotic family to which the ARG confer resistance to

-   Genefam: gene family build by CD-HIT alignment at 90% identity of
    ARG from ResFinder (and manually inspected to give the name
    according to the most abundant gene within each cluster or by the
    first number of the cluster within bla and addA genes)

-   Gene: gene name, removing everything after "\_" from the hit

-   Hit: hit or gene name according to ResFinder database

To obtain the *Family*, *Gene* and *Hit* columns information, we
recommend to download the *phenotypes.txt* and *all.fsa* files from
ResFinder
(<https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/>),
to run CD-HIT and following script:

```         
cdhit -i all.fsa -o cdhit.txt -c 0.9
ruby extract_cluster_table.rb
```

You must curate (and add) the genefam column manually according to
cluster and gene columns resulted.



The obtained counts matrix was processed by R-scripts to calculate the
counts per million reads (CPM) adding a bacterial marker modification
according to the formula:

CPM = (target_genes_reads \* 10\^6 \*
Bacterial_Markers_alignment)/total_reads

where CPM is the total counts per million reads value for each gene;
target_genes_reads are the number of reads that match with the target
genes; total_reads is the total number of reads obtained on each sample;
and Bacterial_Markers_alignment is the value obtained from viromeQC,
which was employed using the script 01.viromeqc.rb by:

```         
ruby 01.viromeqc.rb folder_filtered_reads 
```

The *Bacterial Markers alignment* value indicates the proportion of
bacterial DNA on a metagenomic fastq file, thus applying this parameter
on the equation we remove those reads not assigned to bacterial taxa.

An example of how to calculate this CPM using R-script is located on
"calculate_cpm" folder, which contains *calculate_cpm.R*,
*matrix_ResFinder.txt* and *viromeqc_all.txt* files.
