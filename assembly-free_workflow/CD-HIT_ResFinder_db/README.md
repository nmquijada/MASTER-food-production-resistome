The submitted *cdhit.txt* and *phenotypes_20250421.txt* were not from the same version of ResFinder database (*cdhit.txt* is the *all.fsa* version from the manuscript XXXXXX).

To prove the pipeline we recommend to use the last version for both *phenotypes.txt* and *all.fsa* files (https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) and take care about the possibility of some missing genes on phenotyhpes.txt file. In this example, both *phenotypes_20250421.txt* and *all_20250421.fsa* were download the same day, as it is indicated on their names.

To run *extract_cluster_table.rb* script use: 

```
ruby extract_cluster_table.rb phenotype_file.txt cdhit_clstr_file.clstr
```
where in that case *phenotype_file.txt* is the phenotype.txt file from ResFinder datababase repository (*phenotype_20250421.txt*); and cdhit_clstr_file is the *.clstr* output from cd-hit.

Now, a new column must be manually included (*genefam*) according to *cluster* and *gene* columns. Enjoy.
