The submitted *cdhit.txt* and *phenotypes_20250421.txt* were not from the same version of ResFinder database (*cdhit.txt* is the *all.fsa* version from the manuscript XXXXXX).

To prove the pipeline we recommend to use the last version for both *phenotypes.txt* and *all.fsa* files (https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) and take care about the possibility of some missing genes on phenotyhpes.txt file. In this example, both *phenotypes_20250421.txt* and *all_20250421.fsa* were download the same day, as it is indicated on their names.

After running the *extract_cluster_table.rb* script, a new column must be manually included (*genefam*) according to *cluster* and *gene* columns. Enjoy.
