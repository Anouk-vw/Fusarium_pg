# nucmer
Snakefile to run nucmer all-vs-all on a set of genomes.
Genomes specifies in genome_names.txt with name and path.
output is a per genome table with the coverage of all other genomes per window

Ref         TR4      x       y
0 1000      1000    800     800
1000 2000   1000    800     1000

# scripts

Scripts for downstream analysis

- analyse_ncumer.py:
  Script to analyse nucmer output and produce all_vs_all dataframe
- AR_similarity:
  Script to go from nucmer_directory -> pairwise similairty analysis
