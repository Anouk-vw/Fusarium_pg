# nucmer
Snakefile to run nucmer all-vs-all on a set of genomes.
Genomes specifies in genome_names.txt with name and path.
output is a per genome table with the coverage of all other genomes per window

Ref         TR4      x       y
0 1000      1000    800     800
1000 2000   1000    800     1000

scripts/analyse_ncumer.py:

script to analyse nucmer output and produce all_vs_all dataframe

