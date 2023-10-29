# nucmer
Snakefile to run nucmer all-vs-all on a set of genomes.
Genomes specifies in genome_names.txt with name and path.
output is a per genome table with the coverage of all other genomes per window

Ref         TR4      x       y
0 1000      1000    800     800
1000 2000   1000    800     1000

----
Notebook:

Get_acc.ipynb 
can create accessory regions into a bed file (make script. And plot an heatmap of mean number of synteny blocks.

Acc_circos 
can use these bed files to create a circos plot with links between accessory regions.
