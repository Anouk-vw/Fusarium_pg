# Genomes -> annotation
This folder contains everything to go from:
genomes 1-> graph 2-> Comparative Annotation 3-> OGs 4-> Pangenome 5-> 
effector profile

1. create graph
2. CAT .. 
	* one transcript per gene: scripts/gff_to_protein.py
3. Brocoli v1.1
	python /path/to/broccoli.py -dir brocoli_in/ -threads 24 
4. scripts/ 
	* get singletons: scripts/singletons_brocCol.py
	python singletons_brocCol.py brocolli_dir_step3/table_protein_names.txt brocoli_in extension > new_protein_singletons.txt
	* Get itteration_graph: scripts/itterate_counts.py
	python itterate_counts.py brocolli_dir_step3/table_protein_counts.txt new_protein_singletons.txt new_protein.fasta itterate output.pdf
	* Get fasta per OGs:  scripts/OGtoFasta.py
	+Run mafft, and raxml 
	python OGtoFasta.py brocolli_dir_step3/table_protein_names.txt brocoli_in path/to/cat/out 1 5
5. exonerate
	* get effector OGs: scripts/OGstoExonerate.py
	python cactus_complete/scripts/OGstoExonerate.py brocoli_results/dir_step3/table_OGs_protein_counts.txt brocoli_results/dir_step3/table_OGs_protein_names.txt Predicted_effectors.txt brocoli_results/unique_genes.txt prot_fasta/
	* get PAV in set of genomes: (not part of analysis)
	for f in $(ls /path/to/fasta_files) ;\
	do exonerate --model protein2genome --query effectors.fasta \
	--target $f --ryo '%ps' > $f-exonerate.txt ; done
