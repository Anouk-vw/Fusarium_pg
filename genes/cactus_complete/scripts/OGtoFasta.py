"""
script to read brocoli OG table and create fasta file per OG
"""
import utils
from sys import argv
from Bio import SeqIO
import os
import subprocess
from Bio.Phylo.PAML import codeml
from Bio.Align.Applications import MafftCommandline
from itertools import groupby

def read_table(table, subset = False, range_OG = [0]):
	"""
	Read OG table and return dict {OG : gene names}
	"""
	OG_dict = utils.OGnames_dict(table)
	if subset == True:
		#for x in range(int(range_OG[0]), int(range_OG[1])):
		#	OG_name = f'OG_{x}'
		#	print(OG_name)
		OG_dict = {f'OG_{i}':OG_dict[f'OG_{i}'] for i in range(int(range_OG[0]), int(range_OG[1]))}
	#subset genomes:
	OG_dict = {k:v for k,v in OG_dict.items() if sum("TR4_II5" in vg for vg in v) == 1 and sum("CR1.1" in vg for vg in v) == 1} 
	return OG_dict

def dict_to_fasta(OG_name, gene_list, cat_dir, ext):
	"""
	Write fasta file with genes / OG
	"""
	print('fetch fasta records')
	print(f'working in {OG_name}')
	OG_records_dict = {}
	#for k,v in zip(OG, gene_list):
	records_list = []
	for genes in gene_list:
		prefix = genes.split('_T')[0]
		#if prefix == '15_KTG06-B':
		#	fasta_loc = f'{cat_dir}/{prefix}2{ext}'
		if prefix == 'TR4_II5' or prefix == 'CR1.1':
			fasta_loc = f'{cat_dir}/{prefix}{ext}'
			with open(fasta_loc) as temp_fas:
				record = [r for r in SeqIO.parse(temp_fas, "fasta") if r.id in genes]
			records_list.extend(record)  # use any record ID
	OG_records_dict[OG_name] = records_list
	return OG_records_dict

def write_fasta(broc_table, CAT_dir, nuc_dir, subset, range_list):
	"""
	From brocoli input to prot and nuc fasta files
	"""
	OGs = []
	print('writing fasta')
	OGgenes_dict = read_table(broc_table, subset, range_OG = range_list)
	for OG_id,genes in OGgenes_dict.items():
		#does prot exist?
		name = f"prot_fasta/{OG_id}_prot.fasta"
		if os.path.isfile(name):
			print(f"{name}_exists")
		else:
		#no fetch records and write
			prot_records_dict = dict_to_fasta(OG_id,genes, CAT_dir, '.new_protein.fasta')
			with open(name, "w+") as output_handle:
				for x in prot_records_dict[OG_id]:
					SeqIO.write(x, output_handle,"fasta")
		#does nuc exists?
		name = f"nuc_fastas/{OG_id}_nuc.fasta"
		if os.path.isfile(name):
			print(f"{name}_exists")
		else:
		#no fetch records and write
			nuc_records_dict = dict_to_fasta(OG_id,genes, nuc_dir, '.consensus.fasta')
			with open(name, "w+") as output_handle:
					for x in nuc_records_dict[OG_id]:
						SeqIO.write(x, output_handle,"fasta")
		OGs.append(OG_id)
	return OGs

def align(fasta, prefix):
	"""align fasta"""
	out = f'{prefix}.aln'
	if os.path.isfile(out):
		print('mafft exists')
	else:
		cmd = ['mafft', '--thread 5', fasta, '>', out]
		print(f'running mafft as {" ".join(cmd)}')
		maff = MafftCommandline(input=fasta)
		stdout, stderr = maff()
		with open(out, "w") as handle:
			handle.write(stdout)
		#subprocess.run(cmd, shell = True)
	return

def RAXML(alignment, out):
	outdir = 'trees'
	if os.path.isfile(f'{outdir}/RAxML_bestTree.{out}'):
		print(f'{out} exists')
	else:
		print('creating tree')
		cmd = ['raxmlHPC', '-f', 'a', '-m', 'PROTGAMMABLOSUM62', '-s', \
		alignment, '-n', out, '-x', '12', '-N', '5', '-T', '10', '-w', \
		'home/anouk/anouk2/pangenome/cactus_complete/genes/trees', \
		'-p', '12']
		print(' '.join(cmd))
		subprocess.run(cmd)
	return

def ral2nal(alignment, nuc_fasta, out):
	print('running pal2nal')
	if os.path.isfile(out):
		print('pal2nal exists')
	else:
		cmd = ['bash', 'pal2nal.sh', alignment, nuc_fasta, out]
		#['perl', 'align/pal2nal.v14/pal2nal.pl', alignment, nuc_fasta, \
		#'-output', 'paml', '-nogap', '>', out]
		print(' '.join(cmd))
		pn = subprocess.Popen(cmd, shell=True)
		print(f'finished pal2nal with {pn.communicate()[1]}')
	return 

def run_codeml(ral2nalout, tree, out):
	if os.path.isfile(out):
		print(f'{out}_exists')
	else:
		print('running codeml')
		cml = codeml.Codeml(
		alignment = ral2nalout,
		tree=tree,
		out_file=out,
		working_dir="./",
		)
		cml.read_ctl_file('align/test.cnt')
		cml.run(verbose = True)
	return

def parse_codeML(codeML_in):
	results = codeml.read(codeML_in)
	omega = results['NSsites'][0]['parameters'] #[v.items() for v in results['pairwise'].values()]])#([k for k,v in results.items()]))
	return omega.values()

def fasta_length(fasta_file):
	"""
	return average length of fasta file accessions
	"""
	tot = 0
	num = 0
	with open(fasta_file) as handle:
		for header, group in groupby(handle, lambda x:x.startswith('>')):
			if not header:
				num += 1
				tot += sum(map(lambda x: len(x.strip()), group))
	result = float(tot)/num
	return result

if __name__ == '__main__':
	#argv 3 = subset range
	if argv[1] == 'stats':
		dnds_dict = {}
		length_dict = {}
		print('stats')
		for OG in [x.split('_prot')[0] for x in os.listdir(argv[2]) if '_prot' in x]:
			#codeml_out = parse_codeML(f'codeml/{OG}.codeml')
			#dnds_dict[OG] = codeml_out
			length_dict[OG] = fasta_length(f'{OG}_prot.fasta')
			print(length_dict)
	else:
		gene_names = argv[1] #broc_table_names
		CAT = argv[2] #new_protein_files/
		nuc_dir = argv[3] #CAT/consensus
		#read_table(gene_names)
		written_OGs = write_fasta(gene_names, CAT, nuc_dir, False, [0,0])#, True, [argv[4], argv[5]])
		#written_OGs contains a list of fasta files in gene_names
		for OG in written_OGs:
			#start with alignments only, convinient will be used in another step
			align(f'prot_fasta/{OG}_prot.fasta', f'align/{OG}')
		for OG in written_OGs:
			RAXML(f'align/{OG}.aln', f'{OG}_tree')
		"""
			#ral2nal(f'align/{OG}.aln', f'{OG}_nuc.fasta', f'align/ral2nal/{OG}.codon')
			ral2nal(f'align/{OG}.aln', f'{OG}_nuc.fasta', f'align/ral2nal/{OG}.codon')
			run_codeml(f'align/ral2nal/{OG}.codon', f'trees/RAxML_bestTree.{OG}_tree', f'codeml/{OG}.codeml')
			codeml_out = parse_codeML(f'codeml/{OG}.codeml')
			dnds_dict[OG] = codeml_out
		"""

