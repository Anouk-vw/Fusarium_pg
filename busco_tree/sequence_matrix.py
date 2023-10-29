#argv[1]='align dir'
#argv[2]='end of file' ('trimmed')

from sys import argv
import datetime
import os

begin_time = datetime.datetime.now()

start = 0
dict_seq = {}
parts = []
seq = []
for aln_f in [x for x in os.listdir(argv[1]) if x.endswith(argv[2])]:
    genome = 0
    with open('{}/{}'.format(argv[1], aln_f)) as alignment:
        for line in alignment:
            #if genome is 0, this is the first line
            if line.startswith('>') and genome == 0:
                #get genome and gene, init sequence list
                genome, gene = line.split('|')
                seq = []
            if line.startswith('>') and genome != 0:
                #if genome not is 0, a sequence has been parsed
                #create sequence
                sequence = ''.join(seq)
                try:
                    #append sequence and genome to dictionary
                    if len(sequence) > 50:
                        dict_seq[genome].append(sequence)
                except KeyError:
                    if len(sequence) > 50:
                        dict_seq[genome] = []
                        dict_seq[genome].append(sequence)
                #get genome and gene to proceed
                genome, gene = line.split('|')
                seq = []
                #if not '>' in line, create sequence
            if not line.startswith('>'):
                seq.append(line.strip())
        #do this also after parsing the file to include the final thing
        sequence = ''.join(seq)
        try:
        #append sequence and genome to dictionary
            if len(sequence) > 50:
                dict_seq[genome].append(sequence)
        except KeyError:
            if len(sequence) > 50:
                dict_seq[genome] = []
                dict_seq[genome].append(sequence)
         
        #determine sequence length and define start and stop of fasta-region
        #all lengths are similair, because seq is aligned
        if len(sequence) > 50:
            stop = len(sequence) + start
            parts.append([aln_f.split('_')[0], '=', str(start + 1), str(stop)])
            start += len(sequence)

print([(k, len(''.join(v))) for k,v in dict_seq.items()])

del dict_seq['>BUSCO_PHYCPC1']

with open('sequence_matrix.fasta', 'w') as seqmat:
    for k,v in dict_seq.items():
        seqmat.write(str(k))
        seqmat.write('\n')
        #print(len(''.join(v)))
        seq = ''.join(v)
        seqmat.write(seq)
        seqmat.write('\n')

with open('partitioning.txt', 'w') as partitioning:
    print('part', parts)
    for p in parts:
        partitioning.write('charset, ')
        partitioning.write(p[0])
        partitioning.write(' = ')
        partitioning.write(p[2])
        partitioning.write('-')
        partitioning.write(p[3])
        partitioning.write('\n')


#print((datetime.datetime.now() - begin_time))
print([(k,len(v)) for k,v in dict_seq.items()])

