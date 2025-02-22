from Bio import SeqIO


input_fasta = 'genome.fa'
output_fasta = 'genome-nochr.fa'

contig_rename_dict = {}
for i in range(1, 23):
    contig_rename_dict[f'{i}'] = f'chr{i}'
for c in ['X', 'Y', 'M']:
    contig_rename_dict[f'{c}'] = f'chr{c}'

with open(input_fasta, 'r') as infile:
    with open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            print(record.id)
            if record.id in contig_rename_dict:
                record.id = contig_rename_dict[record.id]
            SeqIO.write(record, outfile, 'fasta')
