#!/usr/bin/env python

import sys
import argparse

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print(f'Script requires BioPython module to be installed run\npip3 install BioPython --user')
    exit()

parser = argparse.ArgumentParser(usage=f'To run give paths to vcf file, fasta and output\n'
           f'eg: vcf_to_fasta test_files/barcode3_indels.vcf test_files/template_gal_k2.fasta output.fasta\n',
           description=f'Script requires BioPython module to be installed run\npip3 install BioPython --user')

parser.add_argument('vcf_file', help='path to vcf file')
parser.add_argument('fasta_file', help='path to fasta file')
parser.add_argument('output', help='path and file name for output')

args = parser.parse_args()

vcf_file = args.vcf_file
fasta_file = args.fasta_file
output = args.output

fasta_dict = {}
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        fasta_dict[record.id] = record

mut_dict ={}
removed = []
seq_id = None
with open(vcf_file, 'r') as vcf_reader:
    for line in vcf_reader:
        if line.startswith('#'): continue
        splits = line.split()
        if seq_id is None:
            seq_id = splits[0]
        if seq_id != splits[0]:
            print(f'Error expected only a single sequence id in the VCF file\n{repr(seq_id)}\n{repr(splits[0])}')
            exit()
        pos = int(splits[1])
        ref = splits[3]
        alt = splits[4]
        af = float(splits[-1].split(';')[1].split('=')[1])
        if pos in mut_dict.keys():
            if mut_dict[pos][2] < af:
                removed.append((pos, mut_dict[pos], (ref, alt, af)))
                mut_dict[pos] = (ref, alt, af)
        else:
            mut_dict[pos] = (ref, alt, af)
fasta_array = list(fasta_dict[seq_id].seq)
for pos, values in mut_dict.items():
    fasta_array[pos-1] = values[1]
new_str = ''.join(fasta_array)
new_fasta = SeqRecord(Seq(new_str), seq_id)
with open(output, "w") as output_handle:
    SeqIO.write(new_fasta, output_handle, "fasta")

for lost in removed:
    print(f'While running the following positions had multiple mutations, highest freq inserted in fasta:\n'
          f'Pos: {lost[0]} ref: {lost[1][0]} alt1: {lost[1][1]} freq1 {lost[1][2]} '
          f'alt2: {lost[2][1]} freq2 {lost[2][2]}')