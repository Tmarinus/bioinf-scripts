#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description=f"Parse a sam file that contains barcode positions in its description.\n"
                                             f"Please make sure pysam is installed:\n"
                                             f"pip install pysam\n"
                                             f"Please make sure Biopython is installed:\n"
                                             f"pip install biopython\n"
                                             f"Please make sure numpy is installed:\n"
                                             f"pip install numpy")
parser.add_argument('sam_file')
parser.add_argument('ref_fasta')

parser.add_argument('-qt', metavar='--quality-threshold', help='Threshold for counting a nucleotide (used for coverage). Default=0, cause currently only counting total coverage', default=0)
parser.add_argument('-o', metavar='--output', help='Output name for the csv', default='')
args = parser.parse_args()

sam_file_name = args.sam_file.split('/')[-1].split('.')[0]
if args.o:
    csv_name = args.o.strip('.csv')
    if not csv_name.endswith('/'):
        csv_name += '_'

else:
    csv_name = sam_file_name

# Import modules
try:
    import pysam
    from Bio import SeqIO
    from collections import defaultdict
    import numpy as np
    import json
    import pandas as pd
except ImportError as e:
    pip_dict = {
        'Bio': 'biopython',
    }
    mis_package = str(e).split('\'')[-2].strip("'")
    try:
        install_name = pip_dict[mis_package]
    except KeyError:
        install_name = mis_package
    print(f"Package '{mis_package}' was not found, please install by running:\n"
          f"pip install {install_name}")
    exit()

samfile = pysam.AlignmentFile(args.sam_file)

ref_trans = {}
position_dict = {}
ref_sequence = {}
for record in SeqIO.parse(args.ref_fasta, 'fasta'):
    ref_trans[record.id] = np.zeros(len(record))
    ref_sequence[record.id] = str(record.seq)
    position_dict[record.id] = {}

barcode_reads = defaultdict(list)
ins = 0
dels = 0
corr_barcode = 0
wildtype = 0
unmapped = 0
total_reads = 0

def adjust_bar_pos(cigar, bar_pos, ref_start):
    curr_pos = 0
    adjustment = 0
    for cig in cigar:
        if cig[0] == 0: # no mutations
            if curr_pos + cig[1] > bar_pos:
                return ref_start + bar_pos + adjustment
            curr_pos += cig[1]
        elif cig[0] == 1: # ins in ref
            adjustment -= cig[1]
        elif cig[0] == 1: # del in ref
            adjustment += cig[1]
    return ref_start + bar_pos + adjustment



cnt_strange_start = 0
cnt_strange_start_sequence = 0
for read in samfile:
    total_reads += 1
    if read.is_unmapped:
        unmapped += 1
        continue
    q_q = read.query_qualities
    pairs = read.get_aligned_pairs(matches_only=True)
    ref_positions = [pair[1] for pair in pairs if q_q[pair[0]] >= args.qt]
    ref_trans[read.reference_name][ref_positions] += 1
    found_barcodes = read.query_name.split('_')[:-1]
    if read.reference_start != 4838 and read.query_alignment_sequence.startswith('CTATGAAA'):
        cnt_strange_start += 1
        if 'gccaata'.upper() in read.get_tag('MD'):
            cnt_strange_start_sequence += 1
    if not found_barcodes:
        continue
    for barcode in found_barcodes:
        barcode_id, bar_positions = barcode.split(':')
        read_positions = [int(pos) for pos in bar_positions.split(',')]
        mapped_positions = [int(pos)+read.reference_start for pos in bar_positions.split(',')]
        tuples_cigar = [x[0] for x in read.cigartuples]
        if 1 in tuples_cigar:
            ins += 1
        if 2 in tuples_cigar:
            dels += 1
        barcode_reads[barcode_id.split(':')[0]].append((read.query_name.split('_')[-1], read_positions, mapped_positions))
        adjusted_bar_pos = [adjust_bar_pos(read.cigartuples, y, read.reference_start) for y in [int(x) for x in bar_positions.split(',')]]

        try:
            position_dict[read.reference_name][barcode_id.split(':')[0]][adjusted_bar_pos] += 1
        except KeyError:
            position_dict[read.reference_name][barcode_id.split(':')[0]] = np.zeros(len(ref_trans[read.reference_name]))
            position_dict[read.reference_name][barcode_id.split(':')[0]][adjusted_bar_pos] += 1

print(f'total reads: {total_reads}')
print(f'Alignment not starting at 4838 even if seq starts with CTATGAAA: {cnt_strange_start}\n'
      f'These having TCAACCACAATGTGTTGTCCGGCCAATACTTGTTGCAA: {cnt_strange_start_sequence} {cnt_strange_start_sequence/cnt_strange_start*100:0.1f}')
print(f"Totals reads with an insert: {ins}")
print(f"Total reads with an deletion: {dels}")
orf_corr = 0
orf_false = 0
for trans_id, sequence in ref_sequence.items():
    cov = ref_trans[trans_id]
    barcodes = position_dict[trans_id]
    csv_dict = {trans_id: list(sequence), 'tot_coveraeg': list(cov)}
    for barcode, values in barcodes.items():
        for idx, val in enumerate(values):
            if idx < 4869: continue
            if (idx - 4869) % 3 == 0:
                orf_corr += val
            else:
                orf_false += val
        csv_dict[barcode] = list(values)
    df = pd.DataFrame(csv_dict)
    df.to_csv(f"{csv_name}{trans_id}.csv")

print(f"\nORF positions of the barcodes:\n"
      f"Correct: {orf_corr:,.0f} {(orf_corr/(orf_corr+orf_false))*100:.3f}\n"
      f"Incorrect: {orf_false:,.0f} {(orf_false/(orf_corr+orf_false))*100:.3f}")

for barcode, values in barcode_reads.items():
    read_ids, read_positions, mapped_positions = [], [], []
    for v in values:
        read_ids.append(v[0])
        read_positions.append(":".join([str(x) for x in v[1]]))
        mapped_positions.append(":".join([str(x) for x in v[2]]))
    csv_dict = {'read_ids': read_ids, 'read_positions': read_positions, 'mapped_positions': mapped_positions}
    df = pd.DataFrame(csv_dict)
    df.to_csv(f"{csv_name}{barcode}_mapping.csv")


print(f"\nstored as: {csv_name}{trans_id}.csv")