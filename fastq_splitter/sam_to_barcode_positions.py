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
# parser.add_argument('-bf', metavar='--barcode_file', help='Fasta file containing one or more barcodes to split on.')
# parser.add_argument('-r', default=False, action='store_true', help='Removes barcode from the read.')
# parser.add_argument('-kd', default=False, action='store_true',
#                     help='If multiple barcodes are found in a read also put them in the individual fastq files.')
args = parser.parse_args()

sam_file_name = args.sam_file.split('/')[-1].split('.')[0]
if args.o:
    csv_name = args.o.strip('.csv')
else:
    csv_name = sam_file_name

try:
    import pysam
    from Bio import SeqIO
    from collections import defaultdict
    import numpy as np
    import json
    import pandas as pd
except ImportError as e:
    pip_dict = {
        'Bio': 'biopython'
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

for read in samfile:
    if read.is_unmapped: continue
    q_q = read.query_qualities
    pairs = read.get_aligned_pairs(matches_only=True)
    ref_positions = [pair[1] for pair in pairs if q_q[pair[0]] >= args.qt]
    ref_trans[read.reference_name][ref_positions] += 1
    found_barcodes = read.query_name.split('_')[:-1]
    for barcode in found_barcodes:
        barcode_id, bar_positions = barcode.split(':')
        try:
            position_dict[read.reference_name][barcode_id.split(':')[0]][[int(pos)+read.reference_start for pos in bar_positions.split(',')]] += 1
        except KeyError:
            position_dict[read.reference_name][barcode_id.split(':')[0]] = np.zeros(len(ref_trans[read.reference_name]))
            position_dict[read.reference_name][barcode_id.split(':')[0]][[int(pos)+read.reference_start for pos in bar_positions.split(',')]] += 1

for trans_id, sequence in ref_sequence.items():
    cov = ref_trans[trans_id]
    barcodes = position_dict[trans_id]
    csv_dict = {trans_id: list(sequence), 'tot_coveraeg': list(cov)}
    for barcode, values in barcodes.items():
        csv_dict[barcode] = list(values)
    df = pd.DataFrame(csv_dict)
    df.to_csv(f"{csv_name}_{trans_id}.csv")
