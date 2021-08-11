#!/usr/bin/env python
import argparse
import re

parser = argparse.ArgumentParser(description=f"Split fastq files based on read sequences.\n"
                                             f"Please make sure Biopython is installed:\n"
                                             f"pip install biopython")
parser.add_argument('fastq_file')
parser.add_argument('-b', metavar='--barcode', help='Barcode to split on, multiple barcodes can be given by comma separation.')
parser.add_argument('-bf', metavar='--barcode_file', help='Fasta file containing one or more barcodes to split on.')
parser.add_argument('-o', metavar='--output', help='Output directory for the splitted fastq files', default='')
parser.add_argument('-mm', metavar='--miss-match-cutoff', help='Maximum number of miss matches for barcode matching', default=1, type=int)
parser.add_argument('-ms', metavar='--minimum-sequence-length', help='Minimum sequence length when using sequence matching for the barcodes',
                    default=6, type=int)
parser.add_argument('-seq-hit', default=False, action='store_true', help='Instead of having maximum mis matches cutoff, use a minimum sequantial '
                                                                         'hit for identification of the presence of a tag.')

parser.add_argument('-em', default=False, action='store_true', help='Report extended metrics, eg. barcode match error rates'
                                                                                                  'and mismatch locations.')
parser.add_argument('-el', metavar='--edge-length', help='Minimum edge length to qualify for a barcode hit on the edges.',
                    default=0, type=int)
parser.add_argument('-nr', default=False, action='store_true', help='Do not remove barcode from the read.')
parser.add_argument('-kd', default=False, action='store_true',
                    help='If multiple barcodes are found in a read also put them in the individual fastq files.')
args = parser.parse_args()

from Bio import SeqIO
from collections import defaultdict

# Check for barcode file or barcodes given as arguments
barcodes = {}
if args.b:
    barcodes += dict(zip(args.b.split(','), range(0, len(args.b.split(',')))))
if args.bf:
    for record in SeqIO.parse(args.bf, 'fasta'):
        barcodes[record.id] = str(record.seq)
fastq_name = args.fastq_file.split('/')[-1]

# Parse output naming and location
if args.o and not args.o.endswith('/'):
    args.o += '/'
# Open file handlers for the output fastq files
open_files = {
    'nobar': open(args.o+f"no_barcode-{fastq_name}", 'w')
}
for barcode in barcodes.keys():
    open_files[barcode] = open(args.o+f"bar_{barcode}-{fastq_name}", 'w')

open_files['all'] = open(args.o+f"all-{fastq_name}", 'w')


# Append read to fastq file
def append_fastq(seq, path):
    SeqIO.write(seq, path, 'fastq')

def find_tag(barcode, read, phreds, min_seq, metrics=False):
    bar_positions = []
    e_metrics = []
    for i in range(0, len(read)-len(barcode)):
        total_identical = 0
        longest_sequence = 0
        current_sequence = 0
        metric_arr = []
        for char_idx, (charB, charR, phred) in enumerate(zip(barcode, read[i:i+len(barcode)], phreds[i:i+len(barcode)])):
            if charB == charR:
                total_identical += 1
                current_sequence += 1
            else:
                if metrics:
                    metric_arr.append((i+char_idx, charR, charB, phred))
                longest_sequence = max(longest_sequence, current_sequence)
                current_sequence = 0
        if longest_sequence >= min_seq:
            bar_positions.append(i)
            if metrics and metric_arr:
                e_metrics.append(metric_arr)
    if metrics:
        return bar_positions, e_metrics
    return bar_positions


def euclidian_dist(barcode, read, phreds, mm_threshold, metrics=False, edge_size=6):
    bar_positions = []
    e_metrics = []
    for i in range(0, len(read)-len(barcode)):
        euclid_dist = 0
        metric_arr = []
        for char_idx, (charB, charR, phred) in enumerate(zip(barcode, read[i:i+len(barcode)], phreds[i:i+len(barcode)])):
            if charB != charR:
                euclid_dist += 1
                if metrics:
                    metric_arr.append((i+char_idx, charR, charB, phred))
            if euclid_dist > mm_threshold: break
        if euclid_dist <= mm_threshold:
            bar_positions.append(i)
            if metrics and metric_arr:
                e_metrics.append(metric_arr)
    if metrics:
        return bar_positions, e_metrics
    return bar_positions

def check_edges(barcode, read, edge_size, mm_threshold):
    b_start = barcode[:edge_size]
    bar_start_match = list(re.finditer(b_start, read))
    if bar_start_match:
        bar_start_match = bar_start_match[-1].start()
    else:
        return False
    euclid_dist = 0
    for charB, charR in zip(barcode[edge_size:], read[bar_start_match+edge_size:]):
        if charB != charR:
            euclid_dist += 1
        if euclid_dist > mm_threshold: break
    if euclid_dist <= mm_threshold:
        print(list(re.finditer(b_start, read)))
        print(read)
        print(barcode[edge_size:], read[bar_start_match+edge_size:])
        print("yeah!!")
        exit()


tot_barcodes = 0
barcode_metrics = defaultdict(int)
if args.em:
    import csv
    csv_fh = open('csv_test.csv', 'w')
    csv_writer = csv.writer(csv_fh, delimiter=',')
    csv_writer.writerow(['read', 'barcode', 'mut_pos', 'ref', 'mut', 'phred'])

tot_reads = 0
tot_reads_with_barcode = 0
# Loop through records in fastq file
for record in SeqIO.parse(args.fastq_file, 'fastq'):
    tot_reads += 1
    barcode_hits = []
    record_id = record.id
    record_cnt_barcodes = 0
    for bar_id, barcode in barcodes.items():
        if args.em:
            if args.seq_hit:
                pos_list, e_metrics = find_tag(barcode, str(record.seq), record.letter_annotations['phred_quality'], args.ms, args.em)
            else:
                pos_list, e_metrics = euclidian_dist(barcode, str(record.seq), record.letter_annotations['phred_quality'], args.mm, args.em)
            for mut in e_metrics:
                if len(mut) > 1:
                    mut = "|".join([str(m[0]) for m in mut]), "|".join([str(m[1]) for m in mut]), "|".join([str(m[2]) for m in mut]), \
                          "|".join([str(m[3]) for m in mut])
                else:
                    mut = mut[0]
                csv_writer.writerow([record_id, bar_id, mut[0], mut[1], mut[2], mut[3]])
        else:
            pos_list = euclidian_dist(barcode, str(record.seq), record.letter_annotations['phred_quality'], args.mm)
        if args.el:
            edge_hits = check_edges(barcode, str(record.seq), args.el, args.mm)
        tot_barcodes += len(pos_list)
        record_cnt_barcodes += len(pos_list)
        if pos_list:
            barcode_hits.append(bar_id)
            if len(pos_list) > 1:
                barcode_metrics['dup'] += 1
            if not args.nr:
                removed = 0
                for idx, pos in enumerate(pos_list):
                    record = record[:pos-removed] + record[(pos+len(barcode))-removed:]
                    pos_list[idx] = pos_list[idx] - removed
                    removed += len(barcode)
            record.id = f"{bar_id}:{','.join(map(str, pos_list))}_{record.id}"
    append_fastq(record, open_files['all'])
    if len(barcode_hits) == 0:
        append_fastq(record, open_files['nobar'])
    elif len(barcode_hits) == 1:
        tot_reads_with_barcode += 1
        append_fastq(record, open_files[barcode_hits[0]])
    elif len(barcode_hits) > 1:
        tot_reads_with_barcode += 1
        barcode_metrics['mult'] += 1
        if args.kd:
            for bar_id in barcode_hits:
                append_fastq(record, open_files[bar_id])
        bar_ids = tuple(barcode_hits)
        try:
            append_fastq(record, open_files[bar_ids])
        except KeyError:
            open_files[bar_ids] = open(args.o+f"bars_{'_'.join(bar_ids)}-{fastq_name}", 'w')
            append_fastq(record, open_files[bar_ids])
print(f"Total reads: {tot_reads}\nReads with barcode: {tot_reads_with_barcode}\n"
      f"Reads without barcode: {tot_reads-tot_reads_with_barcode}\nTotal barcodes: {tot_barcodes}")
print(barcode_metrics)
for file in open_files.values():
    file.close()
if args.em:
    csv_fh.close()
