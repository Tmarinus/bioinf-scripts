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
parser.add_argument('-seq-hit', default=False, action='store_true', help='Instead of having maximum mis matches cutoff, use a minimum sequential '
                                                                         'hit for identification of the presence of a tag.')

parser.add_argument('-em', default=False, action='store_true', help='Report extended metrics, eg. barcode match error rates'
                                                                                                  'and mismatch locations.')
parser.add_argument('-el', metavar='--edge-length', help='Minimum edge length to qualify for a barcode hit on the edges.',
                    default=0, type=int)
parser.add_argument('-qc', metavar='--quality-cutoff', default=0, help='If phred score is below this score do not count as mismatch. Default off (0)')
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

def extend_match(sub_bar, sub_read, sub_phreds, max_mm):
    euclid_dist = 0
    metric_arr = []
    for char_idx, (charB, charR, phred) in enumerate(zip(sub_bar, sub_read, sub_phreds)):
        if charB != charR:
            euclid_dist += 1
            if args.em:
                metric_arr.append((char_idx, charR, charB, phred))
        if euclid_dist > max_mm:
            break
    if euclid_dist <= max_mm:
        return True, metric_arr
    return False, []


def find_tag(barcode, read, phreds, min_seq, quality_cutoff):
    bar_positions = []
    e_metrics = []
    for i in range(0, len(read)-len(barcode)):
        total_identical = 0
        longest_sequence = 0
        current_sequence = 0
        metric_arr = []
        for char_idx, (charB, charR, phred) in enumerate(zip(barcode, read[i:i+len(barcode)], phreds[i:i+len(barcode)])):
            if charB == charR or phred <= quality_cutoff:
                total_identical += 1
                current_sequence += 1
            else:
                if args.em:
                    metric_arr.append((i+char_idx, charR, charB, phred))
                longest_sequence = max(longest_sequence, current_sequence)
                current_sequence = 0
        if longest_sequence >= min_seq:
            bar_positions.append(i)
            if args.em and metric_arr:
                e_metrics.append(metric_arr)
    return bar_positions, e_metrics


def euclidian_dist(barcode, read, phreds, mm_threshold, quality_cutoff):
    bar_positions = []
    e_metrics = []
    for i in range(0, len(read)-len(barcode)):
        euclid_dist = 0
        metric_arr = []
        for char_idx, (charB, charR, phred) in enumerate(zip(barcode, read[i:i+len(barcode)], phreds[i:i+len(barcode)])):
            if charB != charR and phred > quality_cutoff:
                euclid_dist += 1
                if args.em:
                    metric_arr.append((i+char_idx, charR, charB, phred))
            if euclid_dist > mm_threshold: break
        if euclid_dist <= mm_threshold:
            bar_positions.append(i)
            if metric_arr:
                e_metrics.append(metric_arr)
    if not len(bar_positions) and args.el:
        bar_position, metric_arr = check_edges(barcode, read, phreds, mm_threshold)
        if bar_position is not False:
            bar_positions.append(bar_position)
            if metric_arr:
                e_metrics.append(metric_arr)
    return bar_positions, e_metrics

def check_edges(barcode, read, phreds, mm_threshold):
    if len(barcode) >= len(read): return False, []
    try:
        result = [_.start() for _ in re.finditer(barcode[:args.el], read[-len(barcode):])][-1]
        result += len(read)-len(barcode)
    except:
        result = False
    if result:
        if mm_threshold >= len(read[result+args.el:])-(len(barcode)-args.el):
            euclid_dist = max(len(read[result+args.el:])-(len(barcode)-args.el), 0)
            metric_arr = []
            for char_idx, (charB, charR, phred) in enumerate(zip(barcode[args.el:], read[result+args.el:], phreds[result+args.el:])):
                if charB != charR:
                    euclid_dist += 1
                    if args.em:
                        metric_arr.append((result+char_idx, charR, charB, phred))
                if euclid_dist > mm_threshold: break
            if euclid_dist <= mm_threshold:
                return result, metric_arr
    rev_barcode, rev_read, rev_phreds = barcode[::-1], read[::-1], phreds[::-1]
    try:
        result = [_.start() for _ in re.finditer(rev_barcode[:args.el], rev_read[-len(rev_barcode):])][-1]
        result += len(rev_read)-len(rev_barcode)
    except:
        result = False
    if result:
        if mm_threshold >= len(rev_read[result+args.el:])-(len(rev_barcode)-args.el):
            euclid_dist = max(len(rev_read[result+args.el:])-(len(rev_barcode)-args.el), 0)
            metric_arr = []
            for char_idx, (charB, charR, phred) in enumerate(zip(rev_barcode[args.el:], rev_read[result+args.el:], rev_phreds[result+args.el:])):
                if charB != charR:
                    euclid_dist += 1
                    if args.em:
                        metric_arr.append((char_idx, charR, charB, phred))
                if euclid_dist > mm_threshold: break
            if euclid_dist <= mm_threshold:
                return (0, len(read)-result), metric_arr
    # try:
    #     result = [_.start() for _ in re.finditer(barcode[-args.el:], read[:len(barcode)])][0]
    # except:
    #     result = False
    # if result:
    #     print('second', result)
    #     if mm_threshold >= len(read[result-1::-1])-(len(barcode)-args.el):
    #         euclid_dist = max(len(read[result-1::-1])-(len(barcode)-args.el), 0)
    #         metric_arr = []
    #         for char_idx, (charB, charR, phred) in enumerate(zip(barcode[:args.el:-1], read[:result-args.el:-1], phreds[:result-args.el:-1])):
    #             if charB != charR:
    #                 euclid_dist += 1
    #                 if args.em:
    #                     metric_arr.append((result+char_idx, charR, charB, phred))
    #             if euclid_dist > mm_threshold: break
    #         if euclid_dist <= mm_threshold:
    #             return 0, metric_arr
    return False, []

tot_barcodes = 0
barcode_metrics = defaultdict(int)
csv_writer = None
if args.em:
    import csv
    csv_fh = open(args.o+'extended_metrics.csv', 'w')
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
        if args.seq_hit:
            pos_list, e_metrics = find_tag(barcode, str(record.seq), record.letter_annotations['phred_quality'], args.ms, args.qc)
        else:
            pos_list, e_metrics = euclidian_dist(barcode, str(record.seq), record.letter_annotations['phred_quality'], args.mm, args.qc)
        for mut in e_metrics:
            if len(mut) > 1:
                mut = "|".join([str(m[0]) for m in mut]), "|".join([str(m[1]) for m in mut]), "|".join([str(m[2]) for m in mut]), \
                      "|".join([str(m[3]) for m in mut])
            else:
                mut = mut[0]
            if args.em:
                csv_writer.writerow([record_id, bar_id, mut[0], mut[1], mut[2], mut[3]])
        tot_barcodes += len(pos_list)
        record_cnt_barcodes += len(pos_list)
        if pos_list:
            barcode_hits.append(bar_id)
            if len(pos_list) > 1:
                barcode_metrics['dup'] += 1
            if not args.nr:
                removed = 0
                for idx, pos in enumerate(pos_list):
                    if isinstance(pos, tuple):
                        record = record[pos[1]:]
                        pos_list = [pos[0]]
                    else:
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
      f"Reads without barcode: {tot_reads-tot_reads_with_barcode}\nTotal barcodes: {tot_barcodes}\n\n"
      f"Number of reads with duplicate barcodes: {barcode_metrics['dup']}\nReads with multiple barcodes: {barcode_metrics['mult']}")
for file in open_files.values():
    file.close()
if args.em:
    csv_fh.close()
