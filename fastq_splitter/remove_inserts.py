#!/usr/bin/env python
import argparse
import re

parser = argparse.ArgumentParser(description=f"Removes inserts from reads")
parser.add_argument('fastq_file')
parser.add_argument('-i', metavar='--insert', help='Insert to be removed')
parser.add_argument('-r', metavar='--replacement', help='Sequence to be inserted instead of the insert')
parser.add_argument('-o', metavar='--output', help='Output directory for the splitted fastq files', default='')
parser.add_argument('-mm', metavar='--miss-match-cutoff', help='Maximum number of miss matches for insert matching', default=1, type=int)
args = parser.parse_args()


from Bio import SeqIO

insert = args.i
replacement = args.r
replacement_phred = len(replacement) * [37]
fastq_name = args.fastq_file.split('/')[-1]

# Parse output naming and location
if args.o and not args.o.endswith('/'):
    args.o += '/'
cleaned_fq = open(args.o+f"cleaned-{fastq_name}", 'w')

# Append read to fastq file
def append_fastq(seq, path):
    SeqIO.write(seq, path, 'fastq')


def euclidian_dist(barcode, read, mm_threshold):
    bar_positions = []
    euclid_distances = []
    for i in range(0, len(read)-len(barcode)):
        euclid_dist = 0
        for char_idx, (charB, charR) in enumerate(zip(barcode, read[i:i+len(barcode)])):
            if charB != charR:
                euclid_dist += 1
            if euclid_dist > mm_threshold: break
        if euclid_dist <= mm_threshold:
            bar_positions.append(i)
            euclid_distances.append(euclid_dist)
    return bar_positions, euclid_distances

tot_reads = 0
tot_reads_with_insert = 0
distances = [0]*len(insert)
# Loop through records in fastq file
for record in SeqIO.parse(args.fastq_file, 'fastq'):
    tot_reads += 1
    barcode_hits = []
    record_id = record.id
    record_cnt_barcodes = 0
    pos_list, euclid_distances = euclidian_dist(insert, str(record.seq), args.mm)
    tot_reads_with_insert += len(pos_list)
    record_cnt_barcodes += len(pos_list)
    # print(record)
    #         print(record.letter_annotations['phred_quality'])
    if pos_list:
        for dist in euclid_distances:
            distances[dist] += 1
        removed = 0
        for idx, pos in enumerate(pos_list):
            phred_quality = record.letter_annotations['phred_quality'][:pos-removed] + replacement_phred + \
                                                         record.letter_annotations['phred_quality'][(pos+len(insert))-removed:]
            record.letter_annotations = {}
            record.seq = record.seq[:pos-removed] + replacement + record.seq[(pos+len(insert))-removed:]
            record.letter_annotations = {'phred_quality': phred_quality}
            # record = SeqRecord(tmp_seq, id=record.id, )
            pos_list[idx] = pos_list[idx] - removed
            removed += len(insert)
    append_fastq(record, cleaned_fq)
print(f"Total reads: {tot_reads}\nReplaced {sum(distances)} inserts")
import numpy
for number, dist_cnt in enumerate(distances[:numpy.max(numpy.nonzero(distances))+1]):
    print(f"Mismatches {number} occurrences: {dist_cnt}")

cleaned_fq.close()