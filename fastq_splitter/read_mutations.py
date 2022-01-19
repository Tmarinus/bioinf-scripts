#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description="Calculate mutations from sam file.")
parser.add_argument('sam_file')
args = parser.parse_args()

file_p = args.sam_file
ins_ranges = [range(4867, 4956), range(4957, 5046), range(5047, 5133)]


import pysam
samfile = pysam.AlignmentFile(args.sam_file)

total_reads = 0
unmapped = 0
wildtype = 0
dels_and_ins = 0
dels_or_ins = 0
insertions = []
deletions = []
r_count_insertions = []
r_count_deletions = []
inserted_ending = 0
size_change = 0
for read in samfile:
    total_reads += 1
    if read.is_unmapped:
        unmapped += 1
        continue
    if not (320 < len(read.seq) < 340):
        # print(read.seq)
        size_change += 1
    tmp_ins, tmp_dels = 0, 0
    for cigar in read.cigartuples:
        if 1 == cigar[0]:
            tmp_ins += cigar[1]
        if 2 == cigar[0]:
            tmp_dels += cigar[1]
    if tmp_ins > 0 and tmp_dels > 0:
        dels_and_ins += 1
    if tmp_ins > 0 or tmp_dels > 0:
        dels_or_ins += 1
    if tmp_ins:
        if 'gccaata'.upper() in read.get_tag('MD'):
            inserted_ending += 1

        insertions.append(tmp_ins)
    if tmp_dels: deletions.append(tmp_dels)
    if tmp_ins == tmp_dels == 0 and not read.query_name.split('_')[:-1]: wildtype += 1
print(f"size: {size_change}")
print(f"total reads: {total_reads}\nReads with Insertions: {len(insertions)} {(len(insertions)/total_reads)*100:0.1f}%"
      f" and deletions: {len(deletions)} {(len(deletions)/total_reads)*100:0.1f}%\n"
      f"Both: {dels_and_ins} {(dels_and_ins/total_reads)*100:0.1f}% one or both: {dels_or_ins} {(dels_or_ins/total_reads)*100:0.1f}%\n\n"
      f"Wildtypes: {wildtype} {(wildtype/total_reads)*100:0.1f}% unmapped: {unmapped} {(unmapped/total_reads)*100:0.1f}%")

print(f"inserted ending: {inserted_ending} {(inserted_ending/len(insertions))*100:0.1f}%")