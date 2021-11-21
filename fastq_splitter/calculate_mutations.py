#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description="Calculate mutations from sam file.")
parser.add_argument('sam_file')
parser.add_argument('fasta_file')
args = parser.parse_args()
import subprocess
import pysam

file_p = args.sam_file
ins_ranges = [range(4867, 4956), range(4957, 5046), range(5047, 5133)]


insert_bio = 0
insert_seq = 0
out_bio = 0
out_seq = 0
# file_p = '/home/dwarrel/projects/data_testing/sonja/sonja_HMK_data/sonja_old_dataset/fastq_join/test/variant/mutation_lines.txt'
tot = 0
cmd = subprocess.Popen(f"sam2tsv -R {args.fasta_file} {args.sam_file}", shell=True, stdout=subprocess.PIPE)
cmd.stdout.readline()
current_read = None
curr_bar_positions = []
active_ranges = []
prev_type = 'M'
prev_length = 0
prev_pos = 0
curr_read = [0,0,0]

MID_reads = []
cmd2 = subprocess.Popen(f"samtools view -c {args.sam_file}", shell=True, stdout=subprocess.PIPE)
cmd3 = subprocess.Popen(f"samtools view -c -F 4 {args.sam_file}", shell=True, stdout=subprocess.PIPE)
num_reads = int(cmd2.stdout.readline())
mapped = int(cmd3.stdout.readline())
unmapped = num_reads-mapped


from collections import defaultdict
type_counter = defaultdict(list)
for line in cmd.stdout:
    line = line.decode()
    tot += 1
    splits = line.split()
    type_s = splits[-1]
    if type_s.upper() == '.':
        if prev_type != '.':
            type_counter[prev_type].append(prev_length)
            prev_length = 0
            prev_type = '.'
        continue
    name = splits[0]
    try:
        pos = int(splits[-3])
        prev_pos = pos
    except ValueError:
        pos = prev_pos

    if name != current_read:
        MID_reads.append(curr_read)
        current_read = name
        curr_read = [0,0,0]
        curr_bar_positions = []
        for barcode in name.split('_')[:-1]:
            barcode_id, bar_positions = barcode.split(':')
            for bar_pos in bar_positions.split(','):
                curr_bar_positions.append(int(bar_pos))
        active_ranges = []
        for p in curr_bar_positions:
            for r in ins_ranges:
                if p+pos in r:
                    active_ranges += list(r)
        type_counter[prev_type].append(prev_length)
        prev_length = 0
        prev_type = type_s.upper()
        # active_ranges = [x if ]
    if prev_type != type_s.upper():
        type_counter[prev_type].append(prev_length)
        prev_length = 0
        prev_type = type_s.upper()
    if type_s.upper() == 'M':
        curr_read[0] = 1
    if type_s.upper() == 'I':
        curr_read[1] = 1
    if type_s.upper() == 'D':
        curr_read[2] = 1
    prev_length += 1
    ref = splits[-2]
    read = splits[-5]
    quality = splits[-4]
    if ref.upper() != read.upper():
        if pos in active_ranges:
            if quality == 'F':
                insert_bio += 1
            else:
                insert_seq += 1
        else:
            if quality == 'F':
                out_bio += 1
            else:
                out_seq += 1
        # break
if type_s.upper() == 'M':
    curr_read[0] = 1
if type_s.upper() == 'I':
    curr_read[1] = 1
if type_s.upper() == 'D':
    curr_read[2] = 1
MID_reads.append(curr_read)
MID_reads = MID_reads[1:]
tot_reads = len(MID_reads)
# print(MID_reads)
mut = 0
indels = 0
dels = 0
for x in MID_reads:
    mut += x[0]
    indels += x[1]
    dels += x[2]
print(f"\nTotal reads: {num_reads} unmapped {unmapped} {((unmapped)/num_reads)*100:.1f}%")
help_dict = {}
help_dict['M'] = mut
help_dict['I'] = indels
help_dict['D'] = dels
# print(type_counter)
for k, v in type_counter.items():
    if k == '.': continue
    print(f"{k} occurrences {help_dict[k]} average {sum(v)/len(v)}")
tot_mut = insert_seq+insert_bio+out_bio+out_seq
print(tot_mut)
print(f"\n\nTotal nucleotides: {tot:,}, total mutations {tot_mut:,}\n\n"
      f"Inserted region: {insert_bio+insert_seq:,} {((insert_bio+insert_seq)/tot_mut)*100:.1f}%, "
      f"wild region: {out_seq+out_bio:,} {((out_bio+out_seq)/tot_mut)*100:.1f}%\n\n"
      f"Seq error: {out_seq+insert_seq:,} {((insert_seq+out_seq)/tot_mut)*100:.1f}%\n"
      f"Bio error: {insert_bio+out_bio:,} {((out_bio+insert_bio)/tot_mut)*100:.1f}%\n\n")
if insert_bio+insert_seq > 0:
    print(f"Total mut in inserted regions\nBiomut: {insert_bio:,} {(insert_bio/(insert_bio+insert_seq))*100:.1f}% "
          f"Seq_err: {insert_seq:,} {(insert_seq/(insert_bio+insert_seq))*100:.1f}%\n\n")
print(out_bio, out_seq)
print(f"Total mut in wildtype regions\nBiomut: {out_bio:,} {(out_bio/(out_bio+out_seq))*100:.1f}% "
      f"Seq_err: {out_seq:,} {(out_seq/(out_bio+out_seq))*100:.1f}%")