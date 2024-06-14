#!/usr/bin/env python
import sys
try:
    import os
    import pandas as pd
    import pysam
    import numpy as np
    from Bio import SeqIO
    from alive_progress import alive_bar
    from Bio.SeqRecord import SeqRecord
    from Bio.Data.CodonTable import CodonTable as CT
    from Bio.Seq import Seq
    import argparse
    import matplotlib.pyplot as plt
    from collections import defaultdict
except ImportError:
    print(f'Script requires pysam, pandas, numpy, Biopython, matplotlib modules to be installed run\n'
          f'pip3 install pysam --user\npip3 install BioPython --user\n'
          f'pip3 install matplotlib --user\npip3 install pandas --user')
    exit()

AA_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', 'Z', 'J', 'U', '*']

AA_conv_dict = {
    'A': 'Alanine',
    'C': 'Cysteine',
    'D': 'Aspartic-acid',
    'E': 'Glutamic-acid',
    'F': 'Phenylalanine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'K': 'Lysine',
    'L': 'Leucine',
    'M': 'Methionine',
    'N': 'Asparagine',
    'P': 'Proline',
    'Q': 'Glutamine',
    'R': 'Arginine',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine'
}

parser = argparse.ArgumentParser(usage=f'Param order; sam_path, fasta_path, output_path, AA_start (in fasta), AA_end.\n'
                                       f'Optional use --split_sam to split the sam file based on read conditions. 🐻',
           description=f'Script requires BioPython module to be installed run\npip3 install BioPython --user')

parser.add_argument('sam_file', help='path to sam file')
parser.add_argument('fasta_file', help='path to fasta file used for read mapping')
parser.add_argument('output', help='file name and path for output (will create output.xslx and if split sam enabled a sam folder for the output)')
parser.add_argument('aa_start', help='AA starting position in fasta file (start counting at 1)', type=int)
parser.add_argument('aa_end', help='AA ending position in fasta file (this is the last nucleotide that is still part of the AA seq)', type=int)
parser.add_argument('--split_sam', help='Split sam files, depending on read type, will be stored in sam folder', action='store_true')
parser.add_argument('--mut_quality', help=f'Output mutation quality. Will be split on (expected) target aa mut or non target aa mut'
                                          f'Default calculation is minimum quality of the three nucleotides of the codon', action='store_true')
parser.add_argument('--mut_quality_max', help=f'Set mut_quality to maximum of three codons', action='store_true')
parser.add_argument('--mut_quality_avg', help=f'Set mut_quality to average of three codons', action='store_true')
parser.add_argument('--set_target_aa', nargs=2, help=f'Set the target amino acid for the scan. By default this is Alanine (incase reference is Alanine then Glycine).'+\
    f"Example --set_target_aa A G")

args = parser.parse_args()

target_aa = args.set_target_aa
t_aa = 'A'
m_aa = 'G'
if target_aa:
    t_aa = target_aa[0]
    m_aa = target_aa[1]

fasta_file = args.fasta_file
sam_file = args.sam_file
AA_start = args.aa_start - 1
AA_end = args.aa_end - 1
#TODO hack
# amino_str = args.aa
# ins_site = 'CCTGCAGGG'  # sbf1
output = os.path.abspath(args.output)

if output.endswith('/'):
    sam_output = output+'sam_files'
    if not os.path.exists(output):
        os.makedirs(output)
    output = output+'output'
else:
    sam_output = "/".join(output.split('/')[:-1])+f"/sam_files_{output.split('/')[-1]}"

AA_length = int((AA_end+1 - AA_start) / 3)

fasta_dict = {}
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        fasta_dict[record.id] = record
f_str = str(list(fasta_dict.values())[0].seq)

samfile = pysam.AlignmentFile(sam_file, "rb")

read_passed_cnt, read_mut_cnt, target_mut_cnt, alt_mut_cnt = 0, [0, 0, 0, 0, 0], 0, 0
read_count = np.array([0]*len(f_str))
aa_target_only = np.array([0]*AA_length)
aa_alt_muts = np.array([0]*AA_length)
idx = 0
mut_list = [dict((AA, 0) for AA in AA_list) for x in range(AA_length)]
sam_writers = {}
mut_quality_target = []
mut_quality_other = []
mut_quality = None
tmp_muts = []

mut_qualities = {
    f'singl_{AA_conv_dict[t_aa]}': [],
    f'not_alone_{AA_conv_dict[t_aa]}': [],
    'singl_other': [],
    f'other_besides_{AA_conv_dict[t_aa]}': [],
    f'other_no_{AA_conv_dict[t_aa]}': [],
}

if args.split_sam:
    if not os.path.exists(sam_output):
        os.makedirs(sam_output)
    sam_writers = {
        'H': pysam.AlignmentFile(f"{sam_output}/H.bam", "wb", template=samfile),
        'S': pysam.AlignmentFile(f"{sam_output}/S.bam", "wb", template=samfile),
        'indel': pysam.AlignmentFile(f"{sam_output}/indel.bam", "wb", template=samfile),
        'no_mut': pysam.AlignmentFile(f"{sam_output}/no_mut.bam", "wb", template=samfile),
    }

def write_sams(key, read):
    if not args.split_sam: return
    if key not in sam_writers.keys():
        sam_writers[key] = pysam.AlignmentFile(f"{sam_output}/{key}.bam", "wb", template=samfile)
    sam_writers[key].write(read)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

for idx, read in enumerate(samfile.fetch()):
    # Skipping all reads that are not perfectly mapped!
    if len(read.get_blocks()) != 1 or read.flag:
        write_sams('indel', read)
        continue
    if 'H' in read.cigarstring:
        write_sams('H', read)
        continue
    if 'S' in read.cigarstring:
        write_sams('S', read)
        continue
    forward_str = read.get_forward_sequence()
    start, end = read.reference_start, read.reference_end
    # Ignore the part of the read that is not within the gene region of interest.
    if start < AA_start:
        forward_str = forward_str[AA_start-start:]
        start = AA_start
    if end > AA_end:
        forward_str = forward_str[:AA_end-end]
        end = AA_end
    if end-start != len(forward_str):
        print(f"Significant error found within a read. This can either be that the end and start of region of interest are mixed, or a"
              f"incorrect read. Program will continue to run but check this: \n\n {read}")
        continue
    read_passed_cnt += 1
    read_count[start:end] += 1

    codon_start_offset = ((AA_start-start) % 3)
    read_codon_start = codon_start_offset + start
    codon_end_offset = ((end - read_codon_start) % 3)
    read_codon_end = end - codon_end_offset

    if ((read_codon_start-AA_start)/3).is_integer():
        read_first_aa = int((read_codon_start-AA_start)/3)
    else:
        print(f"Warning bug in code, codon does not match read idx: {idx}. exiting")
        exit()
    if f_str[read_codon_start:read_codon_end] == forward_str[codon_start_offset:len(forward_str)-codon_end_offset]:
        read_mut_cnt[0] += 1
        write_sams('no_mut', read)
        continue

    phred_quality = [ord(x)-33 for x in read.qual][codon_start_offset:len(forward_str)-codon_end_offset]
    f_nucl = f_str[read_codon_start:read_codon_end]
    read_nucl = forward_str[codon_start_offset:len(forward_str)-codon_end_offset]
    f_aa = Seq(f_nucl).translate()
    read_aa = Seq(read_nucl).translate()
    # TODO rename x,y,q into something logical.
    aa_diff = [(pos, pos+read_first_aa, x, y, q) for pos, (x, y, q) in enumerate(zip(f_aa, read_aa, chunks(phred_quality, 3))) if x != y]
    while len(aa_diff) >= len(read_mut_cnt): # Expand array size to be able to count number of mutations
        read_mut_cnt.append(0)
    read_mut_cnt[len(aa_diff)] += 1
    write_sams(len(aa_diff), read)
    target_aa_found = False
    for mut in aa_diff: # iterate over mutations
        if args.mut_quality:
            if args.mut_quality_max:
                mut_quality = max(mut[4])
            elif args.mut_quality_avg:
                mut_quality = sum(mut[4])/3
            else:
                mut_quality = min(mut[4])

        try:
            mut_list[mut[1]][mut[3]] += 1
        except KeyError:
            print(mut)
            print("error was found in above mutation. Contact support.")
            exit()
        # Check if target aa mutation (or guanine)
        if mut[3] == t_aa or (mut[2] == t_aa and mut[3] == m_aa):
            target_mut_cnt += 1
            mut_quality_target.append(mut_quality)
            target_aa_found = True
            # Was target aa the only mutation?
            if len(aa_diff) == 1:
                mut_qualities[f'singl_{AA_conv_dict[t_aa]}'].append(mut_quality)
                aa_target_only[mut[1]] += 1
            else:
                mut_qualities[f'not_alone_{AA_conv_dict[t_aa]}'].append(mut_quality)
        else:
            alt_mut_cnt += 1
            aa_alt_muts[mut[1]] += 1
            mut_quality_other.append(mut_quality)
            if len(aa_diff) == 1:
                mut_qualities['singl_other'].append(mut_quality)
            else:
                tmp_muts.append(mut_quality)
    if tmp_muts:
        if target_aa_found:
            mut_qualities[f'other_besides_{AA_conv_dict[t_aa]}'].extend(tmp_muts)
        else:
            mut_qualities[f'other_no_{AA_conv_dict[t_aa]}'].extend(tmp_muts)
        tmp_muts = []


read_cnt = idx+1
print(f"{AA_conv_dict[t_aa]} mutations {target_mut_cnt} other mutations {alt_mut_cnt}\n"
      f"Reads with only a single {AA_conv_dict[t_aa]} mutation: {sum(aa_target_only)} {(sum(aa_target_only)/read_passed_cnt)*100:0.1f}% of passed reads")
print(f"Total reads {read_cnt} used reads: {read_passed_cnt} {(read_passed_cnt/read_cnt)*100:0.1f}% not used: "
      f"{read_cnt-read_passed_cnt} {((read_cnt-read_passed_cnt)/read_cnt)*100:0.1f}% ")
aa_coverage = []
for codon_start in range(AA_start, AA_end, 3):
    aa_coverage.append(int(sum(read_count[codon_start:codon_start+3])/3))

# Converting data to pandas dataframe
df = pd.DataFrame(mut_list[0], index=[0])
for muts in mut_list[1:]:
    df = df.append(muts, ignore_index=True)
df.insert(0, 'aa_pos', list(range(1, AA_length+1)), True)
# TODO fix the error
print(f"length {len(f_str)}, {AA_start} {AA_end}")
df.insert(1, 'amino_acid', list(Seq(f_str[AA_start:AA_end+2]).translate()), True)
print("B")
df.insert(2, 'coverage', aa_coverage, True)
df.insert(3, f'single_{AA_conv_dict[t_aa]}', aa_target_only, True)

# Writing dataframe to excel file
writer = pd.ExcelWriter(output+'.xlsx', engine='xlsxwriter')
workbook = writer.book
df.to_excel(writer, sheet_name="mutations")
mutations_sheet = writer.sheets["mutations"]
cell_format = workbook.add_format({'bold': False, 'border': False})
cell_format.set_font_size(10)
cell_format.set_font_name('Liberation Sans')
mutations_sheet.set_column(0, 100, 5, cell_format=cell_format)

#Write read mut count
num_mut_df = pd.DataFrame({'#mutations': list(range(len(read_mut_cnt))), 'occurrences': read_mut_cnt})
num_mut_df.to_excel(writer, sheet_name="muts_per_read", index=False)
#Write mutation quality
num_mut_df = pd.DataFrame({f'expected ({AA_conv_dict[t_aa]})': pd.Series(mut_quality_target),
                           'other': pd.Series(mut_quality_other),
                           f'Single {AA_conv_dict[t_aa]}': pd.Series(mut_qualities[f'singl_{AA_conv_dict[t_aa]}']),
                           f'{AA_conv_dict[t_aa]} with others': pd.Series(mut_qualities[f'not_alone_{AA_conv_dict[t_aa]}']),
                           f'Single non {AA_conv_dict[t_aa]}': pd.Series(mut_qualities['singl_other']),
                           f'Other next to {AA_conv_dict[t_aa]}': pd.Series(mut_qualities[f'other_besides_{AA_conv_dict[t_aa]}']),
                           f'Other no {AA_conv_dict[t_aa]}': pd.Series(mut_qualities[f'other_no_{AA_conv_dict[t_aa]}']),
                           })
num_mut_df.to_excel(writer, sheet_name="mutation_qualities", index=False)

writer.save()
print(f"\nsamfile location: {sam_output}/")
print(f"xlsx output file: {output}.xlsx")

for s_writer in sam_writers.values():
    s_writer.close()
