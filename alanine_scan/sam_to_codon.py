#!/usr/bin/env python

import sys
try:
    import pandas as pd
    import pysam
    import numpy as np
    from Bio import SeqIO
    from alive_progress import alive_bar
    from Bio.SeqRecord import SeqRecord
    from Bio.Data.CodonTable import CodonTable as CT
    from Bio.Seq import Seq
    import matplotlib.pyplot as plt
except ImportError:
    print(f'Script requires pysam module to be installed run\npip3 install pysam --user\npip3 install BioPython --user\n'
          f'pip3 install alive_progress --user\npip3 install matplotlib --user\npip3 install pandas --user')
    exit()
debug = False

if len(sys.argv) != 4 and not debug:
    # print(f'To run give paths to sam/bam, nucleotide reference, aa reference,  and output\n'
    #       f'eg: sam_to_codon test_files/template_gal_k2.fasta test_files/template_gal_k2_AA.fasta output.fasta')
    print(f'TMP can only do k2 killer toxin. just provide sam/bam fasta_ref and output path')
    print(f'Make sure the bam/sam files have an index\nsamtools index file.bam')
    exit()


# Codon table:
# bases = "TCAG".upper()
# codons = [a + b + c for a in bases for b in bases for c in bases]
# amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
# codon_table = dict(zip(codons, amino_acids))


if debug:
    fasta_file = '/home/dwarrel/projects/data_testing/k2_toxin/template_gal_k2.fasta'
    sam_file = '/home/dwarrel/projects/data_testing/k2_toxin/bowtie_barcode3.sam'
    AA_start = 2530  # in fasta 2531
    AA_end = 3616
    #TODO hack
    amino_str = 'MKETTTSLMQDELTLGEPATQARMCVRLLRFFIGLTITAFIIAACIIKSATGGSGYSKAVAVRGEADTPSTIVGQLVERGGFQAWAVGAGIYLFAKIAYDTSKVTAAVCNPEALIAITSYVAYAPTLCAGAYVIGAMSGAMSAGLALYAGYKGWQWGGPGGMAEREDVASFYSPLLNNTLYVGGDHTADYDSELATILGSVYNDVVHLGVYYDNSTGIVKRDSRPSMISWTVLHDNMMITSYHRPDQLGAAATAYKAYTTNTTRVGKRQDGEWVSYSVYGENVDYERYPVAHLQEEADACYESLGNMITSQVQPCTQRECYAMDQKVCAAVGFSSDAGVNSAMVGEAYFYAYGGVDGECDSG'
    ins_site = 'CCTGCAGGG'  # sbf1
    output = 'test'
else:
    AA_start = 2530  # in fasta 2531
    AA_end = 3616
    amino_str = 'MKETTTSLMQDELTLGEPATQARMCVRLLRFFIGLTITAFIIAACIIKSATGGSGYSKAVAVRGEADTPSTIVGQLVERGGFQAWAVGAGIYLFAKIAYDTSKVTAAVCNPEALIAITSYVAYAPTLCAGAYVIGAMSGAMSAGLALYAGYKGWQWGGPGGMAEREDVASFYSPLLNNTLYVGGDHTADYDSELATILGSVYNDVVHLGVYYDNSTGIVKRDSRPSMISWTVLHDNMMITSYHRPDQLGAAATAYKAYTTNTTRVGKRQDGEWVSYSVYGENVDYERYPVAHLQEEADACYESLGNMITSQVQPCTQRECYAMDQKVCAAVGFSSDAGVNSAMVGEAYFYAYGGVDGECDSG'
    sam_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output = sys.argv[3]
# import seaborn as sns
fasta_dict = {}
# with open('/home/tycho/tmp_r/test_files/template_gal_k2.fasta', "r") as file:
with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        fasta_dict[record.id] = record
f_str = str(list(fasta_dict.values())[0].seq)

samfile = pysam.AlignmentFile(sam_file, "rb")
read_passed_cnt, read_mut_cnt, alanine_mut_cnt, alt_mut_cnt = 0, [0, 0, 0, 0, 0], 0, 0
read_count = np.array([0]*len(f_str))
aa_alanine_muts = np.array([0]*len(amino_str))
aa_alt_muts = np.array([0]*len(amino_str))
idx = 0
for idx, read in enumerate(samfile.fetch()):
    if len(read.get_blocks()) == 1 and not read.flag:
        # TODO fix clippings
        if 'H' in read.cigarstring or 'S' in read.cigarstring: continue
        read_passed_cnt += 1
        # if idx != 11: continue
        forward_str = read.get_forward_sequence()
        start, end = read.reference_start, read.reference_end
        read_count[start:end] += 1
        if start < AA_start:
            forward_str = forward_str[AA_start-start:]
            start = AA_start
        if end > AA_end:
            forward_str = forward_str[:AA_end-end]
            end = AA_end
        if end-start != len(forward_str):
            print('error, end and start are mixed?')
            exit()

        codon_start_offset = ((AA_start-start) % 3)
        read_codon_start = codon_start_offset + start
        codon_end_offset = ((end - read_codon_start) % 3)
        read_codon_end = end - codon_end_offset
        if ((read_codon_start-AA_start)/3).is_integer():
            read_first_aa = int((read_codon_start-AA_start)/3)
        else:
            print(f"Warning bug in code, codon does not match read idx: {idx}. exiting")
            exit()
        if f_str[read_codon_start:read_codon_end] == forward_str[codon_start_offset:len(forward_str)-codon_end_offset]: continue
        f_nucl = f_str[read_codon_start:read_codon_end]
        read_nucl = forward_str[codon_start_offset:len(forward_str)-codon_end_offset]
        f_aa = Seq(f_nucl).translate()
        read_aa = Seq(read_nucl).translate()
        aa_diff = [(pos, pos+read_first_aa, x, y) for pos, (x, y) in enumerate(zip(f_aa, read_aa)) if x != y]
        # print(aa_diff)
        while len(aa_diff) > len(read_mut_cnt):
            read_mut_cnt.append(0)
        read_mut_cnt[len(aa_diff)-1] += 1
        for mut in aa_diff:
            if mut[3] == 'A' or (mut[2] == 'A' and mut[3] == 'G'):
                alanine_mut_cnt += 1
                aa_alanine_muts[mut[1]] += 1
            else:
                alt_mut_cnt += 1
                aa_alt_muts[mut[1]] += 1
        # , alanine_mut_cnt, alt_mut_cnt
        # print(read_codon_start, (read_codon_start-AA_start)/3, codon_start_offset)
        # if idx > 10: exit()
        continue
read_cnt = idx+1
print(read_mut_cnt)
print(f"alanine mutations {alanine_mut_cnt} alternativemutations {alt_mut_cnt}")
print(f"Readcnt {read_cnt} used reads {read_passed_cnt}")
aa_coverage = []
for codon_start in range(AA_start, AA_end, 3):
    aa_coverage.append(sum(read_count[codon_start:codon_start+3])/3)
result_dict = {
    'amino_acid': list(amino_str),
    'coverage': aa_coverage,
    'alanine_muts':aa_alanine_muts,
    'alt_muts':aa_alt_muts
}
df = pd.DataFrame(result_dict)

df.to_csv(output.split('.')[0]+'.csv')
