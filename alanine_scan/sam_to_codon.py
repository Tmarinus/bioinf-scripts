#!/usr/bin/env python

import sys
try:
    import pandas as pd
    import pysam
    import numpy as np
    from Bio import SeqIO
    from alive_progress import alive_bar
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import matplotlib.pyplot as plt
except ImportError:
    print(f'Script requires pysam module to be installed run\npip3 install pysam --user\npip3 install BioPython --user\n'
          f'pip3 install alive_progress --user\npip3 install matplotlib --user\npip3 install pandas --user')
    exit()


if len(sys.argv) != 4:
    print(f'To run give paths to vcf file, sam/bam and output\n'
          f'eg: sam_to_codon test_files/barcode3_indels.vcf test_files/template_gal_k2.fasta output.fasta')
    print(f'Make sure the bam/sam files have an index\nsamtools index file.bam')
    exit()

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
# if len(sys.argv) != 4:
#     print(f'To run give paths to vcf file, fasta and output\n'
#           f'eg: vcf_to_fasta test_files/barcode3_indels.vcf test_files/template_gal_k2.fasta output.fasta')
#     exit()
#
# vcf_file = sys.argv[1]
# fasta_file = sys.argv[2]
# output = sys.argv[3]
# save = pysam.set_verbosity(0)
# samfile = pysam.AlignmentFile("/home/tycho/tmp_r/barcode3.sorted.bam", "rb")
samfile = pysam.AlignmentFile(sam_file, "rb")
# pysam.set_verbosity(save)

#TODO hack
amino_str = 'MKETTTSLMQDELTLGEPATQARMCVRLLRFFIGLTITAFIIAACIIKSATGGSGYSKAVAVRGEADTPSTIVGQLVERGGFQAWAVGAGIYLFAKIAYDTSKVTAAVCNPEALIAITSYVAYAPTLCAGAYVIGAMSGAMSAGLALYAGYKGWQWGGPGGMAEREDVASFYSPLLNNTLYVGGDHTADYDSELATILGSVYNDVVHLGVYYDNSTGIVKRDSRPSMISWTVLHDNMMITSYHRPDQLGAAATAYKAYTTNTTRVGKRQDGEWVSYSVYGENVDYERYPVAHLQEEADACYESLGNMITSQVQPCTQRECYAMDQKVCAAVGFSSDAGVNSAMVGEAYFYAYGGVDGECDSG'

a_codons = ['GCT', 'GCC', 'GCA', 'GCG']
orf_offset = 0
read_count = np.array([0]*len(f_str))
a_count = np.array([0]*len(f_str))
a_count1 = np.array([0]*len(f_str))
a_count2 = np.array([0]*len(f_str))
cnt = 0
for read in samfile.fetch():
    if len(read.get_blocks()) == 1 and not read.flag:
        # TODO fix clippings
        if 'H' in read.cigarstring or 'S' in read.cigarstring: continue
        cnt+=1
        forward_str = read.get_forward_sequence()
        start, end = read.reference_start, read.reference_end
        if end-start != len(forward_str):
            print('errro')
        read_count[start:end] += 1
        codon_start_offset = (read.reference_start-2527) % 3
        if codon_start_offset: codon_start_offset = 3-codon_start_offset
        #TODO brainfart
        # print(read.reference_start, read.query_length, read.reference_end)
        start_offset, end_offset = 3-read.reference_start % 3, read.reference_end % 3
        for codon in range(codon_start_offset, end-start, 3):
            if len(forward_str[codon:codon+3]) < 3:
                break
            aa = Seq(forward_str[codon:codon+3])
            if str(aa) == 'GCT' and Seq(f_str[start+codon:start+codon+3]).translate() is not 'A':
                a_count[start+codon:start+codon+3] += 1
            if str(aa) == 'GGT' and Seq(f_str[start+codon:start+codon+3]).translate() is 'A':
                a_count[start+codon:start+codon+3] += 1
for z,x in zip(a_count, read_count):
    if x == 0 and z > 0:
        print('error')
normalized = np.divide(a_count, read_count, where=read_count != 0)
print(cnt, samfile.count())
normalized[np.isnan(normalized)] = 0
normalized[np.isinf(normalized)] = 0
amino_list = list(normalized)[2530:2530+1086:3]
print('average', np.average(amino_list))
read_avg = []
coverage_counted_reads = read_count/cnt
for idx in range(2530, 2530+1086, 3):
    read_avg.append(np.average(coverage_counted_reads[idx:idx+3]))

print(len(read_avg), len(amino_list))
result_dict = {
    'amino_acid':list(amino_str),
    'coverage_total': read_avg,
    'alanine_mut':amino_list
}
df = pd.DataFrame(result_dict)

df.to_csv(output.split('.')[0]+'.csv')

fig = plt.figure()
plt.title(output)
plt.plot(list(range(0, len(amino_list))), amino_list)
plt.xlabel('AminoAcid Pos')
plt.ylabel('MutFreq\nNumAlanineMutations/Reads')
plt.ylim(0,0.05)
fig.savefig(output, dpi=200)

