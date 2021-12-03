#!/usr/bin/env python
try:
    import os
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Data.CodonTable import CodonTable as CT
    from Bio.Seq import Seq
    import argparse
except ImportError:
    print(f'Script requires pysam module to be installed run\npip3 install BioPython --user\n')
    exit()

parser = argparse.ArgumentParser(usage=f'Split fastq files based on length of a read.🐻',
           description=f'Script requires BioPython module to be installed run\npip3 install BioPython --user')

parser.add_argument('fastq_file', help='path to fastq file')
parser.add_argument('length', help='Read length to split on. len(read) >= length', type=int)
parser.add_argument('-o', metavar='--output', help='Output folder, default is same folder as input, two files shorter and longer will be created.',
                    default=None)
args = parser.parse_args()

if args.o:
    output_folder = os.path.join(args.o, '')
else:
    output_folder = os.path.join(os.path.dirname(args.fastq_file), '')

fastq_name = os.path.basename(args.fastq_file)
cnt_l = 0
cnt_s = 0
with open(output_folder+f"longer-{fastq_name}", 'w') as longer, open(output_folder+f"shorter-{fastq_name}", 'w') as shorter:
    for record in SeqIO.parse(args.fastq_file, 'fastq'):
        if len(record) >= args.length:
            SeqIO.write(record, longer, 'fastq')
            cnt_l += 1
        else:
            SeqIO.write(record, shorter, 'fastq')
            cnt_s += 1
print(f"Reads longer or equal than {args.length}:\n{cnt_l:,} Longer\n{cnt_s:,} Shorter")