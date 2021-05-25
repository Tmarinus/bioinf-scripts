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
parser.add_argument('-nr', default=False, action='store_true', help='Do not remove barcode from the read.')
parser.add_argument('-kd', default=False, action='store_true',
                    help='If multiple barcodes are found in a read also put them in the individual fastq files.')
args = parser.parse_args()

from Bio import SeqIO

barcodes = {}
if args.b:
    barcodes += dict(zip(args.b.split(','), range(0, len(args.b.split(',')))))
if args.bf:
    for record in SeqIO.parse(args.bf, 'fasta'):
        barcodes[record.id] = str(record.seq)
fastq_name = args.fastq_file.split('/')[-1]

if args.o and not args.o.endswith('/'):
    args.o += '/'
open_files = {
    'nobar': open(args.o+f"no_barcode-{fastq_name}", 'w')
}
for barcode in barcodes.keys():
    open_files[barcode] = open(args.o+f"bar_{barcode}-{fastq_name}", 'w')

open_files['all'] = open(args.o+f"all-{fastq_name}", 'w')


def append_fastq(seq, path):
    SeqIO.write(seq, path, 'fastq')


for record in SeqIO.parse(args.fastq_file, 'fastq'):
    barcode_hits = []
    tmp = len(record)
    old_rec = record
    for bar_id, barcode in barcodes.items():
        pos_list = [i.start() for i in re.finditer(barcode, str(record.seq))]
        if pos_list:
            barcode_hits.append(bar_id)
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
        append_fastq(record, open_files[barcode_hits[0]])
    elif len(barcode_hits) > 1:
        if args.kd:
            for bar_id in barcode_hits:
                append_fastq(record, open_files[bar_id])
        bar_ids = tuple(barcode_hits)
        try:
            append_fastq(record, open_files[bar_ids])
        except KeyError:
            open_files[bar_ids] = open(args.o+f"bars_{'_'.join(bar_ids)}-{fastq_name}", 'w')
            append_fastq(record, open_files[bar_ids])

for file in open_files.values():
    file.close()
