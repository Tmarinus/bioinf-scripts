# bioinf-scripts
Scripts for simple and quick bioinformatics data analysis. This is mainly for myself to keep track of my scripts.

Feel free to use or adjust. Bugs can be reported. Do note that these scripts are mostly quick python scripts made on a Friday evening to support colleagues. They are written without any real optimizations in mind.

## Alanine-scan
#### Work in progress!

Code to analyse NGS data of alanine scanning experiments. Cells were grouped in halo size, DNA library prep was performed per group with an unique bar code. Goal of this script was to analyse the mapped sam/bam file and find all amino acids that are mutated into Alanine.

## fastq_splitter
Two scripts that are used for extracting specific sequences from reads, and reporting where these sequences were found in a reference.
fastq_splitter: Splits fastq reads based on sequence found in the reads. By default also removes the removed barcode and adds the position where it was found to the read ID. (Once bowtie2 2.4.2 is more main stream this can become the read description.).

sam_to_barcode_position, after mapping the splitted and barcode removed fastq file, this script creates a csv file with per base read coverage and reports where the specific barcode sequences were found. (pos 5 means barcode was removed nucleotide after position 5.)