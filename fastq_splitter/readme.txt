Mutation increase, check the difference in barcodes found for more missmatch allowed:
Basically no change in barcodes, when increasing to more mismatches allowed

Dupplications are really in the reads



~/projects/scripts/fastq_splitter/remove_inserts.py KSNGS_join.fastq -i AGATCAACCACAATGTGTTGGGAGATCAACCACAATGTGTTG -r AGATCAACCACAATGTGTTGGG -mm 1

Total reads: 126002
Replaced 24 inserts
Mismatches 0 occurrences: 5
Mismatches 1 occurrences: 19

----> cleaned_KSNGS_join.fastq

time ~/projects/scripts/fastq_splitter/fastq_splitter.py  -bf ../../Streptag_EarItag.fasta cleaned-KSNGS_join.fastq -o test -el 6 -mm 1 -em 

Total reads: 126002
Reads with barcode: 97437 (0.773)
Reads without barcode: 28565 (0.227)
Total barcodes: 112763

Number of reads with duplicate barcodes: 12185 (0.097)
Reads with multiple barcodes: 2474 (0.020)

Barcodes found:
Strep-tag: 59074
EarI-tag: 53689

Total bases of barcodes: 2921068, total mutations found: 3489
Mutation rate in barcodes: 0.001

---> 
all-cleaned-KSNGS_join.fastq
bar_EarI-tag-cleaned-KSNGS_join.fastq
bars_Strep-tag_EarI-tag-cleaned-KSNGS_join.fastq
bar_Strep-tag-cleaned-KSNGS_join.fastq
no_barcode-cleaned-KSNGS_join.fastq
cleaned-KSNGS_join_extended_metrics.csv
--->

Important, use very sensitive bowtie setting

bowtie2 -p 12 -x /home/dwarrel/projects/data_testing/sonja/sonja_HMK_data/pRS423_GAL1p_HMK_ENO1t -U test/all-cleaned-KSNGS_join.fastq --very-sensitive -S test/all-KSNGS.sam

126002 reads; of these:
  126002 (100.00%) were unpaired; of these:
    3565 (2.83%) aligned 0 times
    122437 (97.17%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
97.17% overall alignment rate

---> all-KSNGS.sam

~/projects/scripts/fastq_splitter/sam_to_barcode_positions.py test/all-KSNGS.sam /home/dwarrel/projects/data_testing/sonja/sonja_HMK_data/pRS423_GAL1p_HMK_ENO1t.fasta -o test/

EarI-tag_mapping.csv
Strep-tag_mapping.csv
pSBxx__MoClopRS423type__GAL1p-HMKoriginalSeq-ENO1t_Segmented.csv



picard CreateSequenceDictionary R=../../pRS423_GAL1p_HMK_ENO1t.fasta  O=../../pRS423_GAL1p_HMK_ENO1t.dict