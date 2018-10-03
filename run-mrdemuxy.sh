#!/bin/bash

# this program demultiplexes a set of fastq.gz files using Mr_Demuxy

# it needs files for the forward and reverse barcodes
echo -e "A   GGTAC
B   CAACAC
C   ATCGGTT
D   TCGGTCAA
E   AAGCG
F   GCCACA
G   CTGGATG
H   TGATTGAC" > forward_bcs.txt

echo -e "01  AGGAA
02  GAGTGG
03  CCACGTC
04  TTCTCAGC
05  CTAGG
06  TGCTTA
07  GCGAAGT
08  AATCCTAT
09  ATCTG
10  GAGACT
11  CGATTCC
12  TCTCAATC" > reverse_bcs.txt

#Mr_Demuxy only takes unzipped files
ls *.gz | parallel gunzip {}

# assumes fastq file name is like sample1_blah_blah_R1_001.fastq
for F1 in `ls *.fastq | cut -f1 -d'_' | sort -u`; do
 echo "demultiplexing $F1";
 pe_demuxer.py -r1 $F1*_R1_001.fastq -r2 $F1*_R2_001.fastq -r1_bc forward_bcs.txt -r2_bc reverse_bcs.txt

# rename them from Mr_Demuxy's default structure
 cd pe_demuxer_output
 for D2 in R1 R2; do
   cd $D2
   for FILE in *.fastq; do
    mv $FILE $F1"_"$FILE
   done;
   cd ..
 done;

cd ..

#rezip files
 gzip pe_demuxer_output/R1/*.fastq
 gzip pe_demuxer_output/R2/*.fastq

mv pe_demuxer_output/*/*.gz ./

rm -rf pe_demuxer_output

done;
