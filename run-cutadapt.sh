#!/usr/bin/env bash

# this program trims primers and adapters
# input DIRECTORY as $1

# takes input after initial demultiplexing of outer barcode using Mr_Demuxy

for F1 in `find $1 -name "*fastq.gz" | cut -f1-2 -d'_' | uniq`; do
echo $F1

	DIR=`dirname $F1`
	FILE=`basename $F1`

	for FQ in `find $DIR -name "$FILE*R1*.fastq.gz" | grep -v trimmed`; do

    RQ=${FQ/_R1/_R2} #input reverse fastq

    FQO=${FQ/fastq\.gz/ctrimmed\.fastq\.gz} #final output file
    RQO=${RQ/fastq\.gz/ctrimmed\.fastq\.gz}
    ST=`basename $FQO .fastq.gz`

    echo -e ">515F
    GTGCCAGCMGCCGCGGTAA
    >806RB
    GGACTACNVGGGTWTCTAAT
    >TAReuk454FWD1
    CCAGCASCYGCGGTAATTCC
    >TAReukREV3_modified
    ACTTTCGTTCTTGATYRATGA
    >ITS-F
    GAACCWGCGGARGGATCATTA
    >ITS-R
    GCTGCGTTCTTCATCGATG" > primers-all.fa

    echo -e ">515F
    GTGCCAGCMGCCGCGGTAA
    >806RB
    GGACTACNVGGGTWTCTAAT" > primers-bac.fa

    echo -e ">TAReuk454FWD1
    CCAGCASCYGCGGTAATTCC
    >TAReukREV3_modified
    ACTTTCGTTCTTGATYRATGA" > primers-euk.fa

    echo -e ">demult-ITS-F
    GAACCWGCGGARGGATCA
    >demult-ITS-R
    GCTGCGTTCTTCATCGATG" > primers-its.fa

    PRIMERS=primers-all.fa;

    # we use multiple primers, so naming scheme assumes sample names start with "16S", "18S", or "ITS" to distinguish them
    if [[ "$ST" =~ ^16.* ]];
          then PRIMERS=primers-bac.fa;
          elif [[ "$ST" =~ ^18.* ]];
          then PRIMERS=primers-euk.fa;
          elif [[ "$ST" =~ ^ITS.* ]];
          then PRIMERS=primers-its.fa;
          else PRIMERS=primers-all.fa;
      fi

    # using cutadapt: https://github.com/marcelm/cutadapt

    /work/cryomics/apps/cutadapt -g file:$PRIMERS -G file:$PRIMERS --discard-untrimmed -j 8 -e 0.12 -o $FQO -p $RQO $FQ $RQ

	done;
done;
