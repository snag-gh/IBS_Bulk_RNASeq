#! /bin/bash

FASTQDIR=/data/Catiele/usftp1.novogene.com/X202SC20020886-Z01-F001/rawdata
mkdir $1

salmon quant -i /data/users/sushmanagaraj/ref/mm/vM25/salmon_index/ -l A -1 ${FASTQDIR}/${1}_1.fq.gz -2 ${FASTQDIR}/${1}_2.fq.gz --validateMappings -p 20 --gcBias --seqBias -o ${1}/transcripts_quant 


