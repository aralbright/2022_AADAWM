#!/bin/bash

# flexbar loop to trim NEBNext Single-cell adapters

# run in bin folder to access adapter fastas 


FILES="/Volumes/albright_postdoc/tubulinRNAseq/raw/AWM01/"

OUTPUTDIR="/Volumes/albright_postdoc/tubulinRNAseq/trimmed"

for first in $FILES*R1_001.fastq.gz
do

  output=$(echo $first | cut -c51-55)
  echo "$OUTPUTDIR/$output"
  echo "$first"
  echo "${first/_R1_001.fastq.gz/_R2_001.fastq.gz}"

  flexbar --reads $first \
          --reads2 ${first/_R1_001.fastq.gz/_R2_001.fastq.gz} \
          --stdout-reads \
          --adapters tso_g_wo_hp.fasta \
          --adapter-trim-end LEFT \
          --adapter-revcomp ON \
          --adapter-revcomp-end RIGHT \
          --htrim-left GT \
          --htrim-right CA \
          --htrim-min-length 3 \
          --htrim-max-length 5 \
          --htrim-max-first \
          --htrim-adapter \
          --min-read-length 2 \
          --threads 4 | \

          flexbar \
          --reads - \
          --interleaved \
          --target $OUTPUTDIR/$output \
          --adapters ilmn_20_2_seqs.fasta \
          --adapter-trim-end RIGHT \
          --min-read-length 2 \
          --threads 4
done

##################################################

FILES="/Volumes/albright_postdoc/tubulinRNAseq/raw/AWM02/"

OUTPUTDIR="/Volumes/albright_postdoc/tubulinRNAseq/trimmed"

for first in $FILES*R1_001.fastq.gz
do

  output=$(echo $first | cut -c51-55)
  echo "$OUTPUTDIR/$output"
  echo "$first"
  echo "${first/_R1_001.fastq.gz/_R2_001.fastq.gz}"

  flexbar --reads $first \
          --reads2 ${first/_R1_001.fastq.gz/_R2_001.fastq.gz} \
          --stdout-reads \
          --adapters tso_g_wo_hp.fasta \
          --adapter-trim-end LEFT \
          --adapter-revcomp ON \
          --adapter-revcomp-end RIGHT \
          --htrim-left GT \
          --htrim-right CA \
          --htrim-min-length 3 \
          --htrim-max-length 5 \
          --htrim-max-first \
          --htrim-adapter \
          --min-read-length 2 \
          --threads 4 | \

          flexbar \
          --reads - \
          --interleaved \
          --target $OUTPUTDIR/$output \
          --adapters ilmn_20_2_seqs.fasta \
          --adapter-trim-end RIGHT \
          --min-read-length 2 \
          --threads 4
done
