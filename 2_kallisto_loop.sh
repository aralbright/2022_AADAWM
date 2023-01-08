# Kallisto loop

# Kallisto loop
#!/bin/bash

kallisto index -i scoe-cds.idx S_coeruleus_Nov2016_CDS.fasta

########################

INDEX="/Volumes/albright_postdoc/tubulinRNAseq/genome/scoe-cds.idx"

FILES="/Volumes/albright_postdoc/tubulinRNAseq/trimmed/"

OUTPUTDIR="/Volumes/albright_postdoc/tubulinRNAseq/output_20221118"

for first in $FILES*_1.fastq.gz
do
  output=$(echo $first | cut -c49-53)
  echo "$OUTPUTDIR/$output"
  echo "$first"
  echo "${first/_1.fastq.gz/_2.fastq.gz}"

  kallisto quant -i $INDEX -o $OUTPUTDIR/$output -b30 $first ${first/_1.fastq.gz/_2.fastq.gz}

done
