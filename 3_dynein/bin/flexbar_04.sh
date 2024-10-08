#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 2 GiB of RAM per slot
#$ -l h_rt=24:00:00   # job requires up to x hours of runtime


## 0. In case TMPDIR is not set, e.g. on development nodes, set
##    it to local /scratch, if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

## Use temporary working directory
cd "$TMPDIR"

# flexbar loop to trim NEBNext Single-cell adapters

module load CBI miniconda3/23.5.2-0-py311

conda activate flexbar

FILES="/scratch/aralbright/AAWM04/"

OUTPUTDIR="con_trimmed"

for first in $FILES*R1_001.fastq.gz
do

  output=$(basename $first | cut -c48-52)
  echo "$OUTPUTDIR/$output"
  echo "$first"
  echo "${first/_R1_001.fastq.gz/_R2_001.fastq.gz}"

  flexbar --reads $first \
          --reads2 ${first/_R1_001.fastq.gz/_R2_001.fastq.gz} \
          --stdout-reads \
          --adapters /wynton/home/marshall/aralbright/dynein/bin/tso_g_wo_hp.fasta \
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
          --threads 8 | \

          flexbar \
          --reads - \
          --interleaved \
          --target $OUTPUTDIR/$output \
          --adapters /wynton/home/marshall/aralbright/dynein/bin/ilmn_20_2_seqs.fasta \
          --adapter-trim-end RIGHT \
          --min-read-length 2 \
          --threads 8

   # Check if the flexbar command was successful
  if [ $? -ne 0 ]; then
    echo "Flexbar command failed for $first"
    exit 1
  fi

done

