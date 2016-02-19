#!/bin/bash

# -- SGE options (whose lines must begin with #$)

#$ -S /bin/bash       # Our jobscript is written for the bash shell
#$ -V                 # Inherit environment settings (e.g., from loaded modulefiles)
#$ -cwd               # Run the job in the current directory

# -- the commands to be executed (programs to be run) on a compute node:

module load apps/samtools/0.1.19/gcc-4.4.7

for readfile in "$@"; do
    readfilebase="$(basename $readfile .sam.gz)";
    newfilename="$readfilebase""_sorted.bam";
    echo "`date`: Writing $readfile to $newfilename";
    samtools view -bS "$readfile" | samtools sort - "$newfilename";
done

echo "`date` Finished writing BAM files"
