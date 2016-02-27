#!/bin/bash
wc -l *.txt | sort -nr -k 1 | awk '{ if ($1 < 1) print $2 }' | cut -d '.' -f 1  > /gits/ensembl-assembly-exceptions/core_assemblies_without_patches.txt
