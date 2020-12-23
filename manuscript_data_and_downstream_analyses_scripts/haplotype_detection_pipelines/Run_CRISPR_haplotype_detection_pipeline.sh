#! /bin/bash

## Needs python 3+, and the following modules installed: argparse, re, os, sys, Biopython,  fuzzysearch, collections, pandas, numpy, gzip

python download_CRISPR_data.py -dir ./CRISPR_fastq
python call_CRISPR_haplotypes.py -i ./CRISPR_fastq -o ./CRISPR_haplos