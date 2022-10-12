#!/bin/sh
set -e
fastq-dump *.sra
echo "unzip done"
mkdir SRA ; mv *.sra SRA
gzip *.fastq ; mkdir fastq ; mv *.gz fastq
echo "all done"

