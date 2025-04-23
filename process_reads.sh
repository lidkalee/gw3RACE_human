#!/bin/bash

# Użycie: ./process_reads.sh READ1 READ2 PREFIX1 PREFIX2

# Sprawdzenie liczby argumentów
if [ "$#" -ne 4 ]; then
    echo "Użycie: $0 READ1 READ2 PREFIX1 PREFIX2"
    exit 1
fi

READ1=$1
READ2=$2
PREFIX1=$3
PREFIX2=$4

# Step 1: Align R1
echo "Aligning $READ1 with Bowtie2..."
bowtie2 -x fasta/index/ -U "$READ1" -S "sam/${PREFIX1}_mapped.sam" --no-unal

# Step 2: Extract aligned read names
echo "Extracting aligned read names..."
awk '{if ($1 !~ /^@/) print $1}' "sam/${PREFIX1}_mapped.sam" | sort | uniq > aligned_names.txt

# Step 3: Filter R2 based on aligned R1
echo "Filtering $READ2 based on aligned R1..."
seqtk subseq "$READ2" aligned_names.txt > "fastq/${PREFIX2}_filtered.fastq"

# Step 4: Grep adapter-containing reads
echo "Filtering adapter-containing reads in R2..."
grep -E -A 2 -B 1 --no-group-separator '^[[:alpha:]]{15}GTCAG' "fastq/${PREFIX2}_filtered.fastq" > "fastq/${PREFIX2}_grep_filtered.fastq"

# Step 5: Paired filtering
echo "Filtering paired reads..."
fastq_pair "fastq/${PREFIX2}_grep_filtered.fastq" "$READ1"
rm fastq/*.single.*

# Step 6: fastp filtering
echo "Running fastp filtering..."
fastp -i "fastq/${PREFIX1}.fastq.paired.fq" \
      -I "fastq/${PREFIX2}_grep_filtered.fastq.paired.fq" \
      -o "fastq/${PREFIX1}_fp.fq" \
      -O "fastq/${PREFIX2}_fp.fq" \
      -Q --umi --umi_loc=read2 --umi_len=15 --trim_front2=5 -l 80 \
      -h "fastq/${PREFIX1}_report.html"
