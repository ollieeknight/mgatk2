# NUMT Blacklist Files

## Overview

These BED files contain genomic coordinates of **NUMT regions** (Nuclear MiTochondrial DNA sequences), which are segments of the nuclear genome that are homologous to mitochondrial DNA, which can cause false-positive alignments and contaminate single-cell mitochondrial DNA analysis. They are taken directly from the [mitoblacklist repo](https://github.com/caleblareau/mitoblacklist).

| File | Genome | Build | NUMT Regions | Description |
|------|--------|-------|--------------|-------------|
| `hg38_numts.bed` | Human | GRCh38/hg38 | 820 | Most recent human reference |
| `hg19_numts.bed` | Human | GRCh37/hg19 | 902 | Previous human reference |
| `mm10_numts.bed` | Mouse | GRCm38/mm10 | 523 | Most recent mouse reference |
| `mm9_numts.bed` | Mouse | GRCm37/mm9 | 3,415 | Previous mouse reference |

### How these files were generated

The blacklists were created by simulating mtDNA reads and identifying where they mis-align to nuclear chromosomes:

```bash
# Define reasonable parameters
THREADS=4
MERGEDIST=10
NREADS=10000000

# Set paths
REF=/path/to/bowtie2/reference
ENCODE_BL_FILE=/path/to/encode/blacklist

# Consider these carefully
MITOCHR="chrM" # update depending on your organism / reference genome
genome="hg19" # update depending on your organism / reference genome; 
# genome defines the prefix of output files; assumes "fasta/${genome}_mitochromosome.fa" exists

## DEPENDENCIES
# bowtie2
# macs2
# bedtools

# Step 1: Simulate Reads from mitochondrial genome 
./art_illumina -ss NS50 -c $NREADS -i "fasta/${genome}_mitochromosome.fa" -l 20 -o "fastq/${genome}_dat"
rm "fastq/${genome}_dat.aln"
gzip "fastq/${genome}_dat.fq"

# Step 2: Align to the full reference genome with bowtie2
bowtie2 -p $THREADS -x $REF -U "fastq/${genome}_dat.fq.gz" | \
  samtools view -bS - | \
  samtools sort -@ $THREADS - -o "bam/${genome}.all.sorted.bam"
samtools index "bam/${genome}.all.sorted.bam"

# Step 3: Extract non-mitochondrial chromosomes
chrs=`samtools view -H "bam/${genome}.all.sorted.bam" | grep SQ | cut -f2 | sed 's/SN://g' | grep -v $MITOCHR`

# Step 4: Generate clean .bam file (only nuclear alignments)
samtools view -b "bam/${genome}.all.sorted.bam" -o "bam/${genome}.nomito.bam" `echo $chrs`
samtools index "bam/${genome}.nomito.bam"

# Step 5: Call peaks on the mis-aligned reads
macs2 callpeak -t"bam/${genome}.nomito.bam" --nomodel --nolambda --keep-dup all -n "peaks/${genome}"

# Step 6: Create the full blacklist file
cat $ENCODE_BL_FILE "peaks/${genome}_peaks.narrowPeak" | \
  awk '{print $1"\t"$2"\t"$3}' | \
  sortBed | \
  mergeBed -i stdin -d $MERGEDIST > "combinedBlacklist/${genome}.full.blacklist.bed"
```

If you use these blacklists in your research, please cite:

**Original methodology:**
- Lareau CA, et al. (2023). mitoblacklist: Filtering mitochondrial DNA false-positive variants. GitHub repository: https://github.com/caleblareau/mitoblacklist