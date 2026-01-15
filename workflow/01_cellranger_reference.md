# Cell Ranger reference building with hard-masked genomes

This guide shows how to build CellRanger-ARC reference genomes with NUMT-masked mitochondrial regions using mgatk2. We use CellRanger-ARC because it can function as both an ATAC and RNA reference.

## Prerequisites

```bash
# Required software
- mgatk2 (this package)
- cellranger-arc >= 2.0.2

# Install mgatk2
pip install git+https://github.com/ollieeknight/mgatk2.git
```

## Human (GRCh38) Reference

### Quick Start

```bash
#!/bin/bash

# Genome metadata
genome="GRCh38"
version="2024-A"

# Set up source and build directories
build="${genome}-${version}-build"
source="${genome}-${version}-reference-sources"
mkdir -p "$build" "$source"

# Download source files
# NOTE: Using Ensembl release 109 instead of 110 to avoid PAR issues
# Release 110 moved from GRCh38.p13 to GRCh38.p14, which unmasked the 
# pseudo-autosomal region causing ambiguous mappings to PAR locus genes
fasta_url="http://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v44.primary_assembly.annotation.gtf"
motifs_url="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
motifs_in="${source}/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi
if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi

# Modify FASTA headers to match GENCODE format
# Input:  >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
# Output: >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /; s/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffixes from gene/transcript/exon IDs
# Input:  gene_id "ENSG00000223972.5";
# Output: gene_id "ENSG00000223972"; gene_version "5";
gtf_modified="$build/$(basename "$gtf_in").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define biotype patterns for filtering
BIOTYPE_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Create gene allowlist (filter by biotype and exclude readthrough transcripts)
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

# Filter GTF based on gene allowlist and remove PAR_Y genes
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    | awk -F "\t" '$1 != "chrY" || $1 == "chrY" && $4 >= 2752083 && $4 < 56887903 && !/ENSG00000290840/' \
    >> "$gtf_filtered"

# Format motifs (NAME_ID format)
motifs_modified="$build/$(basename "$motifs_in").modified"
awk 'substr($1, 1, 1) == ">" { print ">" $2 "_" substr($1,2) } substr($1, 1, 1) != ">" { print }' \
    "$motifs_in" > "$motifs_modified"

# Hard-mask NUMT regions with mgatk2
fasta_hardmasked="${build}/$(basename "$fasta_in").hardmasked"
mgatk2 hardmask-fasta \
    -i "$fasta_modified" \
    -o "$fasta_hardmasked" \
    -g hg38

# Create CellRanger config
config_in="${build}/config"
echo """{
    organism: \"Homo_sapiens\"
    genome: [\"${genome}\"]
    input_fasta: [\"${fasta_hardmasked}\"]
    input_gtf: [\"${gtf_filtered}\"]
    input_motifs: \"${motifs_modified}\"
    non_nuclear_contigs: [\"chrM\"]
}""" > "$config_in"

# Build CellRanger reference
cellranger-arc mkref --ref-version="$version" \
    --config="$config_in" --nthreads=16

rm -rf ${build} ${source}
```


## Mouse (mm10/GRCm38) Reference

### Quick Start

```bash
#!/bin/bash

# Genome metadata
genome="GRCm38"
version="2024-A"

# Set up source and build directories
build="${genome}-${version}-build"
source="${genome}-${version}-reference-sources"
mkdir -p "$build" "$source"

# Download source files
# Using Ensembl release 98 for GRCm38 (mm10)
fasta_url="https://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Mus_musculus.GRCm38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.vM25.primary_assembly.annotation.gtf"
motifs_url="https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
motifs_in="${source}/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi
if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi

# Modify FASTA headers to match GENCODE format
# Input:  >1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF
# Output: >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /; s/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffixes from gene/transcript/exon IDs
gtf_modified="$build/$(basename "$gtf_in").modified"
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define biotype patterns for filtering
BIOTYPE_PATTERN=\
"(protein_coding|protein_coding_LoF|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Create gene allowlist (filter by biotype and exclude readthrough transcripts)
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

# Filter GTF based on gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
grep -Ff "${build}/gene_allowlist" "$gtf_modified" >> "$gtf_filtered"

# Format motifs (NAME_ID format)
motifs_modified="$build/$(basename "$motifs_in").modified"
awk 'substr($1, 1, 1) == ">" { print ">" $2 "_" substr($1,2) } substr($1, 1, 1) != ">" { print }' \
    "$motifs_in" > "$motifs_modified"

# Hard-mask NUMT regions with mgatk2
fasta_hardmasked="${build}/$(basename "$fasta_in").hardmasked"
mgatk2 hardmask-fasta \
    -i "$fasta_modified" \
    -o "$fasta_hardmasked" \
    -g mm10

# Create CellRanger config
config_in="${build}/config"
echo """{
    organism: \"Mus_musculus\"
    genome: [\"${genome}\"]
    input_fasta: [\"${fasta_hardmasked}\"]
    input_gtf: [\"${gtf_filtered}\"]
    input_motifs: \"${motifs_modified}\"
    non_nuclear_contigs: [\"chrM\"]
}""" > "$config_in"

# Build CellRanger reference
cellranger-arc mkref --ref-version="$version" \
    --config="$config_in" --nthreads=16

rm -rf ${build} ${source}
```