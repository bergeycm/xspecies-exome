#!/bin/sh

# ------------------------------------------------------------------------------
# --- Run sampe to generate SAM file
# ------------------------------------------------------------------------------

# Check that genome FASTA and genome code were passed as parameters
USAGE="$0 genome.fasta genome_code";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

# Strip ending from $1 (fasta)
GENOME_PATH=$(echo $1 | sed 's/.[^.]*$//g')

# Infer genome code from genome FASTA filename
GENOME_CODE=$2

$BWA/bwa sampe $GENOME_PATH results/read1.bwa.${GENOME_CODE}.sai results/read2.bwa.${GENOME_CODE}.sai $READ1 $READ2 > results/bwa.${GENOME_CODE}.sam
echo results/bwa.${GENOME_CODE}.sam

exit;