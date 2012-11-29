#!/bin/sh

# ------------------------------------------------------------------------------
# --- Convert SAM to BAM
# ------------------------------------------------------------------------------

# Check that genome FASTA and genome code were passed as parameters
USAGE="$0 genome.index.fai genome_code";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

# Infer genome code from genome FASTA filename
GENOME_CODE=$2

$SAMTOOLS/samtools view -b -t ${1} -o results/bwa.${GENOME_CODE}.sam.bam results/bwa.${GENOME_CODE}.sam
echo results/bwa.${GENOME_CODE}.sam.bam

exit;