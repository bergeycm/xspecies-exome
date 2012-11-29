#!/bin/sh

# ------------------------------------------------------------------------------
# --- Align reads to genome with BWA
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

$BWA/bwa aln $BWA_ALN_PARAM $GENOME_PATH $READ1 > results/read1.bwa.${GENOME_CODE}.sai
echo results/read1.bwa.${GENOME_CODE}.sai;

$BWA/bwa aln $BWA_ALN_PARAM $GENOME_PATH $READ2 > results/read2.bwa.${GENOME_CODE}.sai
echo results/read2.bwa.${GENOME_CODE}.sai;

exit;