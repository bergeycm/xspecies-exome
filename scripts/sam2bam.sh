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

INDEX_FAI=$1
GENOME_CODE=$2

$SAMTOOLS/samtools view \
	-b \
	-t ${INDEX_FAI} \
	-o results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam

echo results/${IND_ID}.bwa.${GENOME_CODE}.sam

exit;