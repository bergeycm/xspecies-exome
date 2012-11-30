#!/bin/sh

# ------------------------------------------------------------------------------
# --- Sort and Index BAM
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

# Sort BAM
# Output file is *.sorted.bam

$SAMTOOLS/samtools sort \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted

# Index BAM

$SAMTOOLS/samtools index \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam

echo results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam;

exit;