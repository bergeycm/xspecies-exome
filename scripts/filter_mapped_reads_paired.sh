#!/bin/sh

# ------------------------------------------------------------------------------
# --- Filter mapped reads for mapping, pairing, and proper paired
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

${BAMTOOLS}/bamtools filter \
	-isMapped true \
	-isPaired true \
	-isProperPair true \
	-in results/${IND_ID}.bwa.${GENOME_CODE}.fixed.bam \
	-out results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.bam

exit;