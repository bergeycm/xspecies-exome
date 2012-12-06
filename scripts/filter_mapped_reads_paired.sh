#!/bin/sh

# ------------------------------------------------------------------------------
# --- Filter mapped reads for mapping and pairing
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

# Removed "-isProperPair true"

${BAMTOOLS}/bamtools filter \
	-isMapped true \
	-isPaired true \
	-in results/${IND_ID}.bwa.${GENOME_CODE}.fixed.bam \
	-out results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.bam

exit;