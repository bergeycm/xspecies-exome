#!/bin/sh

# ------------------------------------------------------------------------------
# --- Filter mapped reads for quality
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

${BAMTOOLS}/bamtools filter \
	-mapQuality ">=20" \
	-in results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.nodup.RG.bam \
	-out results/${IND_ID}.bwa.${GENOME_CODE}.passed.bam

exit;