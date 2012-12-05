#!/bin/sh

# ------------------------------------------------------------------------------
# --- Get statistics on hybrid selection
# ------------------------------------------------------------------------------

# Check that input BAM and genome name were passed as parameters
USAGE="$0 input.bam genome_name";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BAM=$1
GENOME_NAME=$2

echo "java -Xmx10g -jar ${PICARD}/CalculateHsMetrics.jar \
	BAIT_INTERVALS=${IN_BAM}.picard.baits.bed \
	TARGET_INTERVALS=${IN_BAM}.picard.ccds.bed \
	INPUT=${IN_BAM} \
	OUTPUT=reports/${IND_ID}.bwa.${GENOME_NAME}.hsmetrics.txt";

java -Xmx10g -jar ${PICARD}/CalculateHsMetrics.jar \
	BAIT_INTERVALS=${IN_BAM}.picard.baits.bed \
	TARGET_INTERVALS=${IN_BAM}.picard.ccds.bed \
	INPUT=${IN_BAM} \
	OUTPUT=reports/${IND_ID}.bwa.${GENOME_NAME}.hsmetrics.txt

exit;




