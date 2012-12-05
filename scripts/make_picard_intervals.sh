#!/bin/sh

# ------------------------------------------------------------------------------
# --- Convert BED to Picard intervals file
# ------------------------------------------------------------------------------

# Check that input BAM, baits BED file, and CCDS BED file were passed as parameters
USAGE="$0 input.bam targets.bed ccds.bed";
if [ -z "$3" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BAM=$1
TARGETS=$2
CCDS=$3

# Isolate headers from BAM file
echo "${SAMTOOLS}/samtools view \
	-H ${IN_BAM} \
	> ${IN_BAM}.header";

${SAMTOOLS}/samtools view \
	-H ${IN_BAM} \
	> ${IN_BAM}.header

# Copy header to new baits and CCDS intervals files
echo "cp ${IN_BAM}.header ${IN_BAM}.picard.baits.bed";
cp ${IN_BAM}.header ${IN_BAM}.picard.baits.bed
echo "cp ${IN_BAM}.header ${IN_BAM}.picard.ccds.bed";
cp ${IN_BAM}.header ${IN_BAM}.picard.ccds.bed

# Remove header
echo "rm ${IN_BAM}.header"
rm ${IN_BAM}.header

# Convert baits and CCDS intervals files
echo "gawk 'BEGIN { OFS=\"\t\" } { print $1,$2,$3,$6,$4 }' ${TARGETS} \
	>> ${IN_BAM}.picard.baits.bed";
gawk 'BEGIN { OFS="\t" } { print $1,$2,$3,$6,$4 }' ${TARGETS} \
	>> ${IN_BAM}.picard.baits.bed
echo "gawk 'BEGIN { OFS=\"\t\" } { print $1,$2,$3,$6,$4 }' ${CCDS} \
	>> ${IN_BAM}.picard.ccds.bed"
gawk 'BEGIN { OFS="\t" } { print $1,$2,$3,$6,$4 }' ${CCDS} \
	>> ${IN_BAM}.picard.ccds.bed

exit;
