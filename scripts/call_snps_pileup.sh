#!/bin/sh

# ------------------------------------------------------------------------------
# --- Call SNPs with samtools (not m!) pileup
# ------------------------------------------------------------------------------

# Check that input BAM and genome were passed as parameters
USAGE="$0 input.bam genome.fa";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BAM=$1
OUT_PILEUP=$(echo $IN_BAM | sed 's/\.bam/.pileup/')
GENOME_FA=$2

# pileup is deprecated in the latest version of SAMtools, but BSNP requires its output
# The last version of SAMtools to include pileup is 0.1.16, so I have it installed
# and its location is listed in config.mk as $OLD_SAMTOOLS

echo "${OLD_SAMTOOLS}/samtools pileup \
	-csf ${GENOME_FA} \
	${IN_BAM} \
	> ${OUT_PILEUP}";

${OLD_SAMTOOLS}/samtools pileup \
	-csf ${GENOME_FA} \
	${IN_BAM} \
	> ${OUT_PILEUP}

exit;

