#!/bin/sh

# ------------------------------------------------------------------------------
# --- Call BSNP
# ------------------------------------------------------------------------------

# Check that input pileup file was passed as parameter
USAGE="$0 input.pileup";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_PILEUP=$1
OUT_SNP=$(echo $IN_PILEUP | sed 's/\.pileup/.bsnp.snp.out/')
OUT_AUX=$(echo $IN_PILEUP | sed 's/\.pileup/.bsnp.aux.out/')
OUT_SUM=$(echo $IN_PILEUP | sed 's/\.pileup/.bsnp.sum.out/')

echo "${BSNP}/BSNP \
	-i $IN_PILEUP \
	-o $OUT_SNP \
	-on $OUT_AUX \
	-os $OUT_SUM \
	-ph 0.5 -ka 1.0 -p0 0.001 -mq 1 -mp 0.000 \
	-th 0.85 -tm 0.00 -si 1 -ns 1 -ig 1 -v 1 -pb 1";

${BSNP}/BSNP \
	-i $IN_PILEUP \
	-o $OUT_SNP \
	-on $OUT_AUX \
	-os $OUT_SUM \
	-ph 0.5 -ka 1.0 -p0 0.001 -mq 1 -mp 0.000 \
	-th 0.85 -tm 0.00 -si 1 -ns 1 -ig 1 -v 1 -pb 1

exit;

