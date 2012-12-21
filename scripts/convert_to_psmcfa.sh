#!/bin/sh

# ------------------------------------------------------------------------------
# --- Convert FASTQ to PSMCFA
# ------------------------------------------------------------------------------

# Check that consensus FASTQ was passed as parameter
USAGE="$0 consensus.fq.gz";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_FQ=$1
OUT_PSMCFA=$(echo $IN_FQ | sed -e 's/\.consensus\.fq\.gz/.diploid.psmcfa/')

# fq2psmcfa has been modified to decrease. Line 101:
# if ((double)n_good_bases / len >= 0.2 && n_good_bases >= n_min_good) {
# has been changed to 
# if ((double)n_good_bases / len >= 0.002 && n_good_bases >= n_min_good) {

echo "CMD: ${PSMC}/utils/fq2psmcfa -q2 \
	${IN_FQ} \
	> $OUT_PSMCFA";

${PSMC}/utils/fq2psmcfa -q2 \
	${IN_FQ} \
	> $OUT_PSMCFA;

exit;
