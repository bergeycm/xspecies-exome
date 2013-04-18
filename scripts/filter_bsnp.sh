#!/bin/sh

# ------------------------------------------------------------------------------
# --- Filter BSNP output to basepairs with > 4 reads
# ------------------------------------------------------------------------------

# Check that input pileup file was passed as parameter
USAGE="$0 bsnp.out";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BSNP=$1

awk '{ if ($31 > 4) print }' $IN_BSNP;


exit;

