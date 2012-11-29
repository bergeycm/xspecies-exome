#!/bin/sh

# ------------------------------------------------------------------------------
# --- Merge overlapping intervals in BED file
# ------------------------------------------------------------------------------

# Check that BED file was passed as parameter
USAGE="$0 intervals.bed";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

echo "---${BEDTOOLS}/mergeBed -i $1  -nms > ${1}_MERGED;"
	
exit;