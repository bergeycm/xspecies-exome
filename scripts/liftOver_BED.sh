#!/bin/sh

# ------------------------------------------------------------------------------
# --- LiftOver BED file coordinates
# ------------------------------------------------------------------------------

# Check that BED file was passed as parameter
USAGE="$0 intervals.bed";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

INPUT_BED=$1

# Secondary genome coordinates are liftOver'ed version of the human BED files.
# Chain file downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/*.over.chain.gz

# LiftOver commands:
# Syntax: liftOver oldFile chain.file newFile unMapped

# Add human line number to end of human targets, so we know which ones mapped
cat $INPUT_BED | sed s'/$/\t/' | awk '{ print $0,NR }' > ${INPUT_BED}_w_ids

# Human target BED to secondary genome
$LIFTOVER/liftOver -minMatch=0.1 \
	${INPUT_BED}_w_ids \
	$CHAINFILE \
	${INPUT_BED}_2nd_liftover.bed \
	${INPUT_BED}_2nd_liftover.unmapped.bed

# Output lines in original human targets that found matches in secondary genome
# In other words, output humans that have the same ID as is found in secondary genome
gawk 'FNR==NR{f1[$5];next}$5 in f1' ${INPUT_BED}_2nd_liftover.bed ${INPUT_BED}_w_ids > ${INPUT_BED}_hg_liftover.bed

# Figue out total number of regions that mapped
wc -l ${1}_hg_liftover.bed > results/liftOver_output.txt
wc -l ${1}_2nd_liftover.bed >> results/liftOver_output.txt

exit

