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

# Secondary genome coordinates are liftOver'ed version of the human BED files.
# Chain file downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/*.over.chain.gz

# LiftOver commands:
# Syntax: liftOver oldFile chain.file newFile unMapped

echo "1 is $1\n";

# Add human line number to end of human targets, so we know which ones mapped
cat $1 | sed s'/$/\t/' | awk '{ print $0,NR }' > ${1}_w_ids

# Human target BED to secondary genome
#$LIFTOVER/liftOver -minMatch=0.1 \
#	${1}_w_ids \
#	$CHAINFILE \
#	${1}_2nd_liftover.bed \
#	${1}_2nd_liftover.unmapped.bed

# Output lines in original human targets that found matches in secondary genome
# In other words, output humans that have the same ID as is found in secondary genome
gawk 'FNR==NR{f1[$5];next}$5 in f1' ${1}_2nd_liftover.bed ${1}_w_ids > ${1}_hg_liftover.bed

# Figue out total number of regions that mapped
wc -l ${1}_hg_liftover.bed > results/liftOver_output.txt
wc -l ${1}_2nd_liftover.bed >> results/liftOver_output.txt

exit

