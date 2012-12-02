#!/bin/sh

# ------------------------------------------------------------------------------
# --- Fix Mate Pair Info
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

# Make temp folder
TMP_DIR=tmp/$RANDOM
mkdir -p $TMP_DIR

# Then fix mate pair info with Picard:
# Also shorten the filename
java -Djava.io.tmpdir=${TMP_DIR} \
	-jar ${PICARD}/FixMateInformation.jar \
	INPUT= results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam \
	OUTPUT=results/${IND_ID}.bwa.${GENOME_CODE}.fixed.bam \
	SO=coordinate \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=true

# Delete temp folder
rm -r $TMP_DIR

exit;