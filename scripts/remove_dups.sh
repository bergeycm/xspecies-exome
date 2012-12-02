#!/bin/sh

# ------------------------------------------------------------------------------
# --- Remove Duplicates
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

# Then remove duplicates with Picard:

echo "CMD: java -Djava.io.tmpdir=${TMP_DIR} \
	-jar ${PICARD}/MarkDuplicates.jar \
	INPUT=results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.bam \
	OUTPUT=results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.nodup.bam \
	M=reports/duplicate_report.txt \
	VALIDATION_STRINGENCY=SILENT \
	REMOVE_DUPLICATES=true";

java -Djava.io.tmpdir=${TMP_DIR} \
	-jar ${PICARD}/MarkDuplicates.jar \
	INPUT=results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.bam \
	OUTPUT=results/${IND_ID}.bwa.${GENOME_CODE}.fixed.filtered.nodup.bam \
	M=reports/duplicate_report.txt \
	VALIDATION_STRINGENCY=SILENT \
	REMOVE_DUPLICATES=true

# Delete temp folder
rm -r $TMP_DIR

exit;