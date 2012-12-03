#!/bin/sh

# ------------------------------------------------------------------------------
# --- Local realignment, step 1: ID realign targets
# ------------------------------------------------------------------------------

# Check that genome code and genome were passed as parameters
USAGE="$0 genome_code genome.fa";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1
GENOME_FA=$2

echo "CMD: java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${GENOME_FA} \
	-o results/${IND_ID}.bwa.${GENOME_CODE}.passed.bam.list \
	-I results/${IND_ID}.bwa.${GENOME_CODE}.passed.bam";

java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R ${GENOME_FA} \
	-o results/${IND_ID}.bwa.${GENOME_CODE}.passed.bam.list \
	-I results/${IND_ID}.bwa.${GENOME_CODE}.passed.bam

exit;