#!/bin/sh

# ------------------------------------------------------------------------------
# --- Sort and Index BAM
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

# Sort BAM
# Output file is *.sorted.bam

$SAMTOOLS/samtools sort results/bwa.${GENOME_CODE}.sam.bam results/bwa.${GENOME_CODE}.sam.bam.sorted

# Index BAM

$SAMTOOLS/samtools index results/bwa.${GENOME_CODE}.sam.bam.sorted.bam

echo results/bwa.${GENOME_CODE}.sam.bam.sorted.bam;

exit;