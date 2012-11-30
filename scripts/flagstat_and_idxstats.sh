#!/bin/sh

# ------------------------------------------------------------------------------
# --- Analyze alignment output with flagstat and idxstats
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1

# Run flagstat
$SAMTOOLS/samtools flagstat \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam \
	> reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.txt

# Run idxstats
$SAMTOOLS/samtools idxstats \
	results/bwa.${GENOME_CODE}.sam.bam.sorted.bam \
	>> reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.txt

# Run bamtools stats


echo results/flagstat_and_idxstats_output.txt;

exit;