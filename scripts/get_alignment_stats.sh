#!/bin/sh

# ------------------------------------------------------------------------------
# --- Analyze alignment output with flagstat, idxstats, and bamtools stats
# ------------------------------------------------------------------------------

# Check that genome code was passed as parameter
USAGE="$0 genome_code";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

GENOME_CODE=$1
REPORT=reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.txt	

# TO DO: Filename needs to be dynamic so it can handle the later reports

echo "Alignment statistics for results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam" > $REPORT

# Run flagstat
echo "samtools flagstat:" >> $REPORT
$SAMTOOLS/samtools flagstat \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam \
	>> reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.txt
echo "" >> $REPORT

# Run idxstats
echo "samtools idxstats:" >> $REPORT
$SAMTOOLS/samtools idxstats \
	results/${IND_ID}.bwa.${GENOME_CODE}.sam.bam.sorted.bam \
	>> $REPORT
echo "" >> $REPORT

# Run bamtools stats
echo "bamtools stats:" >> $REPORT
echo "" >> $REPORT

echo $REPORT;

exit;