#!/bin/sh

# ------------------------------------------------------------------------------
# --- Call consensus sequence
# ------------------------------------------------------------------------------

# Check that input BAM, genome FASTA, genome name, and interval BED were passed as parameters
USAGE="$0 input.bam genome.fa genome_name intervals.bed";
if [ -z "$4" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BAM=$1
GENOME_FA=$2
GENOME_NAME=$3
INTERVALS=$4

echo "${SAMTOOLS}/samtools mpileup \
	-l ${INTERVALS} \
	-Auf ${GENOME_FA} \
	${IN_BAM} | \
	${BCFTOOLS}/bcftools view -c - | \
	${BCFTOOLS}/vcfutils.pl vcf2fq -d ${SNP_MIN_COV} -D ${SNP_MAX_COV} | \
	gzip > results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq.gz";

# -d = minimum read depth
# -D = maximum read depth
# -d should be one third of avg depth
# -D should be twice avg depth

${SAMTOOLS}/samtools mpileup \
	-l ${INTERVALS} \
	-Auf ${GENOME_FA} \
	${IN_BAM} | \
	${BCFTOOLS}/bcftools view -c - | \
	${BCFTOOLS}/vcfutils.pl vcf2fq -d ${SNP_MIN_COV} -D ${SNP_MAX_COV} \
	> results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq;

echo "gzip -cf results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq > results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq.gz";
gzip -cf results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq > results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq.gz

###	# Reduce consensus down to just targeted region
###	
###	echo "${BEDTOOLS}/fastaFromBed \
###		-fi results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq.gz \
###		-bed ${INTERVALS} \
###		-fo results/${IND_ID}.bwa.${GENOME_NAME}.consensus.targets.fq.gz"
###	
###	${BEDTOOLS}/fastaFromBed \
###		-fi results/${IND_ID}.bwa.${GENOME_NAME}.consensus.fq \
###		-bed ${INTERVALS} \
###		-fo results/${IND_ID}.bwa.${GENOME_NAME}.consensus.targets.fq

exit;

