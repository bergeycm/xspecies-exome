#!/bin/sh

# ------------------------------------------------------------------------------
# --- Filter SNPs for quality
# ------------------------------------------------------------------------------

# Check that input raw BCF was passed as parameter
USAGE="$0 raw_snps.bcf";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_BCF=$1
OUT_VCF=$(echo $IN_BCF | sed 's/\.raw\.bcf/.flt.vcf/')

echo "${BCFTOOLS}/bcftools view ${IN_BCF} | \
	${BCFTOOLS}/vcfutils.pl varFilter -d 2 -D 100 \
	> ${OUT_VCF};";

${BCFTOOLS}/bcftools view ${IN_BCF} | \
	${BCFTOOLS}/vcfutils.pl varFilter -d 2 -D 100 \
	> ${OUT_VCF};

exit;

