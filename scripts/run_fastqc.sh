#!/bin/sh

# ------------------------------------------------------------------------------
# --- Analyze FASTQ file of reads with FASTQC
# ------------------------------------------------------------------------------

# Check that reads file and output directory name were passed as parameters
USAGE="$0 reads.fq out_dir_name";
if [ -z "$2" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

READS_FQ=$1
OUT_FILE=$2

${FASTQC}/fastqc ${READS_FQ}

FQC_OUT=$(echo ${READS_FQ} | sed 's/\.fastq/\_fastqc/')

# Move output files into reports directory and rename them
mv ${FQC_OUT} reports/${OUT_FILE}
mv ${FQC_OUT}.zip reports/${OUT_FILE}.zip

exit;