#!/bin/sh

# ------------------------------------------------------------------------------
# --- Call PSMC
# ------------------------------------------------------------------------------

# Check that PSMCFA file was passed as parameter
USAGE="$0 consensus.psmcfa";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_PSMCFA=$1
OUT_PSMC=$(echo $IN_PSMCFA | sed -e 's/\.psmcfa/.psmc/')

echo "CMD: ${PSMC}/psmc \
	-N6 -t15 -r5 \
	-p \"4+25*2+4+6\" \
	-o ${OUT_PSMC} \
	${IN_PSMCFA}";

${PSMC}/psmc \
	-N6 -t15 -r5 \
	-p "4+25*2+4+6" \
	-o ${OUT_PSMC} \
	${IN_PSMCFA}


exit;
