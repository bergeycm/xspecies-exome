#!/bin/sh

# ------------------------------------------------------------------------------
# --- Generate the ms command line that simulates the history inferred by PSMC
# --- Also plot the result
# ------------------------------------------------------------------------------

# Check that PSMC file was passed as parameter
USAGE="$0 consensus.psmc";
if [ -z "$1" ]; then
	echo "ERROR: $USAGE";
	exit 1;
fi

IN_PSMC=$1
OUT_MS=$(echo $IN_PSMC | sed -e 's/\.psmc/.ms_cmd.sh/')
OUT_PLOT=$(echo $IN_PSMC | sed -e 's/\.psmc/.plot/')


echo "CMD: ${PSMC}/utils/psmc2history.pl \
	${IN_PSMC} \
	| ${PSMC}/utils/history2ms.pl \
	> ${OUT_MS}";

${PSMC}/utils/psmc2history.pl \
	${IN_PSMC} \
	| utils/history2ms.pl \
	> ${OUT_MS}

echo "CMD: ${PSMC}/utils/psmc_plot.pl ${OUT_PLOT} ${IN_PSMC}";

${PSMC}/utils/psmc_plot.pl ${OUT_PLOT} ${IN_PSMC}

exit;



