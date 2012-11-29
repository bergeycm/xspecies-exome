# -------------------------------------------------------------------------------------- #
# --- Makefile to run exome pipeline. 
# --- Called by the executable shell script, full_analysis
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

SHELL = /bin/bash

HUMAN_GENOME_DIR=$(dir ${HUMAN_GENOME_FA})
2ND_GENOME_DIR=$(dir ${2ND_GENOME_FA})

HUMAN_GENOME_CODE=
2ND_GENOME_CODE=

# Output files of BWA index
BWA_INDEX_ENDINGS = .amb .ann .bwt .pac .sa
PROTO_HUMAN_BWA_INDEX = $(addprefix ${HUMAN_GENOME_FA}, ${BWA_INDEX_ENDINGS})
HUMAN_BWA_INDEX = $(subst .fa,,${PROTO_HUMAN_BWA_INDEX})
PROTO_2ND_BWA_INDEX = $(addprefix ${2ND_GENOME_FA}, ${BWA_INDEX_ENDINGS})
2ND_BWA_INDEX = $(subst .fa,,${PROTO_2ND_BWA_INDEX})

index_genome : ${HUMAN_GENOME_FA}i ${HUMAN_BWA_INDEX} ${2ND_GENOME_FA}i ${2ND_BWA_INDEX}
merge_beds : ${TARGETS}_MERGED ${CCDS}_MERGED
liftover_beds : ${TARGETS} ${CCDS} ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed ${TARGETS}_2nd_liftover.unmapped.bed ${CCDS}_2nd_liftover.unmapped.bed results/liftOver_output.txt ${TARGETS}_2nd_liftover.bed_MERGED ${CCDS}_2nd_liftover.bed_MERGED

all : ${HUMAN_GENOME_FA}i ${HUMAN_BWA_INDEX} ${2ND_GENOME_FA}i ${2ND_BWA_INDEX} ${TARGETS}_MERGED ${CCDS}_MERGED ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed ${TARGETS}_2nd_liftover.unmapped.bed ${CCDS}_2nd_liftover.unmapped.bed results/liftOver_output.txt ${TARGETS}_2nd_liftover.bed_MERGED ${CCDS}_2nd_liftover.bed_MERGED

# Hack to be able to export Make variables to child scripts
MAKE_ENV := $(shell echo '$(.VARIABLES)' | awk -v RS=' ' '/^[a-zA-Z0-9]+$$/')
SHELL_EXPORT := $(foreach v,$(MAKE_ENV),$(v)='$($(v))')

# -------------------------------------------------------------------------------------- #
# --- Index genome
# -------------------------------------------------------------------------------------- #

# The .fai output of samtools depends on the genome, BWA, samtools, & index_genome.sh
${HUMAN_GENOME_FA}i : ${HUMAN_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* scripts/index_genome.sh
	@echo "# === Indexing human genome === #";
	${SHELL_EXPORT} ./scripts/index_genome.sh ${HUMAN_GENOME_FA};
	@sleep 2
	@touch ${HUMAN_GENOME_FA}i ${HUMAN_BWA_INDEX}
${2ND_GENOME_FA}i : ${2ND_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* scripts/index_genome.sh
	@echo "# === Indexing secondary genome === #";
	${SHELL_EXPORT} ./scripts/index_genome.sh ${2ND_GENOME_FA};
	@sleep 2
	@touch ${2ND_GENOME_FA}i ${2ND_BWA_INDEX}

# The output files of bwa depend on the output of samtools.
# A hack to deal with the problem make has with multiple targets dependent on one rule
# See for details:
# http://www.cmcrossroads.com/ask-mr-make/12908-rules-with-multiple-outputs-in-gnu-make
${HUMAN_BWA_INDEX} : ${HUMAN_GENOME_FA}i
${2ND_BWA_INDEX} : ${2ND_GENOME_FA}i
	
# -------------------------------------------------------------------------------------- #
# --- Merge overlapping intervals in BED files of targets and CCDS
# -------------------------------------------------------------------------------------- #

# The BED file of merged target intervals depends on the original targets BED, bedtools, & merge_bed.sh
${TARGETS}_MERGED : ${BEDTOOLS}/* ${TARGETS} scripts/merge_bed.sh
	@echo "# === Merging targets BED === #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS};

${CCDS}_MERGED : ${BEDTOOLS}/* ${CCDS} scripts/merge_bed.sh
	@echo "# === Merging CCDS BED === #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${CCDS};

# -------------------------------------------------------------------------------------- #
# --- LiftOver BED files of targets and CCDS into coordinates of the secondary genome
# -------------------------------------------------------------------------------------- #

# The LiftOver'd BED files depend on LiftOver, the chain file, the original BED, & liftOver_BED.sh
${TARGETS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${TARGETS} scripts/liftOver_BED.sh
	@echo "# === LiftingOver targets BED === #";
	${SHELL_EXPORT} ./scripts/liftOver_BED.sh ${TARGETS};
${CCDS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${CCDS} scripts/liftOver_BED.sh
	@echo "# === LiftingOver CCDS BED === #";
	${SHELL_EXPORT} ./scripts/liftOver_BED.sh ${CCDS};

# Other output files from LiftingOvering depend on the above mentioned output BED file
${TARGETS}_2nd_liftover.unmapped.bed ${TARGETS}_hg_liftover.bed : ${TARGETS}_2nd_liftover.bed
${CCDS}_2nd_liftover.unmapped.bed ${CCDS}_hg_liftover.bed : ${CCDS}_2nd_liftover.bed

# As does the LiftOver output file
results/liftOver_output.txt : ${TARGETS}_2nd_liftover.bed
results/liftOver_output.txt : ${CCDS}_2nd_liftover.bed

# Final merged versions depend on original BED files, plus bedtools and merge_bed.sh
${TARGETS}_2nd_liftover.bed_MERGED : ${BEDTOOLS}/* ${TARGETS}_2nd_liftover.bed scripts/merge_bed.sh
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS}_2nd_liftover.bed;
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS}_hg_liftover.bed;

${CCDS}_2nd_liftover.bed_MERGED : ${BEDTOOLS}/* ${CCDS}_2nd_liftover.bed scripts/merge_bed.sh
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${CCDS}_2nd_liftover.bed;
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${CCDS}_hg_liftover.bed;

# -------------------------------------------------------------------------------------- #
# --- Filter and trim reads
# -------------------------------------------------------------------------------------- #


# -------------------------------------------------------------------------------------- #
# --- Align reads to genome with BWA
# -------------------------------------------------------------------------------------- #

# Alignment output (*.sai) depends on bwa, the reads FASTAs, the genome (index), and align.sh
results/read1.bwa.human.sai : ${BWA}/* ${READS1} ${READS2} ${HUMAN_GENOME_FA}i scripts/align.sh
	@echo "# === Aligning reads to human genome === #";
	${SHELL_EXPORT} ./scripts/align.sh ${HUMAN_GENOME_FA} human;
results/read1.bwa.other.sai : ${BWA}/* ${READS1} ${READS2} ${2ND_GENOME_FA}i scripts/align.sh
	@echo "# === Aligning reads to human genome === #";
	${SHELL_EXPORT} ./scripts/align.sh ${2ND_GENOME_FA} other;

# -------------------------------------------------------------------------------------- #
# --- sampe
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Convert SAM file to BAM file
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Sort and index BAM
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Get info with flagstat and idxstats
# -------------------------------------------------------------------------------------- #

