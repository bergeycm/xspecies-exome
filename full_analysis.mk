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
align : results/read*.bwa.*.sai
sampe : results/bwa.*.sam
sam2bam : results/bwa.*.sam.bam
sort_and_index_bam : results/bwa.*.sam.bam.sorted.bam
flagstat_idxstats : results/flagstat_and_idxstats_output.txt

all : ${HUMAN_GENOME_FA}i ${HUMAN_BWA_INDEX} ${2ND_GENOME_FA}i ${2ND_BWA_INDEX} ${TARGETS}_MERGED ${CCDS}_MERGED ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed ${TARGETS}_2nd_liftover.unmapped.bed ${CCDS}_2nd_liftover.unmapped.bed results/liftOver_output.txt ${TARGETS}_2nd_liftover.bed_MERGED ${CCDS}_2nd_liftover.bed_MERGED results/read*.bwa.*.sai results/bwa.*.sam results/bwa.*.sam.bam results/bwa.*.sam.bam.sorted.bam results/flagstat_and_idxstats_output.txt


# Hack to be able to export Make variables to child scripts
MAKE_ENV := $(shell echo '$(.VARIABLES)' | awk -v RS=' ' '/^[a-zA-Z0-9]+$$/')
SHELL_EXPORT := $(foreach v,$(MAKE_ENV),$(v)='$($(v))')

# -------------------------------------------------------------------------------------- #
# --- Index genome
# -------------------------------------------------------------------------------------- #

# The .fai output of samtools depends on the genome, BWA, samtools, & index_genome.sh
${HUMAN_GENOME_FA}i : ${HUMAN_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* scripts/index_genome.sh
	@echo "# === Indexing human genome =================================================== #";
	${SHELL_EXPORT} ./scripts/index_genome.sh ${HUMAN_GENOME_FA};
	@sleep 2
	@touch ${HUMAN_GENOME_FA}i ${HUMAN_BWA_INDEX}
${2ND_GENOME_FA}i : ${2ND_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* scripts/index_genome.sh
	@echo "# === Indexing secondary genome =============================================== #";
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
	@echo "# === Merging targets BED ===================================================== #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS};

${CCDS}_MERGED : ${BEDTOOLS}/* ${CCDS} scripts/merge_bed.sh
	@echo "# === Merging CCDS BED ======================================================== #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${CCDS};

# -------------------------------------------------------------------------------------- #
# --- LiftOver BED files of targets and CCDS into coordinates of the secondary genome
# -------------------------------------------------------------------------------------- #

# The LiftOver'd BED files depend on LiftOver, the chain file, the original BED, & liftOver_BED.sh
${TARGETS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${TARGETS} scripts/liftOver_BED.sh
	@echo "# === LiftingOver targets BED ================================================= #";
	${SHELL_EXPORT} ./scripts/liftOver_BED.sh ${TARGETS};
${CCDS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${CCDS} scripts/liftOver_BED.sh
	@echo "# === LiftingOver CCDS BED ==================================================== #";
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
# --- Get intersection of CCDS and targets
# -------------------------------------------------------------------------------------- #

# tr_2

# -------------------------------------------------------------------------------------- #
# --- Compute summary stats on targeted intervals, CCDS, and their intersection
# -------------------------------------------------------------------------------------- #

# tr_4

# -------------------------------------------------------------------------------------- #
# --- Analyze reads with FastQC. Total sequence bp, Maximum possible sequence depth
# -------------------------------------------------------------------------------------- #

# bep_1_1

# -------------------------------------------------------------------------------------- #
# --- Filter and trim reads
# -------------------------------------------------------------------------------------- #

# bep_1_2

# -------------------------------------------------------------------------------------- #
# --- Remove pairs whose partners were filtered, using cdbfasta & cdbyank
# -------------------------------------------------------------------------------------- #

# bep_1_3

# -------------------------------------------------------------------------------------- #
# --- [Optional] randomly subsample reads, if say you want to compare two different runs
# -------------------------------------------------------------------------------------- #

# bep_3

# -------------------------------------------------------------------------------------- #
# --- Align reads to genome with BWA
# -------------------------------------------------------------------------------------- #

# Alignment output (*.sai) depends on bwa, the reads FASTAs, the genome (index), and align.sh
results/read1.bwa.human.sai : ${BWA}/* ${READS1} ${READS2} ${HUMAN_GENOME_FA}i scripts/align.sh
	@echo "# === Aligning reads to human genome ========================================== #";
	${SHELL_EXPORT} ./scripts/align.sh ${HUMAN_GENOME_FA} human;
results/read1.bwa.other.sai : ${BWA}/* ${READS1} ${READS2} ${2ND_GENOME_FA}i scripts/align.sh
	@echo "# === Aligning reads to human genome ========================================== #";
	${SHELL_EXPORT} ./scripts/align.sh ${2ND_GENOME_FA} other;

# -------------------------------------------------------------------------------------- #
# --- Run sampe to generate SAM files
# -------------------------------------------------------------------------------------- #

# sampe output (*.sam) depends on *.sai files and sampe.sh
results/bwa.human.sam : results/read*.bwa.human.sai scripts/sampe.sh
	@echo "# === Combining reads to make SAM file for human genome ======================= #";
	${SHELL_EXPORT} ./scripts/sampe.sh ${HUMAN_GENOME_FA} human;
results/bwa.other.sam : results/read*.bwa.other.sai scripts/sampe.sh
	@echo "# === Combining reads to make SAM file for other genome ======================= #";
	${SHELL_EXPORT} ./scripts/sampe.sh ${2ND_GENOME_FA} other;

# -------------------------------------------------------------------------------------- #
# --- Convert SAM file to BAM file
# -------------------------------------------------------------------------------------- #

# BAM file depends on SAM file, samtools, genome .fai index, and scripts/sam2bam.sh
results/bwa.human.sam.bam : results/bwa.human.sam ${SAMTOOLS}/* ${HUMAN_GENOME_FA}i scripts/sam2bam.sh
	@echo "# === Converting SAM file to BAM file for human genome ======================== #";
	${SHELL_EXPORT} ./scripts/sam2bam.sh ${HUMAN_GENOME_FA}i human;
results/bwa.other.sam.bam : results/bwa.other.sam ${SAMTOOLS}/* ${2ND_GENOME_FA}i scripts/sam2bam.sh
	@echo "# === Converting SAM file to BAM file for other genome ======================== #";
	${SHELL_EXPORT} ./scripts/sam2bam.sh ${2ND_GENOME_FA}i other;

# -------------------------------------------------------------------------------------- #
# --- Sort and index BAM
# -------------------------------------------------------------------------------------- #

# Sorted BAM file depends on unsorted BAM file and scripts/sort_and_index_bam.sh
results/bwa.human.sam.bam.sorted.bam : results/bwa.human.sam.bam scripts/sort_and_index_bam.sh
	@echo "# === Sorting and Indexing BAM file for human genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_and_index_bam.sh human;
results/bwa.other.sam.bam.sorted.bam : results/bwa.other.sam.bam scripts/sort_and_index_bam.sh
	@echo "# === Sorting and Indexing BAM file for other genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_and_index_bam.sh other;

# -------------------------------------------------------------------------------------- #
# --- Analyze alignment output with flagstat and idxstats
# -------------------------------------------------------------------------------------- #

results/flagstat_and_idxstats_output.txt : results/bwa.human.sam.bam.sorted.bam results/bwa.other.sam.bam.sorted.bam scripts/flagstat_idxstats.sh
	@echo "# === Analyzing alignment output for human genome ============================= #";
	${SHELL_EXPORT} ./scripts/flagstat_idxstats.sh human;
	@echo "# === Analyzing alignment output for other genome ============================= #";
	${SHELL_EXPORT} ./scripts/flagstat_idxstats.sh other;

# -------------------------------------------------------------------------------------- #
# --- Calculate coverage of targeted regions
# -------------------------------------------------------------------------------------- #

# cov_0_1

# -------------------------------------------------------------------------------------- #
# --- Extract Number of bases with greater than 1x coverage, target region lengths, and 
# --- % coverages from coverageBed output
# -------------------------------------------------------------------------------------- #

# cov_0_1

# -------------------------------------------------------------------------------------- #
# --- Run coverageBed with histogram option
# -------------------------------------------------------------------------------------- #

# cov_0_2

# -------------------------------------------------------------------------------------- #
# --- Get simple percentage of target covered by at least one read
# -------------------------------------------------------------------------------------- #

# cov_1

# -------------------------------------------------------------------------------------- #
# --- Calculate avg coverage of intervals in CCDS
# -------------------------------------------------------------------------------------- #

# cov_2

# -------------------------------------------------------------------------------------- #
# --- Total target BP
# --- Report percentage and count of target bp covered by at least N reads
# --- Print bp count and percentage of target covered by N reads
# --- Get average percent coverage of each interval
# --- Average coverage of basepair in exome intervals
# -------------------------------------------------------------------------------------- #

# cov_2_1/summarize_covBed.R

# -------------------------------------------------------------------------------------- #
# For each sample, calculate number of targeted exons covered in their entirety
# at a miminum of 5x 10x and 20x coverage. Make plot and table.
# -------------------------------------------------------------------------------------- #

# cov_3

# -------------------------------------------------------------------------------------- #
# --- Get FASTA of sequences from reference genome that overlap targets
# -------------------------------------------------------------------------------------- #

# get_targeted_seqs_from_genomes

# -------------------------------------------------------------------------------------- #
# --- Align with Muscle, get pairwise distance with PAUP, get GC%, get number of indels
# -------------------------------------------------------------------------------------- #

# compare_hg_rhe_seqs.pl

# -------------------------------------------------------------------------------------- #
# --- Build linear model (avg. coverage as response variable and pairwise distance, 
# --- indels, exon length, GC\% as the predictors.
# -------------------------------------------------------------------------------------- #

# fac_6

# -------------------------------------------------------------------------------------- #
# --- George ABySS methods (not itemized yet)
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- SNP calling methods (not itemized yet)
# -------------------------------------------------------------------------------------- #



