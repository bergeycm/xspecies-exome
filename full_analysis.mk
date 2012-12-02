# -------------------------------------------------------------------------------------- #
# --- Makefile to run exome pipeline. 
# --- Called by the executable shell script, full_analysis
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

HUMAN_GENOME_DIR=$(dir ${HUMAN_GENOME_FA})
SECOND_GENOME_DIR=$(dir ${SECOND_GENOME_FA})

# Output files of BWA index. None of these variables are exported since they begin with "_"
_BWA_INDEX_ENDINGS = .amb .ann .bwt .pac .sa
_PROTO_HUMAN_BWA_INDEX = $(addprefix ${HUMAN_GENOME_FA}, ${BWA_INDEX_ENDINGS})
_HUMAN_BWA_INDEX = $(subst .fa,,${PROTO_HUMAN_BWA_INDEX})
_PROTO_SECOND_BWA_INDEX = $(addprefix ${SECOND_GENOME_FA}, ${BWA_INDEX_ENDINGS})
_SECOND_BWA_INDEX = $(subst .fa,,${PROTO_SECOND_BWA_INDEX})

# Steps. Can be called one-by-one with something like, make index_genome
# --- preliminary_steps
index_genome : ${HUMAN_GENOME_FA}i ${_HUMAN_BWA_INDEX} ${SECOND_GENOME_FA}i ${_SECOND_BWA_INDEX}
merge_beds : ${TARGETS}_MERGED ${CCDS}_MERGED
liftover_beds : ${TARGETS} ${CCDS} ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed ${TARGETS}_2nd_liftover.unmapped.bed ${CCDS}_2nd_liftover.unmapped.bed results/liftOver_output.txt ${TARGETS}_2nd_liftover.bed_MERGED ${CCDS}_2nd_liftover.bed_MERGED 
# --- pre_aln_filtering_steps:
fastqc_raw : reports/${IND_ID}.read1.raw.stats.zip reports/${IND_ID}.read2.raw.stats.zip
filter_reads : ${READ1}.filtered.fastq ${READ2}.filtered.fastq
fastqc_filtered : reports/${IND_ID}.read1.filtered.stats.zip reports/${IND_ID}.read2.filtered.stats.zip
# --- alignment_steps
align : results/${IND_ID}.read1.bwa.human.sai results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai
sampe : results/${IND_ID}.bwa.human.sam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam
sam2bam : results/${IND_ID}.bwa.human.sam.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam
sort_and_index_bam : results/${IND_ID}.bwa.human.sam.bam.sorted.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam
get_alignment_stats : reports/${IND_ID}.bwa.human.aln_stats.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt
# --- post_alignment_filtering_steps
fix_mate_pairs : results/${IND_ID}.bwa.human.fixed.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam
filter_unmapped : results/${IND_ID}.bwa.human.fixed.filtered.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam
remove_dups : results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam
add_read_groups : results/${IND_ID}.bwa.human.fixed.filtered.nodup.RG.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.RG.bam
filter_bad_qual : results/${IND_ID}.bwa.human.passed.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam

# Group steps together
preliminary_steps : index_genome merge_beds liftover_beds
pre_aln_filtering_steps : fastqc_raw filter_reads fastqc_filtered
alignment_steps : align sampe sam2bam sort_and_index_bam get_alignment_stats
post_alignment_filtering_steps : fix_mate_pairs filter_unmapped remove_dups add_read_groups filter_bad_qual

all : preliminary_steps pre_aln_filtering_steps alignment_steps post_alignment_filtering_steps

# Hack to be able to export Make variables to child scripts
# Don't export variables from make that begin with non-alphanumeric character
# After that, underscores are OK
MAKE_ENV := $(shell echo '$(.VARIABLES)' | awk -v RS=' ' '/^[a-zA-Z0-9][a-zA-Z0-9_]+$$/')
SHELL_EXPORT := $(foreach v,$(MAKE_ENV),$(v)='$($(v))')

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Preliminary Steps
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Index genome
# -------------------------------------------------------------------------------------- #

# The .fai output of samtools depends on the genome, BWA, samtools, & index_genome.sh
${HUMAN_GENOME_FA}i : ${HUMAN_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* #scripts/index_genome.sh
	@echo "# === Indexing human genome =================================================== #";
	${SHELL_EXPORT} ./scripts/index_genome.sh ${HUMAN_GENOME_FA};
	@sleep 2
	@touch ${HUMAN_GENOME_FA}i ${_HUMAN_BWA_INDEX}
${SECOND_GENOME_FA}i : ${SECOND_GENOME_FA} ${BWA}/* ${SAMTOOLS}/* #scripts/index_genome.sh
	@echo "# === Indexing secondary genome =============================================== #";
	${SHELL_EXPORT} ./scripts/index_genome.sh ${SECOND_GENOME_FA};
	@sleep 2
	@touch ${SECOND_GENOME_FA}i ${_SECOND_BWA_INDEX}

# The output files of bwa depend on the output of samtools.
# A hack to deal with the problem make has with multiple targets dependent on one rule
# See for details:
# http://www.cmcrossroads.com/ask-mr-make/12908-rules-with-multiple-outputs-in-gnu-make
${_HUMAN_BWA_INDEX} : ${HUMAN_GENOME_FA}i
${_SECOND_BWA_INDEX} : ${SECOND_GENOME_FA}i
	
# -------------------------------------------------------------------------------------- #
# --- Merge overlapping intervals in BED files of targets and CCDS
# -------------------------------------------------------------------------------------- #

# The BED file of merged target intervals depends on the original targets BED, bedtools, & merge_bed.sh
${TARGETS}_MERGED : ${BEDTOOLS}/* ${TARGETS} #scripts/merge_bed.sh
	@echo "# === Merging targets BED ===================================================== #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS};

${CCDS}_MERGED : ${BEDTOOLS}/* ${CCDS} #scripts/merge_bed.sh
	@echo "# === Merging CCDS BED ======================================================== #";
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${CCDS};

# -------------------------------------------------------------------------------------- #
# --- LiftOver BED files of targets and CCDS into coordinates of the secondary genome
# -------------------------------------------------------------------------------------- #

# The LiftOver'd BED files depend on LiftOver, the chain file, the original BED, & liftOver_BED.sh
${TARGETS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${TARGETS} #scripts/liftOver_BED.sh
	@echo "# === LiftingOver targets BED ================================================= #";
	${SHELL_EXPORT} ./scripts/liftOver_BED.sh ${TARGETS};
${CCDS}_2nd_liftover.bed : ${LIFTOVER} ${CHAIN} ${CCDS} #scripts/liftOver_BED.sh
	@echo "# === LiftingOver CCDS BED ==================================================== #";
	${SHELL_EXPORT} ./scripts/liftOver_BED.sh ${CCDS};

# Other output files from LiftingOvering depend on the above mentioned output BED file
${TARGETS}_2nd_liftover.unmapped.bed ${TARGETS}_hg_liftover.bed : ${TARGETS}_2nd_liftover.bed
${CCDS}_2nd_liftover.unmapped.bed ${CCDS}_hg_liftover.bed : ${CCDS}_2nd_liftover.bed

# As does the LiftOver output file
results/liftOver_output.txt : ${TARGETS}_2nd_liftover.bed
results/liftOver_output.txt : ${CCDS}_2nd_liftover.bed

# Final merged versions depend on original BED files, plus bedtools and merge_bed.sh
${TARGETS}_2nd_liftover.bed_MERGED : ${BEDTOOLS}/* ${TARGETS}_2nd_liftover.bed #scripts/merge_bed.sh
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS}_2nd_liftover.bed;
	${SHELL_EXPORT} ./scripts/merge_bed.sh ${TARGETS}_hg_liftover.bed;

${CCDS}_2nd_liftover.bed_MERGED : ${BEDTOOLS}/* ${CCDS}_2nd_liftover.bed #scripts/merge_bed.sh
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

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Read filtering
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Analyze reads with FastQC. Total sequence bp, Maximum possible sequence depth
# -------------------------------------------------------------------------------------- #

# FastQC reports depend on read files, FastQC, and run_fastqc.sh
reports/${IND_ID}.read1.raw.stats.zip : ${READ1} ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (1st pair) before filtering ================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ1} ${IND_ID}.read1.raw.stats;
reports/${IND_ID}.read2.raw.stats.zip : ${READ2} ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (2nd pair) before filtering ================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ2} ${IND_ID}.read2.raw.stats;
	
# -------------------------------------------------------------------------------------- #
# --- Filter and trim reads
# -------------------------------------------------------------------------------------- #

# Filtered reads FASTQs depends on read files, FastX, and filter_reads.sh
${READ1}.filtered.fastq : ${READ1} ${FASTX}/* #scripts/filter_reads.sh
	@echo "# === Filtering reads (1st pair) with FastX =================================== #";
	${SHELL_EXPORT} ./scripts/filter_reads.sh ${READ1};
${READ2}.filtered.fastq : ${READ2} ${FASTX}/* #scripts/filter_reads.sh
	@echo "# === Filtering reads (2nd pair) with FastX =================================== #";
	${SHELL_EXPORT} ./scripts/filter_reads.sh ${READ2};

# Call FastQC again
reports/${IND_ID}.read1.filtered.stats.zip : ${READ1}.filtered.fastq ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (1st pair) after filtering =================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ1}.filtered.fastq ${IND_ID}.read1.filtered.stats;
reports/${IND_ID}.read2.filtered.stats.zip : ${READ2}.filtered.fastq ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (2nd pair) after filtering =================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ2}.filtered.fastq ${IND_ID}.read2.filtered.stats;

# -------------------------------------------------------------------------------------- #
# --- Remove pairs whose partners were filtered, using cdbfasta & cdbyank
# -------------------------------------------------------------------------------------- #

# Is this really necessary? We can get rid of them after mapping with bamtools filter

# bep_1_3

# Call fastqc again

# -------------------------------------------------------------------------------------- #
# --- [Optional] randomly subsample reads, if say you want to compare two different runs
# -------------------------------------------------------------------------------------- #

# bep_3
# Call fastqc again

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mapping to reference genomes
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Align reads to genome with BWA
# -------------------------------------------------------------------------------------- #

# Alignment output (*.sai) depends on bwa, the filtered reads FASTAs, the genome (index), and align.sh
# Using the first read as a stand in for the both 
results/${IND_ID}.read1.bwa.human.sai : ${BWA}/* ${READ1}.filtered.fastq ${READ2}.filtered.fastq ${HUMAN_GENOME_FA}i #scripts/align.sh
	@echo "# === Aligning reads to human genome ========================================== #";
	${SHELL_EXPORT} ./scripts/align.sh ${HUMAN_GENOME_FA} human;
results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai : ${BWA}/* ${READ1}.filtered.fastq ${READ2}.filtered.fastq ${SECOND_GENOME_FA}i #scripts/align.sh
	@echo "# === Aligning reads to other genome ========================================== #";
	${SHELL_EXPORT} ./scripts/align.sh ${SECOND_GENOME_FA} ${SECOND_GENOME_NAME};

# Read 2 depends on read 1
results/${IND_ID}.read2.bwa.human.sai : results/${IND_ID}.read1.bwa.human.sai
results/${IND_ID}.read2.bwa.${SECOND_GENOME_NAME}.sai : results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai

# -------------------------------------------------------------------------------------- #
# --- Run sampe to generate SAM files
# -------------------------------------------------------------------------------------- #

# sampe output (*.sam) depends on *.sai files and sampe.sh
# Using the first read as a stand in for the both
results/${IND_ID}.bwa.human.sam : results/${IND_ID}.read1.bwa.human.sai #scripts/sampe.sh
	@echo "# === Combining reads to make SAM file for human genome ======================= #";
	${SHELL_EXPORT} ./scripts/sampe.sh ${HUMAN_GENOME_FA} human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam : results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai #scripts/sampe.sh
	@echo "# === Combining reads to make SAM file for other genome ======================= #";
	${SHELL_EXPORT} ./scripts/sampe.sh ${SECOND_GENOME_FA} ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Convert SAM file to BAM file
# -------------------------------------------------------------------------------------- #

# BAM file depends on SAM file, samtools, genome .fai index, and scripts/sam2bam.sh
results/${IND_ID}.bwa.human.sam.bam : results/${IND_ID}.bwa.human.sam ${SAMTOOLS}/* ${HUMAN_GENOME_FA}i #scripts/sam2bam.sh
	@echo "# === Converting SAM file to BAM file for human genome ======================== #";
	${SHELL_EXPORT} ./scripts/sam2bam.sh ${HUMAN_GENOME_FA}i human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam ${SAMTOOLS}/* ${SECOND_GENOME_FA}i #scripts/sam2bam.sh
	@echo "# === Converting SAM file to BAM file for other genome ======================== #";
	${SHELL_EXPORT} ./scripts/sam2bam.sh ${SECOND_GENOME_FA}i ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Sort and index BAM
# -------------------------------------------------------------------------------------- #

# Sorted BAM file depends on unsorted BAM file, scripts/sort_bam, and scripts/index_bam.sh
results/${IND_ID}.bwa.human.sam.bam.sorted.bam : results/${IND_ID}.bwa.human.sam.bam #scripts/sort_bam scripts/index_bam.sh
	@echo "# === Sorting and Indexing BAM file for human genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_bam.sh results/${IND_ID}.bwa.human.sam.bam;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.sam.bam.sorted.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam #scripts/sort_bam scripts/index_bam.sh
	@echo "# === Sorting and Indexing BAM file for other genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam;

# -------------------------------------------------------------------------------------- #
# --- Analyze alignment output with flagstat, idxstats, and stats
# -------------------------------------------------------------------------------------- #

# Align stats report depends on the sorted BAM and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.txt : results/${IND_ID}.bwa.human.sam.bam.sorted.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome ============================= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh human;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome ============================= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh ${SECOND_GENOME_NAME};

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Post-alignment filtering steps
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Fix mate pairs info
# -------------------------------------------------------------------------------------- #

# BAM with fixed mate pair info depends on output BAM from sort_and_index.sh, Picard, and scripts/fix_mate_pairs.sh
results/${IND_ID}.bwa.human.fixed.bam : results/${IND_ID}.bwa.human.sam.bam.sorted.bam ${PICARD}/* # scripts/fix_mate_pairs.sh
	@echo "# === Fixing mate pair info for human genome alignment ======================== #";
	${SHELL_EXPORT} ./scripts/fix_mate_pairs.sh human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam ${PICARD}/* # scripts/fix_mate_pairs.sh
	@echo "# === Fixing mate pair info for other genome alignment ======================== #";
	${SHELL_EXPORT} ./scripts/fix_mate_pairs.sh ${SECOND_GENOME_NAME};

# Run flagstat, idxstats, bedtools stats. reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.pairsfix.txt

# -------------------------------------------------------------------------------------- #
# --- Filtering for mapping, pairing, and proper paired
# -------------------------------------------------------------------------------------- #

# Filtered BAM depends on output BAM from fix_mate_pairs.sh, BAMtools, and scripts/filter_mapped_reads_paired.sh
results/${IND_ID}.bwa.human.fixed.filtered.bam : results/${IND_ID}.bwa.human.fixed.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_paired.sh
	@echo "# === Filtering unpaired reads mapped to human genome ========================= #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_paired.sh human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_paired.sh
	@echo "# === Filtering reads mapped to other genome ================================== #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_paired.sh ${SECOND_GENOME_NAME};

# Run flagstat, idxstats, bedtools stats. reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.pairsfix.fltr.txt

# -------------------------------------------------------------------------------------- #
# --- Remove duplicates
# -------------------------------------------------------------------------------------- #

# BAM sans dups depends on output BAM from filter_mapped_reads_paired.sh, Picard, and scripts/remove_dups.sh
results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam : results/${IND_ID}.bwa.human.fixed.filtered.bam ${PICARD}/* # scripts/remove_dups.sh
	@echo "# === Removing duplicate reads mapped to human genome ================================== #";
	${SHELL_EXPORT} ./scripts/remove_dups.sh human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam ${PICARD}/* # scripts/remove_dups.sh
	@echo "# === Removing duplicate reads mapped to other genome ================================== #";
	${SHELL_EXPORT} ./scripts/remove_dups.sh ${SECOND_GENOME_NAME};

# Run flagstat, idxstats, bedtools stats. reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.pairsfix.fltr.nodups.txt

# -------------------------------------------------------------------------------------- #
# --- Add read groups
# -------------------------------------------------------------------------------------- #

# BAM without RGs depends on output BAM from remove_dups.sh, Picard, and scripts/add_read_groups.sh
results/${IND_ID}.bwa.human.fixed.filtered.nodup.RG.bam : results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam ${PICARD}/* # scripts/add_read_groups.sh
	@echo "# === Adding read groups for reads mapped to human genome ===================== #";
	${SHELL_EXPORT} ./scripts/add_read_groups.sh human;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.RG.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam ${PICARD}/* # scripts/add_read_groups.sh
	@echo "# === Adding read groups for reads mapped to other genome ===================== #";
	${SHELL_EXPORT} ./scripts/add_read_groups.sh ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Remove reads with low mapping quality
# -------------------------------------------------------------------------------------- #

# Filtered BAM depends on output BAM from add_read_groups.sh, BAMtools, and scripts/filter_mapped_reads_quality.sh
results/${IND_ID}.bwa.human.passed.bam : results/${IND_ID}.bwa.human.fixed.filtered.nodup.RG.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_quality.sh
	@echo "# === Filtering low quality reads mapped to human genome ====================== #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_quality.sh human;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.passed.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.RG.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_quality.sh
	@echo "# === Filtering low quality reads mapped to other genome ====================== #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_quality.sh ${SECOND_GENOME_NAME};
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam;


# Index the bam

# Run flagstat, idxstats, bedtools stats. reports/${IND_ID}.bwa.${GENOME_CODE}.aln_stats.pairsfix.fltr.nodups.highqual.txt

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Coverage calculations
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

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

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Explore factors affecting coverage
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

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

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- George ABySS methods (not itemized yet)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- SNP calling methods (not itemized yet)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #



