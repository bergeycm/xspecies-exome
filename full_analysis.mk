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
# --- pre_aln_analysis_steps:
fastqc : reports/${IND_ID}.read1.stats.zip reports/${IND_ID}.read2.stats.zip
# --- alignment_steps
align : results/${IND_ID}.read1.bwa.human.sai results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai
sampe : results/${IND_ID}.bwa.human.sam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam
sam2bam : results/${IND_ID}.bwa.human.sam.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam
sort_and_index_bam : results/${IND_ID}.bwa.human.sam.bam.sorted.bam.bai results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam.bai
get_alignment_stats : reports/${IND_ID}.bwa.human.aln_stats.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt
# --- post_alignment_filtering_steps
fix_mate_pairs : results/${IND_ID}.bwa.human.fixed.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam reports/${IND_ID}.bwa.human.aln_stats.pairsfix.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.txt
filter_unmapped : results/${IND_ID}.bwa.human.fixed.filtered.bam.bai results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam.bai reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.txt
remove_dups : results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam.bai results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam.bai reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.nodup.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.nodup.txt
add_read_groups : results/${IND_ID}.bwa.human.fixed.filtered.nodup.RG.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.RG.bam
filter_bad_qual : results/${IND_ID}.bwa.human.passed.bam.bai results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam.bai reports/${IND_ID}.bwa.human.aln_stats.passed.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.txt
# --- abyss_steps

# --- snp_calling_steps
local_realign_targets : results/${IND_ID}.bwa.human.passed.bam.list results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam.list
local_realign : results/${IND_ID}.bwa.human.passed.realn.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam reports/${IND_ID}.bwa.human.aln_stats.passed.realn.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.realn.txt
call_snps : results/${IND_ID}.bwa.human.passed.realn.raw.bcf results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.raw.bcf
filter_snps : results/${IND_ID}.bwa.human.passed.realn.flt.vcf results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf
get_snp_stats : reports/${IND_ID}.bwa.human.passed.realn.flt.vcf.stats.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf.stats.txt
call_consensus : results/${IND_ID}.bwa.human.consensus.fq.gz results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.consensus.fq.gz
# --- coverage_calc_steps
make_picard_intervals : results/${IND_ID}.bwa.human.passed.realn.bam.picard.baits.bed results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.baits.bed
get_hsmetrics : reports/${IND_ID}.bwa.human.hsmetrics.txt reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.hsmetrics.txt
# --- psmc_steps
#fastq_to_psmcfa : results/${IND_ID}.bwa.human.diploid.psmcfa results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmcfa
#psmc : results/${IND_ID}.bwa.human.diploid.psmc results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmc
#psmc_ms_plot : results/${IND_ID}.bwa.human.diploid.plot results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.plot
# --- demog_steps
index_snps : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup
call_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out
filter_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4
bed_from_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed
# --- annotate_steps
convert_annovar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annonvar
annovar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annonvar.exonic_variant_function

# Group steps together
preliminary_steps : index_genome merge_beds liftover_beds
pre_aln_analysis_steps : fastqc
alignment_steps : align sampe sam2bam sort_and_index_bam get_alignment_stats
post_alignment_filtering_steps : fix_mate_pairs filter_unmapped remove_dups add_read_groups filter_bad_qual
snp_calling_steps : local_realign_targets local_realign call_snps filter_snps get_snp_stats call_consensus
coverage_calc_steps : make_picard_intervals get_hsmetrics
#psmc_steps : fastq_to_psmcfa psmc psmc_ms_plot
pre_demog_steps : index_snps call_bsnp filter_bsnp bed_from_bsnp
annotate_steps : convert_annovar annovar

all : preliminary_steps pre_aln_analysis_steps alignment_steps post_alignment_filtering_steps snp_calling_steps coverage_calc_steps pre_demog_steps annotate_steps

# Hack to be able to export Make variables to child scripts
# Don't export variables from make that begin with non-alphanumeric character
# After that, underscores are OK
# Also get rid of newlines.
#MAKE_ENV := $(shell echo '$(.VARIABLES)' | awk -v RS=' ' '/^[a-zA-Z0-9][a-zA-Z0-9_]+$$/')
#SHELL_EXPORT := $(foreach v,$(MAKE_ENV),$(v)='$($(v))')
# Also get rid of newlines and module=() {  eval `/opt/Modules/bin/modulecmd bash `}
#SHELL_EXPORT := $(shell echo ${SHELL_EXPORT} | tr '\n' ' ' | sed -e 's/module=.*\}//g' | sed -e 's/rm -f/"rm -f"/g')
SHELL_EXPORT := 

# Export Make variables to child scripts
.EXPORT_ALL_VARIABLES :

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
# --- Analyze reads
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Analyze reads with FastQC. Total sequence bp, Maximum possible sequence depth
# -------------------------------------------------------------------------------------- #

# FastQC reports depend on read files, FastQC, and run_fastqc.sh
reports/${IND_ID}.read1.stats.zip : ${READ1} ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (1st pair) before filtering ================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ1} ${IND_ID}.read1.stats;
reports/${IND_ID}.read2.stats.zip : ${READ2} ${FASTQC}/* #scripts/run_fastqc.sh
	@echo "# === Analyzing quality of reads (2nd pair) before filtering ================== #";
	${SHELL_EXPORT} ./scripts/run_fastqc.sh ${READ2} ${IND_ID}.read2.stats;
	
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
results/${IND_ID}.read1.bwa.human.sai : ${BWA}/* ${READ1} ${READ2} ${HUMAN_GENOME_FA}i #scripts/align.sh
	@echo "# === Aligning reads to human genome ========================================== #";
	${SHELL_EXPORT} ./scripts/align.sh ${HUMAN_GENOME_FA} human;
results/${IND_ID}.read1.bwa.${SECOND_GENOME_NAME}.sai : ${BWA}/* ${READ1} ${READ2} ${SECOND_GENOME_FA}i #scripts/align.sh
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

# Sorted BAM file index depends on unsorted BAM file, scripts/sort_bam, and scripts/index_bam.sh
results/${IND_ID}.bwa.human.sam.bam.sorted.bam.bai : results/${IND_ID}.bwa.human.sam.bam #scripts/sort_bam scripts/index_bam.sh
	@echo "# === Sorting and Indexing BAM file for human genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_bam.sh results/${IND_ID}.bwa.human.sam.bam;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.sam.bam.sorted.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam.bai : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam #scripts/sort_bam scripts/index_bam.sh
	@echo "# === Sorting and Indexing BAM file for other genome ========================== #";
	${SHELL_EXPORT} ./scripts/sort_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam;

# -------------------------------------------------------------------------------------- #
# --- Analyze alignment output with flagstat, idxstats, and stats
# -------------------------------------------------------------------------------------- #

# Align stats report depends on the sorted BAM and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.txt : results/${IND_ID}.bwa.human.sam.bam.sorted.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome ============================= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.sam.bam.sorted.bam reports/${IND_ID}.bwa.human.aln_stats.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome ============================= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt	

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

# Align stats report depends on the BAM with fixed mate pair info and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.pairsfix.txt : results/${IND_ID}.bwa.human.fixed.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome (post mate pair fix) ======== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.fixed.bam reports/${IND_ID}.bwa.human.aln_stats.pairsfix.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome (post mate pair fix) ======== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.txt;

# -------------------------------------------------------------------------------------- #
# --- Filtering for mapped and paired
# -------------------------------------------------------------------------------------- #

# Filtered BAM [index file] depends on output BAM from fix_mate_pairs.sh, BAMtools, and scripts/filter_mapped_reads_paired.sh
results/${IND_ID}.bwa.human.fixed.filtered.bam.bai : results/${IND_ID}.bwa.human.fixed.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_paired.sh
	@echo "# === Filtering unpaired reads mapped to human genome ========================= #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_paired.sh human;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.fixed.filtered.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam.bai : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_paired.sh
	@echo "# === Filtering unpaired reads mapped to other genome ========================= #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_paired.sh ${SECOND_GENOME_NAME};
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam;

# Align stats report depends on filtered BAM and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.txt : results/${IND_ID}.bwa.human.fixed.filtered.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome (filtered for paired) ======= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.fixed.filtered.bam reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome (filtered for paired) ======= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.txt;

# -------------------------------------------------------------------------------------- #
# --- Remove duplicates
# -------------------------------------------------------------------------------------- #

# BAM sans dups [index file] depends on output BAM from filter_mapped_reads_paired.sh, Picard, and scripts/remove_dups.sh
results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam.bai : results/${IND_ID}.bwa.human.fixed.filtered.bam ${PICARD}/* # scripts/remove_dups.sh
	@echo "# === Removing duplicate reads mapped to human genome ========================= #";
	${SHELL_EXPORT} ./scripts/remove_dups.sh human;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam.bai : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam ${PICARD}/* # scripts/remove_dups.sh
	@echo "# === Removing duplicate reads mapped to other genome ========================= #";
	${SHELL_EXPORT} ./scripts/remove_dups.sh ${SECOND_GENOME_NAME};
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam;

# Align stats report depends on BAM sans dups and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.nodup.txt : results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome (duplicates removed) ======== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.nodup.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.nodup.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome (duplicates removed) ======== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.nodup.txt;

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
results/${IND_ID}.bwa.human.passed.bam.bai : results/${IND_ID}.bwa.human.fixed.filtered.nodup.RG.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_quality.sh
	@echo "# === Filtering low quality reads mapped to human genome ====================== #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_quality.sh human;
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.human.passed.bam;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam.bai : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.RG.bam ${BEDTOOLS}/* # scripts/filter_mapped_reads_quality.sh
	@echo "# === Filtering low quality reads mapped to other genome ====================== #";
	${SHELL_EXPORT} ./scripts/filter_mapped_reads_quality.sh ${SECOND_GENOME_NAME};
	${SHELL_EXPORT} ./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam;

# Align stats report depends on quality-filtered BAM and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.passed.txt : results/${IND_ID}.bwa.human.passed.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome (after qual filtering) ====== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.passed.bam reports/${IND_ID}.bwa.human.aln_stats.passed.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome (after qual filtering) ====== #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.txt;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- George ABySS methods (not itemized yet)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #


# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- SNP calling methods
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Local realignment, step 1: ID realign targets
# -------------------------------------------------------------------------------------- #

# List of intervals to realign depends on BAM of reads that passed filtering, GATK, and scripts/local_realign_get_targets.sh
results/${IND_ID}.bwa.human.passed.bam.list : results/${IND_ID}.bwa.human.passed.bam ${GATK}/* #scripts/local_realign_get_targets.sh
	@echo "# === Identifying intervals in need or local realignment for human genome ===== #";
	${SHELL_EXPORT} ./scripts/local_realign_get_targets.sh human ${HUMAN_GENOME_FA};
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam.list : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam ${GATK}/* #scripts/local_realign.sh
	@echo "# === Identifying intervals in need or local realignment for other genome ===== #";
	${SHELL_EXPORT} ./scripts/local_realign_get_targets.sh ${SECOND_GENOME_NAME} ${SECOND_GENOME_FA};

# -------------------------------------------------------------------------------------- #
# --- Local realignment, step 2: realign around indels
# -------------------------------------------------------------------------------------- #

# Realigned BAM depends on list of realign targets, BAM of reads that passed filtering, GATK, and scripts/local_realign.sh
results/${IND_ID}.bwa.human.passed.realn.bam : results/${IND_ID}.bwa.human.passed.bam.list results/${IND_ID}.bwa.human.passed.bam ${GATK}/* #scripts/local_realign.sh
	@echo "# === Doing local realignment for human genome ================================ #";
	${SHELL_EXPORT} ./scripts/local_realign.sh human ${HUMAN_GENOME_FA};
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam.list results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam ${GATK}/* #scripts/local_realign.sh
	@echo "# === Doing local realignment for other genome ================================ #";
	${SHELL_EXPORT} ./scripts/local_realign.sh ${SECOND_GENOME_NAME} ${SECOND_GENOME_FA};

# Align stats report depends on realigned BAM and scripts/get_alignment_stats.sh
reports/${IND_ID}.bwa.human.aln_stats.passed.realn.txt : results/${IND_ID}.bwa.human.passed.realn.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for human genome (locally realigned) ========= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.passed.realn.bam reports/${IND_ID}.bwa.human.aln_stats.passed.realn.txt;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.realn.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam #scripts/get_alignment_stats.sh
	@echo "# === Analyzing alignment output for other genome (locally realigned) ========= #";
	${SHELL_EXPORT} ./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.realn.txt;

# -------------------------------------------------------------------------------------- #
# --- Call SNPs
# -------------------------------------------------------------------------------------- #

# Raw SNPs file depends on realigned BAM, VCFtools, BCFtools, and scripts/call_snps.sh
results/${IND_ID}.bwa.human.passed.realn.raw.bcf : results/${IND_ID}.bwa.human.passed.realn.bam #${VCFTOOLS}/* ${BCFTOOLS}/* #scripts/call_snps.sh
	@echo "# === Calling raw SNPs relative to human genome =============================== #";
	${SHELL_EXPORT} ./scripts/call_snps.sh results/${IND_ID}.bwa.human.passed.realn.bam ${HUMAN_GENOME_FA};
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.raw.bcf : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam #${VCFTOOLS}/* ${BCFTOOLS}/* #scripts/call_snps.sh
	@echo "# === Calling raw SNPs relative to other genome =============================== #";
	${SHELL_EXPORT} ./scripts/call_snps.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_FA};
	
# -------------------------------------------------------------------------------------- #
# --- Filter SNPs for quality
# -------------------------------------------------------------------------------------- #

# Filtered SNP file depends on raw SNP file, BCFtools, and scripts/filter_snps.sh
results/${IND_ID}.bwa.human.passed.realn.flt.vcf : results/${IND_ID}.bwa.human.passed.realn.raw.bcf #${BCFTOOLS}/* #scripts/filter_snps.sh
	@echo "# === Filtering raw SNPs relative to human genome ============================= #";
	${SHELL_EXPORT} ./scripts/filter_snps.sh results/${IND_ID}.bwa.human.passed.realn.raw.bcf;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.raw.bcf #${BCFTOOLS}/* #scripts/filter_snps.sh
	@echo "# === Filtering raw SNPs relative to other genome ============================= #";
	${SHELL_EXPORT} ./scripts/filter_snps.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.raw.bcf;

# -------------------------------------------------------------------------------------- #
# --- Get basic stats on SNPs
# -------------------------------------------------------------------------------------- #

# File of SNP stats depends on VCF file, VCFtools, and scripts/get_snp_stats.sh
reports/${IND_ID}.bwa.human.passed.realn.flt.vcf.stats.txt : results/${IND_ID}.bwa.human.passed.realn.flt.vcf #${VCFTOOLS}/* #scripts/get_snp_stats.sh
	@echo "# === Getting basic SNPs stats for human genome =============================== #";
	${SHELL_EXPORT} ./scripts/get_snp_stats.sh results/${IND_ID}.bwa.human.passed.realn.flt.vcf;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf.stats.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf #${VCFTOOLS}/* #scripts/get_snp_stats.sh
	@echo "# === Getting basic SNPs stats for other genome =============================== #";
	${SHELL_EXPORT} ./scripts/get_snp_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf;

# -------------------------------------------------------------------------------------- #
# --- Call consensus sequence
# -------------------------------------------------------------------------------------- #

# Consensus sequence depends on realigned BAM, SAMtools, BCFtools, and scripts/call_consensus.sh
results/${IND_ID}.bwa.human.consensus.fq.gz : results/${IND_ID}.bwa.human.passed.realn.bam ${SAMTOOLS}/* ${BCFTOOLS}/* #scripts/call_consensus.sh
	@echo "# === Calling consensus sequence relative to human genome ===================== #";
	${SHELL_EXPORT} ./scripts/call_consensus.sh results/${IND_ID}.bwa.human.passed.realn.bam ${HUMAN_GENOME_FA} human ${TARGETS}_MERGED;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.consensus.fq.gz : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SAMTOOLS}/* ${BCFTOOLS}/* #scripts/call_consensus.sh
	@echo "# === Calling consensus sequence relative to other genome ===================== #";
	${SHELL_EXPORT} ./scripts/call_consensus.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_FA} ${SECOND_GENOME_NAME} ${TARGETS}_2nd_liftover.bed_MERGED;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Coverage calculations
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Make Picard-friendly intervals lists
# -------------------------------------------------------------------------------------- #

# Picard intervals lists depends on realigned BAM, target BED file, CCDS BED file, SAMtools and scripts/make_picard_intervals.sh
results/${IND_ID}.bwa.human.passed.realn.bam.picard.baits.bed : results/${IND_ID}.bwa.human.passed.realn.bam ${TARGETS} ${CCDS} ${SAMTOOLS}/* #scripts/make_picard_intervals.sh
	@echo "# === Making Picard-friendly intervals files for human genome ================= #";
	${SHELL_EXPORT} ./scripts/make_picard_intervals.sh results/${IND_ID}.bwa.human.passed.realn.bam ${TARGETS} ${CCDS};
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.baits.bed : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed ${SAMTOOLS}/* #scripts/make_picard_intervals.sh
	@echo "# === Making Picard-friendly intervals files for other genome ================= #";
	${SHELL_EXPORT} ./scripts/make_picard_intervals.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed;

# -------------------------------------------------------------------------------------- #
# --- Calculate coverage of targeted regions
# -------------------------------------------------------------------------------------- #

#java -Xmx10g -jar ${PICARD}/CalculateHsMetrics.jar \
#	BAIT_INTERVALS=${IN_BAM}.picard.baits.bed \
#	TARGET_INTERVALS=${IN_BAM}.picard.ccds.bed \
#	INPUT=${IN_BAM} \
#	OUTPUT=reports/${IND_ID}.bwa.${GENOME_NAME}.hsmetrics.txt

# Picards HsMetrics output depends on realigned BAM, Picard-formatted BED files for targets and CCDS, Picard, and scripts/get_hsmetrics.sh
reports/${IND_ID}.bwa.human.hsmetrics.txt : results/${IND_ID}.bwa.human.passed.realn.bam results/${IND_ID}.bwa.human.passed.realn.bam.picard.baits.bed results/${IND_ID}.bwa.human.passed.realn.bam.picard.ccds.bed ${PICARD}/* #scripts/get_hsmetrics.sh
	@echo "# === Calculating coverage statistics with Picard HsMetrics for human genome == #";
	${SHELL_EXPORT} ./scripts/get_hsmetrics.sh results/${IND_ID}.bwa.human.passed.realn.bam human;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.hsmetrics.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.baits.bed results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.ccds.bed ${PICARD}/* #scripts/get_hsmetrics.sh
	@echo "# === Calculating coverage statistics with Picard HsMetrics for other genome == #";
	${SHELL_EXPORT} ./scripts/get_hsmetrics.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_NAME};


# Old way:
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

	###	# ====================================================================================== #
	###	# -------------------------------------------------------------------------------------- #
	###	# --- Run PSMC to infer demographic history
	###	# -------------------------------------------------------------------------------------- #
	###	# ====================================================================================== #
	###	
	###	# -------------------------------------------------------------------------------------- #
	###	# --- Convert FASTQ to PSMCFA
	###	# -------------------------------------------------------------------------------------- #
	###	
	###	# PSMCFA file depends on consensus FASTQ, PSMC, and scripts/convert_to_psmcfa.sh
	###	results/${IND_ID}.bwa.human.diploid.psmcfa : results/${IND_ID}.bwa.human.consensus.fq.gz ${PSMC}/* #scripts/convert_to_psmcfa.sh
	###		@echo "# === Converting FASTQ to PSMCFA for human genome ============================= #";
	###		${SHELL_EXPORT} ./scripts/convert_to_psmcfa.sh results/${IND_ID}.bwa.human.consensus.fq.gz;
	###	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmcfa : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.consensus.fq.gz ${PSMC}/* #scripts/convert_to_psmcfa.sh
	###		@echo "# === Converting FASTQ to PSMCFA for other genome ============================= #";
	###		${SHELL_EXPORT} ./scripts/convert_to_psmcfa.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.consensus.fq.gz;
	###	
	###	# -------------------------------------------------------------------------------------- #
	###	# --- Run PSMC
	###	# -------------------------------------------------------------------------------------- #
	###	
	###	# PSMC output file depends on PSMCFA file, PSMC, and scripts/run_psmc.sh
	###	results/${IND_ID}.bwa.human.diploid.psmc : results/${IND_ID}.bwa.human.diploid.psmcfa ${PSMC}/* #scripts/run_psmc.sh
	###		@echo "# === Running PSMC for consensus from human genome ============================ #";
	###		${SHELL_EXPORT} ./scripts/run_psmc.sh results/${IND_ID}.bwa.human.diploid.psmcfa;
	###	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmc : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmcfa ${PSMC}/* #scripts/run_psmc.sh
	###		@echo "# === Running PSMC for consensus from other genome ============================ #";
	###		${SHELL_EXPORT} ./scripts/run_psmc.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmcfa;
	###	
	###	# -------------------------------------------------------------------------------------- #
	###	# --- Call psmc2history.pl and Plot PSMC results
	###	# -------------------------------------------------------------------------------------- #
	###	
	###	# PSMC plot file depends on PSMC file, PSMC, and scripts/psmc_to_ms_and_plot.sh
	###	results/${IND_ID}.bwa.human.diploid.plot : results/${IND_ID}.bwa.human.diploid.psmc ${PSMC}/* #scripts/psmc_to_ms_and_plot.sh
	###		@echo "# === Generating ms command and plot from PSMC for human genome =============== #";
	###		${SHELL_EXPORT} ./scripts/psmc_to_ms_and_plot.sh results/${IND_ID}.bwa.human.diploid.psmc;
	###	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.plot : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmc ${PSMC}/* #scripts/psmc_to_ms_and_plot.sh
	###		@echo "# === Generating ms command and plot from PSMC for other genome =============== #";
	###		${SHELL_EXPORT} ./scripts/psmc_to_ms_and_plot.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.diploid.psmc;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Do individual processing in preparation to infer demographic history
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Index SNPs with (not m!) pileup
# -------------------------------------------------------------------------------------- #

# SNPs pileup file depends on realigned BAM, SAMtools, and scripts/call_snps_pileup.sh
results/${IND_ID}.bwa.human.passed.realn.pileup : results/${IND_ID}.bwa.human.passed.realn.bam #${SAMTOOLS}/* #scripts/call_snps_pileup.sh
	@echo "# === Indexing raw SNPs relative to human genome ============================== #";
	${SHELL_EXPORT} ./scripts/call_snps_pileup.sh results/${IND_ID}.bwa.human.passed.realn.bam ${HUMAN_GENOME_FA};
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam #${SAMTOOLS}/* #call_snps_pileup/call_snps.sh
	@echo "# === Indexing raw SNPs relative to other genome ============================== #";
	${SHELL_EXPORT} ./scripts/call_snps_pileup.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_FA};

# -------------------------------------------------------------------------------------- #
# --- Call BSNP
# -------------------------------------------------------------------------------------- #

# BSNP main output file depends on pileup output, BSNP, and scripts/call_bsnp.sh
results/${IND_ID}.bwa.human.passed.realn.bsnp.snp.out : results/${IND_ID}.bwa.human.passed.realn.pileup #${BSNP}/* #scripts/call_bsnp.sh
	@echo "# === Calling BSNP for SNPs relative to human genome ========================== #";
	${SHELL_EXPORT} ./scripts/call_bsnp.sh results/${IND_ID}.bwa.human.passed.realn.pileup;
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup #${BSNP}/* #scripts/call_bsnp.sh
	@echo "# === Calling BSNP for SNPs relative to other genome ========================== #";
	${SHELL_EXPORT} ./scripts/call_bsnp.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup;

# -------------------------------------------------------------------------------------- #
# --- Filter BSNP-called SNPs with < 5 reads (second genome only)
# -------------------------------------------------------------------------------------- #

# Filtered BSNP file depends on original BSNP file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out ./scripts/filter_bsnp.sh
	@echo "# === Filtering BSNP-called SNPs with < 5 reads in 2nd genome only ============ #";
	${SHELL_EXPORT} ./scripts/filter_bsnp.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out > results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4;

# -------------------------------------------------------------------------------------- #
# --- Write BED of contigs covered by 5 or more reads (second genome only)
# -------------------------------------------------------------------------------------- #

# BED of regions to analyze depends on the filtered BSNP file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 #scripts/bed_from_bsnp.pl
	@echo "# === Writing BED of contigs covered 5+ reads in 2nd genome only ============ #";
	${SHELL_EXPORT} perl scripts/bed_from_bsnp.pl results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 > results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Annotate SNPs with ANNOVAR
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert SNPs to ANNOVAR format
# -------------------------------------------------------------------------------------- #

# ANNOVAR formatted file depends on pileup formatted file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annonvar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup #${ANNOVAR}/convert2annovar.pl
	@echo "# === Converting SNPs to ANNOVAR format in 2nd genome only ============ #";
	${SHELL_EXPORT} ${ANNOVAR}/convert2annovar.pl results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup;

# -------------------------------------------------------------------------------------- #
# --- Run ANNOVAR to annotate SNPs
# -------------------------------------------------------------------------------------- #

# ANNOVAR's exonic_variant_function output depends on ANNOVAR formatted file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annonvar.exonic_variant_function : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annonvar #${ANNOVAR}/annotate_variation.pl
	@echo "# === Running ANNOVAR to annotate SNPs in 2nd genome only ============ #";
	${ANNOVAR}/annotate_variation.pl --buildver ${ANNOVAR_BUILDVER} results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar ${ANNOVAR_DB_PATH};

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Look for Runs of Homozygosity with plink
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# First convert to plink's PED format
#~/exome_macaque/bin/vcftools_0.1.9/bin/vcftools --vcf results/george.bwa.rhesus.passed.realn.flt.vcf      --plink --out results/george.bwa.rhesus.passed.realn.flt

# Then convert the PED to a binary PED file, and make FAM files, etc.
# ~/exome_macaque/bin/plink-1.07-x86_64/plink --file results/george.bwa.rhesus.passed.realn.flt      --make-bed --out results/george.bwa.rhesus.passed.realn.flt

# Then do ROH analysis in plink
# Should be own script?
#~/exome_macaque/bin/plink-1.07-x86_64/plink --bfile results/vallender.bwa.rhesus.passed.realn.flt   --homozyg-window-kb 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --homozyg-snp 5 --homozyg-kb 1 --allow-no-sex --out results/vallender.ROH