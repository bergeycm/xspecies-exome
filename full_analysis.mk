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
# --- demog_steps
index_snps : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup
call_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out
filter_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4
bed_from_bsnp : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed
# --- annotate_steps
convert_annovar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar
annovar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar.exonic_variant_function
# --- roh_steps
vcf_to_ped : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.ped
binary_ped : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.fam
plink_roh : results/${IND_ID}.ROH.hom
# --- pre_gphocs
intersect_indiv_beds : results/all.bsnp.snp.out.gt4.large.bed
make_gphocs_seq : results/all.combined.gphocs.seq
convert_to_nexus : results/all.combined.gphocs.nex
nj_tree : results/all.combined.gphocs.tre
# --- run_gphocs
run_gphocs : results/all.combined.gphocs.trace
# --- untr_pre_gphocs
get_untr_bed : results/all.combined.gphocs.seq.untranscribed.bed
make_gphocs_untr_seq : results/all.combined.gphocs.untranscribed.seq
convert_to_nexus_untr : results/all.combined.gphocs.untranscribed.nex
nj_tree_untr : results/all.combined.gphocs.untranscribed.tre
# --- run_gphocs_untr
run_gphocs_untr : results/all.combined.gphocs.untranscribed.trace
# --- filter_gphocs
calc_tajima_d : results/tajimas_d_full.txt
filter_1_5 : results/macaque_FULL.d1.5.seq
filter_2_0 : results/macaque_FULL.d2.seq
filter_2_5 : results/macaque_FULL.d2.5.seq
filter_3_0 : results/macaque_FULL.d3.seq
filter_5_0 : results/macaque_FULL.d5.seq
run_gphocs_d1_5 : results/macaque_FULL.d1.5.trace.log
run_gphocs_d2_0 : results/macaque_FULL.d2.trace.log
run_gphocs_d2_5 : results/macaque_FULL.d2.5.trace.log
run_gphocs_d3_0 : results/macaque_FULL.d3.trace.log
run_gphocs_d5_0 : results/macaque_FULL.d5.trace.log

# Group steps together

# Individual steps
preliminary_steps : index_genome merge_beds liftover_beds
pre_aln_analysis_steps : fastqc
alignment_steps : align sampe sam2bam sort_and_index_bam get_alignment_stats
post_alignment_filtering_steps : fix_mate_pairs filter_unmapped remove_dups add_read_groups filter_bad_qual
snp_calling_steps : local_realign_targets local_realign call_snps filter_snps get_snp_stats call_consensus
coverage_calc_steps : make_picard_intervals get_hsmetrics
pre_demog_steps : index_snps call_bsnp filter_bsnp bed_from_bsnp
annotate_steps : convert_annovar annovar
roh_steps : vcf_to_ped binary_ped plink_roh

# Comparative steps
pre_gphocs : intersect_indiv_beds make_gphocs_seq convert_to_nexus nj_tree
gphocs : run_gphocs 
untr_pre_gphocs : get_untr_bed make_gphocs_untr_seq convert_to_nexus_untr nj_tree_untr
gphocs_untr : run_gphocs_untr
filter_for_gphocs : calc_tajima_d filter_1_5 filter_2_0 filter_2_5 filter_3_0 filter_5_0 
gphocs_filtered : run_gphocs_d1_5 run_gphocs_d2_0 run_gphocs_d2_5 run_gphocs_d3_0 run_gphocs_d5_0

# Steps for individuals
indiv : preliminary_steps pre_aln_analysis_steps alignment_steps post_alignment_filtering_steps snp_calling_steps coverage_calc_steps pre_demog_steps annotate_steps roh_steps

# Steps for group
# Note! G-PhoCS on the full, untranscribed, and filtered datasets is not included in "group"
# This is only because they take four days to run.
# Command to include them would be:
compare_w_gphocs : pre_gphocs untr_pre_gphocs gphocs gphocs_untr filter_for_gphocs gphocs_filtered
compare : pre_gphocs untr_pre_gphocs filter_for_gphocs

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

# Picards HsMetrics output depends on realigned BAM, Picard-formatted BED files for targets and CCDS, Picard, and scripts/get_hsmetrics.sh
reports/${IND_ID}.bwa.human.hsmetrics.txt : results/${IND_ID}.bwa.human.passed.realn.bam results/${IND_ID}.bwa.human.passed.realn.bam.picard.baits.bed results/${IND_ID}.bwa.human.passed.realn.bam.picard.ccds.bed ${PICARD}/* #scripts/get_hsmetrics.sh
	@echo "# === Calculating coverage statistics with Picard HsMetrics for human genome == #";
	${SHELL_EXPORT} ./scripts/get_hsmetrics.sh results/${IND_ID}.bwa.human.passed.realn.bam human;
reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.hsmetrics.txt : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.baits.bed results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam.picard.ccds.bed ${PICARD}/* #scripts/get_hsmetrics.sh
	@echo "# === Calculating coverage statistics with Picard HsMetrics for other genome == #";
	${SHELL_EXPORT} ./scripts/get_hsmetrics.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_NAME};

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
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out #./scripts/filter_bsnp.sh
	@echo "# === Filtering BSNP-called SNPs with < 5 reads in 2nd genome only ============ #";
	touch results/standin.bsnp.gt4;
	${SHELL_EXPORT} ./scripts/filter_bsnp.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out > results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4;

# -------------------------------------------------------------------------------------- #
# --- Write BED of contigs covered by 5 or more reads (second genome only)
# -------------------------------------------------------------------------------------- #

# BED of regions to analyze depends on the filtered BSNP file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 #scripts/bed_from_bsnp.pl
	@echo "# === Writing BED of contigs covered 5+ reads in 2nd genome only ============ #";
	touch results/standin.bsnp.gt4.bed;
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
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup #${ANNOVAR}/convert2annovar.pl
	@echo "# === Converting SNPs to ANNOVAR format in 2nd genome only ==================== #";
	${SHELL_EXPORT} ${ANNOVAR}/convert2annovar.pl results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup > results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar;

# -------------------------------------------------------------------------------------- #
# --- Run ANNOVAR to annotate SNPs
# -------------------------------------------------------------------------------------- #

# ANNOVAR's exonic_variant_function output depends on ANNOVAR formatted file
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar.exonic_variant_function : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar #${ANNOVAR}/annotate_variation.pl
	@echo "# === Running ANNOVAR to annotate SNPs in 2nd genome only ===================== #";
	${ANNOVAR}/annotate_variation.pl --buildver ${ANNOVAR_BUILDVER} results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar ${ANNOVAR_DB_PATH};

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Look for Runs of Homozygosity with plink
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- First convert to plink's PED format
# -------------------------------------------------------------------------------------- #

# plink PED file depends on VCF formated file and VCFtools 
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.ped : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf #${VCFTOOLS}/*
	@echo "# === Converting to plink format in 2nd genome only =========================== #";
	${VCFTOOLS}/vcftools --vcf results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf --plink --out results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt;

# -------------------------------------------------------------------------------------- #
# --- Then convert the PED to a binary PED file, and make FAM files, etc.
# -------------------------------------------------------------------------------------- #

# plink FAM file depends on plink PED file and plink 
results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.fam : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.ped #${PLINK}/*
	@echo "# === Converting to binary plink format in 2nd genome only ==================== #";
	${PLINK}/plink --file results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt --make-bed --out results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt;

# -------------------------------------------------------------------------------------- #
# --- Then do ROH analysis in plink
# -------------------------------------------------------------------------------------- #

# Should be own script, since there are a ton of parameters?
# ROH output depends on plink FAM binary file and plink 
results/${IND_ID}.ROH.hom : results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.fam #${PLINK}/*
	@echo "# === Performing ROH analysis in 2nd genome only ============================== #";
	${PLINK}/plink --bfile results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt --homozyg-window-kb 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --homozyg-snp 5 --homozyg-kb 1 --allow-no-sex --out results/${IND_ID}.ROH

# ====================================================================================== #
# ====================================================================================== #
# ============================= Comparative analyses below ============================= #
# ====================================================================================== #
# ====================================================================================== #

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Intersect BEDs to get regions to include in G-PhoCS analysis (autosomal, >500 bp)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# BED of regions to run depends on individual BED files' stand-in, Perl script, and BEDtools
results/all.bsnp.snp.out.gt4.large.bed : results/standin.bsnp.gt4.bed #scripts/intersect_bsnp_beds.pl ${BEDTOOLS}/*
	@echo "# === Intersecting individual BEDs to get targets of G-PhoCS analysis ========= #";
	perl scripts/intersect_bsnp_beds.pl;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Make FASTA of the seqs for G-PhoCS and combine them into G-PhoCS sequence file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# G-PhoCS sequence file depends on big BED, BSNP files' stand-in, and Perl scripts
results/all.combined.gphocs.seq : results/all.bsnp.snp.out.gt4.large.bed results/standin.bsnp.gt4 #scripts/make_gphocs_seq_file.pl scripts/reduce_BSNP_via_BED scripts/bsnp_fastas_to_gphocs_seq_file
	@echo "# === Making individual FASTAs and then combining into G-PhoCS sequence file == #";
	perl scripts/make_gphocs_seq_file.pl;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Convert G-PhoCS sequence file into a NEXUS file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# NEXUS file depends on G-PhoCS sequence file and Perl script
results/all.combined.gphocs.nex : results/all.combined.gphocs.seq #scripts/gphocs_to_nexus.pl
	@echo "# === Converting G-PhoCS sequence file into a NEXUS file ====================== #";
	perl scripts/gphocs_to_nexus.pl results/all.combined.gphocs.seq > results/all.combined.gphocs.nex

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Infer NJ tree (results/all.combined.gphocs.tre)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# NJ tree depends on NEXUS file and PAUP
results/all.combined.gphocs.tre : results/all.combined.gphocs.nex #${PAUP}/*
	@echo "# === Inferring NJ tree ======================================================= #";
	${PAUP}/paup results/all.combined.gphocs.nex;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Call G-PhoCS after making sure seq filename, etc. is right in the control file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# G-PhoCS trace output depends on G-PhoCS seq file and control file and G-PhoCS
results/all.combined.gphocs.trace : results/all.combined.gphocs.seq ${GPHOCS_CTL_FULL} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on full dataset ========================================= #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_FULL}

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Get BED of untranscribed regions, results/all.combined.gphocs.seq.untranscribed.bed
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Note, this does A LOT more than just make a BED of the untranscribed stuff.
# It also makes BEDs of transcribed stuff and generates alignments of exons

# BED of untranscribed regions depends on G-PhoCS sequence file and refGene file
results/all.combined.gphocs.seq.untranscribed.bed : results/all.combined.gphocs.seq ${REFGENE};
	@echo "# === Making BED of untranscribed regions ===================================== #";
	perl scripts/compare_seqs_to_refgene.pl

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Make FASTA of untranscribed seqs and combine them into G-PhoCS sequence file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Untranscribed G-PhoCS sequence file depends on big (untr) BED, BSNP files' stand-in, and Perl scripts
results/all.combined.gphocs.untranscribed.seq : results/all.combined.gphocs.seq.untranscribed.bed results/standin.bsnp.gt4 #scripts/make_gphocs_seq_file.pl scripts/reduce_BSNP_via_BED scripts/bsnp_fastas_to_gphocs_seq_file
	@echo "# === Making indiv. untranscribed FASTAs and combining into G-PhoCS seq file == #";
	perl scripts/make_gphocs_seq_file.pl untranscribed;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Convert untranscribed G-PhoCS sequence file into a NEXUS file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# NEXUS file depends on G-PhoCS sequence file and Perl script
results/all.combined.gphocs.untranscribed.nex : results/all.combined.gphocs.untranscribed.seq #scripts/gphocs_to_nexus.pl
	@echo "# === Converting untranscribed G-PhoCS sequence file into a NEXUS file ======== #";
	perl scripts/gphocs_to_nexus.pl results/all.combined.gphocs.untranscribed.seq > results/all.combined.gphocs.untranscribed.nex

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Infer untranscribed regions NJ tree (results/all.combined.gphocs.untranscribed.tre)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# NJ tree depends on NEXUS file and PAUP
results/all.combined.gphocs.untranscribed.tre : results/all.combined.gphocs.untranscribed.nex #${PAUP}/*
	@echo "# === Inferring NJ tree for untranscribed regions ============================= #";
	${PAUP}/paup results/all.combined.gphocs.untranscribed.nex;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Call G-PhoCS on untranscribed after making sure all is right in the control file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Untr. G-PhoCS trace output depends on untr. G-PhoCS seq file and untr. control file and G-PhoCS
results/all.combined.gphocs.untranscribed.trace : results/all.combined.gphocs.untranscribed.seq ${GPHOCS_CTL_UNTR} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on untranscribed dataset ================================ #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_UNTR}

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Find Tajima's D for each G-PhoCS sequence
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# File of Tajima's D values depends on full G-PhoCS seq file, Perl script, and R script
results/tajimas_d_full.txt : results/all.combined.gphocs.seq # scripts/find_selected_loci.pl scripts/calc_tajima_d.R
	@echo "# === Calculating Tajima's D for all loci ===================================== #";
	perl scripts/find_selected_loci.pl results/all.combined.gphocs.seq > results/tajimas_d_full.txt

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Create filtered datasets by removing sequences with extreme Tajima's D values
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Filtered sequence file depends on file of Tajima's D values, full G-PhoCS seq file, and Perl script
# Filter: |D| < 1.5
results/macaque_FULL.d1.5.seq : results/all.combined.gphocs.seq results/tajimas_d_full.txt # scripts/remove_tajima_d_outlier_seqs.pl
	@echo "# === Filtering dataset |Tajima's D| < 1.5 ==================================== #";
	perl scripts/remove_tajima_d_outlier_seqs.pl 1.5 > results/macaque_FULL.d1.5.seq
	grep -c "chr" results/macaque_FULL.d1.5.seq > results/macaque_FULL.d1.5.seq.tmp
	echo >> results/macaque_FULL.d1.5.seq.tmp
	cat results/macaque_FULL.d1.5.seq >> results/macaque_FULL.d1.5.seq.tmp
	cp results/macaque_FULL.d1.5.seq.tmp results/macaque_FULL.d1.5.seq
	rm results/macaque_FULL.d1.5.seq.tmp
# Filter: |D| < 2
results/macaque_FULL.d2.seq : results/all.combined.gphocs.seq results/tajimas_d_full.txt # scripts/remove_tajima_d_outlier_seqs.pl
	@echo "# === Filtering dataset |Tajima's D| < 2.0 ==================================== #";
	perl scripts/remove_tajima_d_outlier_seqs.pl 2.0 > results/macaque_FULL.d2.seq
	grep -c "chr" results/macaque_FULL.d2.seq > results/macaque_FULL.d2.seq.tmp
	echo >> results/macaque_FULL.d2.seq.tmp
	cat results/macaque_FULL.d2.seq >> results/macaque_FULL.d2.seq.tmp
	cp results/macaque_FULL.d2.seq.tmp results/macaque_FULL.d2.seq
	rm results/macaque_FULL.d2.seq.tmp
# Filter: |D| < 2.5
results/macaque_FULL.d2.5.seq : results/all.combined.gphocs.seq results/tajimas_d_full.txt # scripts/remove_tajima_d_outlier_seqs.pl
	@echo "# === Filtering dataset |Tajima's D| < 2.5 ==================================== #";
	perl scripts/remove_tajima_d_outlier_seqs.pl 2.5 > results/macaque_FULL.d2.5.seq
	grep -c "chr" results/macaque_FULL.d2.5.seq > results/macaque_FULL.d2.5.seq.tmp
	echo >> results/macaque_FULL.d2.5.seq.tmp
	cat results/macaque_FULL.d2.5.seq >> results/macaque_FULL.d2.5.seq.tmp
	cp results/macaque_FULL.d2.5.seq.tmp results/macaque_FULL.d2.5.seq
	rm results/macaque_FULL.d2.5.seq.tmp
# Filter: |D| < 3
results/macaque_FULL.d3.seq : results/all.combined.gphocs.seq results/tajimas_d_full.txt # scripts/remove_tajima_d_outlier_seqs.pl
	@echo "# === Filtering dataset |Tajima's D| < 3.0 ==================================== #";
	perl scripts/remove_tajima_d_outlier_seqs.pl 3.0 > results/macaque_FULL.d3.seq
	grep -c "chr" results/macaque_FULL.d3.seq > results/macaque_FULL.d3.seq.tmp
	echo >> results/macaque_FULL.d3.seq.tmp
	cat results/macaque_FULL.d3.seq >> results/macaque_FULL.d3.seq.tmp
	cp results/macaque_FULL.d3.seq.tmp results/macaque_FULL.d3.seq
	rm results/macaque_FULL.d3.seq.tmp
# Filter: |D| < 5
results/macaque_FULL.d5.seq : results/all.combined.gphocs.seq results/tajimas_d_full.txt # scripts/remove_tajima_d_outlier_seqs.pl
	@echo "# === Filtering dataset |Tajima's D| < 5.0 ==================================== #";
	perl scripts/remove_tajima_d_outlier_seqs.pl 5.0 > results/macaque_FULL.d5.seq
	grep -c "chr" results/macaque_FULL.d5.seq > results/macaque_FULL.d5.seq.tmp
	echo >> results/macaque_FULL.d5.seq.tmp
	cat results/macaque_FULL.d5.seq >> results/macaque_FULL.d5.seq.tmp
	cp results/macaque_FULL.d5.seq.tmp results/macaque_FULL.d5.seq
	rm results/macaque_FULL.d5.seq.tmp

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Call G-PhoCS on these filtered datasets
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# G-PhoCS trace output depends on filtered G-PhoCS seq file and control file and G-PhoCS
results/macaque_FULL.d1.5.trace.log : results/macaque_FULL.d1.5.seq ${GPHOCS_CTL_1_5} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on dataset |Tajima's D| < 1.5 =========================== #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_1_5}
results/macaque_FULL.d2.trace.log : results/macaque_FULL.d2.seq ${GPHOCS_CTL_2_0} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on dataset |Tajima's D| < 2.0 =========================== #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_2_0}
results/macaque_FULL.d2.5.trace.log : results/macaque_FULL.d2.5.seq ${GPHOCS_CTL_2_5} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on dataset |Tajima's D| < 2.5 =========================== #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_2_5}
results/macaque_FULL.d3.trace.log : results/macaque_FULL.d3.seq ${GPHOCS_CTL_3_0} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on dataset |Tajima's D| < 3.0 =========================== #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_3_0}
results/macaque_FULL.d5.trace.log : results/macaque_FULL.d5.seq ${GPHOCS_CTL_5_0} #${GPHOCS}/*
	@echo "# === Calling G-PhoCS on dataset |Tajima's D| < 5.0 =========================== #";
	${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_5_0}