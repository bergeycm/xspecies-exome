#!/bin/sh

# Pipeline for the analysis of macaque exome data.
# These steps are contained in a Makefile, full_analysis.mk.
# The individual steps can be run by calling the shell script, full_analysis.
# The comparative steps can be run by calling the shell script, compare_analysis.
# User-specific parameters can be set with the configuration file, config.mk.

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Preliminary Steps
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Index genome
# -------------------------------------------------------------------------------------- #

echo "# === Indexing human genome =================================================== #";
./scripts/index_genome.sh ${HUMAN_GENOME_FA};
echo "# === Indexing secondary genome =============================================== #";
./scripts/index_genome.sh ${SECOND_GENOME_FA};
	
# -------------------------------------------------------------------------------------- #
# --- Merge overlapping intervals in BED files of targets and CCDS
# -------------------------------------------------------------------------------------- #

echo "# === Merging targets BED ===================================================== #";
./scripts/merge_bed.sh ${TARGETS};

echo "# === Merging CCDS BED ======================================================== #";
./scripts/merge_bed.sh ${CCDS};

# -------------------------------------------------------------------------------------- #
# --- LiftOver BED files of targets and CCDS into coordinates of the secondary genome
# -------------------------------------------------------------------------------------- #

echo "# === LiftingOver targets BED ================================================= #";
./scripts/liftOver_BED.sh ${TARGETS};
echo "# === LiftingOver CCDS BED ==================================================== #";
./scripts/liftOver_BED.sh ${CCDS};

./scripts/merge_bed.sh ${TARGETS}_2nd_liftover.bed;
./scripts/merge_bed.sh ${TARGETS}_hg_liftover.bed;

./scripts/merge_bed.sh ${CCDS}_2nd_liftover.bed;
./scripts/merge_bed.sh ${CCDS}_hg_liftover.bed;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Analyze reads
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Analyze reads with FastQC. Total sequence bp, Maximum possible sequence depth
# -------------------------------------------------------------------------------------- #

echo "# === Analyzing quality of reads (1st pair) =================================== #";
./scripts/run_fastqc.sh ${READ1} ${IND_ID}.read1.stats;
echo "# === Analyzing quality of reads (2nd pair) =================================== #";
./scripts/run_fastqc.sh ${READ2} ${IND_ID}.read2.stats;
	
# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mapping to reference genomes
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Align reads to genome with BWA
# -------------------------------------------------------------------------------------- #

echo "# === Aligning reads to human genome ========================================== #";
./scripts/align.sh ${HUMAN_GENOME_FA} human;
echo "# === Aligning reads to other genome ========================================== #";
./scripts/align.sh ${SECOND_GENOME_FA} ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Run sampe to generate SAM files
# -------------------------------------------------------------------------------------- #

echo "# === Combining reads to make SAM file for human genome ======================= #";
./scripts/sampe.sh ${HUMAN_GENOME_FA} human;
echo "# === Combining reads to make SAM file for other genome ======================= #";
./scripts/sampe.sh ${SECOND_GENOME_FA} ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Convert SAM file to BAM file
# -------------------------------------------------------------------------------------- #

echo "# === Converting SAM file to BAM file for human genome ======================== #";
./scripts/sam2bam.sh ${HUMAN_GENOME_FA}i human;
echo "# === Converting SAM file to BAM file for other genome ======================== #";
./scripts/sam2bam.sh ${SECOND_GENOME_FA}i ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Sort and index BAM
# -------------------------------------------------------------------------------------- #

echo "# === Sorting and Indexing BAM file for human genome ========================== #";
./scripts/sort_bam.sh results/${IND_ID}.bwa.human.sam.bam;
./scripts/index_bam.sh results/${IND_ID}.bwa.human.sam.bam.sorted.bam;
echo "# === Sorting and Indexing BAM file for other genome ========================== #";
./scripts/sort_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam;
./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam;

# -------------------------------------------------------------------------------------- #
# --- Analyze alignment output with flagstat, idxstats, and stats
# -------------------------------------------------------------------------------------- #

echo "# === Analyzing alignment output for human genome ============================= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.human.sam.bam.sorted.bam \
	reports/${IND_ID}.bwa.human.aln_stats.txt;
echo "# === Analyzing alignment output for other genome ============================= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.sam.bam.sorted.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.txt	

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Post-alignment filtering steps
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Fix mate pairs info
# -------------------------------------------------------------------------------------- #

echo "# === Fixing mate pair info for human genome alignment ======================== #";
./scripts/fix_mate_pairs.sh human;
echo "# === Fixing mate pair info for other genome alignment ======================== #";
./scripts/fix_mate_pairs.sh ${SECOND_GENOME_NAME};

echo "# === Analyzing alignment output for human genome (post mate pair fix) ======== #";
./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.fixed.bam \
	reports/${IND_ID}.bwa.human.aln_stats.pairsfix.txt;
echo "# === Analyzing alignment output for other genome (post mate pair fix) ======== #";
./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.txt;

# -------------------------------------------------------------------------------------- #
# --- Filtering for mapped and paired
# -------------------------------------------------------------------------------------- #

echo "# === Filtering unpaired reads mapped to human genome ========================= #";
./scripts/filter_mapped_reads_paired.sh human;
./scripts/index_bam.sh results/${IND_ID}.bwa.human.fixed.filtered.bam;
echo "# === Filtering unpaired reads mapped to other genome ========================= #";
./scripts/filter_mapped_reads_paired.sh ${SECOND_GENOME_NAME};
./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam;

echo "# === Analyzing alignment output for human genome (filtered for paired) ======= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.human.fixed.filtered.bam \
	reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.txt;
echo "# === Analyzing alignment output for other genome (filtered for paired) ======= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.txt;

# -------------------------------------------------------------------------------------- #
# --- Remove duplicates
# -------------------------------------------------------------------------------------- #

echo "# === Removing duplicate reads mapped to human genome ========================= #";
./scripts/remove_dups.sh human;
./scripts/index_bam.sh results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam;
echo "# === Removing duplicate reads mapped to other genome ========================= #";
./scripts/remove_dups.sh ${SECOND_GENOME_NAME};
./scripts/index_bam.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam;

echo "# === Analyzing alignment output for human genome (duplicates removed) ======== #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.human.fixed.filtered.nodup.bam \
	reports/${IND_ID}.bwa.human.aln_stats.pairsfix.flt.nodup.txt;
echo "# === Analyzing alignment output for other genome (duplicates removed) ======== #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.fixed.filtered.nodup.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.pairsfix.flt.nodup.txt;

# -------------------------------------------------------------------------------------- #
# --- Add read groups
# -------------------------------------------------------------------------------------- #

echo "# === Adding read groups for reads mapped to human genome ===================== #";
./scripts/add_read_groups.sh human;
echo "# === Adding read groups for reads mapped to other genome ===================== #";
./scripts/add_read_groups.sh ${SECOND_GENOME_NAME};

# -------------------------------------------------------------------------------------- #
# --- Remove reads with low mapping quality
# -------------------------------------------------------------------------------------- #

echo "# === Filtering low quality reads mapped to human genome ====================== #";
./scripts/filter_mapped_reads_quality.sh human;
./scripts/index_bam.sh results/${IND_ID}.bwa.human.passed.bam;
echo "# === Filtering low quality reads mapped to other genome ====================== #";
./scripts/filter_mapped_reads_quality.sh ${SECOND_GENOME_NAME};
./scripts/index_bam.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam;

echo "# === Analyzing alignment output for human genome (after qual filtering) ====== #";
./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.human.passed.bam \
	reports/${IND_ID}.bwa.human.aln_stats.passed.txt;
echo "# === Analyzing alignment output for other genome (after qual filtering) ====== #";
./scripts/get_alignment_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.txt;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- SNP calling methods
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Local realignment, step 1: ID realign targets
# -------------------------------------------------------------------------------------- #

echo "# === Identifying intervals in need or local realignment for human genome ===== #";
./scripts/local_realign_get_targets.sh human ${HUMAN_GENOME_FA};
echo "# === Identifying intervals in need or local realignment for other genome ===== #";
./scripts/local_realign_get_targets.sh ${SECOND_GENOME_NAME} ${SECOND_GENOME_FA};

# -------------------------------------------------------------------------------------- #
# --- Local realignment, step 2: realign around indels
# -------------------------------------------------------------------------------------- #

echo "# === Doing local realignment for human genome ================================ #";
./scripts/local_realign.sh human ${HUMAN_GENOME_FA};
echo "# === Doing local realignment for other genome ================================ #";
./scripts/local_realign.sh ${SECOND_GENOME_NAME} ${SECOND_GENOME_FA};

echo "# === Analyzing alignment output for human genome (locally realigned) ========= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.human.passed.realn.bam \
	reports/${IND_ID}.bwa.human.aln_stats.passed.realn.txt;
echo "# === Analyzing alignment output for other genome (locally realigned) ========= #";
./scripts/get_alignment_stats.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam \
	reports/${IND_ID}.bwa.${SECOND_GENOME_NAME}.aln_stats.passed.realn.txt;

# -------------------------------------------------------------------------------------- #
# --- Call SNPs
# -------------------------------------------------------------------------------------- #

echo "# === Calling raw SNPs relative to human genome =============================== #";
./scripts/call_snps.sh results/${IND_ID}.bwa.human.passed.realn.bam \
	${HUMAN_GENOME_FA};
echo "# === Calling raw SNPs relative to other genome =============================== #";
./scripts/call_snps.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam \
	${SECOND_GENOME_FA};
	
# -------------------------------------------------------------------------------------- #
# --- Filter SNPs for quality
# -------------------------------------------------------------------------------------- #

echo "# === Filtering raw SNPs relative to human genome ============================= #";
./scripts/filter_snps.sh results/${IND_ID}.bwa.human.passed.realn.raw.bcf;
echo "# === Filtering raw SNPs relative to other genome ============================= #";
./scripts/filter_snps.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.raw.bcf;

# -------------------------------------------------------------------------------------- #
# --- Get basic stats on SNPs
# -------------------------------------------------------------------------------------- #

echo "# === Getting basic SNPs stats for human genome =============================== #";
./scripts/get_snp_stats.sh results/${IND_ID}.bwa.human.passed.realn.flt.vcf;
echo "# === Getting basic SNPs stats for other genome =============================== #";
./scripts/get_snp_stats.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf;

# -------------------------------------------------------------------------------------- #
# --- Call consensus sequence
# -------------------------------------------------------------------------------------- #

echo "# === Calling consensus sequence relative to human genome ===================== #";
./scripts/call_consensus.sh \
	results/${IND_ID}.bwa.human.passed.realn.bam \
	${HUMAN_GENOME_FA} human ${TARGETS}_MERGED;
echo "# === Calling consensus sequence relative to other genome ===================== #";
./scripts/call_consensus.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam \
	${SECOND_GENOME_FA} ${SECOND_GENOME_NAME} ${TARGETS}_2nd_liftover.bed_MERGED;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Coverage calculations
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Make Picard-friendly intervals lists
# -------------------------------------------------------------------------------------- #

echo "# === Making Picard-friendly intervals files for human genome ================= #";
./scripts/make_picard_intervals.sh \
	results/${IND_ID}.bwa.human.passed.realn.bam ${TARGETS} ${CCDS};
echo "# === Making Picard-friendly intervals files for other genome ================= #";
./scripts/make_picard_intervals.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam \
	${TARGETS}_2nd_liftover.bed ${CCDS}_2nd_liftover.bed;

# -------------------------------------------------------------------------------------- #
# --- Calculate coverage of targeted regions
# -------------------------------------------------------------------------------------- #

echo "# === Calculating coverage statistics with Picard HsMetrics for human genome == #";
./scripts/get_hsmetrics.sh \
	results/${IND_ID}.bwa.human.passed.realn.bam human;
echo "# === Calculating coverage statistics with Picard HsMetrics for other genome == #";
./scripts/get_hsmetrics.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_NAME};

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Do individual processing in preparation to infer demographic history
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Index SNPs with (not m!) pileup
# -------------------------------------------------------------------------------------- #

echo "# === Indexing raw SNPs relative to human genome ============================== #";
./scripts/call_snps_pileup.sh \
	results/${IND_ID}.bwa.human.passed.realn.bam ${HUMAN_GENOME_FA};
echo "# === Indexing raw SNPs relative to other genome ============================== #";
./scripts/call_snps_pileup.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bam ${SECOND_GENOME_FA};

# -------------------------------------------------------------------------------------- #
# --- Call BSNP
# -------------------------------------------------------------------------------------- #

echo "# === Calling BSNP for SNPs relative to human genome ========================== #";
./scripts/call_bsnp.sh results/${IND_ID}.bwa.human.passed.realn.pileup;
echo "# === Calling BSNP for SNPs relative to other genome ========================== #";
./scripts/call_bsnp.sh results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup;

# -------------------------------------------------------------------------------------- #
# --- Filter BSNP-called SNPs with < 5 reads (second genome only)
# -------------------------------------------------------------------------------------- #

echo "# === Filtering BSNP-called SNPs with < 5 reads in 2nd genome only ============ #";
./scripts/filter_bsnp.sh \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out \
	> results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4;

# -------------------------------------------------------------------------------------- #
# --- Write BED of contigs covered by 5 or more reads (second genome only)
# -------------------------------------------------------------------------------------- #

echo "# === Writing BED of contigs covered 5+ reads in 2nd genome only ============ #";
perl scripts/bed_from_bsnp.pl \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4 \
	> results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.bsnp.snp.out.gt4.bed;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Annotate SNPs with ANNOVAR
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert SNPs to ANNOVAR format
# -------------------------------------------------------------------------------------- #

echo "# === Converting SNPs to ANNOVAR format in 2nd genome only ==================== #";
${ANNOVAR}/convert2annovar.pl \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.pileup \
	> results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar;

# -------------------------------------------------------------------------------------- #
# --- Run ANNOVAR to annotate SNPs
# -------------------------------------------------------------------------------------- #

echo "# === Running ANNOVAR to annotate SNPs in 2nd genome only ===================== #";
${ANNOVAR}/annotate_variation.pl --buildver ${ANNOVAR_BUILDVER} \
	results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.annovar ${ANNOVAR_DB_PATH};

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Look for Runs of Homozygosity with plink
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- First convert to plink's PED format
# -------------------------------------------------------------------------------------- #

echo "# === Converting to plink format in 2nd genome only =========================== #";
${VCFTOOLS}/vcftools \
	--vcf results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt.vcf --plink \
	--out results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt;

# -------------------------------------------------------------------------------------- #
# --- Then convert the PED to a binary PED file, and make FAM files, etc.
# -------------------------------------------------------------------------------------- #

echo "# === Converting to binary plink format in 2nd genome only ==================== #";
${PLINK}/plink --file results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt \
	--make-bed --out results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt;

# -------------------------------------------------------------------------------------- #
# --- Then do ROH analysis in plink
# -------------------------------------------------------------------------------------- #

# Should be own script, since there are a ton of parameters?
echo "# === Performing ROH analysis in 2nd genome only ============================== #";
${PLINK}/plink --bfile results/${IND_ID}.bwa.${SECOND_GENOME_NAME}.passed.realn.flt \
	--homozyg-window-kb 1000 --homozyg-window-snp 50 --homozyg-window-het 1 \
	--homozyg-window-missing 5 --homozyg-window-threshold 0.05 --homozyg-snp 5 \
	--homozyg-kb 1 --allow-no-sex --out results/${IND_ID}.ROH;

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

echo "# === Intersecting individual BEDs to get targets of G-PhoCS analysis ========= #";
perl scripts/intersect_bsnp_beds.pl;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Make FASTA of the seqs for G-PhoCS and combine them into G-PhoCS sequence file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Making individual FASTAs and then combining into G-PhoCS sequence file == #";
perl scripts/make_gphocs_seq_file.pl;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Convert G-PhoCS sequence file into a NEXUS file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Converting G-PhoCS sequence file into a NEXUS file ====================== #";
perl scripts/gphocs_to_nexus.pl results/all.combined.gphocs.seq \
	> results/all.combined.gphocs.nex

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Infer NJ tree (results/all.combined.gphocs.tre)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Inferring NJ tree ======================================================= #";
${PAUP}/paup results/all.combined.gphocs.nex;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Call G-PhoCS after making sure seq filename, etc. is right in the control file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Calling G-PhoCS on full dataset ========================================= #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_FULL}

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Get BED of untranscribed regions, results/all.combined.gphocs.seq.untranscribed.bed
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Note, this does A LOT more than just make a BED of the untranscribed stuff.
# It also makes BEDs of transcribed stuff and generates alignments of exons

echo "# === Making BED of untranscribed regions ===================================== #";
perl scripts/compare_seqs_to_refgene.pl

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Make FASTA of untranscribed seqs and combine them into G-PhoCS sequence file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Making indiv. untranscribed FASTAs and combining into G-PhoCS seq file == #";
perl scripts/make_gphocs_seq_file.pl untranscribed;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Convert untranscribed G-PhoCS sequence file into a NEXUS file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Converting untranscribed G-PhoCS sequence file into a NEXUS file ======== #";
perl scripts/gphocs_to_nexus.pl results/all.combined.gphocs.untranscribed.seq \
	> results/all.combined.gphocs.untranscribed.nex

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Infer untranscribed regions NJ tree (results/all.combined.gphocs.untranscribed.tre)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Inferring NJ tree for untranscribed regions ============================= #";
${PAUP}/paup results/all.combined.gphocs.untranscribed.nex;

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Call G-PhoCS on untranscribed after making sure all is right in the control file
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Calling G-PhoCS on untranscribed dataset ================================ #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_UNTR}

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Find Tajima's D for each G-PhoCS sequence
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Calculating Tajima's D for all loci ===================================== #";
perl scripts/find_selected_loci.pl results/all.combined.gphocs.seq > results/tajimas_d_full.txt

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Create filtered datasets by removing sequences with extreme Tajima's D values
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

echo "# === Filtering dataset |Tajima's D| < 1.5 ==================================== #";
perl scripts/remove_tajima_d_outlier_seqs.pl 1.5 > results/macaque_FULL.d1.5.seq
grep -c "chr" results/macaque_FULL.d1.5.seq > results/macaque_FULL.d1.5.seq.tmp
echo >> results/macaque_FULL.d1.5.seq.tmp
cat results/macaque_FULL.d1.5.seq >> results/macaque_FULL.d1.5.seq.tmp
cp results/macaque_FULL.d1.5.seq.tmp results/macaque_FULL.d1.5.seq
rm results/macaque_FULL.d1.5.seq.tmp

echo "# === Filtering dataset |Tajima's D| < 2.0 ==================================== #";
perl scripts/remove_tajima_d_outlier_seqs.pl 2.0 > results/macaque_FULL.d2.seq
grep -c "chr" results/macaque_FULL.d2.seq > results/macaque_FULL.d2.seq.tmp
echo >> results/macaque_FULL.d2.seq.tmp
cat results/macaque_FULL.d2.seq >> results/macaque_FULL.d2.seq.tmp
cp results/macaque_FULL.d2.seq.tmp results/macaque_FULL.d2.seq
rm results/macaque_FULL.d2.seq.tmp

echo "# === Filtering dataset |Tajima's D| < 2.5 ==================================== #";
perl scripts/remove_tajima_d_outlier_seqs.pl 2.5 > results/macaque_FULL.d2.5.seq
grep -c "chr" results/macaque_FULL.d2.5.seq > results/macaque_FULL.d2.5.seq.tmp
echo >> results/macaque_FULL.d2.5.seq.tmp
cat results/macaque_FULL.d2.5.seq >> results/macaque_FULL.d2.5.seq.tmp
cp results/macaque_FULL.d2.5.seq.tmp results/macaque_FULL.d2.5.seq
rm results/macaque_FULL.d2.5.seq.tmp

echo "# === Filtering dataset |Tajima's D| < 3.0 ==================================== #";
perl scripts/remove_tajima_d_outlier_seqs.pl 3.0 > results/macaque_FULL.d3.seq
grep -c "chr" results/macaque_FULL.d3.seq > results/macaque_FULL.d3.seq.tmp
echo >> results/macaque_FULL.d3.seq.tmp
cat results/macaque_FULL.d3.seq >> results/macaque_FULL.d3.seq.tmp
cp results/macaque_FULL.d3.seq.tmp results/macaque_FULL.d3.seq
rm results/macaque_FULL.d3.seq.tmp

echo "# === Filtering dataset |Tajima's D| < 5.0 ==================================== #";
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

echo "# === Calling G-PhoCS on dataset |Tajima's D| < 1.5 =========================== #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_1_5}

echo "# === Calling G-PhoCS on dataset |Tajima's D| < 2.0 =========================== #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_2_0}

echo "# === Calling G-PhoCS on dataset |Tajima's D| < 2.5 =========================== #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_2_5}

echo "# === Calling G-PhoCS on dataset |Tajima's D| < 3.0 =========================== #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_3_0}

echo "# === Calling G-PhoCS on dataset |Tajima's D| < 5.0 =========================== #";
${GPHOCS}/G-PhoCS-1-2-1 ${GPHOCS_CTL_5_0}
