# -------------------------------------------------------------------------------------- #
# --- Makefile to generate and analyze fake human exomes
# --- Called by the executable shell script, analyze_fake_human_exomes
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

HUMAN_GENOME_DIR=$(dir ${HUMAN_GENOME_FA})

# Steps. Can be called one-by-one with something like, make liftover_bed
# --- prelim_grab_genome:
liftover_bed_full : data/gphocs_regions_human_full.bed
liftover_bed_untr : data/gphocs_regions_human_untr.bed
clean_bed_full : data/gphocs_regions_human_full.noRandom.bed
clean_bed_untr : data/gphocs_regions_human_untr.noRandom.bed
merge_bed_full : data/gphocs_regions_human_full.noRandom.merged.bed
merge_bed_untr : data/gphocs_regions_human_untr.noRandom.merged.bed
grab_genome_seq_full : data/hg18_gphocs_targets.full.fa
grab_genome_seq_untr : data/hg18_gphocs_targets.untr.fa
# --- prelim_grab_chimp:
liftover_to_chimp_full : data/gphocs_regions_human_full.noRandom.merged.panTro2.bed
liftover_to_chimp_untr : data/gphocs_regions_human_untr.noRandom.merged.panTro2.bed
clean_chimp_bed_full : data/gphocs_regions_human_full.noRandom.merged.panTro2.noRandom.bed
clean_chimp_bed_untr : data/gphocs_regions_human_untr.noRandom.merged.panTro2.noRandom.bed
grab_chimp_seq_full : data/panTro2_gphocs_targets.full.fa
grab_chimp_seq_untr : data/panTro2_gphocs_targets.untr.fa
# --- generate_exomes
make_exomes_full : fake_human_exomes_full.seq 
make_exomes_untr : fake_human_exomes_untr.seq

# Group steps together

prelim_grab_genome : liftover_bed_full liftover_bed_untr clean_bed_full clean_bed_untr merge_bed_full merge_bed_untr grab_genome_seq_full grab_genome_seq_untr 
prelim_grab_chimp : liftover_to_chimp_full liftover_to_chimp_untr clean_chimp_bed_full clean_chimp_bed_untr grab_chimp_seq_full grab_chimp_seq_untr
generate_exomes : make_exomes_full make_exomes_untr

all : prelim_grab_genome prelim_grab_chimp generate_exomes

SHELL_EXPORT := 

# Export Make variables to child scripts
.EXPORT_ALL_VARIABLES :

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Preliminary steps to grab human genome sequences
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert coordinates from rhesus exomes to human coordinates - Full dataset
# -------------------------------------------------------------------------------------- #

# Liftover'd human BED file depends on LiftOver, macaque BED file, and chain file
data/gphocs_regions_human_full.bed : ${LIFTOVER}/liftOver data/all.bsnp.snp.out.gt4.large.bed
	@echo "# === Converting rhesus coordinates to human - Full dataset =================== #";
	${LIFTOVER}/liftOver data/all.bsnp.snp.out.gt4.large.bed ${TO_HUMAN_CHAINFILE} data/gphocs_regions_human_full.bed data/gphocs_regions_human_full.unMapped.bed;

# -------------------------------------------------------------------------------------- #
# --- Convert coordinates from rhesus exomes to human coordinates - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Liftover'd human BED file depends on LiftOver, macaque BED file, and chain file
data/gphocs_regions_human_untr.bed : ${LIFTOVER}/liftOver data/all.combined.gphocs.seq.untranscribed.bed
	@echo "# === Converting rhesus coordinates to human - Untranscribed dataset ========== #";
	${LIFTOVER}/liftOver data/all.combined.gphocs.seq.untranscribed.bed ${TO_HUMAN_CHAINFILE} data/gphocs_regions_human_untr.bed data/gphocs_regions_human_untr.unMapped.bed;

# -------------------------------------------------------------------------------------- #
# --- Remove *_random and *_hap* entries in BED - Full dataset
# -------------------------------------------------------------------------------------- #

# Cleaned BED file depends on original BED file
data/gphocs_regions_human_full.noRandom.bed : data/gphocs_regions_human_full.bed
	@echo "# === Cleaning BED file - Full dataset ======================================== #";
	grep -v "_random" data/gphocs_regions_human_full.bed | grep -v "_hap" > data/gphocs_regions_human_full.noRandom.bed

# -------------------------------------------------------------------------------------- #
# --- Remove *_random and *_hap* entries in BED - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Cleaned BED file depends on original BED file
data/gphocs_regions_human_untr.noRandom.bed : data/gphocs_regions_human_untr.bed
	@echo "# === Cleaning BED file - Untranscribed dataset =============================== #";
	grep -v "_random" data/gphocs_regions_human_untr.bed | grep -v "_hap" > data/gphocs_regions_human_untr.noRandom.bed

# -------------------------------------------------------------------------------------- #
# --- Merge overlapping regions in BED file - Full dataset
# -------------------------------------------------------------------------------------- #

# Merged BED file depends on cleaned BED file
data/gphocs_regions_human_full.noRandom.merged.bed : data/gphocs_regions_human_full.noRandom.bed
	@echo "# === Merging overlapping regions in BED file - Full dataset ================== #";
	${BEDTOOLS}/mergeBed -i data/gphocs_regions_human_full.noRandom.bed > data/gphocs_regions_human_full.noRandom.merged.bed

# -------------------------------------------------------------------------------------- #
# --- Merge overlapping regions in BED file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Merged BED file depends on cleaned BED file and BEDtools
data/gphocs_regions_human_untr.noRandom.merged.bed : data/gphocs_regions_human_untr.noRandom.bed ${BEDTOOLS}/*
	@echo "# === Merging overlapping regions in BED file - Untranscribed dataset ========= #";
	${BEDTOOLS}/mergeBed -i data/gphocs_regions_human_untr.noRandom.bed > data/gphocs_regions_human_untr.noRandom.merged.bed

# -------------------------------------------------------------------------------------- #
# --- Grab human genome sequences for all regions in BED file - Full dataset
# -------------------------------------------------------------------------------------- #

# FASTA file depends on merged BED file, genome file, and BEDtools
data/hg18_gphocs_targets.full.fa : data/gphocs_regions_human_full.noRandom.merged.bed ${HUMAN_GENOME_FA} ${BEDTOOLS}/*
	@echo "# === Making FASTA of genome sequences - Full dataset ========================= #";
	${BEDTOOLS}/fastaFromBed -fi ${HUMAN_GENOME_FA} -bed data/gphocs_regions_human_full.noRandom.merged.bed -fo data/hg18_gphocs_targets.full.fa

# -------------------------------------------------------------------------------------- #
# --- Grab human genome sequences for all regions in BED file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# FASTA file depends on merged BED file, genome file, and BEDtools
data/hg18_gphocs_targets.untr.fa : data/gphocs_regions_human_untr.noRandom.merged.bed ${HUMAN_GENOME_FA} ${BEDTOOLS}/*
	@echo "# === Making FASTA of genome sequences - Untranscribed dataset ================ #";
	${BEDTOOLS}/fastaFromBed -fi ${HUMAN_GENOME_FA} -bed data/gphocs_regions_human_untr.noRandom.merged.bed -fo data/hg18_gphocs_targets.untr.fa

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Preliminary steps to grab chimpanzee genome sequences
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert human coordinates to chimp - Full dataset
# -------------------------------------------------------------------------------------- #

# Liftover'd chimp coordinates depend on human merged BED file, LiftOver, and chain file
data/gphocs_regions_human_full.noRandom.merged.panTro2.bed : data/gphocs_regions_human_full.noRandom.merged.bed ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE}
	@echo "# === Making FASTA of genome sequences - Untranscribed dataset ================ #";
	${LIFTOVER}/liftOver data/gphocs_regions_human_full.noRandom.merged.bed ${TO_CHIMP_CHAINFILE} data/gphocs_regions_human_full.noRandom.merged.panTro2.bed data/gphocs_regions_human_full.noRandom.merged.panTro2.unmapped.bed

# -------------------------------------------------------------------------------------- #
# --- Convert human coordinates to chimp - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Liftover'd chimp coordinates depend on human merged BED file, LiftOver, and chain file
data/gphocs_regions_human_untr.noRandom.merged.panTro2.bed : data/gphocs_regions_human_untr.noRandom.merged.bed ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE}
	@echo "# === Making FASTA of genome sequences - Untranscribed dataset ================ #";
	${LIFTOVER}/liftOver data/gphocs_regions_human_untr.noRandom.merged.bed ${TO_CHIMP_CHAINFILE} data/gphocs_regions_human_untr.noRandom.merged.panTro2.bed data/gphocs_regions_human_untr.noRandom.merged.panTro2.unmapped.bed

# -------------------------------------------------------------------------------------- #
# --- Remove *_random and *_hap* entries in chimp BED - Full dataset
# -------------------------------------------------------------------------------------- #

# Cleaned chimp BED file depends on original chimp BED file
data/gphocs_regions_human_full.noRandom.merged.panTro2.noRandom.bed : data/gphocs_regions_human_full.noRandom.merged.panTro2.bed
	@echo "# === Cleaning chimp BED file - Full dataset ================================== #";
	grep -v "_random" data/gphocs_regions_human_full.noRandom.merged.panTro2.bed | grep -v "_hap" > data/gphocs_regions_human_full.noRandom.merged.panTro2.noRandom.bed

# -------------------------------------------------------------------------------------- #
# --- Remove *_random and *_hap* entries in chimp BED - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Cleaned chimp BED file depends on original chimp BED file
data/gphocs_regions_human_untr.noRandom.merged.panTro2.noRandom.bed : data/gphocs_regions_human_untr.noRandom.merged.panTro2.bed
	@echo "# === Cleaning chimp BED file - Untranscribed dataset ========================= #";
	grep -v "_random" data/gphocs_regions_human_untr.noRandom.merged.panTro2.bed | grep -v "_hap" > data/gphocs_regions_human_untr.noRandom.merged.panTro2.noRandom.bed

# -------------------------------------------------------------------------------------- #
# --- Grab chimp genome sequences for all regions in BED file - Full dataset
# -------------------------------------------------------------------------------------- #

# FASTA file depends on merged BED file, genome file, and BEDtools
data/panTro2_gphocs_targets.full.fa : data/gphocs_regions_human_full.noRandom.merged.panTro2.noRandom.bed ${CHIMP_GENOME_FA} ${BEDTOOLS}/*
	@echo "# === Making FASTA of chimp genome sequences - Full dataset =================== #";
	${BEDTOOLS}/fastaFromBed -fi ${CHIMP_GENOME_FA} -bed data/gphocs_regions_human_full.noRandom.merged.panTro2.noRandom.bed -fo data/panTro2_gphocs_targets.full.fa

# -------------------------------------------------------------------------------------- #
# --- Grab chimp genome sequences for all regions in BED file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# FASTA file depends on merged BED file, genome file, and BEDtools
data/panTro2_gphocs_targets.untr.fa : data/gphocs_regions_human_untr.noRandom.merged.panTro2.noRandom.bed ${CHIMP_GENOME_FA} ${BEDTOOLS}/*
	@echo "# === Making FASTA of chimp genome sequences - Untranscribed dataset ========== #";
	${BEDTOOLS}/fastaFromBed -fi ${CHIMP_GENOME_FA} -bed data/gphocs_regions_human_untr.noRandom.merged.panTro2.noRandom.bed -fo data/panTro2_gphocs_targets.untr.fa

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Generate fake human exome alignments
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Make fake human exome seq file - Full dataset
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
#  qsub -t 1-5017:100 pbs/make_exomes_full_parallel.pbs
# Then recombine output with:
#  cat `ls -v *full_parallel*.o*` > fake_human_exomes_full.seq

# Fake exome sequence file depends on human sequence FASTA, LiftOver, chain file, and chimp sequence FASTA
fake_human_exomes_full.seq : data/hg18_gphocs_targets.full.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.full.fa
	@echo "# === Generating fake exomes - Full dataset =================================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.full.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.full.fa > fake_human_exomes_full.seq

# -------------------------------------------------------------------------------------- #
# --- Make fake human exome seq file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3698:100 pbs/make_exomes_untr_parallel.pbs
# Then recombine output with:
# cat `ls -v *untr_parallel*.o*` > fake_human_exomes_untr.seq

# Fake exome sequence file depends on human sequence FASTA, LiftOver, chain file, and chimp sequence FASTA
fake_human_exomes_untr.seq : data/hg18_gphocs_targets.untr.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.untr.fa
	@echo "# === Generating fake exomes - Untranscribed dataset ========================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.untr.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.untr.fa > fake_human_exomes_untr.seq

# -------------------------------------------------------------------------------------- #
# --- Reduce to autosomes - Full dataset
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Reduce to autosomes - Untranscribed only
# -------------------------------------------------------------------------------------- #


# Manually get rid of chrY and chrX

# -------------------------------------------------------------------------------------- #
# --- Remove aberrant sequences (and one San) - Full dataset
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Remove aberrant sequences (and one San) - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Remove funky seqs, one of San
# perl remove_missing_sequences.pl fake_human_exomes_FULL.seq > fake_human_exomes_FULL.filtered.seq

# -------------------------------------------------------------------------------------- #
# --- Mask CpG regions in sequences - Full dataset
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Mask CpG regions in sequences - Untranscribed only
# -------------------------------------------------------------------------------------- #


# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Compute various pop gen statistics
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Compute pop gen stats for filtration
# -------------------------------------------------------------------------------------- #

# Just full dataset
# perl scripts/get_aln_stats_all_loci_human.pl *.seq
# Automatically makes *.stats.txt

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Filter sequences
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - NoNA
# -------------------------------------------------------------------------------------- #

# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt tajima 100 > filtered/human.CpG.noNA.seq

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Tajima's D
# -------------------------------------------------------------------------------------- #

# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt tajima 2 > filtered/human.CpG.taj.2.0.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt tajima 2.5 > filtered/human.CpG.taj.2.5.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt tajima 3 > filtered/human.CpG.taj.3.0.seq

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Fu & Li's D
# -------------------------------------------------------------------------------------- #

# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli 2 > filtered/human.CpG.fuli.2.0.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli 2.5 > filtered/human.CpG.fuli.2.5.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli 3 > filtered/human.CpG.fuli.3.0.seq

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Fu & Li's D*
# -------------------------------------------------------------------------------------- #
 
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli_star 2 > filtered/human.CpG.fuli_star.2.0.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli_star 2.5 > filtered/human.CpG.fuli_star.2.5.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt fuli_star 3 > filtered/human.CpG.fuli_star.3.0.seq

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - GC percentage
# -------------------------------------------------------------------------------------- #

# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt gc 0.5 > filtered/human.CpG.gc.50.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt gc 0.55 > filtered/human.CpG.gc.55.seq
# perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_FULL.filtered.CpGmasked.seq fake_human_exomes_FULL.filtered.CpGmasked.stats.txt gc 0.6 > filtered/human.CpG.gc.60.seq

# -------------------------------------------------------------------------------------- #
# --- Fix header (with loci number) of all filtered seq files
# -------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------- #
# --- Compute pop gen stats on filtered datasets
# -------------------------------------------------------------------------------------- #

# perl scripts/get_aln_stats_all_loci_human.pl *.seq
# Automatically makes *.stats.txt

# -------------------------------------------------------------------------------------- #
# --- Compute mu for all seq files
# -------------------------------------------------------------------------------------- #

# Estimate mu as pi per site 
# From each *.stats.txt, using estimate_mu_from_pi.R

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mask out CCDS regions
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# perl filter_gphocs_seq.pl fake_human_exomes_FULL.filtered.CpGmasked.seq ../xspecies-exome/targets/ccdsGene.hg19.4apr12.bed > fake_human_exomes_UNTR.filtered.CpGmasked.seq

# perl remove_missing_sequences.pl fake_human_exomes_UNTR.CpGmasked.seq > fake_human_exomes_UNTR.CpGmasked.filtered.seq

# Fix header

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mask by codon
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Find codon locations with RefGene and get_codons_from_refgene.R

# Get rid of anything with only one column:
# awk -F'\t' '$2 != ""' hg18_refGene_codon1.bed > hg18_refGene_codon1.fixed.bed

# Combine codons to get [1st and 2nd] and [1st, 2nd, and NA]

# Sort with BEDtools’ sortBed

# Mask by codon

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Randomly subset neutral dataset
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Use filtering_scripts/randomly_subset_seqs.pl

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Run G-PhoCS (many, many times)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Make ctl file human_exome_FULL.ctl (and for untr)
# and PBS file call_gphocs_FULL.pbs (and for untr)

# Run G-PhoCS, twice for each

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Analyze G-PhoCS output
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# Extract results

# Scale results

# Grab numbers for paper