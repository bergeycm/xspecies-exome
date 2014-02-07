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

# Fake exome sequence file depends on human sequence FASTA, LiftOver, chain file, and chimp sequence FASTA
fake_human_exomes_full.seq : data/hg18_gphocs_targets.full.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.full.fa
	@echo "# === Generating fake exomes - Full dataset =================================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.full.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.full.fa > fake_human_exomes_full.seq

# -------------------------------------------------------------------------------------- #
# --- Make fake human exome seq file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Fake exome sequence file depends on human sequence FASTA, LiftOver, chain file, and chimp sequence FASTA
fake_human_exomes_untr.seq : data/hg18_gphocs_targets.untr.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.untr.fa
	@echo "# === Generating fake exomes - Untranscribed dataset ========================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.untr.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.untr.fa > fake_human_exomes_untr.seq
