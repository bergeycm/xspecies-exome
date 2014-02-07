# -------------------------------------------------------------------------------------- #
# --- Makefile to generate and analyze fake human exomes
# --- Called by the executable shell script, analyze_fake_human_exomes
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

HUMAN_GENOME_DIR=$(dir ${HUMAN_GENOME_FA})

# Steps. Can be called one-by-one with something like, make liftover_bed
# --- prepare_to_make_exomes:
liftover_bed_full : data/gphocs_regions_human_full.bed
liftover_bed_untr : data/gphocs_regions_human_untr.bed
clean_bed_full : data/gphocs_regions_human_full.noRandom.bed
clean_bed_untr : data/gphocs_regions_human_untr.noRandom.bed
merge_bed_full : data/gphocs_regions_human_full.noRandom.merged.bed
merge_bed_untr : data/gphocs_regions_human_untr.noRandom.merged.bed
grab_genome_seq_full : data/hg18_gphocs_targets.full.fa
grab_genome_seq_untr : data/hg18_gphocs_targets.untr.fa
# --- 

# Group steps together

generate_exomes : liftover_bed_full liftover_bed_untr clean_bed_full clean_bed_untr merge_bed_full merge_bed_untr grab_genome_seq_full grab_genome_seq_untr

all : generate_exomes

SHELL_EXPORT := 

# Export Make variables to child scripts
.EXPORT_ALL_VARIABLES :

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Generate fake human exome alignments
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Convert coordinates from rhesus exomes to human coordinates - Full dataset
# -------------------------------------------------------------------------------------- #

# Liftover'd human BED file depends on LiftOver, macaque BED file, and chain file
data/gphocs_regions_human_full.bed : ${LIFTOVER}/* data/all.bsnp.snp.out.gt4.large.bed
	@echo "# === Converting rhesus coordinates to human - Full dataset =================== #";
	${LIFTOVER}/liftOver data/all.bsnp.snp.out.gt4.large.bed ${TO_HUMAN_CHAINFILE} data/gphocs_regions_human_full.bed data/gphocs_regions_human_full.unMapped.bed;

# -------------------------------------------------------------------------------------- #
# --- Convert coordinates from rhesus exomes to human coordinates - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Liftover'd human BED file depends on LiftOver, macaque BED file, and chain file
data/gphocs_regions_human_untr.bed : ${LIFTOVER}/* data/all.combined.gphocs.seq.untranscribed.bed
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



