# -------------------------------------------------------------------------------------- #
# --- Makefile to generate and analyze fake human exomes
# --- Called by the executable shell script, analyze_fake_human_exomes
# -------------------------------------------------------------------------------------- #

# Get user editable variables
include config.mk

HUMAN_GENOME_DIR=$(dir ${HUMAN_GENOME_FA})

# List of seq files for G-PhoCS analyses
SEQ_ORIG=fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_untr.filtered.CpGmasked.seq
SEQ_NONA=filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.seq 
SEQ_FILTERED=filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.seq
SEQ_MASKED=fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.seq fake_human_exomes_full.filtered.CpGmasked.noCodon2.seq fake_human_exomes_full.filtered.CpGmasked.noCodon12.seq fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.seq 
SEQ_SUBSETS_FULL=subsets/full_subset1.seq subsets/full_subset2.seq subsets/full_subset3.seq subsets/full_subset4.seq subsets/full_subset5.seq 
SEQ_SUBSETS_UNTR=subsets/untr_subset1.seq subsets/untr_subset2.seq subsets/untr_subset3.seq subsets/untr_subset4.seq subsets/untr_subset5.seq 
SEQ_SUBSETS_NONA=subsets/noNA_subset1.seq subsets/noNA_subset2.seq subsets/noNA_subset3.seq subsets/noNA_subset4.seq subsets/noNA_subset5.seq 
SEQ_SUBSETS_TAJ2=subsets/taj.2.0_subset1.seq subsets/taj.2.0_subset2.seq subsets/taj.2.0_subset3.seq subsets/taj.2.0_subset4.seq subsets/taj.2.0_subset5.seq
SEQ_SUBSETS_TAJ3=subsets/taj.3.0_subset1.seq subsets/taj.3.0_subset2.seq subsets/taj.3.0_subset3.seq subsets/taj.3.0_subset4.seq subsets/taj.3.0_subset5.seq 
SEQ_SUBSETS_FULI2=subsets/fuli.2.0_subset1.seq subsets/fuli.2.0_subset2.seq subsets/fuli.2.0_subset3.seq subsets/fuli.2.0_subset4.seq subsets/fuli.2.0_subset5.seq 
SEQ_SUBSETS_FULI3=subsets/fuli.3.0_subset1.seq subsets/fuli.3.0_subset2.seq subsets/fuli.3.0_subset3.seq subsets/fuli.3.0_subset4.seq subsets/fuli.3.0_subset5.seq 
SEQ_ALL=${SEQ_ORIG} ${SEQ_NONA} ${SEQ_FILTERED} ${SEQ_MASKED} ${SEQ_SUBSETS_FULL} ${SEQ_SUBSETS_UNTR} ${SEQ_SUBSETS_NONA} ${SEQ_SUBSETS_TAJ2} ${SEQ_SUBSETS_TAJ3} ${SEQ_SUBSETS_FULI2} ${SEQ_SUBSETS_FULI3}

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
make_exomes_full : fake_human_exomes_full_all.seq 
make_exomes_untr : fake_human_exomes_untr_all.seq
autosomes_only_full : fake_human_exomes_full.seq
autosomes_only_untr : fake_human_exomes_untr.seq
filter_seq_full : fake_human_exomes_full.filtered.seq
filter_seq_untr : fake_human_exomes_untr.filtered.seq
mask_CpG_full : fake_human_exomes_full.filtered.CpGmasked.seq
mask_CpG_untr : fake_human_exomes_untr.filtered.CpGmasked.seq
# --- filter_seqs
get_stats_pre_filter : fake_human_exomes_full.filtered.CpGmasked.stats.txt
filter_noNA : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.seq
filter_taj : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.seq
filter_fuli : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.seq
filter_fuli_star : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.seq
filter_gc : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.seq filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.seq
get_filtered_seq_stats : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.stats.txt filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.stats.txt
estimate_mu : filtered_seqs/mu_estimates.txt
# --- make_ccds_masked
mask_ccds : fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq
filter_ccds_masked : fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.seq
# --- mask_by_codon
find_codons : data/hg18_refGene_codon2.bed
combine_codons : data/hg18_refGene_codon12.bed data/hg18_refGene_codon12NA.bed
codon_masking : fake_human_exomes_full.filtered.CpGmasked.noCodon2.seq fake_human_exomes_full.filtered.CpGmasked.noCodon12.seq fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.seq
# --- randomly_subset
get_loci_list : data/neutralLoci-7genomes_locus_list.txt
subset_full_untr : subsets/full_subset5.seq subsets/untr_subset5.seq
subset_filtered : subsets/noNA_subset5.seq subsets/taj.2.0_subset5.seq subsets/taj.3.0_subset5.seq subsets/fuli.2.0_subset5.seq subsets/fuli.3.0_subset5.seq
# --- prep_to_run_gphocs
make_ctl_files_orig : fake_human_exomes_full.filtered.CpGmasked.iter1.ctl fake_human_exomes_untr.filtered.CpGmasked.iter1.ctl 
make_ctl_files_noNA : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.iter1.ctl
make_ctl_files_taj : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.iter1.ctl
make_ctl_files_fuli : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.iter1.ctl
make_ctl_files_fuli_star : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.iter1.ctl
make_ctl_files_gc : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.iter1.ctl filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.iter1.ctl
make_ctl_files_ccds : fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.iter1.ctl
make_ctl_files_codon : fake_human_exomes_full.filtered.CpGmasked.noCodon2.iter1.ctl fake_human_exomes_full.filtered.CpGmasked.noCodon12.iter1.ctl fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.iter1.ctl
make_ctl_files_subsets_full : subsets/full_subset1.iter1.ctl subsets/full_subset2.iter1.ctl subsets/full_subset3.iter1.ctl subsets/full_subset4.iter1.ctl subsets/full_subset5.iter1.ctl
make_ctl_files_subsets_untr : subsets/untr_subset1.iter1.ctl subsets/untr_subset2.iter1.ctl subsets/untr_subset3.iter1.ctl subsets/untr_subset4.iter1.ctl subsets/untr_subset5.iter1.ctl
make_ctl_files_subsets_noNA : subsets/noNA_subset1.iter1.ctl subsets/noNA_subset2.iter1.ctl subsets/noNA_subset3.iter1.ctl subsets/noNA_subset4.iter1.ctl subsets/noNA_subset5.iter1.ctl 
make_ctl_files_subsets_taj.2.0 : subsets/taj.2.0_subset1.iter1.ctl subsets/taj.2.0_subset2.iter1.ctl subsets/taj.2.0_subset3.iter1.ctl subsets/taj.2.0_subset4.iter1.ctl subsets/taj.2.0_subset5.iter1.ctl
make_ctl_files_subsets_taj.3.0 : subsets/taj.3.0_subset1.iter1.ctl subsets/taj.3.0_subset2.iter1.ctl subsets/taj.3.0_subset3.iter1.ctl subsets/taj.3.0_subset4.iter1.ctl subsets/taj.3.0_subset5.iter1.ctl 
make_ctl_files_subsets_fuli.2.0 : subsets/fuli.2.0_subset1.iter1.ctl subsets/fuli.2.0_subset2.iter1.ctl subsets/fuli.2.0_subset3.iter1.ctl subsets/fuli.2.0_subset4.iter1.ctl subsets/fuli.2.0_subset5.iter1.ctl 
make_ctl_files_subsets_fuli.3.0 : subsets/fuli.3.0_subset1.iter1.ctl subsets/fuli.3.0_subset2.iter1.ctl subsets/fuli.3.0_subset3.iter1.ctl subsets/fuli.3.0_subset4.iter1.ctl subsets/fuli.3.0_subset5.iter1.ctl 
make_ctl_all_filtered : make_ctl_files_noNA make_ctl_files_taj make_ctl_files_fuli make_ctl_files_fuli_star make_ctl_files_gc
make_ctl_files_all_subsets : make_ctl_files_subsets_full make_ctl_files_subsets_untr make_ctl_files_subsets_noNA make_ctl_files_subsets_taj.2.0 make_ctl_files_subsets_taj.3.0 make_ctl_files_subsets_fuli.2.0 make_ctl_files_subsets_fuli.3.0
# --- run_gphocs_iter1
call_gphocs_iter1_orig : fake_human_exomes_full.filtered.CpGmasked.iter1.trace.log fake_human_exomes_untr.filtered.CpGmasked.iter1.trace.log 
call_gphocs_iter1_noNA : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.iter1.trace.log
call_gphocs_iter1_taj : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.iter1.trace.log
call_gphocs_iter1_fuli : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.iter1.trace.log
call_gphocs_iter1_fuli_star : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.iter1.trace.log
call_gphocs_iter1_gc : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.iter1.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.iter1.trace.log
call_gphocs_iter1_ccds : fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.iter1.trace.log
call_gphocs_iter1_codon : fake_human_exomes_full.filtered.CpGmasked.noCodon2.iter1.trace.log fake_human_exomes_full.filtered.CpGmasked.noCodon12.iter1.trace.log fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.iter1.trace.log
call_gphocs_iter1_subsets_full : subsets/full_subset1.iter1.trace.log subsets/full_subset2.iter1.trace.log subsets/full_subset3.iter1.trace.log subsets/full_subset4.iter1.trace.log subsets/full_subset5.iter1.trace.log
call_gphocs_iter1_subsets_untr : subsets/untr_subset1.iter1.trace.log subsets/untr_subset2.iter1.trace.log subsets/untr_subset3.iter1.trace.log subsets/untr_subset4.iter1.trace.log subsets/untr_subset5.iter1.trace.log
call_gphocs_iter1_subsets_noNA : subsets/noNA_subset1.iter1.trace.log subsets/noNA_subset2.iter1.trace.log subsets/noNA_subset3.iter1.trace.log subsets/noNA_subset4.iter1.trace.log subsets/noNA_subset5.iter1.trace.log 
call_gphocs_iter1_subsets_taj.2.0 : subsets/taj.2.0_subset1.iter1.trace.log subsets/taj.2.0_subset2.iter1.trace.log subsets/taj.2.0_subset3.iter1.trace.log subsets/taj.2.0_subset4.iter1.trace.log subsets/taj.2.0_subset5.iter1.trace.log
call_gphocs_iter1_subsets_taj.3.0 : subsets/taj.3.0_subset1.iter1.trace.log subsets/taj.3.0_subset2.iter1.trace.log subsets/taj.3.0_subset3.iter1.trace.log subsets/taj.3.0_subset4.iter1.trace.log subsets/taj.3.0_subset5.iter1.trace.log 
call_gphocs_iter1_subsets_fuli.2.0 : subsets/fuli.2.0_subset1.iter1.trace.log subsets/fuli.2.0_subset2.iter1.trace.log subsets/fuli.2.0_subset3.iter1.trace.log subsets/fuli.2.0_subset4.iter1.trace.log subsets/fuli.2.0_subset5.iter1.trace.log 
call_gphocs_iter1_subsets_fuli.3.0 : subsets/fuli.3.0_subset1.iter1.trace.log subsets/fuli.3.0_subset2.iter1.trace.log subsets/fuli.3.0_subset3.iter1.trace.log subsets/fuli.3.0_subset4.iter1.trace.log subsets/fuli.3.0_subset5.iter1.trace.log 
call_gphocs_all_iter1_filtered : call_gphocs_iter1_noNA call_gphocs_iter1_taj call_gphocs_iter1_fuli call_gphocs_iter1_fuli_star call_gphocs_iter1_gc
call_gphocs_all_iter1_subsets : call_gphocs_iter1_subsets_full call_gphocs_iter1_subsets_untr call_gphocs_iter1_subsets_noNA call_gphocs_iter1_subsets_taj.2.0 call_gphocs_iter1_subsets_taj.3.0 call_gphocs_iter1_subsets_fuli.2.0 call_gphocs_iter1_subsets_fuli.3.0
# --- run_gphocs_iter1
call_gphocs_iter2_orig : fake_human_exomes_full.filtered.CpGmasked.iter2.trace.log fake_human_exomes_untr.filtered.CpGmasked.iter2.trace.log 
call_gphocs_iter2_noNA : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.iter2.trace.log
call_gphocs_iter2_taj : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.iter2.trace.log
call_gphocs_iter2_fuli : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.iter2.trace.log
call_gphocs_iter2_fuli_star : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.iter2.trace.log
call_gphocs_iter2_gc : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.iter2.trace.log filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.iter2.trace.log
call_gphocs_iter2_ccds : fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.iter2.trace.log
call_gphocs_iter2_codon : fake_human_exomes_full.filtered.CpGmasked.noCodon2.iter2.trace.log fake_human_exomes_full.filtered.CpGmasked.noCodon12.iter2.trace.log fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.iter2.trace.log
call_gphocs_iter2_subsets_full : subsets/full_subset1.iter2.trace.log subsets/full_subset2.iter2.trace.log subsets/full_subset3.iter2.trace.log subsets/full_subset4.iter2.trace.log subsets/full_subset5.iter2.trace.log
call_gphocs_iter2_subsets_untr : subsets/untr_subset1.iter2.trace.log subsets/untr_subset2.iter2.trace.log subsets/untr_subset3.iter2.trace.log subsets/untr_subset4.iter2.trace.log subsets/untr_subset5.iter2.trace.log
call_gphocs_iter2_subsets_noNA : subsets/noNA_subset1.iter2.trace.log subsets/noNA_subset2.iter2.trace.log subsets/noNA_subset3.iter2.trace.log subsets/noNA_subset4.iter2.trace.log subsets/noNA_subset5.iter2.trace.log 
call_gphocs_iter2_subsets_taj.2.0 : subsets/taj.2.0_subset1.iter2.trace.log subsets/taj.2.0_subset2.iter2.trace.log subsets/taj.2.0_subset3.iter2.trace.log subsets/taj.2.0_subset4.iter2.trace.log subsets/taj.2.0_subset5.iter2.trace.log
call_gphocs_iter2_subsets_taj.3.0 : subsets/taj.3.0_subset1.iter2.trace.log subsets/taj.3.0_subset2.iter2.trace.log subsets/taj.3.0_subset3.iter2.trace.log subsets/taj.3.0_subset4.iter2.trace.log subsets/taj.3.0_subset5.iter2.trace.log 
call_gphocs_iter2_subsets_fuli.2.0 : subsets/fuli.2.0_subset1.iter2.trace.log subsets/fuli.2.0_subset2.iter2.trace.log subsets/fuli.2.0_subset3.iter2.trace.log subsets/fuli.2.0_subset4.iter2.trace.log subsets/fuli.2.0_subset5.iter2.trace.log 
call_gphocs_iter2_subsets_fuli.3.0 : subsets/fuli.3.0_subset1.iter2.trace.log subsets/fuli.3.0_subset2.iter2.trace.log subsets/fuli.3.0_subset3.iter2.trace.log subsets/fuli.3.0_subset4.iter2.trace.log subsets/fuli.3.0_subset5.iter2.trace.log 
call_gphocs_all_iter2_filtered : call_gphocs_iter2_noNA call_gphocs_iter2_taj call_gphocs_iter2_fuli call_gphocs_iter2_fuli_star call_gphocs_iter2_gc
call_gphocs_all_iter2_subsets : call_gphocs_iter2_subsets_full call_gphocs_iter2_subsets_untr call_gphocs_iter2_subsets_noNA call_gphocs_iter2_subsets_taj.2.0 call_gphocs_iter2_subsets_taj.3.0 call_gphocs_iter2_subsets_fuli.2.0 call_gphocs_iter2_subsets_fuli.3.0
# --- summarize_gphocs
RESULTS_ITER_1=$(SEQ_ALL:.seq=.iter1.results)
RESULTS_ITER_2=$(SEQ_ALL:.seq=.iter2.results)
summarize_gphocs_iter_1 : ${RESULTS_ITER_1}
summarize_gphocs_iter_2 : ${RESULTS_ITER_2}

# Group steps together

prelim_grab_genome : liftover_bed_full liftover_bed_untr clean_bed_full clean_bed_untr merge_bed_full merge_bed_untr grab_genome_seq_full grab_genome_seq_untr 
prelim_grab_chimp : liftover_to_chimp_full liftover_to_chimp_untr clean_chimp_bed_full clean_chimp_bed_untr grab_chimp_seq_full grab_chimp_seq_untr
generate_exomes : make_exomes_full make_exomes_untr autosomes_only_full autosomes_only_untr filter_seq_full filter_seq_untr mask_CpG_full mask_CpG_untr
filter_seqs : get_stats_pre_filter filter_noNA filter_taj filter_fuli filter_fuli_star filter_gc get_filtered_seq_stats estimate_mu
make_ccds_masked : mask_ccds filter_ccds_masked
mask_by_codon : find_codons combine_codons codon_masking
randomly_subset : get_loci_list subset_full_untr subset_filtered
prep_to_run_gphocs : make_ctl_files_orig make_ctl_all_filtered make_ctl_files_ccds make_ctl_files_codon make_ctl_files_all_subsets
run_gphocs_iter1 : call_gphocs_iter1_orig call_gphocs_all_iter1_filtered call_gphocs_iter1_ccds call_gphocs_iter1_codon call_gphocs_all_iter1_subsets
run_gphocs_iter2 : call_gphocs_iter2_orig call_gphocs_all_iter2_filtered call_gphocs_iter2_ccds call_gphocs_iter2_codon call_gphocs_all_iter2_subsets
summarize_gphocs : summarize_gphocs_iter_1 summarize_gphocs_iter_2

all : prelim_grab_genome prelim_grab_chimp generate_exomes filter_seqs make_ccds_masked mask_by_codon randomly_subset prep_to_run_gphocs run_gphocs_iter1 run_gphocs_iter2 summarize_gphocs

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
fake_human_exomes_full_all.seq : data/hg18_gphocs_targets.full.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.full.fa
	@echo "# === Generating fake exomes - Full dataset =================================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.full.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.full.fa > fake_human_exomes_full_all.seq

# -------------------------------------------------------------------------------------- #
# --- Make fake human exome seq file - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3698:100 pbs/make_exomes_untr_parallel.pbs
# Then recombine output with:
# cat `ls -v *untr_parallel*.o*` > fake_human_exomes_untr.seq

# Fake exome sequence file depends on human sequence FASTA, LiftOver, chain file, and chimp sequence FASTA
fake_human_exomes_untr_all.seq : data/hg18_gphocs_targets.untr.fa ${LIFTOVER}/liftOver ${TO_CHIMP_CHAINFILE} data/panTro2_gphocs_targets.untr.fa
	@echo "# === Generating fake exomes - Untranscribed dataset ========================== #";
	perl scripts/make_fake_human_exomes.pl --target_fasta data/hg18_gphocs_targets.untr.fa --liftover_path ${LIFTOVER} --hg_to_pan_chain ${TO_CHIMP_CHAINFILE} --chimp_seqs_fasta data/panTro2_gphocs_targets.untr.fa > fake_human_exomes_untr_all.seq

# -------------------------------------------------------------------------------------- #
# --- Reduce to autosomes - Full dataset
# -------------------------------------------------------------------------------------- #

# Careful! This assumes there are 8 individuals per locus 
# and thus nine lines between headers. This might not always be the case.

# Reduced fake exome sequence file depends on seq file with autosomes
fake_human_exomes_full.seq : fake_human_exomes_full_all.seq
	@echo "# === Reduce to autosomes - Full dataset ====================================== #";
	sed -e '/chrX/,+9d' -e '/chrY/,+9d' < fake_human_exomes_full_all.seq > fake_human_exomes_full_tmp.seq
	grep -c "chr" fake_human_exomes_full_tmp.seq > fake_human_exomes_full_tmp.seq.count
	cat fake_human_exomes_full_tmp.seq.count fake_human_exomes_full_tmp.seq > fake_human_exomes_full.seq
	rm fake_human_exomes_full_tmp.seq
	rm fake_human_exomes_full_tmp.seq.count

# -------------------------------------------------------------------------------------- #
# --- Reduce to autosomes - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Careful! This assumes there are 8 individuals per locus 
# and thus nine lines between headers. This might not always be the case.

# Reduced fake exome sequence file depends on seq file with autosomes
fake_human_exomes_untr.seq : fake_human_exomes_untr_all.seq
	@echo "# === Reduce to autosomes - Untranscribed dataset ============================= #";
	sed -e '/chrX/,+9d' -e '/chrY/,+9d' < fake_human_exomes_untr_all.seq > fake_human_exomes_untr_tmp.seq
	grep -c "chr" fake_human_exomes_untr_tmp.seq > fake_human_exomes_untr_tmp.seq.count
	cat fake_human_exomes_untr_tmp.seq.count fake_human_exomes_untr_tmp.seq > fake_human_exomes_untr.seq
	rm fake_human_exomes_untr_tmp.seq
	rm fake_human_exomes_untr_tmp.seq.count

# -------------------------------------------------------------------------------------- #
# --- Remove aberrant sequences (and one San) - Full dataset
# -------------------------------------------------------------------------------------- #

# Filtered sequence file depends on unfiltered seq file
fake_human_exomes_full.filtered.seq : fake_human_exomes_full.seq
	@echo "# === Filter seq file - Full dataset ========================================== #";
	perl scripts/filter_seq_file.pl fake_human_exomes_full.seq > fake_human_exomes_full.filtered_tmp.seq
	grep -c "chr" fake_human_exomes_full.filtered_tmp.seq > fake_human_exomes_full.filtered_tmp.seq.count
	cat fake_human_exomes_full.filtered_tmp.seq.count fake_human_exomes_full.filtered_tmp.seq > fake_human_exomes_full.filtered.seq
	rm fake_human_exomes_full.filtered_tmp.seq
	rm fake_human_exomes_full.filtered_tmp.seq.count

# -------------------------------------------------------------------------------------- #
# --- Remove aberrant sequences (and one San) - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Filtered sequence file depends on unfiltered seq file
fake_human_exomes_untr.filtered.seq : fake_human_exomes_untr.seq
	@echo "# === Filter seq file - Untranscribed dataset ================================= #";
	perl scripts/filter_seq_file.pl fake_human_exomes_untr.seq > fake_human_exomes_untr.filtered_tmp.seq
	grep -c "chr" fake_human_exomes_untr.filtered_tmp.seq > fake_human_exomes_untr.filtered_tmp.seq.count
	cat fake_human_exomes_untr.filtered_tmp.seq.count fake_human_exomes_untr.filtered_tmp.seq > fake_human_exomes_untr.filtered.seq
	rm fake_human_exomes_untr.filtered_tmp.seq
	rm fake_human_exomes_untr.filtered_tmp.seq.count

# -------------------------------------------------------------------------------------- #
# --- Mask CpG regions in sequences - Full dataset
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3349:100 pbs/mask_CpG_full_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_CpG_full_parallel*.o*` > fake_human_exomes_full.filtered.CpGmasked.seq

# CpG masked sequence file depends on filtered sequence file
fake_human_exomes_full.filtered.CpGmasked.seq : fake_human_exomes_full.filtered.seq
	@echo "# === Masking CpG regions - Full dataset ====================================== #";
	${SHELL_EXPORT} perl scripts/mask_gphocs_seq.pl fake_human_exomes_full.filtered.seq data/filter_CpG.bed > fake_human_exomes_full.filtered.CpGmasked.seq

# -------------------------------------------------------------------------------------- #
# --- Mask CpG regions in sequences - Untranscribed only
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-2401:100 pbs/mask_CpG_untr_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_CpG_untr_parallel*.o*` > fake_human_exomes_untr.filtered.CpGmasked.seq

# CpG masked sequence file depends on filtered sequence file
fake_human_exomes_untr.filtered.CpGmasked.seq : fake_human_exomes_untr.filtered.seq
	@echo "# === Masking CpG regions - Untranscribed dataset ============================= #";
	${SHELL_EXPORT} perl scripts/mask_gphocs_seq.pl fake_human_exomes_untr.filtered.seq data/filter_CpG.bed > fake_human_exomes_untr.filtered.CpGmasked.seq

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Filter sequences
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Compute pop gen stats for filtration - Full dataset
# -------------------------------------------------------------------------------------- #

# Stats file depends on CpG-masked sequence file
fake_human_exomes_full.filtered.CpGmasked.stats.txt : fake_human_exomes_full.filtered.CpGmasked.seq
	@echo "# === Computing pop gen stats - Full dataset ================================== #";
	perl scripts/get_aln_stats_all_loci_human.pl fake_human_exomes_full.filtered.CpGmasked.seq

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - NoNA
# -------------------------------------------------------------------------------------- #

# noNA filtered seq file depends on original CpG-masked sequence file and stats file
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing NA sequences ================================== #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt tajima 100 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.tmp.seq*

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Tajima's D
# -------------------------------------------------------------------------------------- #

# Tajima's D filtered seq files depend on original CpG-masked sequence file and stats file
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Tajima's D| > 2.0 sequences ================== #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt tajima 2.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Tajima's D| > 2.5 sequences ================== #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt tajima 2.5 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.5.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Tajima's D| > 3.0 sequences ================== #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt tajima 3.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.tmp.seq*

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Fu & Li's D
# -------------------------------------------------------------------------------------- #

# Fu & Li's D filtered seq files depend on original CpG-masked sequence file and stats file
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D| > 2.0 sequences ================= #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli 2.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.tmp.seq*
	
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D| > 2.5 sequences ================= #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli 2.5 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.5.tmp.seq*
	
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D| > 3.0 sequences ================= #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli 3.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.tmp.seq*

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - Fu & Li's D*
# -------------------------------------------------------------------------------------- #
 
# Fu & Li's D* filtered seq files depend on original CpG-masked sequence file and stats file
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D*| > 2.0 sequences ================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli_star 2.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.0.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D*| > 2.5 sequences ================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli_star 2.5 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.2.5.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing |Fu & Li's D*| > 3.0 sequences ================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt fuli_star 3.0 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli_star.3.0.tmp.seq*

# -------------------------------------------------------------------------------------- #
# --- Filter sequences - GC percentage
# -------------------------------------------------------------------------------------- #

# GC filtered seq files depend on original CpG-masked sequence file and stats file
filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing GC > 50% sequences ============================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt gc 0.50 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.50.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing GC > 55% sequences ============================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt gc 0.55 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.55.tmp.seq*

filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.seq : fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt
	@echo "# === Filtering loci - Removing GC > 60% sequences ============================ #";
	perl scripts/remove_various_outlier_seqs.pl fake_human_exomes_full.filtered.CpGmasked.seq fake_human_exomes_full.filtered.CpGmasked.stats.txt gc 0.60 > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq
	grep -c "chr" filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq.count
	cat filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq.count filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq > filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.seq
	rm filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.gc.60.tmp.seq*

# -------------------------------------------------------------------------------------- #
# --- Compute pop gen stats on filtered datasets
# -------------------------------------------------------------------------------------- #

# Skip running perl scripts/get_aln_stats_all_loci_human.pl since it takes awhile
# Instead just grab the results from the original unfiltered seq file's *.stats.txt file

# Filtered stats files depend on filtered seq files
filtered_seqs/%.stats.txt : filtered_seqs/%.seq
	@echo "# === Getting stats on filtered loci... ======================================= #";
	head -n1 fake_human_exomes_full.filtered.CpGmasked.stats.txt > $@
	grep "chr" $^ | cut -d' ' -f 1 > $^.tmp.chr
	awk 'FNR==NR{a[$$1]=$$1;next}{if (a[$$1]) { print $$0 }}' $^.tmp.chr fake_human_exomes_full.filtered.CpGmasked.stats.txt >> $@
	rm $^.tmp.chr

# -------------------------------------------------------------------------------------- #
# --- Compute mu (as pi per site) for all seq files
# -------------------------------------------------------------------------------------- #

# File of mu estimates depends on all stats files
filtered_seqs/mu_estimates.txt : fake_human_exomes_full.filtered.CpGmasked.stats.txt filtered_seqs/*.stats.txt
	@echo "# === Estimating mu from sequences ============================================ #";
	Rscript scripts/estimate_mu_from_pi.R > filtered_seqs/mu_estimates.txt

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mask out CCDS regions
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Mask out CCDS regions of full dataset
# -------------------------------------------------------------------------------------- #

# Essentially this accomplishes the same goal as the UNTR dataset, 
# but in this case the regions are masked, not removed.

# Parallel version exists too. Call with:
# qsub -t 1-3349:100 pbs/mask_ccds_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_ccds_parallel*.o*` > fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq

# CCDS masked sequence file depends on CpG-masked sequence file
fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq : fake_human_exomes_full.filtered.CpGmasked.seq
	@echo "# === Masking CCDS regions - Full dataset ===================================== #";
	${SHELL_EXPORT} perl scripts/mask_gphocs_seq.pl fake_human_exomes_full.filtered.CpGmasked.seq ${CCDS} > fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq

# -------------------------------------------------------------------------------------- #
# --- Remove aberrant sequences from CCDS-masked seq file
# -------------------------------------------------------------------------------------- #

# Filtered CCDS-masked sequence file depends on CCDS-masked seq file
fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.seq : fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq
	@echo "# === Filter CCDS-masked seq file ============================================= #";
	perl scripts/filter_seq_file.pl fake_human_exomes_full.filtered.CpGmasked.noCCDS.seq > fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq
	grep -c "chr" fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq > fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq.count
	cat fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq.count fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq > fake_human_exomes_full.filtered.CpGmasked.noCCDS.filtered.seq
	rm fake_human_exomes_full.filtered.CpGmasked.noCCDS_tmp.seq*

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Mask by codon
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Find codon locations from RefGene file
# -------------------------------------------------------------------------------------- #

# Last BED file of codon positions depends on data/refGene.txt file
data/hg18_refGene_codon2.bed : data/hg18_refGene.txt
	@echo "# === Find codon positions from RefGene ======================================= #";
	Rscript scripts/get_codons_from_refgene.R

# -------------------------------------------------------------------------------------- #
# --- Combine codons to get [1st and 2nd] and [1st, 2nd, and NA] and sort BED files
# -------------------------------------------------------------------------------------- #

# BED file of 1st and 2nd codons depends on BED files of 1st and of 2nd and on BEDtools
data/hg18_refGene_codon12.bed : data/hg18_refGene_codon1.bed data/hg18_refGene_codon2.bed ${BEDTOOLS}/*
	@echo "# === Combining and sorting codon BED files [1st and 2nd] ===================== #";
	cat data/hg18_refGene_codon1.bed data/hg18_refGene_codon2.bed | ${BEDTOOLS}/sortBed -i stdin > data/hg18_refGene_codon12.bed

# BED file of 1st, 2nd, and unknown codons depends on BED files of 1st, of 2nd, and of NA and on BEDtools
data/hg18_refGene_codon12NA.bed : data/hg18_refGene_codon1.bed data/hg18_refGene_codon2.bed data/hg18_refGene_unknown.bed ${BEDTOOLS}/*
	@echo "# === Combining and sorting codon BED files [1st, 2nd, and NA] ================ #";
	cat data/hg18_refGene_codon1.bed data/hg18_refGene_codon2.bed data/hg18_refGene_unknown.bed | ${BEDTOOLS}/sortBed -i stdin > data/hg18_refGene_codon12NA.bed

# -------------------------------------------------------------------------------------- #
# --- Mask full seq file by codon - Mask [2nd]
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3349:100 pbs/mask_codon2_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_codon2_parallel*.o*` > fake_human_exomes_full.filtered.CpGmasked.noCodon2.seq

# Masked seq file depends on original seq file and codon BED
fake_human_exomes_full.filtered.CpGmasked.noCodon2.seq : fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon2.bed
	@echo "# === Masking out 2nd codons ================================================== #";
	perl scripts/mask_gphocs_seq.pl fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon2.bed > fake_human_exomes_full.filtered.CpGmasked.noCodon2.seq

# -------------------------------------------------------------------------------------- #
# --- Mask full seq file by codon - Mask [1st and 2nd]
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3349:100 pbs/mask_codon12_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_codon12_parallel*.o*` > fake_human_exomes_full.filtered.CpGmasked.noCodon12.seq

# Masked seq file depends on original seq file and codon BED
fake_human_exomes_full.filtered.CpGmasked.noCodon12.seq : fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon12.bed
	@echo "# === Masking out 1st and 2nd codons ========================================== #";
	perl scripts/mask_gphocs_seq.pl fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon12.bed > fake_human_exomes_full.filtered.CpGmasked.noCodon12.seq

# -------------------------------------------------------------------------------------- #
# --- Mask full seq file by codon - Mask [1st, 2nd, and NA]
# -------------------------------------------------------------------------------------- #

# Parallel version exists too. Call with:
# qsub -t 1-3349:100 pbs/mask_codon12NA_parallel.pbs
# Then recombine output with:
# cat `ls -v mask_codon12NA_parallel*.o*` > fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.seq

# Masked seq file depends on original seq file and codon BED
fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.seq : fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon12NA.bed
	@echo "# === Masking out 1st, 2nd, and unknown codons ================================ #";
	perl scripts/mask_gphocs_seq.pl fake_human_exomes_full.filtered.CpGmasked.seq data/hg18_refGene_codon12NA.bed > fake_human_exomes_full.filtered.CpGmasked.noCodon12NA.seq

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Randomly subset neutral dataset
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Get list of loci used in Gronau et al. 2011's whole genome dataset
# -------------------------------------------------------------------------------------- #

# Loci list depends on whole genome seq file
data/neutralLoci-7genomes_locus_list.txt : data/neutralLoci-7genomes.txt
	@echo "# === Getting list of loci in whole genome dataset ============================ #";
	grep "^chr" data/neutralLoci-7genomes.txt | cut -f1 > data/neutralLoci-7genomes_locus_list.txt

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [full dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/full_subset5.seq : fake_human_exomes_full.filtered.CpGmasked.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [full dataset] ============ #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/full_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/full_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [untranscribed dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/untr_subset5.seq : fake_human_exomes_untr.filtered.CpGmasked.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [untranscribed dataset] === #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/untr_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/untr_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [no NA dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/noNA_subset5.seq : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.noNA.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [no NA dataset] =========== #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/noNA_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/noNA_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [Tajima's D < 2.0 dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/taj.2.0_subset5.seq : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.2.0.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [Tajima 2.0 dataset] ====== #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/taj.2.0_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/taj.2.0_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [Tajima's D < 3.0 dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/taj.3.0_subset5.seq : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.taj.3.0.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [Tajima 3.0 dataset] ====== #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/taj.3.0_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/taj.3.0_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [Fu & Li's D < 2.0 dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/fuli.2.0_subset5.seq : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.2.0.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [Fu & Li 2.0 dataset] ===== #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/fuli.2.0_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/fuli.2.0_subset$${number}.seq; \
	done

# -------------------------------------------------------------------------------------- #
# --- Make subset of neutral dataset that is the same size as [Fu & Li's D < 3.0 dataset]
# -------------------------------------------------------------------------------------- #

# Subset seq file depends on original sequence file (for size), neutral loci list file, and neutral loci seq file
subsets/fuli.3.0_subset5.seq : filtered_seqs/fake_human_exomes_full.filtered.CpGmasked.fuli.3.0.seq data/neutralLoci-7genomes_locus_list.txt data/neutralLoci-7genomes.txt
	@echo "# === Making subset of neutral dataset same size as [Fu & Li 3.0 dataset] ===== #";
	for number in 1 2 3 4 5; do \
		head -n1 $< > subsets/fuli.3.0_subset$${number}.seq; \
		perl scripts/randomly_subset_seqs.pl `head -n1 $<` >> subsets/fuli.3.0_subset$${number}.seq; \
	done

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Run G-PhoCS (many, many times)
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Replace SEQ_FILE_HERE in human_exome_template.ctl to make new control files
# -------------------------------------------------------------------------------------- #

# Control file for this run depends on template control file
%.iter1.ctl : human_exome_template.ctl
	@echo "# === Making control file... ================================================== #";
	$(eval filename = $*.seq)
	sed -e "s:SEQ_FILE_HERE:$${filename}:g" $< | sed -e "s/.seq.trace.log/.iter1.trace.log/" > $@
	sed -e "s/iter1/iter2/" $@ > $(subst iter1,iter2,$@)

# -------------------------------------------------------------------------------------- #
# --- Run G-PhoCS, twice for each control file
# -------------------------------------------------------------------------------------- #

# Better to call G-PhoCS from outside this Makefile, since it takes days to run
# Do so by passing the control file to the PBS script, like this:
# qsub -v CONTROL_FILE="fake_human_exomes_full.filtered.CpGmasked.iter1.ctl" pbs/call_gphocs.pbs

# 1st iteration trace files depend on seq file
%.iter1.trace.log : %.seq
	@echo "# === Calling G-PhoCS on iteration 1... ======================================= #";
	${GPHOCS}/G-PhoCS-1-2-1 $(subst trace.log,ctl,$@)

# ...as do 2nd iteration trace files
%.iter2.trace.log : %.seq
	@echo "# === Calling G-PhoCS on iteration 2... ======================================= #";
	${GPHOCS}/G-PhoCS-1-2-1 $(subst trace.log,ctl,$@)

# ====================================================================================== #
# -------------------------------------------------------------------------------------- #
# --- Analyze G-PhoCS output
# -------------------------------------------------------------------------------------- #
# ====================================================================================== #

# -------------------------------------------------------------------------------------- #
# --- Extract results from G-PhoCS output
# -------------------------------------------------------------------------------------- #

# Averaged results depend on trace files
%.iter1.results : %.iter1.trace.log
	@echo "# === Summarizing G-PhoCS output, iteration 1 ================================= #";
	${GPHOCS}/readTrace -d `echo $$(( $$(tail -n1 $< | cut -f1) / 1000))` $< > $@

%.iter2.results : %.iter2.trace.log
	@echo "# === Summarizing G-PhoCS output, iteration 2 ================================= #";
	${GPHOCS}/readTrace -d `echo $$(( $$(tail -n1 $< | cut -f1) / 1000))` $< > $@


# Scale results

# Grab numbers for paper