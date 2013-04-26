#!/usr/bin/perl

# ------------------------------------------------------------------------------
# --- Make FASTA of the region to be analyzed in G-PhoCS
# --- and then combine them into G-PhoCS sequence file.
# ------------------------------------------------------------------------------

# This script calls reduce_BSNP_via_BED.pl on all filtered BSNP files
# to reduce the BSNP output to a FASTA file
# then it calls bsnp_fastas_to_gphocs_seq_file.pl
# to combine this FASTAs into a G-PhoCS sequence file.

# Optionally pass "untranscribed" to only parse untranscribed regions
my $only_untr = shift;
if (!defined $only_untr) {
	$only_untr = '';
}
chomp $only_untr;

use strict;
use warnings;

my $bed_file = "results/all.bsnp.snp.out.gt4.large.bed";

# If we're aiming for untranscibed regions, BED file changes to results/all.combined.gphocs.seq.untranscribed.bed
if ($only_untr eq "untranscribed") {
	my $bed_file = "results/all.combined.gphocs.seq.untranscribed.bed";
}

my @bsnp_files = <results/*.bwa.*.bsnp.snp.out.gt4>;

# Reduce BSNP output to FASTA of just the BPs in BED file intervals

foreach my $bsnp (@bsnp_files) {

	my $perl_cmd = "perl scripts/reduce_BSNP_via_BED.pl $bed_file $bsnp";
	
	# If we're aiming for untranscibed regions, add optional output file,
	# results/*.bwa.*.passed.realn.bsnp.snp.out.gt4.untranscribed
	if ($only_untr eq "untranscribed") {
		my $untr_output = $bsnp . '.untranscribed';
		$perl_cmd .= " $untr_output";
	}
	
	$perl_cmd .= ';';
	
	print STDERR $perl_cmd . "\n";
	system ($perl_cmd);
}

# Combine these FASTAs into a G-PhoCS seq file
my $combine_cmd = "perl scripts/bsnp_fastas_to_gphocs_seq_file.pl ";

# If we're aiming for untranscibed regions, pass optional suffix, 
# so the *.untranscribed file is read, not the default *.results
if ($only_untr eq "untranscribed") {
	$combine_cmd .= 'untranscribed ';
}

$combine_cmd .= "> results/all.combined.gphocs.seq;";

system ($combine_cmd);

exit;

