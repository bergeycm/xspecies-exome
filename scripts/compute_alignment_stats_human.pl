#!/usr/bin/perl

use strict;
use warnings;

use lib '/home/cmb433/local_perl/';

use Bio::PopGen::Utilities;
use Bio::PopGen::Statistics;
use Bio::AlignIO;

my $fasta_in = shift;
chomp $fasta_in;

# ========================================================================================
# --- Read in sequences
# ========================================================================================

my $in = Bio::AlignIO->new(	-file   => $fasta_in,
							-format => 'fasta');

my $aln = $in->next_aln;

my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln,
													-include_monomorphic => 1,
													-site_model => 'all');
#my $happop = $pop->haploid_population;

my @all_inds = $pop->get_Individuals();
my @names = $pop->get_marker_names;

my $locus_length =  scalar @names;

# ========================================================================================
# --- Make three arrays: one of a human and chimp...
# ========================================================================================

my @pi_dyad;

foreach (@all_inds) {
	if ($_->unique_id eq 'kb1') {
		push (@pi_dyad, $_);
	} elsif ($_->unique_id eq 'chimp') {
		push (@pi_dyad, $_);
	}	
}

# ========================================================================================
# --- ...one of all humans, and one of just chimp
# ========================================================================================

my @ingroup;
my @outgroup;

foreach (@all_inds) {
	if ($_->unique_id ne 'chimp') {
		push (@ingroup, $_);
	} else {
		push (@outgroup, $_);
	}	
}

# ========================================================================================
# --- Calculate GC percentage
# ========================================================================================

my $gc_count = 0;
my $len = 0;

foreach my $seq ($aln->each_seq) {

    my @bps = split('', $seq->seq);
    $len += scalar @bps;
    
	foreach my $base (@bps) {
	
		$gc_count++ if $base =~ /[GCS]/i;

		# Don't count ambiguous codes
		$len-- if $base =~ /[RYKMN]/i;
		$len-- if $base =~ /[BV]/i;
		$len-- if $base =~ /[DH]/i;
		$len-- if $base eq '-';
	}	
	
}

my $gc_perc = $gc_count / $len;

# ========================================================================================
# --- Calculate and output statistics
# ========================================================================================

print "locus_length:\t" . $locus_length . "\n";

print "gc_perc:\t" . $gc_perc . "\n";

# --- 

my $all_pi = Bio::PopGen::Statistics->pi($pop);
print "all_pi:\t" . $all_pi . "\n";

my $all_theta = Bio::PopGen::Statistics->theta($pop);
print "all_theta:\t" . $all_theta . "\n";

my $human_chimp_pi = Bio::PopGen::Statistics->pi(\@pi_dyad);
print "human_chimp_pi:\t" . $human_chimp_pi . "\n";

my $human_chimp_pi_per_site = Bio::PopGen::Statistics->pi(\@pi_dyad, $locus_length);
print "pi_per_site:\t" . $human_chimp_pi_per_site . "\n";

my $ingroup_pi = Bio::PopGen::Statistics->pi(\@ingroup);
print "ingroup_pi:\t" . $ingroup_pi . "\n";

my $ingroup_theta = Bio::PopGen::Statistics->theta(\@ingroup);
print "ingroup_theta:\t" . $ingroup_theta . "\n";


# ---

my $FuLi_D;
eval {
	$FuLi_D = Bio::PopGen::Statistics->fu_and_li_D(\@ingroup,\@outgroup);
	1;
} or warn $@;
$FuLi_D = 'NA' if (! $FuLi_D);
print "Fu_and_Li_D:\t" . $FuLi_D . "\n";

my $FuLi_D_star;
eval {
	$FuLi_D_star = Bio::PopGen::Statistics->fu_and_li_D_star(\@ingroup);
	1;
} or warn $@;
$FuLi_D_star = 'NA' if (! $FuLi_D_star);
print "Fu_and_Li_D*:\t" . $FuLi_D_star . "\n";

my $fu_and_li_F;
eval {
	$fu_and_li_F = Bio::PopGen::Statistics->fu_and_li_F(\@ingroup, \@outgroup);
	1;
} or warn $@;
$fu_and_li_F = 'NA' if (! $fu_and_li_F);
print "Fu_and_Li_F:\t" . $fu_and_li_F . "\n";

my $fu_and_li_F_star;
eval {
	$fu_and_li_F_star = Bio::PopGen::Statistics->fu_and_li_F_star(\@ingroup);
	1;
} or warn $@;
$fu_and_li_F_star = 'NA' if (! $fu_and_li_F_star);
print "Fu_and_Li_F*:\t" . $fu_and_li_F_star . "\n";

my $tajima_D = Bio::PopGen::Statistics->tajima_D(\@ingroup);
$tajima_D = 'NA' if (! $tajima_D);
print "Tajima_D:\t" . $tajima_D . "\n";

# ---

my $singletons = Bio::PopGen::Statistics->singleton_count(\@ingroup);
print "Singletons:\t" . $singletons . "\n";

my $segsites = Bio::PopGen::Statistics->segregating_sites_count(\@ingroup);
print "Segregating_sites:\t" . $segsites . "\n";

my @ext = Bio::PopGen::Statistics->derived_mutations(\@ingroup,\@outgroup);
$ext[0] = 'NA' if (! $ext[0]);
$ext[1] = 'NA' if (! $ext[1]);

print "Ancestral_mutations:\t" . $ext[0] . "\n";
print "Derived_mutations:\t" . $ext[1] . "\n";

exit;