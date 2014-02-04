#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

# Macaques:
#my $stats_file = "results/tajimas_d_macaques.txt"; NO!
#my $seq_file = "results/all.combined.gphocs.seq";

# Humans:
#	Tajima's D is in 14th column
#	Fu and Li's D is in 10th column
#	Fu and Li's D* is in 11th column
#	GC percentage is in 3rd column
#my $stats_file = "../human_exome_demog/fake_human_exomes_FULL.filtered.CpGmasked.stats.txt";
#my $seq_file = "../human_exome_demog/fake_human_exomes_FULL.filtered.seq";

my $usage = "perl remove_various_outlier_seqs.pl <sequence_file.seq> ";
$usage .= "<summary_stats.txt> <'tajima'|'fuli'|'fuli_star'|'gc'> <cutoff value>";

my $seq_file = shift;
my $stats_file = shift;
my $stat_to_filter_on = shift;
my $stat_cutoff = shift;
defined ($stat_cutoff) 
	or die "ERROR: Pass correct inputs.\n$usage\n";
chomp $stat_cutoff;

if ($stat_to_filter_on ne 'tajima' && $stat_to_filter_on ne 'fuli' && $stat_to_filter_on ne 'fuli_star' && $stat_to_filter_on ne 'gc') {
	die "ERROR: Pass correct inputs.\n$usage\n";
}

# Read summary stat values into hashes

my %taj_d_all;
my %fuli_d_all;
my %fuli_d_star_all;
my %gc_all;

open(STATFILE, "<$stats_file")
	or die "ERROR: Could not open file of summary statistic data.\n";

while (<STATFILE>) {

	# chr1_6561574_6562260	-1.57597	...
	
	#	Tajima's D is in 14th column
	#	Fu and Li's D is in 10th column
	#	Fu and Li's D* is in 11th column
	#	GC percentage is in 3rd column
	
	my @summary_stat_data = split;
	my $locus = $summary_stat_data[0];
	my $taj_d_val = $summary_stat_data[13];
	my $fuli_d_val = $summary_stat_data[9];
	my $fuli_d_star_val = $summary_stat_data[10];
	my $gc_val = $summary_stat_data[2];
	
	$taj_d_all{$locus} = $taj_d_val;
	$fuli_d_all{$locus} = $fuli_d_val;
	$fuli_d_star_all{$locus} = $fuli_d_star_val;
	$gc_all{$locus} = $gc_val;

}

close STATFILE;

open(SEQ, "<$seq_file")
	or die "ERROR: Could not open sequence file.\n";

while (my $seq_line = <SEQ>) {

	if ($seq_line =~ /^(chr\S+) \d+ \d+/) {

		my $seq_locus = $1;
		
		# Get rid of ones for which Fu & Li's D calculation was impossible
		# Probably due to lack of segregating sites
		
		if (defined($fuli_d_all{$seq_locus})) {
		
			my $this_stat;
			if ($stat_to_filter_on eq 'tajima') {
				$this_stat  = $taj_d_all{$seq_locus};
			} elsif ($stat_to_filter_on eq 'fuli') {
				$this_stat  = $fuli_d_all{$seq_locus};
			} elsif ($stat_to_filter_on eq 'fuli_star') {
				$this_stat  = $fuli_d_star_all{$seq_locus};
			} else {
				$this_stat  = $gc_all{$seq_locus};
			}
			
			if (abs($this_stat) < abs($stat_cutoff)) {
				#while ($seq_line !~ /^vallender/) {
				while ($seq_line !~ /^chimp/) {
					print $seq_line;
					$seq_line = <SEQ>;
				}
				print $seq_line;
			}
		}
	}
}

exit;