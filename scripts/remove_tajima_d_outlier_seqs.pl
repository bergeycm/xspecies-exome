#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

# Macaques:
#my $tajima_fuli_d_file = "results/tajimas_d_macaques.txt";
#my $seq_file = "results/all.combined.gphocs.seq";

# Humans:
#my $tajima_fuli_d_file = "../human_exome_demog/tajima_d_human_nochimp.txt";
#my $seq_file = "../human_exome_demog/fake_human_exomes_FULL.filtered.seq";

my $usage = "perl remove_tajima_d_outlier_seqs.pl <sequence_file.seq> ";
$usage .= "<summary_stats.txt>, <'tajima'|'fuli'> <D cutoff value>";

my $seq_file = shift;
my $tajima_fuli_d_file = shift;
my $taj_or_fuli = shift;
my $tajima_fuli_cutoff = shift;
defined ($tajima_fuli_cutoff) 
	or die "ERROR: Pass correct inputs.\n$usage\n";
chomp $tajima_fuli_cutoff;

if ($taj_or_fuli ne 'tajima' && $taj_or_fuli ne 'fuli') {
	die "ERROR: Pass correct inputs.\n$usage\n";
}

# Read D values into hashes

my %taj_d_all;
my %fuli_d_all;

open(TAJFULI, "<$tajima_fuli_d_file")
	or die "ERROR: Could not open file of Tajima's or Fu & Li's D data.\n";

while (<TAJFULI>) {

	# chr1_6561574_6562260	-1.57597	...
	my @taj_fuli_data = split;
	my $locus = $taj_fuli_data[0];
	my $taj_d_val = $taj_fuli_data[1];
	my $fuli_d_val = $taj_fuli_data[2];
	
	$taj_d_all{$locus} = $taj_d_val;
	$fuli_d_all{$locus} = $fuli_d_val;

}

close TAJFULI;

open(SEQ, "<$seq_file")
	or die "ERROR: Could not open sequence file.\n";

while (my $seq_line = <SEQ>) {

	if ($seq_line =~ /^(chr\S+) \d+ \d+/) {

		my $seq_locus = $1;
		
		# Get rid of ones for which Fu & Li's D calculation was impossible
		# Probably due to lack of segregating sites
		
		if (defined($fuli_d_all{$seq_locus})) {
		
			my $this_d;
			if ($taj_or_fuli eq 'tajima') {
				$this_d  = $taj_d_all{$seq_locus};
			} else {
				$this_d  = $fuli_d_all{$seq_locus};
			}
			
			if (abs($this_d) < abs($tajima_fuli_cutoff)) {
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