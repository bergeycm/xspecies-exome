#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

# On to macaques:
my $tajima_d_file = "results/tajimas_d_macaques.txt";
my $seq_file = "results/all.combined.gphocs.seq";

my $tajima_cutoff = shift;
defined ($tajima_cutoff) 
	or die "ERROR: Must pass Tajima's D cutoff value.\n";
chomp $tajima_cutoff;

# Read Tajima's D values into a hash

my %taj_d_all;

open(TAJ, "<$tajima_d_file")
	or die "ERROR: Could not open file of Tajima's D data.\n";

while (<TAJ>) {

	# chr1_6561574_6562260	-1.57597
	my @taj_data = split;
	my $locus = $taj_data[0];
	my $taj_d_val_all = $taj_data[1];
	
	$taj_d_all{$locus} = $taj_d_val_all;

}

close TAJ;

open(SEQ, "<$seq_file")
	or die "ERROR: Could not open sequence file.\n";

while (my $seq_line = <SEQ>) {

	if ($seq_line =~ /^(chr\S+) \d+ \d+/) {

		my $seq_locus = $1;
		
		if (exists($taj_d_all{$seq_locus})) {
		
			my $d_all = $taj_d_all{$seq_locus};
			
			if (abs($d_all) < abs($tajima_cutoff)) {		# d3
							
				while ($seq_line !~ /^vallender/) {
					print $seq_line;
					$seq_line = <SEQ>;
				}
				print $seq_line;
			}
		}
	}
}

exit;