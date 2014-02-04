#!/usr/bin/perl

use strict;
use warnings;

#my $subset_n = 3197;	# Number of loci in full exome
#my $subset_n = 3157;	# Number of loci with filter Fu & Li's |D| < 2.0
#my $subset_n = 3102;	# Number of loci with filter Fu & Li's |D| < 2.5
#my $subset_n = 3095;	# Number of loci with filter Fu & Li's |D| < 3.0

my $subset_n = shift;
chomp $subset_n;

my $locus_list_file = "gronau_locus_list.txt";
my $seq_file = "neutralLoci-7genomes.txt";

# Read locus names into an array

my @loci_all;

open(LOCI, "<$locus_list_file")
	or die "ERROR: Could not open file of loci in sequence file.\n";

while (<LOCI>) {

	my $locus = $_;
	chomp $locus;
	
	push @loci_all, $locus;

}

close LOCI;

# ========================================================================================
# Take a random subset of the sequences
# ========================================================================================

my %locus_in_subset;

use List::Util qw/shuffle/;  
my @sample = (shuffle(@loci_all))[0..$subset_n-1];
foreach (@sample) {
	$locus_in_subset{$_} = 1;
}

# ========================================================================================

open(SEQ, "<$seq_file")
	or die "ERROR: Could not open sequence file.\n";

while (my $seq_line = <SEQ>) {

	if ($seq_line =~ /^(chr\S+)\s\d+\s\d+/) {

		my $seq_locus = $1;
		
		if (exists($locus_in_subset{$seq_locus})) {
		
			while ($seq_line !~ /^chimp/) {
				print $seq_line;
				$seq_line = <SEQ>;
			}
			print $seq_line;
		}
	}
}


exit;