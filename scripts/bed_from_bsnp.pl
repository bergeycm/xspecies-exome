#!/usr/bin/perl

use strict;
use warnings;

# ------------------------------------------------------------------------------
# --- Get BED file of contiguous intervals covered in a BSNP output file -------
# ------------------------------------------------------------------------------

my $in_bsnp = shift;
chomp $in_bsnp;

open (BSNP, "<$in_bsnp")
	or die "ERROR: Could not open input BSNP file, $in_bsnp.\n";

# Ignore header line
my $header = <BSNP>;

my $contig_chrom = '';
my $contig_start_bp = -1;

my $last_chrom = '';
my $last_bp = -1;

while (my $line = <BSNP>) {

	my @info = split /\s+/, $line;
	
	# If this one is contiguous to last one
	if ($info[1] == $last_bp + 1 && $info[0] eq $last_chrom) {
		
		$last_bp = $info[1];
		$last_chrom = $info[0];

	} else {
		if ($contig_start_bp != -1) {
			print $contig_chrom . "\t" . $contig_start_bp . "\t" . $last_bp . "\n";
		}
		
		$last_bp = $info[1];
		$last_chrom = $info[0];
		
		$contig_start_bp = $info[1];
		$contig_chrom = $info[0];
	
	}
	
}

exit;