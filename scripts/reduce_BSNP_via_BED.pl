#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

# ------------------------------------------------------------------------------
# --- Reduce BSNP output to FASTA of just the BPs in BED file intervals
# ------------------------------------------------------------------------------

my $in_bed = shift;
chomp $in_bed;

my $bsnp_file = shift;
chomp $bsnp_file;

my $reduced_bsnp = shift;
if (! defined ($reduced_bsnp)) {
	$reduced_bsnp = $bsnp_file . ".reduced";
}

print STDERR "Reduced BSNP output is $reduced_bsnp\n";

open (BED, "<$in_bed")
	or die "ERROR: Could not open input BED file, $in_bed.\n";

# Get sample ID from BSNP output path
my $shortname = $bsnp_file;
$shortname =~ s,results/([^\.]+)\..*,$1,g;

open (BSNP, "<$bsnp_file")
	or die "ERROR: Could not open input BSNP file, $bsnp_file.\n";

open (BSNP_LITE, ">$reduced_bsnp")
	or die "ERROR: Could not open output BSNP file, $reduced_bsnp.\n";


BED_LOOP: foreach my $bed (<BED>) {

	print STDERR "Processing BED interval: $bed";
	
	my @bed_info = split /\t/, $bed;

	my $bed_chrom = $bed_info[0];
	my $bed_start = $bed_info[1];
	my $bed_end   = $bed_info[2];
	chomp $bed_end;

	print BSNP_LITE '>' . $bed_chrom . ':' . $bed_start . '-' . $bed_end . "\n";
	
	BSNP_LOOP: while (my $bsnp = <BSNP>) {
	
		my @bsnp_info = split /\t/, $bsnp;

		if ($bsnp_info[0] ne $bed_chrom) {
			# Close and reopen BSNP file if BED intervals 
			# are still waiting to be processed
			if (eof BSNP) {
				print STDERR "Reopening BSNP\n";
				close BSNP;
				open (BSNP, "<$bsnp_file")
					or die "ERROR: Could not open input BSNP file, $bsnp.\n";
			}
			next BSNP_LOOP;
		
		# If we haven't reached the start of the target interval yet, 
		# keep going through BSNP output
		} elsif ($bsnp_info[1] < $bed_start) {
			next BSNP_LOOP;
		
		# If we've hit the end, print the basepair and go on to next interval.
		} elsif ($bsnp_info[1] == $bed_end) {
			print BSNP_LITE $bsnp_info[4];
			print BSNP_LITE "\n";
			next BED_LOOP;
			
		# Otherwise we must be in the interval. Print the basepair
		# and keep going through BSNP output
		} else {
			print BSNP_LITE $bsnp_info[4];
			next BSNP_LOOP;
		}
		
	}

}

close BSNP_LITE;
close BSNP;
close BED;

exit;