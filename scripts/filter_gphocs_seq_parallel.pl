#!/usr/bin/perl

use strict;
use warnings;

# Script to mask out regions, such as CpG islands, from a G-PhoCS sequence file.
# Input is seq file and BED of regions to be masked.

# This hasn't really been tested. Use with caution!.

my $in_seq = shift;
my $in_bed = shift;
my $start_locus = shift;	# 1 is first locus
my $num_locus = shift;		# Number of loci to process
chomp $num_locus;

print STDERR "Filtering out regions from $in_bed from $in_seq\n";
print STDERR "Locus $start_locus and $num_locus loci after it.\n\n";

open (SEQ, '<' . $in_seq)
	or die "ERROR: Could not open sequence file.\n";

my $chr;
my $start;
my $end;

my $locus_count = 0;

while (my $line = <SEQ>) {

	if ($line =~ /^\d+$/) {
		print $line;
		
	} elsif ($line =~ /^chr/) {
		
		$locus_count++;
		if ($locus_count < $start_locus) {
			next;
		} elsif ($locus_count >= ($start_locus + $num_locus)) {
			exit;
		}
		
		if ($line =~ /chr(\d+)_(\d+)_(\d+)/) {
			
			$chr = $1;
			$start = $2;
			$end = $3;
		
			print STDERR "Searching for regions to mask for locus ";
			print STDERR "chr${chr}:${start}-$end\n";
			
			print $line;

			# Find BED intervals that overlap
			my $bedtools_cmd = "echo \"chr$chr\t$start\t$end\" | ";
			$bedtools_cmd .= "~/exome_macaque/bin/BEDTools-Version-2.13.4/bin/intersectBed ";
			$bedtools_cmd .= "-sorted -a $in_bed -b stdin";
			
			my @mask_regions = `$bedtools_cmd`;
			
			my @local_mask_starts;
			my @local_mask_ends;
			
			# Convert mask regions' coordinates into "local" coordinates
			# relative to start of this locus
			foreach (@mask_regions) {
				my @mask_data = split;
				$mask_data[1] -= $start;
				$mask_data[2] -= $start;
				push @local_mask_starts, $mask_data[1];
				push @local_mask_ends, $mask_data[2];
			}
			
			while ($line = <SEQ>) {
				my @seq_data = split /\s/, $line;
				
				# Mask out regions
				for (my $i = 0; $i < scalar @local_mask_starts; $i++) { 
				
					my $this_start = $local_mask_starts[$i];
					my $this_end   = $local_mask_ends[$i];
				
					my $mask_length = $this_end - $this_start;
					my $mask = 'N' x $mask_length;

					substr ($seq_data[1], $this_start, $mask_length) = $mask;
				
				}
				
				print $seq_data[0] . "\t" . $seq_data[1] . "\n";
				
				if ($line =~ /^chimp/) {
					print "\n";
					last;
				}
			}
			
		}
	}

}

exit;