#!/usr/bin/perl

use strict;
use warnings;

# ------------------------------------------------------------------------------
# --- Combine BSNP FASTAs (created with reduce_BSNP_via_BED.pl) to generate a
# --- G-PhoCS sequence input file.
# ------------------------------------------------------------------------------

# Optional suffix for naming output file. Default is "reduced"
my $suffix = shift;
if (! defined ($suffix)) {
	$suffix = "reduced";
}
chomp $suffix;

# Get all individual BSNP output files
my @bsnp_files = `ls results/*.bsnp.snp.out.gt4.$suffix`;
my $bsnp_files_str = join (' ', @bsnp_files);
$bsnp_files_str =~ s/\n//g;

# Interleave BSNP files using paste
my $paste_cmd = 'paste -d \'\n\' ';
$paste_cmd .= $bsnp_files_str;
$paste_cmd .= ' > results/all.bsnp.seqs.combined.tmp';

print STDERR "CMD: $paste_cmd\n";
system($paste_cmd);

# Output of paste is weird merged file with N identical FASTA headers then N sequences,
# where N = number of input BSNP files.

# Get short names of BSNP files for labeling sequences in G-PhoCS file
my @short_names = @bsnp_files;
for (my $n = 0; $n < scalar @short_names; $n++) {
	my @tmp = split (/[\.\/]/, $short_names[$n]);
	$short_names[$n] = $tmp[1];
}

# Get number of individuals
my $num_ind = scalar @short_names;

# Loop through combined BSNP file and write it to STDOUT in G-PhoCS sequence format
# Output should look something like this:
#	locus1 2 12
#	ind_A	GCCCACCTGTTCCATGAC
#	ind_B	GCCCACCTGTTCCATGAC
#	ind_C	GCCCACCTGTTCCATGAC

open (BSNP, "<results/all.bsnp.seqs.combined.tmp")
	or die "ERROR: Couldn't open combined BSNP file. $!\n";

# First thing to do is to find and print the number of loci
my $total_loci = `grep -c '>' $bsnp_files[0]`;
chomp $total_loci;
print $total_loci . "\n\n";

# Keep track of whether header was output
my $header_written = 0;

# Keep track of which individual's sequence we've written
my $current_ind = 0;

while (<BSNP>) {

	# Header line
	if (/^>/) {
		if (!$header_written) {
			# Remove > to get locus name
			my $locus_name = $_;
			$locus_name =~ s/>//;
			chomp $locus_name;
			print $locus_name . ' ';

			# Print number of individuals
			print $num_ind . ' ';
			
			# Infer and print locus length
			#if ($locus_name =~ /^[.+]:(\d+)-(\d+)$/) {
			if ($locus_name =~ /^.*:(\d+)-(\d+)$/) {
				my $loci_length = ($2 - $1) + 1;
				print $loci_length;
			}
			print "\n";
			
			$header_written = 1;
		}
	} else {
	
		# Sequence line
		# Print individual's label then sequence
		# then increment the individual tracker
		
		print $short_names[$current_ind] . "\t";
		$current_ind++;
		
		print $_;
		
		# Test to see if $current_ind is out of bounds
		# That means we're on the last individual's sequence
		if ($current_ind == scalar @short_names) {
			$current_ind = 0;
			$header_written = 0;
			# Print a blank line between loci
			print "\n";
		}
	}
}

# Clean up
system ("rm results/all.bsnp.seqs.combined.tmp");

exit;