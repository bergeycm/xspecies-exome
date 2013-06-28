#!/usr/bin/perl

use strict;
use warnings;

# First load correct version of R
# module load r/intel/2.15.0

my $gphocs_seq = shift;
chomp $gphocs_seq;

open (SEQ, "<$gphocs_seq")
	or die "ERROR: Could not open sequence file, $gphocs_seq. $!\n";

my $seq_hunk = '';

my $num_seqs = <SEQ>;
chomp $num_seqs;

my $total_sequence_count = 0;

foreach (<SEQ>) {

	if (/\S/) {
		
		my @this_line = split /\s/;
	
		if ($this_line[0] =~ /^chr/) {
		
			# Hit a new sequence. Print the old one as a temporary FASTA file.
			
			# Only do if this is a real sequence
			if ($seq_hunk =~ /^chr/) {

				### ----------------------------------------------------------------------
				# Save locus header
				my $locus_id;
				if ($seq_hunk =~ /(chr\S*) \d.*\n/) {
					$locus_id = $1;	
				}
				
				# Convert to FASTA format
				$seq_hunk =~ s/chr.*\n//g;		# Get rid of locus header line
				$seq_hunk =~ s/^/>/g;			# Add > before first sequence header
				$seq_hunk =~ s/\n(\w)/\n>$1/g;	# Add > before other sequence headers
				$seq_hunk =~ s/\t/\n/g;			# Put space between header and sequence
				
				# Print to temporary FASTA
				my @fasta_seqs = split />/, $seq_hunk;
				
				open (TMP_FASTA, ">tmp.fa")
					or die "ERROR: Could not open temporary FASTA file. $!\n";
				
				foreach my $sq (@fasta_seqs) {
					my ($header, $sequence) = split /\n/, $sq;
					if ($sequence) {
						print TMP_FASTA ">" . $header . "\n" . $sequence . "\n";
					}
				}
				
				# Calculate Tajima's D with R (pegas package)
				my $r_cmd = "Rscript --vanilla calc_tajima_d_MACAQUES.R 2> /dev/null";
				my $r_out = `$r_cmd`;
				
				my @r_info = split /\s/, $r_out;
				
				my $tajima_d_full  = $r_info[0];
								
				# Print results
				print "$locus_id\t$tajima_d_full\n";
			
				system ("rm tmp.fa");
				$total_sequence_count++;
				### ----------------------------------------------------------------------
			}
			
			$seq_hunk = $_;
			
		} else {
			if ($this_line[1] =~ /[^N\s]/g) {
				$seq_hunk .= $_;
			}
		} 
	}
}

$total_sequence_count++;

print STDERR "Total sequences output: $total_sequence_count\n";

exit;

