#!/usr/bin/perl

use strict;
use warnings;

my $gphocs_seq = shift;
chomp $gphocs_seq;

open (SEQ, "<$gphocs_seq")
	or die "ERROR: Could not open sequence file, $gphocs_seq. $!\n";

my $hunk_to_print = '';
my $full_hunk = '';

my $num_seqs = <SEQ>;
chomp $num_seqs;

my $total_sequence_count = 0;

foreach (<SEQ>) {

	if (/\S/) {
		
		my @this_line = split /\s/;
	
		if ($this_line[0] =~ /^chr/) {
		
			my $perc_diff;
		
			if ($hunk_to_print =~ /\n.+/) {
			
				# Don't print if chimp-human distance is too long
				# indicating spurious alignment or funky chimp seq
				
				if ($hunk_to_print =~ /chimp\s(\w+)\n/) {
					my $chimp_seq = $1;
					
					if ($hunk_to_print =~ /kb1\s(\w+)\n/) {
						my $kb1_seq = $1;
					
						my @mism_pos;
						for my $i (0 .. length($chimp_seq)) {
							my $chimp_base = substr($chimp_seq,$i,1);
							my $kb1_base   = substr($kb1_seq,$i,1);

							if ($chimp_base ne $kb1_base) {
								push @mism_pos, $i;
							}
						}

						print STDERR "Num. mismatch:\t";
						print STDERR scalar @mism_pos;
						print STDERR "\n";
						
						my $seq_len = length($chimp_seq);
						$perc_diff = (scalar @mism_pos) / $seq_len;
						
						print STDERR "Perc. mismatch:\t$perc_diff\n";
					}
				}

				if ($perc_diff < 0.20) {				
					print $full_hunk;
					$total_sequence_count++;
				}
			}
			# Change the number of individuals
			s/ 8 / 7 /g;
			$hunk_to_print = $_;
			$full_hunk = $_;
		} else {
			if ($this_line[0] ne "abt") {
				s/-/N/g;
				if ($this_line[1] =~ /[^N\s]/g) {
					$hunk_to_print .= $_;
				}
				$full_hunk .= $_;
			}
		} 
	}
}

# Print last sequence
if ($hunk_to_print =~ /\n.+/) {
	print $full_hunk;
	$total_sequence_count++;
}

print STDERR "Total sequences output: $total_sequence_count\n";

exit;