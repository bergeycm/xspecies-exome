#!/usr/bin/perl

use strict;
use warnings;

my $gphocs_seq = shift;
chomp $gphocs_seq;

open (SEQ, "<$gphocs_seq")
	or die "ERROR: Could not open sequence file, $gphocs_seq. $!\n";

my $seq_hunk = '';

my $num_seqs = <SEQ>;
chomp $num_seqs;

my $total_sequence_count = 0;

# Print header
print "locus_id\tlocus_length\tgc_perc\t";
print "all_pi\tall_theta\thuman_chimp_pi\tpi_per_site\t";
print "ingroup_pi\tingroup_theta\tfu_li_d\tfu_li_d_star\t";
print "fu_li_f\tfu_li_f_star\ttajima_d_human\tsingletons\t";
print "segregating_sites\tancestral_mutations\tderived_mutations\n";
	 
foreach (<SEQ>) {

	if (/\S/) {
		
		my @this_line = split /\s/;
	
		if ($this_line[0] =~ /^chr/) {
		
			# Hit a new sequence. Print the old one as a temporary FASTA file.
			
			# Only do if this is a real sequence
			if ($seq_hunk =~ /^chr/) {

				# Save locus header
				my $locus_id;
				if ($seq_hunk =~ /(chr\S*).\d.*\n/) {
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
				
				# ========================================================================
				# --- Calculate alignment statistics with BioPerl script
				# ========================================================================
								
				my $perl_cmd = "perl scripts/compute_alignment_stats.pl ";
				my $perl_out_str = `$perl_cmd`;
				
				# Locus length
				my $locus_length;
				if ($perl_out_str =~ /locus_length:\t(.*)\n/) {
					$locus_length = $1;
				}
				
				# GC percentage
				my $gc_perc;
				if ($perl_out_str =~ /gc_perc:\t(.*)\n/) {
					$gc_perc = $1;
				}
				
				# ------------------------------------------------------------------------

				# All individuals (human+chimp) pi
				my $all_pi;
				if ($perl_out_str =~ /all_pi:\t(.*)\n/) {
					$all_pi = $1;
				}
				
				# All individuals (human+chimp) theta
				my $all_theta;
				if ($perl_out_str =~ /all_theta:\t(.*)\n/) {
					$all_theta = $1;
				}
				
				# Human-chimp dyad pi
				my $human_chimp_pi;
				if ($perl_out_str =~ /human_chimp_pi:\t(.*)\n/) {
					$human_chimp_pi = $1;
				}
				
				# pi per site
				my $pi_per_site;
				if ($perl_out_str =~ /pi_per_site:\t(.*)\n/) {
					$pi_per_site = $1;
				}
				
				# Human only pi
				my $ingroup_pi;
				if ($perl_out_str =~ /ingroup_pi:\t(.*)\n/) {
					$ingroup_pi = $1;
				}
				
				# Human only theta
				my $ingroup_theta;
				if ($perl_out_str =~ /ingroup_theta:\t(.*)\n/) {
					$ingroup_theta = $1;
				}
				
				# ------------------------------------------------------------------------				
								
				# Fu & Li's D
				my $fu_li_d;
				if ($perl_out_str =~ /Fu_and_Li_D:\t(.*)\n/) {
					$fu_li_d = $1;
				}

				# Fu & Li's D*
				my $fu_li_d_star;
				if ($perl_out_str =~ /Fu_and_Li_D\*:\t(.*)\n/) {
					$fu_li_d_star = $1;
				}	
				
				# Fu & Li's F
				my $fu_li_f;
				if ($perl_out_str =~ /Fu_and_Li_F:\t(.*)\n/) {
					$fu_li_f = $1;
				}	
				
				# Fu & Li's F*
				my $fu_li_f_star;
				if ($perl_out_str =~ /Fu_and_Li_F\*:\t(.*)\n/) {
					$fu_li_f_star = $1;
				}	
				
				# Tajima's D
				my $tajima_d_human;
				if ($perl_out_str =~ /Tajima_D:\t(.*)\n/) {
					$tajima_d_human = $1;
				}
				
				# ------------------------------------------------------------------------

				# Singletons
				my $singletons;
				if ($perl_out_str =~ /Singletons:\t(.*)\n/) {
					$singletons = $1;
				}
				
				# Segregating sites
				my $segregating_sites;
				if ($perl_out_str =~ /Segregating_sites:\t(.*)\n/) {
					$segregating_sites = $1;
				}
				
				# Ancestral mutations
				my $ancestral_mutations;
				if ($perl_out_str =~ /Ancestral_mutations:\t(.*)\n/) {
					$ancestral_mutations = $1;
				}
				
				# Derived mutations
				my $derived_mutations;
				if ($perl_out_str =~ /Derived_mutations:\t(.*)\n/) {
					$derived_mutations = $1;
				}
				
				# ------------------------------------------------------------------------
											
				# Print results
				print $locus_id;
				print "\t";
				print $locus_length;
				print "\t";
				print $gc_perc;
				print "\t";

				print $all_pi;
				print "\t";
				print $all_theta;
				print "\t";
				print $human_chimp_pi;
				print "\t";
				print $pi_per_site;
				print "\t";
				print $ingroup_pi;
				print "\t";
				print $ingroup_theta;
				print "\t";

				print $fu_li_d;
				print "\t";
				print $fu_li_d_star;
				print "\t";
				print $fu_li_f;
				print "\t";
				print $fu_li_f_star;
				print "\t";
				print $tajima_d_human;
				print "\t";

				print $singletons;
				print "\t";
				print $segregating_sites;
				print "\t";
				print $ancestral_mutations;
				print "\t";
				print $derived_mutations;
				print "\n";
				
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

