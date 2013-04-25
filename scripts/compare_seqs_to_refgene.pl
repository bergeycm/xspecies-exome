#!/usr/bin/perl

# Program to compare our random genomic hunks that may contain exons to known genes
# given a sequence alignment and a refGene file

use strict;
use warnings;

use Data::Dumper;

my $seq_file = "results/all.combined.gphocs.seq";
#my $seq_file = "tmp.gphocs.seq";
my $ref_gene = "../rhesus_db/rheMac2_refGene.txt";

open (SEQ, "<$seq_file")
	or die "ERROR: Could not open input G-PhoCS sequence file, $seq_file.\n";

# Output BED files
my $out_transcribed = $seq_file . ".transcribed.bed";
open (TRAN, ">$out_transcribed")
	or die "ERROR: Could not open output transcribed BED file, $out_transcribed.\n";
my $out_untranscribed = $seq_file . ".untranscribed.bed";
open (UNTRAN, ">$out_untranscribed")
	or die "ERROR: Could not open output untranscribed BED file, $out_untranscribed.\n";

# Output folder for FASTA alignments
mkdir ("results/loci_alignments");
mkdir ("results/loci_alignments/exons");
mkdir ("results/loci_alignments/untr");
# Clean out old FASTA files
system ("rm results/loci_alignments/exons/*.fasta");
system ("rm results/loci_alignments/untr/*.fasta");

# ----------------------------------------------------------------------------------------
# --- Loop through sequences, comparing coordinates to refGene genes.
# --- Write BED file of untranscribed and transcribed regions.
# --- Also write FASTAs alignments of these interesting regions.
# ----------------------------------------------------------------------------------------

# Store the info on the exons that overlap in hashes
my %all_exon_starts;	# $all_exon_starts{'chr1:626660-627450'} = [100,500,700]
my %all_exon_ends;
my %all_exon_offsets;

# For each sequence in G-PhoCS sequence file...
while (<SEQ>) {

	if (/^chr/) {
		print "Searching for $_\n";
		if (/^(chr\d+):(\d+)-(\d+)/) {
		
			my $chr = $1;
			my $start = $2;
			my $end = $3;
			
			my $hunk_length = $end - $start;
			print "LENGTH:\t$hunk_length\n";
	
			# See if our region overlaps any possible transcribed regions in refGene
			# Test if the possible transcript's start ($5) is <= our region's end
			# and the possible transcript's end ($6) is >= our region's start
			my $awk_cmd = "awk ' { if (\$3 == \"$chr\" && \$5 <= $end && \$6 >= $start ) ";
			$awk_cmd .= "print }' < $ref_gene";
			
	
			my $awk_result = `$awk_cmd`;
			
			#print "RG:\t" . $awk_result . "\n";
			
			if ($awk_result) {
			
				# Overlaps with a transcribed region!
				
				my @ref_gene_data = split /\t/, $awk_result;
				my ($id, $bin, $name, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, 
					$exonCount, $exonStarts, $exonEnds, $score, $name2, $cdsStartStat, 
					$cdsEndStat, $exonFrames) = @ref_gene_data;
				
				my @exon_starts_arr = split ',', $exonStarts;
				my @exon_ends_arr   = split ',', $exonEnds;
				my @exon_frames_arr = split ',', $exonFrames;

				# Translate the coordinates in the refGene interval
				# into coordinates based on our sequence
				
				my $txStart_tr  = $txStart  - $start;
				my $txEnd_tr    = $txEnd    - $start;
				my $cdsStart_tr = $cdsStart - $start;
				my $cdsEnd_tr   = $cdsEnd   - $start;
				
				my @exon_starts_arr_tr;
				foreach (@exon_starts_arr) {
					push @exon_starts_arr_tr, (1 + $_ - $start);
				}
				
				my @exon_ends_arr_tr;
				foreach (@exon_ends_arr) {
					push @exon_ends_arr_tr, ($_ - $start);
				}
				
				
			
				print STDERR "$chr:$start-$end\tMAYBE_TRANSCRIBED\n";
				
				#print "bin\t$bin\n";
				#print "txStart\t$txStart\n";
					#print "\t$txStart_tr\n";
				#print "txEnd\t$txEnd\n";
					#print "\t$txEnd_tr\n";
				#print "cdsStart\t$cdsStart\n";
					#print "\t$cdsStart_tr\n";
				#print "cdsEnd\t$cdsEnd\n";
					#print "\t$cdsEnd_tr\n";
				print "exonCount\t$exonCount\n";
				print "exonStarts\t$exonStarts\n";
					#foreach (@exon_starts_arr_tr) {
					#	print "\t$_\n";
					#}
				print "exonEnds\t$exonEnds\n";
					#foreach (@exon_ends_arr_tr) {
					#	print "\t$_\n";
					#}
				print "exonFrames\t$exonFrames\n";
				
				print "\n\n";
				
				# Loop through exons, looking for ones that overlap with our hunk
				
				# Save the info on the overlapping ones in temporary arrays
				my @tmp_exon_starts;	
				my @tmp_exon_ends;
				my @tmp_exon_offsets;

				
				for (my $i = 0; $i < scalar @exon_starts_arr_tr; $i++) {
				
					my $this_start = $exon_starts_arr_tr[$i];
					my $this_end   = $exon_ends_arr_tr[$i];
					
					my $this_exon_frame = $exon_frames_arr[$i];
					
					#print "\tEXON from $this_start to $this_end\n";
					#print "\tFrame is $this_exon_frame\n";
					
					# Does it overlap?
					# Is the exon start < our end (hunk length)?
					# Is the exon end   > our start (zero)?
					if ($this_start <= $hunk_length && $this_end >= 0) {
						#print "\t\tOVERLAP!\n";
						
						# Figure out local coordinates of overlapping region
						# by adjusting $this_start and $this_end if they
						# are out of bounds
						if ($this_start < 0) {
							# If start is changed, adjust the exonFrame as needed
							# (difference between $this_start and zero) % 3
							
							my $added_bases = 0 - $this_start;
							#print "\tAdded bases: $added_bases\n";
							$this_exon_frame += $added_bases % 3;
							#print "\tFrame is now $this_exon_frame\n";
							
							# Adjust so exonFrame is between -1 and 2
							if ($this_exon_frame > 2) {
								$this_exon_frame -= 3;
							}
							
							$this_start = 0;
						}
						if ($this_end > $hunk_length) {
							$this_end = $hunk_length;
						}
						
						print "\t\tOverlapping region: $this_start to $this_end\n";
						
						push @tmp_exon_starts,  $this_start;
						push @tmp_exon_ends,    $this_end;
						push @tmp_exon_offsets, $this_exon_frame;

					}
					
				}			

				# Add the temporary arrays of exon info to the hashes
				my $locus = "$chr:$start-$end";
				$all_exon_starts {$locus} = \@tmp_exon_starts;	
				$all_exon_ends   {$locus} = \@tmp_exon_ends;
				$all_exon_offsets{$locus} = \@tmp_exon_offsets;
												
				print TRAN "$chr\t$start\t$end\n";
			} else {
				print STDERR "$chr:$start-$end\tNOT_TRANSCRIBED\n";
				print UNTRAN "$chr\t$start\t$end\n";
				
				if (/^(chr\d+):(\d+)-(\d+) (\d+)/) {
					my $chr = $1;
					my $start = $2;
					my $end = $3;
					my $num_taxa = $4;
			
					my $locus = "$chr:$start-$end";
					my $locus_file_name = $locus;
					$locus_file_name =~ s/:/_/g;
				
					my $hunk_length = $end - $start;
					print "LENGTH:\t$hunk_length\n";
					
					# For each taxon
					for (my $t = 0; $t < $num_taxa; $t++) {
						my $taxon_seq_line = <SEQ>;
						# Grab taxon label
						my $taxon_name;
						if ($taxon_seq_line =~ /(\S+)\t/) {
							$taxon_name = $1;
						}
						$taxon_seq_line =~ s/\S+\t(\S+)\n/$1/;
					
						# Infer FASTA output file name
						my $untr_out = "results/loci_alignments/untr/";
						$untr_out .= $locus_file_name . '.fasta';
						
						open (UNTR, ">>$untr_out") 
							or die "ERROR: Could not open output FASTA file $untr_out\n";
						
						print UNTR '>' . $taxon_name . '_' . $locus_file_name . "\n";
						print UNTR $taxon_seq_line . "\n";
						
						print "\n";
						
						close UNTR;
					}
				}
			}
			
			#my $temp = <>;
		}
	}
}

close SEQ;

print "\n\n\n\n";
print "----------------------------------------------------------------------\n";
print "\n\n\n\n";

foreach (keys %all_exon_starts) {

	my @tmp = @{ $all_exon_starts{$_} };
	print "$_ :\n";
	foreach (@tmp) {
		print "\t$_\n";
	}

}

# ----------------------------------------------------------------------------------------
# --- Loop through sequences again, this time pulling out exons, adding Ns,
# --- and writing the new alignments to a NEXUS file
# ----------------------------------------------------------------------------------------

print "\n ---- SECOND LOOP ----\n\n";

open (SEQ, "<$seq_file")
	or die "ERROR: Could not open input G-PhoCS sequence file, $seq_file.\n";

while (<SEQ>) {

	# Search for the header lines
	if (/^chr/) {
		print "Searching for $_\n";
		if (/^(chr\d+):(\d+)-(\d+) (\d+)/) {
			print "mini match\n";
		
			my $chr = $1;
			my $start = $2;
			my $end = $3;
			my $num_taxa = $4;
			
			my $locus = "$chr:$start-$end";
			my $locus_file_name = $locus;
			$locus_file_name =~ s/:/_/g;
			
			print "NUM TAXA\t$num_taxa\n";
			
			my $hunk_length = $end - $start;
			print "LENGTH:\t$hunk_length\n";
			
			# Only do the following if there are exons to be pulled
			if ($all_exon_starts{$locus}) {
				my @these_exon_starts  = @{ $all_exon_starts{$locus} };
				my @these_exon_ends    = @{ $all_exon_ends{$locus} };
				my @these_exon_offsets = @{ $all_exon_offsets{$locus} };
			
				print "\nStarts:\n";
				print @these_exon_starts;
				print "\nEnds:\n";
				print @these_exon_ends;
				print "\n";
				
				# For each taxon
				for (my $t = 0; $t < $num_taxa; $t++) {
					my $taxon_seq_line = <SEQ>;
					# Grab taxon label
					my $taxon_name;
					if ($taxon_seq_line =~ /(\S+)\t/) {
						$taxon_name = $1;
					}
					$taxon_seq_line =~ s/\S+\t(\S+)\n/$1/;
					# For each exon
					for (my $e = 0; $e < scalar @these_exon_starts; $e++) {
					
						# Infer FASTA output file name
						my $exon_out = "results/loci_alignments/exons/";
						$exon_out .= $locus_file_name . '_EXON' . ($e+1) . '.fasta';
						
						open (EXON, ">>$exon_out") 
							or die "ERROR: Could not open output FASTA file $exon_out\n";
						
						print "(Taxon $t, Exon $e:)\nFull: [$taxon_seq_line]\n";
						my $taxon_exon_line = substr $taxon_seq_line, $these_exon_starts[$e], 1 + $these_exon_ends[$e] - $these_exon_starts[$e];
						print "Sub:  [$taxon_exon_line]\n\n";
						
						print "NUMBER OF Ns: $these_exon_offsets[$e]\n";
						
						print EXON '>' . $taxon_name . '_' . $locus_file_name;
						print EXON '_EXON' . ($e+1) . '_';
						print EXON '[' . ($start + $these_exon_starts[$e]) . '-';
						print EXON ($start + $these_exon_ends[$e]) . "]\n";
						print EXON 'N' x $these_exon_offsets[$e];
						print EXON $taxon_exon_line . "\n";
						
						print "\n";
						
						close EXON;
					}
				}
			}
		}		
	}
}

# Now make a file of all the untranscribed sequences, concatenated together
my $paste_cmd = "paste -d'' results/loci_alignments/untr/* ";
$paste_cmd .= "> results/loci_alignments/all_untr.fasta";
my $sed_cmd = "sed -i -e 's/^\(>[^_]*\)_.*/\1/g' results/loci_alignments/all_untr.fasta";

system ($paste_cmd);
system ($sed_cmd);

exit 0;



#945	NM_001190958	chr13	-	47244002	47268268	47244002	47268268	5	

#47244002,47244274,47264066,47265967,47268187,	
#47244072,47244339,47264324,47266027,47268268,	0	C13H2orf61	cmpl	cmpl	2,0,0,0,0,
#      70       65      258      60        81