#!/usr/bin/perl

use strict;
use warnings;

use lib '/home/cmb433/local_perl/';

use Bio::SeqIO;
use File::Temp;

use Getopt::Long;

my $verbose;
my $target_fasta;
my $liftover_path;
my $hg_to_pan_chain;
my $chimp_seqs_fasta;

my $result = GetOptions (	"target_fasta=s"		=> \$target_fasta,
							"liftover_path=s"		=> \$liftover_path,
							"hg_to_pan_chain=s"		=> \$hg_to_pan_chain,
							"chimp_seqs_fasta=s"	=> \$chimp_seqs_fasta,
							"verbose"				=> \$verbose )
	or die "ERROR: Invalid commmand line options.\n";

# BED of locations without chimp synteny
my $chimp_mask_bed = "data/panTro2Mask.bed";

my @pops = ('abt', 'hanChinese', 'kb1', 'na12891', 'na18507', 'sjk', 'venter');

my $inseq = Bio::SeqIO->new(
							-file   => "<$target_fasta",
							-format => "fasta",
							);

my $tmp_bed = File::Temp->new(
							TEMPLATE => 'tmp_hg18_XXXXX',
							DIR      => '.',
							SUFFIX   => '.bed',
							);

LOCUS_LOOP: while (my $seq = $inseq->next_seq) {
	my $region_info = $seq->display_id;
	my ($chr, $start, $end);
	if ($region_info =~ /(chr[\w\d]+):(\d+)-(\d+)/) {
		$chr = $1;
		$start = $2;
		$end = $3;
	}
	
	print STDERR "Analyzing region $chr:${start}-${end}:\n" if ($verbose);
	
	my $locus_length = $end - $start;
	
	# Figure out regions that lack chimp synteny
	print STDERR "\tFinding regions that lack chimp synteny\n" if ($verbose);
	
	my $chimp_mask_cmd = "awk '{ if ( \$1 == \"$chr\" && \$2 >= $start && \$3 <= $end ) ";
	$chimp_mask_cmd .= "print \$0 }' < $chimp_mask_bed";
	
	my $chimp_mask_lines = `$chimp_mask_cmd`;
	my @chimp_mask_data = split /\n/, $chimp_mask_lines;
	
	print STDERR "\t\t" . length(@chimp_mask_data) . " regions found.\n" if ($verbose);
	
	# Grab chimpanzee sequence
	
	# Write locus coordinates (hg18) to temporary BED file
	
	my $tmp_bed_name = $tmp_bed->filename;

	print STDERR "\tWriting hg18 coords to temp BED file [$tmp_bed_name]\n" if ($verbose);

	open TMP_BED, ">", $tmp_bed_name
		or die "ERROR: Could not open temporary BED file, $tmp_bed_name: $!";
	print TMP_BED "$chr\t$start\t$end\n";
	close TMP_BED;
	
	# Call LiftOver on these human coordinates to convert them to panTro2

	if ($verbose) {
		print STDERR "\tConverting to panTro2 coords [${tmp_bed_name}_panTro2.bed]\n";
	}
	
	my $lo_cmd = "${liftover_path}/liftOver $tmp_bed_name ";
	$lo_cmd .= "$hg_to_pan_chain ";
	$lo_cmd .= "${tmp_bed_name}_panTro2.bed ${tmp_bed_name}_panTro2.unMapped.bed";

	system($lo_cmd);
	
	# Read in chimp coordinates
	my $chimp_bed_length = `wc -l ${tmp_bed_name}_panTro2.bed | cut -f 1 -d' '`;
	if ($chimp_bed_length < 1) {
		print STDERR "Nothing in chimp bed, ${tmp_bed_name}_panTro2.bed.\n";
		next LOCUS_LOOP;
	}
	
	my $chimp_coord = `tail -n1 ${tmp_bed_name}_panTro2.bed`;
	my @chimp_coords = split /\s+/, $chimp_coord;
	
	print STDERR "\tChimp coordinates are [@chimp_coords]\n" if ($verbose);
	
	# Grab this region from the set of chimp sequences, in $chimp_seqs_fasta
	my $grep_cmd = "grep -A1 '>" . $chimp_coords[0] . ":" . $chimp_coords[1] . "-";
	$grep_cmd .= $chimp_coords[2] . "' $chimp_seqs_fasta | tail -n 1";
		
	my $chimp_seq = `$grep_cmd`;
	
	if (! defined $chimp_seq) {
		print STDERR "No chimp sequence found. Looping.\n";
		next LOCUS_LOOP;
	}
	chomp $chimp_seq;
	
	if ($locus_length != length($chimp_seq)) {
		print STDERR "Lengths of chimp and human seqs don't match. Skipping.\n";
		next LOCUS_LOOP;
	}
	
	# Remove temporary BED files
	system("rm ${tmp_bed_name} ${tmp_bed_name}_panTro2.bed ${tmp_bed_name}_panTro2.unMapped.bed");
	
	# Print locus header
	my $num_pops = scalar (@pops) + 1;
	
	print $chr . '_' . $start . '_' . $end . ' ';
	print $num_pops . ' ' . $locus_length . "\n";
	
	my $ref_seq = $seq->seq;
	
	if ($verbose && $chimp_mask_lines) {
		print STDERR "\tChimp Mask lines: [$chimp_mask_lines]\n";
	}
	
	foreach my $indiv (@pops) {

		my $snp_file = "data/snps/" . $indiv . "-snps.bed";
	
		my $awk_cmd = "awk '{ if ( \$1 == \"$chr\" && \$2 >= $start && \$3 <= $end ) ";
		$awk_cmd .= "print \$0 }' < $snp_file";

		my $snp_match_lines = `$awk_cmd`;
		
		my @snp_data = split /\n/, $snp_match_lines;
		
		my $this_pop_seq = $ref_seq;
	
		foreach my $snp_line (@snp_data) {
	
			my @snp_info = split /\s/, $snp_line;
			my $snp_location = $snp_info[1] - $start;
		
			if ($verbose) {
				print STDERR "\tSNP: $indiv has a $snp_info[4] at bp $snp_location\n";
			}
					
			substr $this_pop_seq, $snp_location, 1, $snp_info[4];
		}
		
		# Mask regions without chimp synteny
		foreach my $chimp_mask (@chimp_mask_data) {
			my @chimp_mask_info = split /\s/, $chimp_mask;
			
			my $mask_location = $chimp_mask_info[1] - $start;
			
			if ($verbose) {
				print STDERR "\tMasking out location $mask_location ";
				print STDERR "due to lack of chimp synteny.\n";
			}
				
			substr $this_pop_seq, $mask_location, 1, 'N';
		}
		
		print $indiv . "\t" . $this_pop_seq . "\n";
		
	}
	
	# Print chimp sequence
	print "chimp\t" . $chimp_seq . "\n\n";

}

exit;