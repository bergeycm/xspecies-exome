#!/usr/bin/perl

# ------------------------------------------------------------------------------
# --- Intersect all BEDs to get regions we'll include in G-PhoCS analysis
# --- Keep only autosomal regions of a certain size (500 bp or more)
# ------------------------------------------------------------------------------

use strict;
use warnings;

my @bsnp_beds = <results/*.bwa.*.bsnp.snp.out.gt4.bed>;
my $bedtools_path =  $ENV{'BEDTOOLS'};

# Intersect BEDs
my $cmd_1 = "${bedtools_path}/intersectBed ";
$cmd_1 .= "-a $bsnp_beds[0] -b $bsnp_beds[1] ";

for (my $i = 2; $i < scalar @bsnp_beds; $i++) {

	$cmd_1 .= " | ${bedtools_path}/intersectBed -a stdin ";
	$cmd_1 .= "-b $bsnp_beds[$i] ";
		
}		
		
$cmd_1 .= "> results/all.bsnp.snp.out.gt4.bed;";

print STDERR $cmd_1 . "\n";
system ($cmd_1);

# Reduce to intervals of a certain size and get rid of sex chroms

my $cmd_2 = "awk '{ if (\$3 - \$2 > 499) print \$0 }' results/all.bsnp.snp.out.gt4.bed ";
$cmd_2 .= "| grep \"chr[0-9]\" ";
$cmd_2 .= "| ${bedtools_path}/sortBed -i stdin ";
$cmd_2 .= "> results/all.bsnp.snp.out.gt4.large.bed;";

print STDERR $cmd_2 . "\n";
system ($cmd_2);

exit;

