#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

# ------------------------------------------------------------------------------
# --- Convert G-PhoCS sequence file into a concatenated NEXUS file (so an 
# --- NJ tree can be inferred in PAUP).
# ------------------------------------------------------------------------------

my $gphocs_input = shift;
chomp $gphocs_input;

# Loop through once to get loci count and nucleotide count

open (GPHOCS, "<$gphocs_input")
	or die "ERROR: Couldn't open G-PhoCS sequence file. $!\n";

my $num_taxa;
my $num_bp;

while (<GPHOCS>) {
	
	# Test for header line
	if (/^chr.+:(\d+)-(\d+) (\d+) \d+$/) {
		$num_taxa = $3;
		$num_bp += ($2 - $1) + 1;
	}
}

close GPHOCS;

# Print NEXUS header
print "#NEXUS\n";
print "Begin data;\n";
print "Dimensions ntax=" . $num_taxa . " nchar=" . $num_bp . ";\n";
print "Format datatype=dna symbols=\"ACGTRYSWKMBDHVN\" missing=? gap=- interleave;\n";
print "Matrix\n";

# Loop through again to output charset

open (GPHOCS, "<$gphocs_input")
	or die "ERROR: Couldn't open G-PhoCS sequence file. $!\n";

while (<GPHOCS>) {

	if (/^[\w\d]+\s+\w+$/) {
		print $_;
	} else {
		print "\n";
	}	
}

close GPHOCS;

# Print matrix footer
print ";\n";
print "End;\n";

# Print PAUP block
print "begin paup;\n";
print "nj;\n";
#print "savetrees file=results/all.combined.gphocs.tre format=altnex brlens=yes replace;\n";
print "savetrees format=altnex brlens=yes replace;\n";
print "quit WARNTSAVE=No;\n";
print "END;\n";

exit;