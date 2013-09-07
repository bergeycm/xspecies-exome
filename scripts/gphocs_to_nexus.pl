#!/usr/bin/perl

# Script to convert the G-PhoCS input file into a series of NEXUS files,
# one per locus, for input into *BEAST

my $gphocs_in = "results/all.bsnp.snp.out.gt4.gphocs.seq.tmp";

open (GPHOCS, "<$gphocs_in") 
	or die "ERROR: Could not open G-PhoCS file for reading [$gphocs_in]\n";

open (OUT, ">/dev/null");

foreach (<GPHOCS>) {

	if (/^chr/) {
		my @info = split;
		
		# Print NEXUS footer
		print OUT ";\nEND;\n";
		
		open (OUT, '>nexus_files/' . $info[0] . ".nex")
			or die "ERROR: Could not open NEXUS file for writing " . $info[0] . ".nex\n";
		
		# Print NEXUS header
		print OUT "#NEXUS\n";
		print OUT "[" . $info[0] . "]\n";
		print OUT "BEGIN DATA;\n";
		print OUT "DIMENSIONS NTAX =" . $info[1] . " NCHAR=" . $info[2] . ";\n";
		print OUT "FORMAT DATATYPE = DNA GAP = - MISSING = ?;\n";
		print OUT "MATRIX\n";
	} else {
		print OUT $_;
	}
}

# Print NEXUS footer for last NEXUS file
print OUT ";\nEND;\n";

close (GPHOCS);

exit;

