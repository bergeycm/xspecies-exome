#!/usr/bin/env Rscript

# Download RefGene file from:
# http://www.openbioinformatics.org/annovar/download/hg18_refGene.txt.gz

refgene = read.table("hg18_refGene.txt")

get.codon.from.gene = function (rg, output = c("first", "second", "third", "unknown")) {

	# Grab chromosome
	chr    = rg[3]
	strand = rg[4]
	
	# Grab exon info
	exonStarts = as.numeric(unlist(strsplit(as.character(rg[10]), "[,]")))
	exonEnds   = as.numeric(unlist(strsplit(as.character(rg[11]), "[,]")))
	exonFrames = as.numeric(unlist(strsplit(as.character(rg[16]), "[,]")))

	# Change zeros to 3s for easier computation
	exonFrames[exonFrames == 0] = 3
	
	# Loop through exons...	
	for (x in 1:length(exonFrames)) {
		
		thisStart = exonStarts[x]
		thisEnd   = exonEnds[x]
		thisFrame = exonFrames[x]
		
		# Output all bases in exon if user wants "unknown" codons			
		if (thisFrame == -1) {
			
			if (output == "unknown") {
				# Print BED line covering all bases in exon
				cat(paste(chr, thisStart, thisEnd, sep="\t"), sep="\n")
			}
		
		# Otherwise, output "first", "second", or "third" codons depending on user input
		} else {
			
			# Computation depends on strand
			if (strand == "+") {
				firsts = seq(from=thisStart + (3-thisFrame), to=thisEnd, by=3)
				seconds = firsts + 1
				seconds = seconds[seconds < thisEnd]
				thirds = firsts + 2
				thirds = thirds[thirds < thisEnd]
			} else {
				firsts = seq(from=thisEnd - (3-thisFrame), to=thisStart, by=-3)
				seconds = firsts - 1
				seconds = seconds[seconds > thisStart]
				thirds = firsts - 2
				thirds = thirds[thirds > thisStart]
			}
			
			# Output coordinates in BED format
			if (output == "first") {
				cat(paste(
 					paste(chr, firsts, sep="\t"), 
 					(firsts+1), sep="\t"), sep="\n")
			} else if (output == "second") {
				cat(paste(
 					paste(chr, seconds, sep="\t"), 
 					(seconds+1), sep="\t"), sep="\n")
			} else if (output == "third") {
				cat(paste(
 					paste(chr, thirds, sep="\t"), 
 					(thirds+1), sep="\t"), sep="\n")
 			}
		}
	}
}

# Run function to generate BED files for all unknown codons...
sink(file="hg18_refGene_unknown.bed")
	tmp = apply(refgene, 1, function(x) get.codon.from.gene(x, "unknown"))
sink()

# ...first codons...
sink(file="hg18_refGene_codon1.bed")
	tmp = apply(refgene, 1, function(x) get.codon.from.gene(x, "first"))
sink()

# ...and second codons....
sink(file="hg18_refGene_codon2.bed")
	tmp = apply(refgene, 1, function(x) get.codon.from.gene(x, "second"))
sink()

# ...and third codons.
sink(file="hg18_refGene_codon3.bed")
	tmp = apply(refgene, 1, function(x) get.codon.from.gene(x, "third"))
sink()
