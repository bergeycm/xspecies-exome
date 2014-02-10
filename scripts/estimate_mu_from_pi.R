#!/usr/bin/env Rscript 

setwd("filtered_seqs")

files = list.files(".") 

stats_files = txtfiles <- files[grep(glob2rx("*.stats.txt"), files)] 
stats_files = c(stats_files, "../fake_human_exomes_full.filtered.CpGmasked.stats.txt")

estimate.mu.per.site = function(stats_file) {

	this.stats = read.table(stats_file, header=TRUE)
	
	avg.mu = (sum(this.stats$human_chimp_pi) / sum(this.stats$locus_length)) / 6500000

	return (avg.mu)

}

all.mus = sapply(stats_files, estimate.mu.per.site)

# Rename full analysis
stats_files[length(stats_files)] = "full"

short.names = gsub('fake_human_exomes_full.filtered.CpGmasked.', '', stats_files)
short.names = gsub('.stats.txt', '', short.names)

print(data.frame(analysis=short.names, mu=as.vector(all.mus)))
