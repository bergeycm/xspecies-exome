setwd("~/xspecies-exome/final_trace_files/filtered")

files = list.files(".") 

stats_files = txtfiles <- files[grep(glob2rx("*.stats.txt"), files)] 

estimate.mu.per.site = function(stats_file) {

	this.stats = read.table(stats_file, header=TRUE)
	
	avg.mu = (sum(this.stats$human_chimp_pi) / sum(this.stats$locus_length)) / 6500000

	return (c(stats_file, avg.mu))

}

all.mus = lapply(stats_files, estimate.mu.per.site)
