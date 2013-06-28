
require(ape)
library(pegas)

my.dna <- read.dna("tmp.fa", format = "fasta")

full.taj = tajima.test(my.dna)

full.D = full.taj$D

cat(full.D, "\n")
