# -------------------------------------------------------------------------------------- #
# --- Configuration makefile of user-editable variables 
# -------------------------------------------------------------------------------------- #

# All paths must be absolute or relative to the xspecies-exome top directory

# -------------------------------------------------------------------------------------- #
# --- Paths to input files
# -------------------------------------------------------------------------------------- #

# Individual ID (used to name files)
IND_ID=test_monkey

# Paths to input reads files
# Must be in FASTQ format
#READ1=./data/2012-11-28/macIll100.read1.fastq
#READ2=./data/2012-11-28/macIll100.read2.fastq

READ1=./data/2012-12-03/macIll.read1.fastq
READ2=./data/2012-12-03/macIll.read2.fastq

# Paths to genomes files
# Must be in FASTA format
HUMAN_GENOME_FA=./genomes/hg19/hg19.fa
SECOND_GENOME_FA=./genomes/rheMac2/rheMac2.fa

# Common name of secondary genome (used to name files)
SECOND_GENOME_NAME=rhesus

# Path to BED file of exome capture kit targets (relative to human genome)
#Should merge this two, before pipeline starts:
#hg_ill_targets=/scratch/cmb433/macEx/targets/v2/BED/029368_D_BED_20111101.bed
#hg_nimble_targets=/scratch/cmb433/macEx/targets/SeqCap_EZ_Exome_v2.targets.bed
TARGETS=./targets/029368_D_BED_20111101.bed

# Path to BED file of CCDS, human protein coding regions (relative to human genome)
CCDS=./targets/ccdsGene.hg19.4apr12.bed

# Path to LiftOver Chain file for human-to-other genome. Available from:
# http://hgdownload.cse.ucsc.edu/downloads.html
CHAINFILE=/home/cmb433/exome_macaque/bin/liftover/hg19ToRheMac2.over.chain

# -------------------------------------------------------------------------------------- #
# --- Paths to external programs
# -------------------------------------------------------------------------------------- #

FASTQC=/home/cmb433/exome_macaque/bin/FastQC
FASTX=/home/cmb433/exome_macaque/bin/fastx
BWA=/home/cmb433/exome_macaque/bin/bwa-0.6.2
SAMTOOLS=/home/cmb433/exome_macaque/bin/samtools
BEDTOOLS=/home/cmb433/exome_macaque/bin/BEDTools-Version-2.13.4/bin
LIFTOVER=/home/cmb433/exome_macaque/bin/liftover
PICARD=/home/cmb433/exome_macaque/bin/picard-tools-1.77
BAMTOOLS=/home/cmb433/exome_macaque/bin/bamtools/bin
GATK=/home/cmb433/exome_macaque/bin/GATK
BCFTOOLS=/home/cmb433/exome_macaque/bin/samtools/bcftools
VCFTOOLS=/home/cmb433/exome_macaque/bin/vcftools_0.1.9/bin
PSMC=/home/cmb433/exome_macaque/bin/psmc

# -------------------------------------------------------------------------------------- #
# --- Parameters for external programs
# -------------------------------------------------------------------------------------- #

#BWA_ALN_PARAM=-e 63 -i 15 -L -l 31 -t 8 -I 
BWA_ALN_PARAM=-t 8 
FASTX_PARAM=-q 13 -p 80 -v 
SNP_MIN_COV=5
SNP_MAX_COV=1000