# -------------------------------------------------------------------------------------- #
# --- Configuration makefile of user-editable variables 
# -------------------------------------------------------------------------------------- #

# All paths must be absolute or relative to the xspecies-exome top directory

# -------------------------------------------------------------------------------------- #
# --- Paths to input files
# -------------------------------------------------------------------------------------- #

# Individual ID (used to name files)
IND_ID=george

# Paths to input reads files
# Must be in FASTQ format
#READ1=./data/2012-11-28/macIll100.read1.fastq
#READ2=./data/2012-11-28/macIll100.read2.fastq

READ1=./data/2012-12-03/macIll.read1.fastq
READ2=./data/2012-12-03/macIll.read2.fastq

# Paths to genomes files
# Must be in FASTA format
HUMAN_GENOME_FA=./genomes/hg18/hg18.fa
SECOND_GENOME_FA=./genomes/rheMac2/rheMac2.fa
CHIMP_GENOME_FA=./genomes/panTro2/panTro2.fa

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

# Path to LiftOver Chain file for other genome-to-human genome. Available from:
# http://hgdownload.cse.ucsc.edu/downloads.html
TO_HUMAN_CHAINFILE=/home/cmb433/exome_macaque/bin/liftover/rheMac2ToHg18.over.chain

# Path to LiftOver Chain file for human-to-chimp genome. Available from:
# http://hgdownload.cse.ucsc.edu/downloads.html
TO_CHIMP_CHAINFILE=/home/cmb433/exome_macaque/bin/liftover/hg18ToPanTro2.over.chain

# Paths for ANNOVAR annotation
ANNOVAR_BUILDVER=rheMac2
ANNOVAR_DB_PATH=../rhesus_db/

# Paths to G-PhoCS control files, full dataset and only untranscribed
GPHOCS_CTL_FULL=macaque_exome_full_dataset.ctl
GPHOCS_CTL_UNTR=macaque_exome_untranscribed.ctl
# And for filtered datasets
GPHOCS_CTL_1_5=macaque.d1.5.ctl
GPHOCS_CTL_2_0=macaque.d2.ctl
GPHOCS_CTL_2_5=macaque.d2.5.ctl
GPHOCS_CTL_3_0=macaque.d3.ctl
GPHOCS_CTL_5_0=macaque.d5.ctl

# Path to RefGene file
REFGENE=../rhesus_db/rheMac2_refGene.txt

# -------------------------------------------------------------------------------------- #
# --- Paths to external programs
# -------------------------------------------------------------------------------------- #

FASTQC=/home/cmb433/exome_macaque/bin/FastQC
FASTX=/home/cmb433/exome_macaque/bin/fastx
BWA=/home/cmb433/exome_macaque/bin/bwa-0.6.2
SAMTOOLS=/home/cmb433/exome_macaque/bin/samtools
OLD_SAMTOOLS=/home/cmb433/exome_macaque/bin/samtools-0.1.16
BEDTOOLS=/home/cmb433/exome_macaque/bin/BEDTools-Version-2.13.4/bin
LIFTOVER=/home/cmb433/exome_macaque/bin/liftover
PICARD=/home/cmb433/exome_macaque/bin/picard-tools-1.77
BAMTOOLS=/home/cmb433/exome_macaque/bin/bamtools/bin
GATK=/home/cmb433/exome_macaque/bin/GATK
BCFTOOLS=/home/cmb433/exome_macaque/bin/samtools/bcftools
VCFTOOLS=/home/cmb433/exome_macaque/bin/vcftools_0.1.9/bin
PSMC=/home/cmb433/exome_macaque/bin/psmc
BSNP=/home/cmb433/exome_macaque/bin/BSNP/bin
ANNOVAR=/home/cmb433/exome_macaque/bin/annovar/
PLINK=/home/cmb433/exome_macaque/bin/plink-1.07-x86_64/
PAUP=/scratch/disotell/bin/
GPHOCS=/home/cmb433/exome_macaque/bin/G-PhoCS/bin/

# -------------------------------------------------------------------------------------- #
# --- Parameters for external programs
# -------------------------------------------------------------------------------------- #

#BWA_ALN_PARAM=-e 63 -i 15 -L -l 31 -t 8 -I 
BWA_ALN_PARAM=-t 8 
FASTX_PARAM=-q 13 -p 80 -v 
SNP_MIN_COV=5
SNP_MAX_COV=1000