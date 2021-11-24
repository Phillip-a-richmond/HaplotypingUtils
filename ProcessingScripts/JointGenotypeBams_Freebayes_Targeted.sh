#!/bin/bash
#SBATCH --partition=defq
## Change to be your email address

#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage: 64 Gb of RAM for the whole job; using 16 CPUs
#SBATCH --mem=160G

## Using 8 CPUs
#SBATCH --cpus-per-task=20

## Running for a max time of 8 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

# Tools
FREEBAYES=/mnt/common/Precision/Freebayes/freebayes/bin/freebayes-1.3.4-linux-static-AMD64

# Databases
GENOME_FASTA=/mnt/scratch/Precision/Hayden/GRCh37.fa

# Data directory
FULL_BAM_DIR=/mnt/scratch/Precision/Hayden/DataDownload/Processed/20092I_0381_Libo/bam_vcf/
cd $FULL_BAM_DIR

# Step 1 - make a list of bams
BAMLIST=/mnt/scratch/Precision/Hayden/BamList.txt
rm $BAMLIST
for bamdir in $(ls ${FULL_BAM_DIR})
do
	ls $FULL_BAM_DIR${bamdir}/*bam >> $BAMLIST
done

TARGET_BED=/mnt/scratch/Precision/Hayden/ATXN3_Probes_PAR.bed

# Step 2 - Run Freebayes using the bamlist
$FREEBAYES \
	-f $GENOME_FASTA \
	-L $BAMLIST \
	--min-alternate-fraction 0.1 \
	--min-mapping-quality 1 \
	--min-base-quality 15 \
	--min-alternate-count 3 \
	--genotype-qualities \
	-t $TARGET_BED \
	-v /mnt/scratch/Precision/Hayden/ATXN3_Cohort_Freebayes_Targeted.vcf

