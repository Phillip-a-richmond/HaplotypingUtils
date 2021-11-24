#!/bin/bash
#SBATCH --partition=defq
## Change to be your email address

#SBATCH --mail-user=gwright@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage: 64 Gb of RAM for the whole job; using 16 CPUs
#SBATCH --mem=64G

## Using 8 CPUs
#SBATCH --cpus-per-task=8

## Running for a max time of 8 hours
#SBATCH --time=8:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

HAYDEN=/mnt/scratch/Precision/Hayden
multi_VCF=$HAYDEN/ATXN3_Freebayes_GRCh37.merged.vcf
PHASE=$HAYDEN/Phasing
ALLVCF=$HAYDEN/AllVCFs
QC=$HAYDEN/qcVcfs
multi_VCF=$QC/ATXN3_Freebayes_GRCh37.merged.vcf.gz


# Step 0 - Load environment, for bcftools/tabix/bgzip
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Step 1 - Made the directory AllVCFs, and copied all VCFs from the sample set into this directory.
#mkdir -p $ALLVCF
#cp /mnt/scratch/Precision/Hayden/DataDownload/Processed/20092I_0381_Libo/bam_vcf/Sample_20092P0*/*var.vcf.gz $ALLVCF

# Step 2 - Index all VCFs in the directory
#cd $ALLVCF
#for file in *gz ;
#do
#        tabix $file
#done

# Step 3 - Run bcftools merge with the -0 option (uncalled sites get 0/0)

# BCFTools merge
# The -i - fixes the error because summing the DP by default is throwing an error
# Found the fix here:
# https://www.biostars.org/p/431656/
#mkdir -p $QC
#bcftools merge -0 \
#        -O z \
#        -i - \
#        -o ATXN3_Freebayes_GRCh37.merged.vcf.gz \
#        *vcf.gz
#tabix ATXN3_Freebayes_GRCh37.merged.vcf.gz
#mv ATXN3_Freebayes_GRCh37* $QC


# Step 4 - Reformat VCFs and perform QC for phasing
#Remove chr from vcfs
#zcat $multi_VCF | awk '{gsub(/^chr/,""); print}' > $PHASE/ATXN3_Freebayes_GRCh37.merged_chr.vcf
cd $PHASE
#bgzip $PHASE/ATXN3_Freebayes_GRCh37.merged_chr.vcf && tabix -p vcf $PHASE/ATXN3_Freebayes_GRCh37.merged_chr.vcf.gz

#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Run eagle
EAGLE=/mnt/common/Precision/Eagle/Eagle_v2.4.1/eagle
TB=/mnt/common/Precision/Eagle/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz

$EAGLE --vcfTarget $PHASE/ATXN3_Freebayes_GRCh37.merged_chr.vcf.gz \
        --allowRefAltSwap --geneticMapFile=$TB --vcfOutFormat=v \
        --vcfRef=$PHASE/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --outPrefix=$PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased

conda deactivate

# Load VCFTools
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate Hayden

# Consider removing low call rate individuals and SNPs
# Split multi-allelic sites (https://alkesgroup.broadinstitute.org/Eagle/)

# use VCFTools
# Remove rare variants on chr 14
#vcftools --vcf $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased.vcf \
#	--maf 0.01 --recode --recode-INFO-all \
#	--chr 14 --out $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased_maf_0.01

# Extract haplotypes with
PLINK=/mnt/common/Precision/PLINK/plink


conda deactivate

ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Export haplotypes
#bcftools convert --hapsample --vcf-ids $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased_maf_0.01.recode.vcf\
#       	-o $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased_maf_0.01_haps

# Make tab variant file
# Use bcftools query instead of vcf-to-tab
#bcftools query -Hf '%CHROM\t%POS\t%REF[\t%TGT]\n' $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased_maf_0.01.recode.vcf  > $PHASE/ATXN3_Freebayes_GRCh37.merged_chr_phased_maf_0.01.tab
