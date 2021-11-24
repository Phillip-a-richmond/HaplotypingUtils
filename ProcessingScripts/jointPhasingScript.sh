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
multi_VCF=$HAYDEN/ATXN3_Cohort_Freebayes_Targeted.vcf
PHASE=$HAYDEN/Phasing
ALLVCF=$HAYDEN/AllVCFs
QC=$HAYDEN/qcVcfs
multi_VCF=$HAYDEN/ATXN3_Cohort_Freebayes_Targeted.vcf

# Step 0 - Load environment, for bcftools/tabix/bgzip
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Step 4 - Reformat VCFs and perform QC for phasing
#Remove chr from vcfs
# Sort based on postion
cat $multi_VCF | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' |\
awk '{gsub(/^chr/,""); print}' > $QC/ATXN3_Cohort_Freebayes_Targeted_chr.vcf

cd $QC
bgzip $QC/ATXN3_Cohort_Freebayes_Targeted_chr.vcf 
tabix -p vcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr.vcf.gz

conda deactivate

# Load VCFTools
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate Hayden

# Consider removing low call rate individuals and SNPs
# Split multi-allelic sites (https://alkesgroup.broadinstitute.org/Eagle/)

# Get call rate stats by individual
vcftools --gzvcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr.vcf.gz \
        --max-missing 0.95 --missing-indv \
#        --out $QC/ATXN3_Cohort_Freebayes_Targeted_stats

# Exclude individuals with less than 95% call rate
# 24 of 448 excluded
awk '$5<0.05 {print $1}' $QC/ATXN3_Cohort_Freebayes_Targeted_stats.imiss \
	> $QC/ATXN3_Cohort_Freebayes_Targeted_stat.ind.keep

head -1 $QC/ATXN3_Cohort_Freebayes_Targeted_stats.imiss \
	 > $QC/ATXN3_Cohort_Freebayes_Targeted_stat.ind.remove

awk '$5>0.05' $QC/ATXN3_Cohort_Freebayes_Targeted_stats.imiss |\
	sort -n -r -k5 | grep -v INDV \
        >> $QC/ATXN3_Cohort_Freebayes_Targeted_stat.ind.remove

vcftools --gzvcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr.vcf.gz \
	--mac 1 --recode --recode-INFO-all \
	--keep $QC/ATXN3_Cohort_Freebayes_Targeted_stat.ind.keep \
	--out $QC/ATXN3_Cohort_Freebayes_Targeted_chr_ind

# Remove rare variants not on chr 14
# Remove low call rate snps variants
vcftools --vcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr_ind.recode.vcf \
	--recode --recode-INFO-all \
	--chr 14 --out $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind

# Remove low call rate snps variants
vcftools --vcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind.recode.vcf \
	--max-missing 0.95 --recode --recode-INFO-all \
	--out $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc

conda deactivate

ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment
bgzip $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.vcf 
tabix -p vcf $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.vcf.gz

cd $PHASE
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Run eagle
EAGLE=/mnt/common/Precision/Eagle/Eagle_v2.4.1/eagle
TB=/mnt/common/Precision/Eagle/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz

$EAGLE --vcfTarget $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.vcf.gz \
        --allowRefAltSwap --geneticMapFile=$TB --vcfOutFormat=v \
        --vcfRef=$PHASE/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        --outPrefix=$PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased


conda deactivate

ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Merge phased dataset with KGP
bgzip $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf
tabix -p vcf $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf.gz

mkdir -p $PHASE/kgp_intersect
# Extract common variants in atxn3 from kgp
bcftools isec $PHASE/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
        $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf.gz \
        -p $PHASE/kgp_intersect -n =2 -w 1

bgzip $PHASE/kgp_intersect/0000.vcf
tabix -p vcf $PHASE/kgp_intersect/0000.vcf.gz

bcftools merge $PHASE/kgp_intersect/0000.vcf.gz \
        $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf.gz > \
	$PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.vcf

bgzip $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.vcf
tabix -p vcf $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.vcf.gz

# Make tab variant file
# Use bcftools query instead of vcf-to-tab
bcftools query -Hf '%CHROM\t%POS\t%REF[\t%TGT]\n' $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf.gz  > $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.tab

bcftools query -Hf '%CHROM\t%POS\t%REF[\t%TGT]\n' $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.vcf.gz > $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.tab 

bcftools query -Hf '%CHROM\t%POS\t%REF[\t%TGT]\n' $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.vcf.gz > $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.tab


# Copy relevant files to one folder for haplotyping
mkdir -p $HAYDEN/forHaplotyping
cp $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.vcf.gz $HAYDEN/forHaplotyping
cp $PHASE/kgp_intersect/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased_kgp.tab $HAYDEN/forHaplotyping
cp $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.vcf.gz $HAYDEN/forHaplotyping
cp $QC/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc.recode.tab $HAYDEN/forHaplotyping
cp $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.vcf.gz $HAYDEN/forHaplotyping
cp $PHASE/ATXN3_Cohort_Freebayes_Targeted_chr14_ind_qc_phased.tab $HAYDEN/forHaplotyping
cp $QC/ATXN3_Cohort_Freebayes_Targeted_stat.ind.remove $HAYDEN/forHaplotyping 
gunzip $HAYDEN/forHaplotyping/*gz

