# Goal is to take the freebayes VCFs and merge them into a multisample VCF, that Galen can then use for PLINK/KING/whatever

# Step 0 - Load environment, for bcftools/tabix/bgzip
ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
source $ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
conda activate $ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

# Step 1 - Made the directory AllVCFs, and copied all VCFs from the sample set into this directory.
mkdir -p AllVCFs/
cp /mnt/scratch/Precision/Hayden/DataDownload/Processed/20092I_0381_Libo/bam_vcf/Sample_20092P0*/*var.vcf.gz AllVCFs/ 

# Step 2 - Index all VCFs in the directory
cd AllVCFs
for file in *gz ;
do
	tabix $file
done

# Step 3 - Run bcftools merge with the -0 option (uncalled sites get 0/0)

# BCFTools merge
# The -i - fixes the error because summing the DP by default is throwing an error
# Found the fix here:
# https://www.biostars.org/p/431656/
bcftools merge -0 \
	-O z \
	-i - \
	-o ATXN3_Freebayes_GRCh37.merged.vcf.gz \
	*vcf.gz
tabix ATXN3_Freebayes_GRCh37.merged.vcf.gz

