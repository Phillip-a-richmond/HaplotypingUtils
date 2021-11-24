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


# To load shapeit4
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate SHAPEIT4


# Call shapeit4
shapeit4 --help




# To unload shapeit4
conda deactivate

# Load VCFTools
source /mnt/common/Precision/Miniconda3/opt/miniconda3/etc/profile.d/conda.sh
conda activate Hayden

# use VCFTools
vcftools --help



# Use PLINK
/mnt/common/Precision/PLINK/plink --help


# Use EAGLE
/mnt/common/Precision/Eagle/Eagle_v2.4.1/eagle  --help



