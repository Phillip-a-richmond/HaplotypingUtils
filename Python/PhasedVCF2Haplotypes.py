import sys,argparse

def GetOptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-V","--vcf",help="Input phased VCF", required=True)
    parser.add_argument("-B","--bedpe",help="Input SURVIVOR BEDPE file")
    parser.add_argument("-O","--outprefix",help="Output prefix for haplotype file, and sample haplotype label file", required=True)
    args = parser.parse_args()

    invcf = open(args.vcf,'r')
    outHaplo = open("%s_haplotypes.txt"%args.outprefix,'w')
    outHaploGenoSplit = open("%s_haplotypes_genotypeSplit.csv"%args.outprefix,'w')
    outSample = open("%s_samples.txt"%args.outprefix,'w')
    outSampleGenoSplit = open("%s_samples_genotypeSplit.csv"%args.outprefix,'w')

    return invcf,outHaplo,outSample,outSampleGenoSplit,outHaploGenoSplit

# This function will parse the VCF, and produce a dictionary that has {Sample}[Haplotype1,Haplotype2]
# It will also produce vectors for Positions, Refs, and Alts, which will be used in the header of the split-genotype file
def ParseVCFMakeHaploDict(invcf):
        
   # Initialize sample haplotypes dict, where we'll store each sample as key, and the 2 haplotype strings as values
   # the haplotype strings will get added to line-by-line
    SampleHaplotypes = {}
    # We'll also keep a SampleColumns dict which labels the column that a sample belongs to as we iterate with i
    SampleColumns = {}
    # Keep a record of each of the positions for printing out the split-genotype file
    Positions = []
    # Keep the Ref and alt alleles
    Refs = []
    Alts = []

        
    for line in invcf:
        # skip the VCF header lines
        if line[0:2] == '##':
            continue
        # keep the VCF header line starting with #CHROM, that has all the samples
        # We'll need to strip out the following: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
        if line[0] == '#':
            headercols = line.strip('\n').split('\t')
            #print(headercols)
            # starting at headercols[9]. since headercols[8] corresponds to FORMAT
            for i in range(9,len(headercols)):
                # add each to the dict
                sampleID = headercols[i]
                # initiate the SampleHaplotypes with an empty vector for each haplotype string
                SampleHaplotypes[sampleID] = ['','']
                SampleColumns[i] = sampleID
            continue
            #print(SampleHaplotypes)
            #print(SampleColumns)


        # Now we'll parse through each variant, which has the format of: 
        # 14	92524334	.	A	G	948081	PASS	NS=2504;AA=a|||;VT=SNP;DPB=57090;EPPR=129.339;GTI=2;MQMR=59.99;NUMALT=1;ODDS=0.90209;PAIREDR=0.987766;PQR=0;PRO=0;QR=843682;RO=23050;RPPR=27931.6;SRF=13996;SRP=2303.86;SRR=9054;DP=73071;AF=0.60373;EAS_AF=0.3224;AMR_AF=0.147;AFR_AF=0.4274;EUR_AF=0.2237;SAS_AF=0.365;AB=0.485281;ABP=77.2912;AO=34019;CIGAR=1X;DPRA=1.33631;EPP=175.319;LEN=1;MEANALT=1.04612;MQM=59.9883;PAIRED=0.988565;PAO=0;PQA=0;QA=1245550;RPL=4455;RPP=40246.1;RPR=29564;RUN=1;SAF=20860;SAP=3788.54;SAR=13159;TYPE=snp;technology.ILLUMINA=1;AN=5856;AC=2084	GT:GQ:DP:AD:RO:QR:AO:QA:GL	0|0:.:.:.:.:.:.:.:.	1|1:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	1|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	1|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	1|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	1|0:.:.:.:.:.:.:.:.	1|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	1|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	1|1:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|0:.:.:.:.:.:.:.:.	0|1:.:.:.:.:.:.:.:.

        # We need to store the A and G, since 0 corresponds to A, and 1 corresponds to G
        variantcols = line.strip('\n').split('\t')
        ref = variantcols[3]
        Refs.append(ref)
        alt = variantcols[4]
        Alts.append(alt)

        # Here I'm going to put the positions into a vector, and then I'll print this vector out as a header line for the split-geno file
        pos = variantcols[1]
        Positions.append(pos)

        # Now we'll loop through columns, exchange the iterator for the sample ID using SampleColumns, 
        # and then add each genotype to haplotypeA and haplotypeB for the SampleHaplotypes dict
        
        for i in range(9,len(variantcols)):
            # get the sample name
            sampleID = SampleColumns[i]
            #print(sampleID)
            # get the genotype by splitting this: 0|0:.:.:.:.:.:.:.:.
            GT1 = variantcols[i].split(':')[0].split('|')[0]
            GT2 = variantcols[i].split(':')[0].split('|')[1]
        
        # translate (goddamn it's been ahwile since I've coded lol)
            if GT1 == '0':
                genotype1 = ref
            elif GT1 == '1':
                genotype1 = alt
            if GT2 == '0':
                genotype2 = ref
            elif GT2 == '1':
                genotype2 = alt

            # add genotype1 and genotype2 to the haplotype strings inside of our SampleHaplotypes dict
            SampleHaplotypes[sampleID][0] += genotype1
            SampleHaplotypes[sampleID][1] += genotype2
    return SampleHaplotypes,Positions,Refs,Alts
    #print(SampleHaplotypes)
    #print(SampleHaplotypes['HG00096'])
            

# This function will take the output dictionary above, which is in format:
# SampleHaplotypes[SampleID][Haplotype1String,Haplotype2String]

# And what we want are 2 output files: OutHaplo, OutSample, and OutSampleGenoSplit.csv
# OutHaplo: Format of haplotypeID (just a #), HaplotypeString, and TimesObserved, tab-separated
# OutSample: Every sample gets 2 lines: SampleID, haplotypeID, haplotypeString
# OutSampleGenoSplit.csv: Same as OutSample, but each genotype is split into columns
def HaplotypeDictToOutput(SampleHaplotypes,OutHaplo,OutSample,OutSampleGenoSplit,OutHaploGenoSplit,Positions,Refs,Alts):

    # Create header on sample-based split-genotype file
    # Write the variant positions as a header to the genotype file
    OutSampleGenoSplit.write(",Position,")
    for i in range(len(Positions)):
        OutSampleGenoSplit.write("%s,"%Positions[i])
    OutSampleGenoSplit.write("\n")

    # Write the variant ref allele as a header to the genotype file
    OutSampleGenoSplit.write(",RefAllele,")
    for i in range(len(Refs)):
        OutSampleGenoSplit.write("%s,"%Refs[i])
    OutSampleGenoSplit.write("\n")
    
    # Write the variant alt allele as a header to the genotype file
    OutSampleGenoSplit.write(",AltAllele,")
    for i in range(len(Alts)):
        OutSampleGenoSplit.write("%s,"%Alts[i])
    OutSampleGenoSplit.write("\n\n")
    
    # Create header on haplotype-based split-genotype file
    # Write the variant positions as a header to the genotype file
    OutHaploGenoSplit.write(",Position,")
    for i in range(len(Positions)):
        OutHaploGenoSplit.write("%s,"%Positions[i])
    OutHaploGenoSplit.write("\n")

    # Write the variant ref allele as a header to the genotype file
    OutHaploGenoSplit.write(",RefAllele,")
    for i in range(len(Refs)):
        OutHaploGenoSplit.write("%s,"%Refs[i])
    OutHaploGenoSplit.write("\n")
    
    # Write the variant alt allele as a header to the genotype file
    OutHaploGenoSplit.write(",AltAllele,")
    for i in range(len(Alts)):
        OutHaploGenoSplit.write("%s,"%Alts[i])
    OutHaploGenoSplit.write("\n\n")
    

    # start new dict which will be: Haplotypes[HaplotypeString]=[haploIDnumber,TimesObserved]
    Haplotypes = {}

    # iterate through SampleHaplotypes, keeping an iterator to number haplotypes as they are observed
    i = 1
    for sample in SampleHaplotypes:
        print(sample)
        hap1 = SampleHaplotypes[sample][0]
        hap2 = SampleHaplotypes[sample][1]

        # Add to count if it exists
        if hap1 in Haplotypes:
            Haplotypes[hap1][1]+=1
        # if it doesn't exist, create a new entry for this haplotype in Haplotypes{} and increase the iterator to label as a new haplotype
        else:
            Haplotypes[hap1] = [i, 1]
            # increase iterator
            i+=1

        # repeat for hap2
        if hap2 in Haplotypes:
            Haplotypes[hap2][1]+=1
        else:
            Haplotypes[hap2] = [i, 1]
            # increase iterator
            i+=1


        # Now that we have an identifier for each haplotype, we can write to our sample file within this loop
        # format here is 2 lines per sample:
        # Sample    HapID1   HapString1
        # Sample    HapID2   HapString2
        OutSample.write("%s\t%d\t%s\n"%(sample,Haplotypes[hap1][0],hap1))
        OutSample.write("%s\t%d\t%s\n"%(sample,Haplotypes[hap2][0],hap2))

        # Write out to the Split Genotype file
        OutSampleGenoSplit.write("%s,%s,"%(sample,Haplotypes[hap1][0]))
        for gt in hap1:
            OutSampleGenoSplit.write("%s,"%gt)
        OutSampleGenoSplit.write("\n")

        OutSampleGenoSplit.write("%s,%s,"%(sample,Haplotypes[hap2][0]))
        for gt in hap2:
            OutSampleGenoSplit.write("%s,"%gt)
        OutSampleGenoSplit.write("\n")



    # Once the loop for all samples finishes, we can print out our haplotypes to OutHaplo
    for haplotype in Haplotypes:
    # print the haplotypes
        OutHaplo.write("%d\t%s\t%d\n"%(Haplotypes[haplotype][0],haplotype,Haplotypes[haplotype][1]))

    # print the haplotypes but split by genotype
        OutHaploGenoSplit.write("%d,%d,"%(Haplotypes[haplotype][0],Haplotypes[haplotype][1]))
        for gt in haplotype:
            OutHaploGenoSplit.write("%s,"%gt)
        OutHaploGenoSplit.write("\n")
        


def Main():
    Invcf,OutHaplo,OutSample,OutSampleGenoSplit,OutHaploGenoSplit = GetOptions()
    SampleHaplotypes,Positions,Refs,Alts = ParseVCFMakeHaploDict(Invcf)
    HaplotypeDictToOutput(SampleHaplotypes,OutHaplo,OutSample,OutSampleGenoSplit,OutHaploGenoSplit,Positions,Refs,Alts) 




if __name__=="__main__":
	Main()


