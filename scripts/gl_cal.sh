# This script is used to calculate genotype likelihoods (GLs)

GLs need to be computed at all target individuals and all variant sites present in the reference panel of haplotypes used for the imputation. This can be done from sequencing data using BCFtools.

# 1. Using the 1000 Genomes Project reference panel from EBI 1000 genomes ftp site and performing basic QC. Using chromosome 22 as an example, 

# download chromosome 22 data 
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz{,.tbi}

# perform QC
CHR=22
bcftools norm -m -any CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz -Ou --threads 4 | bcftools view -m 2 -M 2 -v snps --threads 4 -Ob -o reference_panel/1000GP.chr22.bcf
bcftools index -f reference_panel/1000GP.chr22.bcf --threads 4

# 2. Extracting variable positions in the reference panel

bcftools view -G -m 2 -M 2 -v snps reference_panel/1000GP.chr22.bcf -Oz -o reference_panel/1000GP.chr22.sites.vcf.gz
bcftools index -f reference_panel/1000GP.chr22.sites.vcf.gz

# convert to the TSV format and index the file
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' reference_panel/1000GP.chr22.sites.vcf.gz | bgzip -c > reference_panel/1000GP.chr22.sites.tsv.gz
tabix -s1 -b2 -e2 reference_panel/1000GP.chr22.sites.tsv.gz

# 3. Computing GLs for a single individual at specific positions

BAM=sample1.chr22.sorted.bam
VCF=reference_panel/1000GP.chr22.sites.vcf.gz
TSV=reference_panel/1000GP.chr22.sites.tsv.gz
REFGEN=reference-genome/hs38DH.chr22.fa.gz # the reference genome version used by the 1000 Genomes Project
OUT=sample_vcf/sample1.chr22.vcf.gz
bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r chr22 ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}
