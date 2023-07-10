# This script is used to infer ancestry of target samples using 1000 genomes reference panel

# prepare 1000 genomes reference panel and only export biallelic SNPs with minor allele frequency above 10% to a new VCF file.

N=48
pids=""
for i in {0..21}
do
((i=i%N)); ((i++==0)) && wait
	bcftools view \
	ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
	--exclude-types indels,mnps,ref,bnd,other \
	--min-alleles 2 \
	--max-alleles 2 \
	--min-af 0.1:minor \
	--phased \
	--exclude 'AN!=2*N_SAMPLES' \
	--output-file chr${i}_1kg_maf10.vcf.gz \
	--output-type z &
pids="$pids $!"
done
wait $pids

# merge all the VCFs of the 1000 genome project.

rm input_vcf_list.txt
for i in {1..22}
do
  ls chr${i}_1kg_maf10.vcf.gz >> input_vcf_list.txt
done

bcftools concat \
  --threads $N \
  --file-list input_vcf_list.txt \
  -o 1kg_maf10.vcf.gz \
  -O z

rm chr*.vcf.gz
fi

if [ -f "1kg.geno" ]; then
    echo "1kg.site exists."
    echo "Skipping step C."
else
echo "1kg.site does not exist."

# convert the 1000 genomes dataset to required input format

laser/vcf2geno/vcf2geno \
--inVcf 1kg_maf10.vcf.gz \
--out 1kg
fi

if [ -f "1kg.bed" ]; then
    echo "1kg.bed exists."
    echo "Skipping step D."
else
echo "1kg.bed does not exist."

# generate a bed file
cat 1kg.site | awk '{if (NR > 1) {print "chr"$1, $2-1, $2;}}' > 1kg.bed
fi

# convert bam to pileup

pids=""
COUNTER=0
while IFS= read -r bam_file || [ -n "${bam_file}" ]
do
((COUNTER=COUNTER%N)); ((COUNTER++==0)) && wait $pids &&  pids=""
SAMPID=$(basename "${bam_file::-4}")
${SAMTOOLS_PATH}/samtools mpileup -q 30 -Q 20 -f ${HG38_PATH}/hg38_primary.fa -l 1kg.bed "${bam_file}" > "./pileup/${SAMPID}.pileup" &
pids="${pids} $!"
done < "bam_list.txt"
wait $pids

echo Step F

ls -A1 pileup/*.pileup > pileup_list.txt

# convert pileup to seq
Python pileup2seq.py \
  -f ${HG38_PATH}/hg38_primary.fa \
  -m 1kg.site \
  -o target_samples \
  $(cat pileup_list.txt)
echo Step Gq

# use COORD_FILE to infer ancestry of the target samples

${LASER_PATH}/laser \
  -c laser.RefPC.coord \
  -g 1kg.geno \
  -s target_samples.seq \
  -fmt 0 \
  -k 4 \
  -minc 0.01 \
  -seed 0 \
  -nt 48 \
  -r 5
