# This script is used to split the genome into chunks and impute

# 1. Chunking a chromosome

bin/GLIMPSE_chunk --input reference_panel/1000GP.chr22.sites.vcf.gz --region chr22 --window-size 2000000 --buffer-size 200000 --output chunks.chr22.txt

# 2. Impute and phase a whole chromosome

VCF=sample_vcf/sample1.chr22.vcf.gz
REF=reference_panel/1000GP.chr22.bcf
MAP=maps/genetic_maps.b38/chr22.b38.gmap.gz
while IFS="" read -r LINE || [ -n "$LINE" ];
do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    OUT=GLIMPSE_imputed/sample${ID}.chr22.bcf
    bin/GLIMPSE_phase --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
    bcftools index -f ${OUT}
done < chunks.chr22.txt

# 3. Ligate multiple chunks together

LST=GLIMPSE_ligated/list.chr22.txt
ls GLIMPSE_imputed/sample1.chr22.imputed.*.bcf > ${LST}
OUT=GLIMPSE_ligated/sample1.chr22.merged.bcf
bin/GLIMPSE_ligate --input ${LST} --output $OUT
bcftools index -f ${OUT}
