# PGTA_aneuploidy
This repository contains the code and description for the paper "Identifying genes associated with the risk of maternal aneuploidy using PGT-A ultra-low-coverage whole-genome sequencing data". The aim of the paper is to introduce a generalizable method that can be leveraged for association studies using ultra-low-coverage whole-genome sequencing data (~.01x) with an example in identifying biomarkers related to infertility. 
## Workflow

![alt text][logo]

[logo]: https://github.com/alohasiqi/PGTA_aneuploidy/blob/main/workflow.png "Analysis workflow"

### Mapping
   
```shellscript
   bwa mem GRCh38_full_analysis_set_plus_decoy_hla.fa -t 10 sample.fastq.gz > sample.sam
   samtools view -S --threads 10 -b sample.sam > sample.bam
   samtools sort sample.bam --threads 10 > sample.sorted.bam
 ```

### Genotype likelihood (gl) calculation
   
See scripts/gl_cal.sh
   
This part is based on the step 3 of GLIMPSE tutorial (https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries)

### Association test

- Ancestry inference

See scripts/ances_infer.sh

This part is based on LASER tutorial (https://genome.sph.umich.edu/wiki/LASER)

- Imputation

See scripts/impute.sh
   
This part is based on the steps 4-6 of GLIMPSE tutorial (https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries)

- Generalized linear regression (GLM) test 
   
```r
   association <- glm(data = calculated_gl, formula = rate ~ ., family = "quasibinomial")
   summary(association)
```
- eQTL Analysis through Genotype-Tissue Expression (GTEx) (https://www.gtexportal.org/home/eqtlDashboardPage) 


