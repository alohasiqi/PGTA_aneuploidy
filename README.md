# PGTA_aneuploidy
This repository contains the code and description for the paper "Identifying genes associated with the risk of maternal aneuploidy using PGT-A ultra-low-coverage whole-genome sequencing data". The aim of the paper is to introduce a generalizable method that can be leveraged for association studies using ultra-low-coverage whole-genome sequencing data (~.01x) with an example in identifying biomarkers related to infertility. 
## Workflow
Steps 1-3 are scripts in BASH, step 4 is in R
1. Mapping
   '''bash
   bwa mem GRCh38_full_analysis_set_plus_decoy_hla.fa -t 10 sample.fastq.gz > sample.sam
   samtools view -S --threads 10 -b sample.sam > sample.bam
   samtools sort sample.bam --threads 10 > sample.sorted.bam
   '''
3. Genotype likelihood (gl) calculation
This part is based on the step 3 of GLIMPSE tutorial (https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries)
3.1 Ancestry inference
This part is based on LASER tutorial (https://genome.sph.umich.edu/wiki/LASER)
3.2 Imputation
This part is based on the steps 4-6 of GLIMPSE tutorial (https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html#run_preliminaries)
5. Association test
   '''R
   glm(data = calculated_gl, formula = rate ~ ., family = "quasibinomial")
   '''

![alt text][logo]

[logo]: https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Analysis workflow"
