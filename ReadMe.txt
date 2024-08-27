# Reads alignment

## Process the reference and do alignments
### Index the reference

bwa index Aculy_genome.fa
samtools faidx Aculy_genome.fa
java -jar picard.jar CreateSequenceDictionary R=Aculy_genome.fa O=Aculy_genome.dict

### Align the raw reads

bwa mem -M -t 8 -R "@RG\tID:SampleID\tPL:ILLUMINA\tLB:SampleID\tSM:SampleID" Aculy_genome.fa SampleID_R1.fq.gz SampleID_R2.fq.gz -o SampleID.sam
samtools view -bS SampleID.sam -o SampleID.bam

### Sort and index of cleaned bam

java -jar picard.jar SortSam INPUT=SampleID.bam OUTPUT=SampleID_sort.bam SORT_ORDER=coordinate
java -jar picard.jar MarkDuplicates INPUT=SampleID_sort.bam OUTPUT=SampleID_dedup.bam METRICS_FILE=SampleID_dedup.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

# Variant calling
## Individual call

gatk HaplotypeCaller -R Aculy_genome.fa -I SampleID_dedup.bam -ERC GVCF -O SampleID.g.vcf

## Joint calling of all samples

gatk CombineGVCFs -R Aculy_genome.fa -V Sample1.g.vcf -V Sample2.g.vcf ...... -V Sample37.g.vcf -O allSample.g.vcf.gz
gatk GenotypeGVCFs -R Aculy_genome.fa -V allSample.g.vcf.gz  -O allSample.vcf.gz

### Separate SNPs and Indels

gatk SelectVariants -R Aculy_genome.fa -select-type SNP -V allSample.vcf.gz -O allSample.snp.raw.g.vcf.gz
gatk SelectVariants -R Aculy_genome.fa -select-type INDEL -V allSample.vcf.gz -O allSample.indel.raw.g.vcf.gz

### Hard filters to both SNPs and Indels

gatk VariantFiltration -R Aculy_genome.fa -V allSample.snp.raw.g.vcf.gz -O allSamplefilter.snp.vcf.gz -filter 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0' --filter-name 'ylffilter'
gatk VariantFiltration -R Aculy_genome.fa -V allSample.indel.raw.g.vcf.gz -O allSamplefilter.indel.vcf.gz -filter 'QD < 2.0 || FS > 100.0 || SOR > 5.0' --filter-name 'ylffilter2'

### Merge SNPs and Indels

java -jar picard.jar MergeVcfs -I allSamplefilter.snp.vcf.gz -I allSamplefilter.indel.vcf.gz -O allSample.filter.vcf.gz

### Filter the allSample vcf file
bcftools filter -S . -e 'FMT/DP<3' -O v -o allSample.vcf allSamplefilter.vcf
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 allSample.vcf \
| bcftools view -m2 -M2 -v snps -O v -o allSample.diploid.gap10.vcf
vcftools --vcf allSample.diploid.gap10.vcf --mac 3 --max-missing 0.7 --recode --out allSample.diploid.gap10.mac3.miss0.7
vcftools --vcf allSample.diploid.gap10.mac3.miss0.1.recode.vcf --missing-indv --out indv-miss