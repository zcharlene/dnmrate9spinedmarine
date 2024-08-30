## set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p namelist)

module load biokit
module load gatk

## read mapping; fixmate; sorted; markduplicate; flagstat; indexing
bwa mem -t4 -M -R "@RG\tID:{} tSM:{} tPL:BGI tLB:{} tPU:1" ../../ninespine_v7/NSPV7.GCA_902500615.3.fasta ../Raw/${name}_1.fq.gz ../Raw/${name}_2.fq.gz | samtools view -F4 -hb -o ${name}.fq.bam
chmod o+r ${name}.fq.bam
samtools fixmate -@ 8 -m ${name}.fq.bam ${name}.fq.fm.bam
chmod o+r ${name}.fq.fm.bam
samtools sort -O bam -o ${name}.fq.fm.sorted.bam ${name}.fq.fm.bam
chmod o+r ${name}.fq.fm.sorted.bam
samtools index ${name}.fq.fm.sorted.bam ${name}.fq.fm.sorted.bai
chmod o+r ${name}.fq.fm.sorted.bai
samtools markdup ${name}.fq.fm.sorted.bam ${name}.fq.fm.sorted.md.bam
chmod o+r ${name}.fq.fm.sorted.md.bam
samtools index ${name}.fq.fm.sorted.md.bam ${name}.fq.fm.sorted.md.bai
chmod o+r ${name}.fq.fm.sorted.md.bai
samtools flagstat ${name}.fq.fm.sorted.md.bam > ${name}.fq.fm.sorted.md.bam.flagstat

## rename groups
picard AddOrReplaceReadGroups \
I= ${name}.fq.fm.sorted.md.bam \
O= ${name}.RG.bam \
ID=${name} \
SM=${name} \
LB=${name} \
PL=DNBseq \
PU=BGI \
VALIDATION_STRINGENCY=SILENT

chmod ${name}.RG.bam

## Recalibration and filtering
gatk BaseRecalibrator \
-I ${name}.RG.bam \
-R ../../ninespine_v7/NSPV7.GCA_902500615.3.fasta \
--known-sites ../../ninespine_v7/snps.vcf.gz \
--known-sites ../../ninespine_v7/indels.vcf.gz \
-O ${name}.table

chmod o+r ${name}.table

gatk ApplyBQSR \
-R ../../ninespine_v7/NSPV7.GCA_902500615.3.fasta \
-I ${name}.RG.bam \
--bqsr-recal-file ${name}.table \
--read-filter GoodCigarReadFilter \
--read-filter NotDuplicateReadFilter \
--read-filter PassesVendorQualityCheckReadFilter \
--read-filter MappingQualityReadFilter \
--read-filter MappingQualityAvailableReadFilter \
--read-filter PrimaryLineReadFilter \
--read-filter MappedReadFilter \
-O ${name}.recalibrated.bam


chmod o+r ${name}.recalibrated.bam
chmod o+r ${name}.recalibrated.bai

touch coverage.txt

## Flagstats and coverage
paste <(echo ${name}.recalibrated.bam | sed 's/.recalibrated.bam//g') \
<(samtools depth ${name}.recalibrated.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}') >> coverage.txt

## CGI bams
samtools view -b -h -L ../../../ninespine_v7/CGI.bed ${name}.recalibrated.bam > ${name}.CGI.bam
chmod o+r ${name}.CGI.bam
samtools index ${name}.CGI.bam ${name}.CGI.bai
chmod o+r ${name}.CGI.bai


## HC
gatk --java-options "-Xmx8g" HaplotypeCaller \
         -R ../../ninespine_v7/NSPV7.GCA_902500615.3.fasta \
         -I ${name}.recalibrated.bam \
         -O ../VCFs/${name}.g.vcf.gz \
         -ERC GVCF \
         -bamout ${name}.bamout.bam \
         --tmp-dir ../../tmp/ \
         -A AlleleFraction \
         -A DepthPerSampleHC \
         -A Coverage \
         -A BaseQuality \
         -A AS_BaseQualityRankSumTest \
         -A MappingQuality \
         -A FisherStrand
