# set the input file to process
name=$(sed -n ${SLURM_ARRAY_TASK_ID}p namelist)
mom=Mother_SampleID
dad=Father_SampleID
ref=REF_GENOME_PATH

module load gatk
module load biokit

gatk --java-options "-Xmx15g" CombineGVCFs  \
   -R ${ref} \
   --variant ../VCFs/${mom}.g.vcf.gz \
   --variant ../VCFs/${dad}.g.vcf.gz \
   --variant ../VCFs/${name}.g.vcf.gz \
   -O Trio_${name}.g.vcf.gz \
   --tmp-dir ../../tmp

gatk --java-options "-Xmx15g" GenotypeGVCFs \
   -R ${ref} \
   -V Trio_${name}.g.vcf.gz \
   -O Trio_${name}.vcf.gz \
   --tmp-dir ../../tmp

gatk SelectVariants \
    -V Trio_${name}.vcf.gz \
    -select-type SNP \
    -O Trio_${name}_snp.vcf.gz


gatk VariantFiltration \
    -V Trio_${name}_snp.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
       -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O Trio_${name}_snp_addfilter.vcf.gz


gatk SelectVariants \
    -V Trio_${name}.vcf.gz \
    -select-type INDEL \
    -O Trio_${name}_indel.vcf.gz

bcftools view  -i 'FILTER="PASS"' Trio_${name}_snp_addfilter.vcf.gz --regions LG1,LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10,LG11,LG12:16900000-33585825,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,FIN-PYO-0_mitoSeq | bgzip -c > Trio_${name}_ats_snp_filtered.vcf.gz
tabix -p vcf Trio_${name}_ats_snp_filtered.vcf.gz

zcat Trio_${name}_ats_snp_filtered.vcf.gz | grep -v '\s\./\.' | grep -v '*' |sed 's#|#/#g' > Trio_${name}_snp_nomissing.vcf

POOHA Trio_${name}_snp_nomissing.vcf ${dad} ${mom} ${name} ../BAMs/bam_recal/${name}.recalibrated.bam > ${name}_phase.txt
