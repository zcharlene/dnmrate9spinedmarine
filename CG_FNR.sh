## Callable genome sizes & number of heterozygotes & FNR


while IFS=" " read -r PED ID DP; do
paste <(echo $ID) \
<(samtools depth $ID.CGI.bam  | awk -v DPmean=$DP '{
    if($3 >= 0.5*DPmean && $3 <= 2*DPmean){print $0}
    }' | awk '
    $1 ~ /^LG/&&($1!="LG12"||($1=="LG12"&&$2>=16900000))
    ' | wc -l ) \
<(samtools depth $ID.recalibrated.bam  | awk -v DPmean=$DP '{
    if($3 >= 0.5*DPmean && $3 <= 2*DPmean){print $0}
    }' | awk '
    $1 ~ /^LG/&&($1!="LG12"||($1=="LG12"&&$2>=16900000))
    ' | wc -l ) \
<(bcftools view ../../VCFs/$ID.g.vcf.gz --regions LG12:1-16900000 | bcftools query -f '%CHROM %POS[ %GT]\n' | awk '
    $3 == "0|1" || $3=="0/1"' | wc -l) \
<(zcat ../../$PED/Trio_${ID}_ats_snp_filtered.vcf.gz | bcftools query -f '%CHROM %POS[ %GT]\n' | awk '(($3 == "0/0"||$3 == "0|0")&& ($4 == "1/1"||$4 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))||(($4 == "0/0"||$4 == "0|0")&&($3 == "1/1"||$3 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))' | wc -l) \
<(zcat ../../$PED/Trio_${ID}_ats_snp_filtered.vcf.gz | bcftools query -f '%CHROM %POS[ %GT][ %AD{0} %AD{1}]\n' | awk '(($3 == "0/0"||$3 == "0|0")&& ($4 == "1/1"||$4 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))||(($4 == "0/0"||$4 == "0|0")&&($3 == "1/1"||$3 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))' | awk '{ADsum = $10 + $11;
    if($11 >= 0.3*ADsum && $11 <= 0.7*ADsum){print $0}}' | wc -l) >> CG_FNR.txt
done < DP.txt

