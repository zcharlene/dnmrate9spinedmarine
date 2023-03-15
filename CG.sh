## Callable genome size & number of heterozygotes


while IFS=" " read -r ID DP; do
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
    $3 == "0|1" || $3=="0/1"' | wc -l) >> CG.txt
done < DP.txt

