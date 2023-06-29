module load biokit

####    Pedigree ID
PEDIGREE="TVA_24"
####    Vcflist
VCFLIST="../VCFlist/${PEDIGREE}_F2_Triovcf_ats.list"
DNMLIST="../DNMlist/${PEDIGREE}/"
ls ${PEDIGREE}Trio*ats_filtered.vcf.gz  > ${VCFLIST}
mkdir ../DNMlist/${PEDIGREE}/
####    Extract info from population vcfs

bcftools query -f '%CHROM %POS %FILTER %REF %ALT{0} %AF{0}[ %GT][ %GQ][ %DP][ %AD{0} %AD{1}]\n' ${PEDIGREE}_snps_filtered.vcf.gz | grep "PASS" > ${DNMLIST}${PEDIGREE}_snps_PASS.txt

paste <(echo "Sample_ID") <(echo "Mendelian") <(echo "indelremove") <(echo "GQ80") <(echo "DP") <(echo "AD_0") \
<(echo "AB")  > ${DNMLIST}${PEDIGREE}_Count_ats.txt

## !!the line number should be changed if other pedigree patterns 
## 2-generation pedigree
### offspring (last generation of a pedigree)


for i in {1..10}
do

## Mendelian Violation
zcat $(sed -n ${i}p $VCFLIST) | bcftools query -l ## double check the sample order in vcfs
zcat $(sed -n ${i}p $VCFLIST) | bcftools query -f '%CHROM %POS %FILTER %REF %ALT{0} %AF{0}[ %GT][ %GQ][ %DP][ %AD{0} %AD{1}]\n' | awk '
($7 == "0|0" || $7 == "0/0") && ($8 == "0|0" || $8 == "0/0")&& ($9 == "0/1" || $9 == "0|1")' | awk '
$1!="LG12"||($1=="LG12"&&$2>=16900000)' > ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_0.txt ## 0/0 -> 0/1
zcat $(sed -n ${i}p $VCFLIST) | bcftools query -f '%CHROM %POS %FILTER %REF %ALT{0} %AF{0}[ %GT][ %GQ][ %DP][ %AD{0} %AD{1}]\n' | awk '
($7 == "1|1" || $7 == "1/1") && ($8 == "1|1" || $8 == "1/1")&& ($9 == "0/1" || $9 == "0|1")' | awk '
$1!="LG12"||($1=="LG12"&&$2>=16900000)' > ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_1.txt ## 1/1 -> 0/1

## remove variants 5bp from indels
zcat $(sed -n ${i}p $VCFLIST | sed 's/ats_snp_filtered/indel/g') | bcftools query -f '%CHROM %POS %REF %ALT{0}[ %GT]\n'  > ${DNMLIST}${PEDIGREE}_F2_0${i}_indel.txt

awk '
  FNR==NR{                                               # 0/0 -> 0/1
    y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2                 
    next                                                 
  }
                                                         
  ($1 in y){                                             
    n=split(y[$1], y1)                                   
                                                         
 for (i=1;i<=n;i++)                                  
      if (y1[i]<$2 && $2-y1[i]<=5)
        print $1,y1[i]
      else if (y1[i]>$2 && y1[i]-$2-length($3)<=5)        
        print $1,y1[i]                                  
  }
' ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_0.txt ${DNMLIST}${PEDIGREE}_F2_0${i}_indel.txt > ${DNMLIST}${PEDIGREE}_F2_0${i}_indel_SNP.txt

awk '
  FNR==NR{                                               # the same for 1/1 -> 0/1
    y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2                 
    next                                                 
  }
                                                         
  ($1 in y){                                             
    n=split(y[$1], y1)                                   
                                                         
 for (i=1;i<=n;i++)                                  
      if (y1[i]<$2 && $2-y1[i]<=5)
        print $1,y1[i]
      else if (y1[i]>$2 && y1[i]-$2-length($3)<=5)        
        print $1,y1[i]        
  }
' ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_1.txt ${DNMLIST}${PEDIGREE}_F2_0${i}_indel.txt >> ${DNMLIST}${PEDIGREE}_F2_0${i}_indel_SNP.txt

awk 'NR == FNR {                              
  k[$1, $2]
  next
  }
!(($1, $2) in k)
  ' ${DNMLIST}${PEDIGREE}_F2_0${i}_indel_SNP.txt ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_0.txt > ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt

awk 'NR == FNR {                              
  k[$1, $2]
  next
  }
!(($1, $2) in k)
  ' ${DNMLIST}${PEDIGREE}_F2_0${i}_indel_SNP.txt ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_1.txt > ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt


## individual filters: GQ[20,100], ADparentalt[0],DPtrio,AB[0.3,0.7]
awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '
    $17 <= 0 && $19 <= 0'  | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}' | awk '{ADsum = $20 + $21;
    if($21 >= 0.3*ADsum && $21 <= 0.7*ADsum){print $0}}' > ${DNMLIST}${PEDIGREE}_F2_0${i}_snp0.txt

awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '
    $16 <= 0 && $18 <= 0'  | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}' | awk '{ADsum = $20 + $21;
    if($20 >= 0.3*ADsum && $20 <= 0.7*ADsum){print $0}}' > ${DNMLIST}${PEDIGREE}_F2_0${i}_snp1.txt

wc -l ${DNMLIST}${PEDIGREE}_F2_0${i}_snp0.txt
wc -l ${DNMLIST}${PEDIGREE}_F2_0${i}_snp1.txt

## print numbers of filtered variants
paste <(echo ${PEDIGREE}_F2_0${i}) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_0.txt | wc -l) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt | wc -l) \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt | wc -l) \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}' | wc -l) \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_0.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '
    $17 <= 0 && $19 <= 0'  | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}' | wc -l) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_snp0.txt | wc -l)  >> ${DNMLIST}${PEDIGREE}_Count_ats.txt

paste <(echo ${PEDIGREE}_F2_0${i}) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_intersection_1.txt | wc -l) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt | wc -l) \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt | wc -l) \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}' | wc -l)  \
<(awk '$10 >= 80 && $11 >= 80 && $12 >= 80' ${DNMLIST}${PEDIGREE}_F2_0${i}_SNPbase_1.txt | awk '
    $13 >= 20 && $14 >= 20 && $15 >= 20 && $13 <= 100 && $14 <= 100 && $15 <= 100' | awk '
    $16 <= 0 && $18 <= 0'  | awk '{DPsum = ($13 + $14 + $15)/3;
    if($15 >= 0.5*DPsum && $15 <= 2*DPsum){print $0}}'| wc -l) \
<(cat ${DNMLIST}${PEDIGREE}_F2_0${i}_snp1.txt | wc -l) >> ${DNMLIST}${PEDIGREE}_Count_ats.txt
done

cat */*_snp0.txt | sed  's/_snp0.txt//g' | sed 's#/# #g' >  ${DNMLIST}${PEDIGREE}_DNMbase.txt
cat */*_snp1.txt | sed  's/_snp0.txt//g' | sed 's#/# #g' >>  ${DNMLIST}${PEDIGREE}_DNMbase.txt

sed 's/LG//g' ${DNMLIST}${PEDIGREE}_DNMbase.txt  | awk '
{       if($3 != l1 || $4 - l2 > 100) {
                if(lp) print l3,l4,l1,l2,l5,l6,l7
                lp = 1
        } else  lp = 0
        l1 = $3 
        l2 = $4
        l3 = $1
        l4 = $2
        l5 = $5
        l6 = $6
        l7 = $7
}
END {   if(lp) print l3,l4,l1,l2,l5,l6,l7
}' > ${DNMLIST}${PEDIGREE}_DNMbase_nocluster.txt

