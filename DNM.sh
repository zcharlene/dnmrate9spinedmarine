
## FNR
touch FNR.txt

for i in {1..10}
do
zcat $(sed -n ${i}p F2_Triovcf_ats.list) | bcftools query -l
paste <(echo F2_0${i}) \
<(zcat $(sed -n ${i}p F2_Triovcf_ats.list) | bcftools query -f '%CHROM %POS[ %GT]\n' | awk '(($3 == "0/0"||$3 == "0|0")&& ($4 == "1/1"||$4 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))||(($4 == "0/0"||$4 == "0|0")&&($3 == "1/1"||$3 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))' | wc -l) \
<(zcat $(sed -n ${i}p F2_Triovcf_ats.list) | bcftools query -f '%CHROM %POS[ %GT][ %AD{0} %AD{1}]\n' | awk '(($3 == "0/0"||$3 == "0|0")&& ($4 == "1/1"||$4 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))||(($4 == "0/0"||$4 == "0|0")&&($3 == "1/1"||$3 == "1|1")&& ($5 == "0/1"||$5 == "0|1"))' | awk '{ADsum = $10 + $11;
    if($11 >= 0.3*ADsum && $11 <= 0.7*ADsum){print $0}}' | wc -l) >> FNR.txt
done


##### annotation

### align V7 mutations to V6 by a liftover file
awk -vinverse=1 -f  liftover.awk V7.agp v7variantsID > v7variantsInContigs
awk -f  liftover.awk NSP_V6.agp v7variantsInContigs > v7variantsInV6Coordinates
#cat mutation_2.list | tr '_' ' '  > mutation.list

### shorten the gff file and label them with the V6 mutation positions
awk '
  FNR==NR{                                               # if file2...
    y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2                 # append Y1-values FS-separated to array `y`
    next                                                 # continue with next record
  }
                                                         # if file1...
  ($1 in y){                                             # col4 matches col1 of file2
    n=split(y[$1], y1)                                   # split Y1-values into array `y1`
                                                         # of length n (using FS as separator)
    for (i=1;i<=n;i++)                                   # loop over y1 values
      if ($4<=y1[i] && y1[i]<=$5)         # Y1 in range?
        print $1"_"y1[i], $0          # print result
  }
' mutation.list Maker_Filtered.gff > new_narrowed.gff



### all V6 gene info according to gene name 
## create a gene list with all info included $1:gene name $2:#cds_fractions $3:cds_length $4:similar $5:CHR $6:direction $7:gene_start $8:gene_end $9:TSS $10:TTS
awk '$3 == "gene"' Maker_Filtered.gff | grep -Po 'ID=[^;]*' | sed 's/ID=//g' > gene_v6.list

while read line || [[ $line ]]; do
paste <(echo $line) \
<(grep $line Maker_Filtered.gff | awk '$3 == "CDS"' | wc -l) \
<(grep $line Maker_Filtered.gff | awk '$3 == "CDS"' | awk '
{sum+=($5-$4+1)}END{print sum}') \
<(grep $line Maker_Filtered.gff | awk '
$3 == "gene"' | grep -Po 'Note=[^;]*' | sed 's/Note=//g' | sed 's/ /_/g') \
<(grep $line Maker_Filtered.gff | awk '
{if($3 == "gene")
print $1"\t"$7}') \
<(grep $line Maker_Filtered.gff | awk '
{if($3 == "gene" && $7 == "+")
print $4"\t"$5;
else if ($3 == "gene" && $7 == "-")
print $5"\t"$4}') \
<(grep $line Maker_Filtered.gff| awk '$3 == "CDS"'| awk '{Prod[$9]++
    min[$9]=Prod[$9]==1||min[$9]>$4?$4:min[$9]
    max[$9]=max[$9]<$5?$5:max[$9]}
    END{ for (var in Prod) print min[var]"\t"max[var]}') >> gene_info_v6.list
done < gene_v6.list

## all V6 UTR info
# $1:gene_ID $2:LG $3:direction $4:five_prime_lower $5:five_prime_upper $6:five_prime_length $7:three_prime_lower $8:three_prime_upper $9:three_prime_length
while read line || [[ $line ]]; do
paste <(grep $line Maker_Filtered.gff | awk '$3 == "gene"{print $1}') \
<(echo $line) \
<(grep $line Maker_Filtered.gff | awk '$3 == "gene"{print $7}') \
<(if [[($(grep $line Maker_Filtered.gff | grep "five_prime_UTR") > 0)]]; then
    echo $(grep $line Maker_Filtered.gff| awk '$3 == "five_prime_UTR"'| awk '{Prod[$9]++
            min[$9]=Prod[$9]==1||min[$9]>$4?$4:min[$9]
            max[$9]=max[$9]<$5?$5:max[$9]}
            END{ for (var in Prod) print min[var],max[var]}')
    else echo "0 0"
fi) \
<(if [[($(grep $line Maker_Filtered.gff | grep "five_prime_UTR") > 0)]]; then
    echo $(grep $line Maker_Filtered.gff | grep "five_prime_UTR" | awk '{sum+=($5-$4+1)}END{print sum}')
    else echo 0
fi) \
<(if [[($(grep $line Maker_Filtered.gff | grep "three_prime_UTR") > 0)]]; then
    echo $(grep $line Maker_Filtered.gff| awk '$3 == "three_prime_UTR"'| awk '{Prod[$9]++
            min[$9]=Prod[$9]==1||min[$9]>$4?$4:min[$9]
            max[$9]=max[$9]<$5?$5:max[$9]}
            END{ for (var in Prod) print min[var],max[var]}')
    else echo "0 0"
fi) \
<(if [[($(grep $line Maker_Filtered.gff | grep "three_prime_UTR") > 0)]]; then
    echo $(grep $line Maker_Filtered.gff | grep "three_prime_UTR" | awk '{sum+=($5-$4+1)}END{print sum}')
    else echo 0
fi) |  tr ' ' '\t' >> UTR_gene_postition.list
done < gene_v6.list 


## map the mutations to the genes if any
awk '
  FNR==NR{                                                                            
    y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2                                              
    next                                                                              
  }                
                                                                                      
  ($1 in y){                                                                          
    n=split(y[$1], y1)                                                                
                                                                                      
    for (i=1;i<=n;i++)                                                                
      if ($3 =="gene" && $4<=y1[i] && y1[i]<=$5 && $7 == "-")                         
        print $1"_"y1[i],$4,$5,y1[i]-$4,$5-y1[i],$5-$4+1,$7,substr($9,"ID=",16),$1,y1[i];      # print result
      else if ($3 =="gene" && $4<=y1[i] && y1[i]<=$5  && $7 == "+")
        print $1"_"y1[i],$4,$5,$5-y1[i],y1[i]-$4,$5-$4+1,$7,substr($9,"ID=",16),$1,y1[i]
  }
' mutation.list Maker_Filtered.gff | sed 's/ID=//g' > gene_mutation.list 
# $4:gene_start, $5:gene_end


## for non_gene_positions
## while read line || [[ $line ]]; do
## if [[($(grep -c $line gene_mutation.list) == 0)]]; then
##  echo $line | sed 's/_/ /g'
## fi >> non_gene.list
## done < mutation_2.list


## the length from mutations to each gene's TSS and TTS
awk 'NR == FNR {
  y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2
  next
  }

  ($5 in y){n=split(y[$5], y1)
    for (i=1;i<=n;i++)
    if ($6 == "+")
      print $5"_"y1[i],$1,sqrt(($9-y1[i])*($9-y1[i])),sqrt(($10-y1[i])*($10-y1[i]));
    else if ($6 == "-")
      print $5"_"y1[i],$1,sqrt(($10-y1[i])*($10-y1[i])),sqrt(($9-y1[i])*($9-y1[i]));
  }' mutation.list gene_info_v6.list > mutation_gene_position.list

awk '{Prod[$1]++
    min[$1]=Prod[$1]==1||min[$1]>$3?$3:min[$1]
    min2[$1]=Prod[$1]==1||min2[$1]>$4?$4:min2[$1]}
    END{ for (var in Prod) print var,min[var],min2[var]}' mutation_gene_position.list > mutation_gene_TSS_TTS.list


## UTR mutations
awk 'NR == FNR {
  k[$8]=(k[$8]=="" ? "" : k[$8] FS) $10
  next
  }
($2 in k){n=split(k[$2], k1)    
   for (i=1;i<=n;i++)
    if (k1[i]>=$4 && k1[i]<=$5)
      print $1"_"k1[i],$2,"5UTR";
      else if (k1[i]>=$7 && k1[i]<=$8)
      print $1"_"k1[i],$2,"3UTR";
      else print $1"_"k1[i],$2,"NA"
  }' gene_mutation.list UTR_gene_postition.list > UTR_mutation.list

## final file
#$1:Chr_Pos $2:type(exon,intron,non-cds) $3:repeat? $4:#cds_fractions $5:cds_length $6:similar_to 
#$7:toTSS $8:toTTS $9:gene_length $10:gene_direction $11:gene_ID $12:UTR_or_not
cat mutation.list| tr '\t' '_' > mutation_2.list 

## look seperately into the mutations on two genes
awk '{print $1}' gene_mutation.list | sort | uniq --count --repeated

while read line || [[ $line ]]; do
paste <(echo $line) \
<(if [[($(grep $line new_narrowed.gff | awk '{print $4}' | grep -c "exon") > 0)]]; then
    echo "Exon"
elif [[($(grep $line new_narrowed.gff | awk '{print $4}' | grep -c "gene") > 0) && ($(grep $line new_narrowed.gff | awk '{print $4}' | grep -c "exon") == 0)]]; then
echo "Intron"
else echo "non-gene"
fi) \
<(if [[($(grep $line new_narrowed.gff | awk '{print $3}' | grep -c "repeatmasker") >0)]]; then
    echo "Y"
else echo "N"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $(grep $line gene_mutation.list | awk '{print $8}') gene_info_v6.list | awk '{print $2}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $(grep $line gene_mutation.list | awk '{print $8}') gene_info_v6.list | awk '{print $3}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $(grep $line gene_mutation.list | awk '{print $8}') gene_info_v6.list | awk '{print $4}')
else echo "NA"
fi) \
<(if [[($(grep -c $line mutation_gene_TSS_TTS.list) >0)]]; then
    echo $(grep $line mutation_gene_TSS_TTS.list | awk '{print $2}')
else echo "NA"
fi) \
<(if [[($(grep -c $line mutation_gene_TSS_TTS.list) >0)]]; then
    echo $(grep $line mutation_gene_TSS_TTS.list | awk '{print $3}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $line gene_mutation.list | awk '{print $6}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $line gene_mutation.list | awk '{print $7}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $line gene_mutation.list | awk '{print $8}')
else echo "NA"
fi) \
<(if [[($(grep -c $line gene_mutation.list) >0)]]; then
    echo $(grep $line UTR_mutation.list | awk '{print $3}')
else echo "NA"
fi) >> mutation_info.list
done < mutation_2.list



## CpG island
awk '
  FNR==NR{                                                                            # if file2...
    y[$1]=(y[$1]=="" ? "" : y[$1] FS) $2                                              # append Y1-values FS-separated to array `y`
    next                                                                              # continue with next record
  }                
                                                                                      # if file1...
  ($1 in y){                                                                          # col1 matches col1 of file2
    n=split(y[$1], y1)                                                                # split Y1-values into array `y1`
    for (i=1;i<=n;i++)                                                                # loop over y1 values
      if ($2<=y1[i] && y1[i]<=$3)                         # Y1 in range?
        print $1"_"y1[i],"Y"    # print result
        }' mutation_v7.list cpgIsland.bed > UCSC_CGI_mutation_2.list 

