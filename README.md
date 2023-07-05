# dnmrate9spinedmarine
We estimated pedigree-based germ line mutation rates in two marine populations in nine-spined sticklebacks (Pungitius pungitius) with whole genome sequencing data of 106 parents-offspring trios.


Array_all.sh: Array job for mapping and haplotype calling
CG_FNR.sh: Calculate callable genome size, double check the sample sex by counting het in sex determination region (X), and calculate false negative rates.
DNM.sh: DNM properties.
DNM_R1.R: Analysing DNMs and plotting
DNM_indv_filtering.sh: Individual filtering.
Trio_combine_vcf.sh: Combine vcf in trio units and site filtering.
