library(ggplot2)
library(plyr)
library(stringr)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(wesanderson)
library(ggchicklet)
library(ggmap)



mutations<- read.csv("mutation_TVA_POR.csv", header = T)
DNMposition <- read.csv("mutation_position_TVA_POR.csv", na.strings=c("", "NA"),header = T)
CGIcount <- read.csv("CGIcount.csv")
WScount <- read.csv("WScount.csv")
filter_stats <- read.csv("Filter.csv")
otherDNM <- read.csv("Otherpaper.csv",stringsAsFactors = FALSE)

###### filters #######################
######################################

Value <- c(filter_stats$Indel,filter_stats$GQ80,filter_stats$DP,filter_stats$AD0,filter_stats$AB,filter_stats$Cluster,filter_stats$RR_F,filter_stats$RR_I,filter_stats$IGV)
Filter <-  c(rep("1:Indel adjacent",each=106), rep("2:GQ",each=106),rep("3:DP",each=106),rep("4:AD",each=106),rep("5:AB",each=106),rep("6:Clustered",each=106),rep("7:Repetitive between families",each=106),rep("8:Repetitive within families",each=106),rep("9:IGV",each=106))
Filter <- data.frame(Value,Filter)
Filter$SampleID <- rep(mutations$Sample_label,9)
Filter$Family <- rep(mutations$Pedigree,9)
Filter$ID <- c(c(1:10),c(1:12),rep(c(1:10),4),c(1:24),c(1:10),c(1:10))
Filter$POP <- c(rep("TVA",each=52),rep("POR",each=54))


ggplot(Filter , aes(fill=Filter,x = factor(ID), y = Value)) + 
  geom_bar(position="stack",stat = "identity",width = 1) +
  scale_x_discrete(name = "Family ID") +
  scale_y_continuous(name = "Number of DNM candidates filtered") +
  theme_classic(base_size = 10) +
  facet_grid(. ~ Family, space = 'free_x', scales = 'free_x')+
  scale_fill_brewer(palette = "Set3") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line=element_line(color="grey"),strip.background =element_rect(color="grey")) 

ggplot(subset(Filter,Filter == "6:Clustered" | Filter == "7:Repetitive between families"| Filter == "8:Repetitive within families"| Filter == "9:IGV"), aes(fill=Filter,x = factor(ID), y = Value)) + 
  geom_bar(position="stack",stat = "identity",width = 1) +
  scale_x_discrete(name = "Family ID") +
  scale_y_continuous(name = "Number of DNM candidates filtered") +
  theme_classic(base_size = 10) +
  facet_grid(. ~ Family, space = 'free_x', scales = 'free_x')+
  scale_fill_manual(values=c("#fdb462", "#b3de69","#fccde5","#d9d9d9"))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line=element_line(color="grey"),strip.background =element_rect(color="grey")) 


###### prev. study ###################
######################################

otherDNM$Species_common <- factor(otherDNM$Species_common, levels = otherDNM$Species_common)
ggplot(data = otherDNM, aes(x = Species_common, y = per.generation..E.08., color=Type)) +
  geom_point() +
  geom_errorbar(aes(x = Species_common, ymin = CI_lowerbound..E.08., ymax = CI_higherbound..E.08.),color="#999999",width=0.1,size=0.3) +
  geom_text(aes(label=X.trios),hjust=-0.7, vjust=-0.6,size=2,angle=90)+
  ylab("DNM rate (E-08/bp/generation)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#ff7f00", "#984ea3","#377eb8","#4daf4a","#e41a1c","#a65628"),name=NULL) +
  theme(legend.position="bottom") 




###### mutation rates ################
######################################

mean(mutations$DNM_rate) ## 4.563396e-09
mean(mutations$DNM_rate)/2 ## 2.281698e-09
sd(mutations$DNM_rate) ##2.867872e-09
l2.model <- lm(DNM_rate ~ 1, mutations)
confint(l2.model, level=0.95)  ##4.011078e-09 5.115714e-09
confint(l2.model, level=0.95) /2 ##2.005539e-09 2.557857e-09
se <- sd(mutations$DNM_rate)/sqrt(106) ##2.699346e-10
qt(p=0.05/2, df=105,lower.tail=F) ##1.982815
4.367264e-09 - 1.982815*se ##3.814946e-09



### between generation

ggplot(data = mutations, aes(x = Generation, y = DNM_rate)) +
  geom_boxplot(aes(fill=Generation),color="#403d77",outlier.shape = NA) +
  geom_point(size=0.8) +
  theme_classic(base_size = 10) +
  xlab("Offspring generation") +
  ylab("Rates of germline mutation (/bp/gen)") +
  stat_compare_means(method = "t.test", label.y = 1.5e-08,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.35e-08,size=3.5) +
  scale_fill_manual(values=c( "#ffffff","#81b1d3"))+
  guides(fill="none") 

t.test(DNM_rate~Generation,data=mutations)

### between Sex
ggplot(data = mutations, aes(x = Sex, y = DNM_rate)) +
  geom_boxplot(aes(fill=Sex),color="#403d77",outlier.shape = NA) +
  geom_point(size=0.8) +
  theme_classic(base_size = 10) +
  xlab("Sex") +
  ylab("Rates of germline mutation (/bp/gen)") +
  stat_compare_means(method = "t.test", label.y = 1.5e-08,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.35e-08,size=3.5) +
  scale_fill_manual(values=c( "#ffffff","#81b1d3"))+
  guides(fill="none") 


### between pedigree type

ggplot(data = subset(mutations,Generation=="F2"), aes(x = BreedType, y = DNM_rate)) +
  geom_boxplot(aes(fill=BreedType),color="#403d77",outlier.shape = NA) +
  geom_point(size=0.8) +
  theme_classic(base_size = 10) +
  xlab("Breed Type") +
  ylab("Rates of germline mutation (/bp/gen)") +
  stat_compare_means(method = "t.test", label.y = 1.5e-08,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.35e-08,size=3.5) +
  scale_fill_manual(values=c( "#ffffff","#81b1d3"))+
  guides(fill="none") 


### transmission rate by pedigree type

ggplot(data = subset(DNMposition,!is.na(BreedType)), aes(x = BreedType, y = TransRate)) +
  geom_boxplot(aes(fill=BreedType),color="#403d77",outlier.shape = NA) +
  geom_jitter(size=0.8) +
  theme_classic(base_size = 10) +
  xlab("Breed Type") +
  ylab("Transmission rate") +
  stat_compare_means(method = "t.test", label.y = 1.1,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1,size=3.5) +
  scale_fill_manual(values=c( "#ffffff","#81b1d3"))+
  guides(fill="none") 

t.test(TransRate~Fam_ID, data=subset(DNMposition,!is.na(TransRate)))



### distribution

ggplot(mutations, aes(x=DNM)) +
  geom_histogram(fill="lightblue",alpha=0.5, position="identity",binwidth = 1) +
  geom_density(aes(y=..count..),color="darkblue", alpha=0.5, position="identity") +
  geom_vline(aes(xintercept=mean(DNM)),
             color="gray", linetype="dashed", size=0.5) +
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  theme_classic(base_size = 10) +
  scale_color_manual(values="#56B4E9") +
  xlab("Total number of DNMs per individual") + 
  ylab("Count")


### between populations

my_comparisons_2 <- c("TVA", "POR")

ggplot(data = subset(mutations), aes(x = Population, y = DNM_rate)) +
  geom_boxplot(outlier.shape = NA,fill="lightblue") +
  geom_point(size = 0.5) +
  theme_classic(base_size = 10) +
  xlab("Population") +
  ylab("DNM rate after curation") +
  # Add horizontal line at base mean
  stat_compare_means(method = "t.test", label.y = 1.6e-08,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.45e-08,size=3.5)

t.test(DNM_rate~Population, data=mutations)

mutations %>% group_by(Population) %>%    
  summarise(mean.rdnm = mean(DNM_rate, na.rm = TRUE),
            sd.rdnm = sd(DNM_rate, na.rm = TRUE),
            n.rdnm = n()) %>%
  mutate(se.rdnm = sd.rdnm / sqrt(n.rdnm),
           lower.ci.rdnm = mean.rdnm - qt(1 - (0.05 / 2), n.rdnm - 1) * se.rdnm,
           upper.ci.rdnm = mean.rdnm + qt(1 - (0.05 / 2), n.rdnm - 1) * se.rdnm)

## Population     mean.rdnm       sd.rdnm n.rdnm  se.rdnm lower.ci.rdnm upper.ci.rdnm
#  1 POR        0.00000000483 0.00000000270     54 3.68e-10 0.00000000409 0.00000000556
#  2 TVA        0.00000000429 0.00000000303     52 4.21e-10 0.00000000345 0.00000000513


###### DNM ###########################
######################################

### chr maps

Chr <- c(1:21)
Pos<-c(30106087,24029954,18520836,33625550,15538166,19164645,17739606,20480644,21145378,17012775,17852859,33585825,21992815,16408376,18287517,19652119,21092886,16178624,20450314,21492179,15364080)
Chrlength <- data.frame(Chr,Pos)
ggplot() + 
  geom_chicklet(data = Chrlength,aes(x=Chr, y=Pos),radius = grid::unit(1.5, "mm"),fill = "floralwhite",color="lightgrey") +
  geom_point(data = subset(DNMposition,PMNumber==1 & !is.na(Position)),aes(x=V7_CHR, y=V7_POS,color = Position),size=4,shape = 108) +
  scale_shape_identity() +
  theme_classic(base_size = 10) +
  coord_flip() + 
  ylab("Position") +
  xlab("Chromosome") +
  scale_color_manual(values=c("#8491B4B2", "#91D1C2B2","#E64B35B2","#999999"),name=NULL) 


### mutation position

ggplot(data =subset(DNMposition,PMNumber==1&Position != "UTR"), aes(y = Position,fill=Parental_Mosaic)) +
  geom_bar(stat="count") +
  theme_classic(base_size = 10) +
  ylab("Mutation location") +
  xlab("Count") +
  scale_fill_manual(values=c("#8491B4B2", "#91D1C2B2"),name=NULL) +
  theme(legend.position="bottom")

subset(DNMposition,PMNumber==1)  %>% count(Position) 
gene_type <- c("non-gene","intron","exon","UTR")
count <- c(165,103,28,12)
length <- c(218274546,161327958,36278040,9953065)
type_of_gene <- data.frame(gene_type,count,length)
chisq.test(count,p=length/425833609)
## overlaped gene and CDS
## gene: 207559063
## CDS: 36278040
## five_prime_UTR: 2553627
## three_prime_UTR: 7399438
## total: 425833609
## intron:165077566


### Exon_type

subset(DNMposition,PMNumber==1&!is.na(Exon_type)) %>% group_by(Exon_type) %>% count(Type_5) 
ggplot(data = subset(DNMposition,PMNumber==1&!is.na(Exon_type))) +
  geom_bar(aes(x=Exon_type,fill=Exon_type),color="#403d77") +
  theme_classic(base_size = 10) +
  ylab("Number of DNMs observed")+
  xlab("Subtitution types") +
  scale_fill_manual(values=c("#ffffff","#81b1d3"),name=NULL) +
  guides(fill="none") 


### Ts vs Tv

ggplot(data = subset(DNMposition,PMNumber==1&!is.na(Type_4))) +
  geom_bar(aes(x=Type_4,fill=Type_4),color="#403d77") +
  theme_classic(base_size = 10) +
  ylab("Number of DNMs observed")+
  xlab("Subtitution types") +
  scale_fill_manual(values=c("#ffffff","#81b1d3"),name=NULL) +
  guides(fill="none") 

subset(DNMposition,PMNumber==1) %>% count(Type_4) 
Tsvcount <- c(170,138)
chisq.test(Tsvcount,p=c(1/3,2/3))


### CGI

# check distribution
plotn<-ggplot(data = CGIcount,aes(x=n,fill=CGI)) +
  geom_histogram(alpha=0.5,bins=10) +
  scale_x_continuous(breaks = seq(0, 10, by = 2))+
  theme_classic(base_size = 10) 
plotDNMrate<-ggplot(data = CGIcount,aes(x=DNMrate*1e08,fill=CGI)) +
  geom_histogram(alpha=0.5,bins=20) +
  theme_classic(base_size = 10)
plotDNMrate_CpG <- ggplot(data = CGIcount,aes(x=DNMrate_CpG*1e08,fill=CGI)) +
  geom_histogram(alpha=0.5,bins=20) +
  theme_classic(base_size = 10) 
plotDNMrate_CG_NonCpG <- ggplot(data = CGIcount,aes(x=DNMrate_CG_NonCpG*1e08,fill=CGI)) +
  geom_histogram(alpha=0.5,bins=20) +
  theme_classic(base_size = 10) 
ggarrange(plotn, plotDNMrate, plotDNMrate_CpG, plotDNMrate_CG_NonCpG, 
          labels = c("A", "B", "C","D"),
          font.label = list(size = 12),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend="bottom")

CGIcount  %>% group_by(CGI) %>% 
  summarise_at(vars(CalGenSize, CalGenSize_CpG, CalGenSize_NonCpG), list(sum=sum))
##    CGI    CalGenSize_sum CalGenSize_CpG_sum CalGenSize_NonCpG_sum
##   CGI        4395421989    1260390091      3135031898
##   Non-CGI    34516242009   1353217928      33163024081

DNMposition %>% group_by(CGI) %>% count(CpG)
##   CGI     CpG         n
##   <chr>   <chr>   <int>
##   1 CGI     CpG        16
##   2 CGI     Non-CpG    15
##   3 Non-CGI CpG        50
##   4 Non-CGI Non-CpG   254


Type <- c("CGI_all","CGI_CpG","CGI_NonCpG","Non-CGI_all","Non-CGI_CpG","Non-CGI_NonCpG")
CS <- c(4395421989,1260390091,3135031898,34516242009,1353217928,33163024081)
DNMnumber <- c(31,16,15,304,50,254)
CGIrate <- data.frame(Type,CS, DNMnumber)
CGIrate$DNMrate <- CGIrate$DNMnumber/(2*CGIrate$CS)

chisq.test(c(16,15),p=c(0.2867507,0.7132493)) ### CGI CpG vs non-CpG
## X-squared = 7.9748, df = 1, p-value = 0.004743

for(i in 1:3) {
print(chisq.test(c(CGIrate$DNMnumber[i],CGIrate$DNMnumber[i+3]),p=c(CGIrate$CS[i],CGIrate$CS[i+3])/(sum(c(CGIrate$CS[i],CGIrate$CS[i+3])))))
}
### CpG CGI vs Non-CGI: X-squared = 11.819, df = 1, p-value = 9.658e-05
### non-CpG CGI vs Non-CGI: X-squared = 2.0579, df = 1, p-value = 0.07393
### CGI vs Non-CGI: X-squared = 0.46578, df = 1, p-value = 0.2377


ggplot(data = CGIcount,aes(x=CGI,y=DNMrate,fill=CGI)) +
  geom_violin(color="#403d77")+
  geom_point(size=0.5)+
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank()) +
  ylab("DNM rates") +
  scale_fill_manual(values=c("#81b1d3","#ffffff"),name=NULL) +
  stat_compare_means(method = "kruskal.test", label.y = 3.2e-08,size=3) +
  theme(legend.position="bottom")

kruskal.test(DNMrate ~ CGI, data = CGIcount) ## Kruskal-Wallis chi-squared = 35.701, df = 1, p-value = 2.3e-09
wilcox.test(DNMrate ~ CGI,  data = CGIcount, mu = 0, alternative = "two.sided")

# CGI, CpG and mutation type

CpGCGI <- DNMposition %>% group_by(CGI,CpG) %>% count(Position)
ggplot(data = CpGCGI,aes(x=Position,y=n,fill=CpG)) +
  geom_bar(position = position_dodge(preserve = "single"),stat = "identity",width = 0.8,color="#403d77") +
  theme_classic(base_size = 10) +
  facet_grid(. ~ CGI, space = 'free_x', scales = 'free_x')+
  ylab("Number of DNMs observed")+
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values=c("#81b1d3","#ffffff"),name=NULL) 

CGItype <- CGI %>% group_by(Type) %>%   
  summarise_at(vars(CpG_length), list(mean=mean,sum=sum))
CGItype$Type <- c("Intergenic CGI","Intragenic CGI","TSS CGI", "TTS CGI")
colnames(CGItype) <- c("Type","CGIlength_mean","CGIlength_sum")
CGItype$nDNM <- c(5,9,1,1)
CGItype$DNMrate <- as.numeric(CGItype$nDNM)/(2*106*as.numeric(CGItype$CGIlength_sum))
NonCGI <- c("Non CGI","NA",12766207,50,1.85E-08)
CGItype <- rbind(CGItype, NonCGI)

ggplot(data = CGItype) +
  geom_bar(aes(x=factor(Type, level = c("TSS CGI", "TTS CGI","Intergenic CGI","Intragenic CGI","Non CGI")),y=as.numeric(DNMrate)*1.0e08),position = position_dodge(preserve = "single"),stat = "identity",width = 0.8,color="#403d77",fill="white") +
  theme_classic(base_size = 10) +
  ylab("Rate of DNM (e-08/bp/generation)")+
  theme(axis.title.x = element_blank()) 


### weak strong pairing

ggplot(data = WScount,aes(x=Type_3,y=DNMrate)) +
  geom_point(size=0.5)+
  geom_boxplot(aes(fill=Type_3),outlier.shape = NA,color="#403d77")+
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank()) +
  ylim(0, 1.15e-08)+
  ylab("DNM rates (currated by FNR)") +
  scale_fill_manual(values=c("#E4F1F7", "#C5E1EF","#9EC9E2","#6CB0D6"),name=NULL)+
  stat_compare_means(method = "kruskal.test", label.y = 1.1e-08,size=3) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.",size = 3,label.y = 1.0e-08,hide.ns = T,show.legend = F)  

subset(WScount)  %>% group_by(Type_3) %>% 
  summarise_at(vars(DNMrate), list(median = median,mean=mean))
# Type_3          median          mean
# 1 S>S           0               4.75e-10
# 2 S>W           0.00000000276   2.42e- 9
# 3 W>S           0.00000000139   1.20e- 9
# 4 W>W           0               4.69e-10


## 0-inflated method - pairing type
WS <- subset(WScount)  %>% group_by(Type_3) %>% 
  summarise_at(vars(n,CalGenSize), list(sum = sum,mean=mean))
WS$DNMrate0 <- WS$n_sum/(2*WS$CalGenSize_sum)

# Type_3 n_sum CalGenSize_sum DNMrate0 
# 1 S>S       35    38911663998 4.497366e-10
# 2 S>W      178    38911663998 2.287232e-09
# 3 W>S       88    38911663998 1.130766e-09
# 4 W>W       34    38911663998 4.368870e-10

ggplot() +
  geom_point(data = WS, aes(x=Type_3, y=DNMrate0),color="#56B4E9") +
  theme_classic(base_size = 10) +
  xlab("Pairing Type") +
  ylim(0,2.5e-09) +
  ylab("DNM rate in 0-inflated method")

subset(DNMposition,PMNumber==1) %>% count(Type_3) 
## Type_3   n
## 1    S>S  32
## 2    S>W 165
## 3    W>S  80
## 4    W>W  32
kruskal.test(DNMrate ~ Type_3, data = WScount)
pairwise.wilcox.test(WScount$DNMrate, WScount$Type_3,
                     p.adjust.method = "BH")



### spectrum

ggplot(data = subset(DNMposition,PMNumber==1), aes(x = Type_2)) +
  geom_bar(aes(y=(..count..),fill=CpG)) +
  geom_text(stat="count",size=3,aes(group=Type_2,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)


# not PM
ggplot(data = subset(DNMposition,is.na(Parental_Mosaic)), aes(x = Type_5)) +
  geom_bar(aes(y=(..count..),fill=CpG),color="darkgray") +
  geom_text(stat="count",size=3,aes(group=Type_5,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  ylim(0,60) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#ffffff","#81b1d3"),name=NULL)

# PM
ggplot(data = subset(DNMposition,!is.na(Parental_Mosaic) & PMNumber==1), aes(x = Type_5)) +
  geom_bar(aes(y=(..count..),fill=CpG),color="#226E9C") +
  geom_text(stat="count",size=3,aes(group=Type_5,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  ylim(0,60) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#ffffff","#81b1d3"),name=NULL)

subset(DNMposition,PMNumber==1) %>% group_by(Type_5,CpG) %>% 
  count(Parental_Mosaic) %>% 
  mutate(prop = case_when(is.na(Parental_Mosaic)  ~ n/248,
                          Parental_Mosaic == "PM" ~ n/60)) %>%
  ggplot(aes(Type_5, prop, color = Parental_Mosaic),alpha=0.5) +
  geom_bar(aes(fill=CpG),stat = 'identity', position = 'dodge', size = 0.5,width=0.8) +
  scale_color_manual(values=c( "#226E9C","darkgray"),name=NULL) +
  scale_fill_manual(values=c("#81b1d3","#ffffff"),name=NULL)+
  theme_classic(base_size = 10) + 
  xlab("Mutation types") +
  ylab("Proportion") +
  theme(legend.position="bottom")

subset(DNMposition,PMNumber==1) %>% group_by(Type_5) %>% count(Parental_Mosaic)
PMcount <- c(4,12,8,7,7,10,12)
chisq.test(PMcount,p=c(22,42,24,41,25,45,49)/248)

subset(DNMposition,PMNumber==1) %>% group_by(CGI) %>% count(Parental_Mosaic)
PMcount_CGI <- c(2,58)
chisq.test(PMcount_CGI,p=c(31,277)/308)

### Parental differences in mutation spectrum

subset(DNMposition,!is.na(Phase) & PMNumber==1) %>% group_by(Type_5,CpG) %>% 
  count(Phase) %>% 
  ggplot(aes(Type_5, n, color = Phase),alpha=0.5) +
  geom_bar(aes(fill=CpG),stat = 'identity', position = 'dodge', size = 0.5,width=0.8) +
  geom_text(aes(Type_5, n,label = n), 
            position = position_dodge(width = 0.8), size=3,vjust=-1) +
  scale_color_manual(values=c("#3690c0", "#016450"),name=NULL) +
  scale_fill_manual(values=c("#81b1d3","#ffffff"),name=NULL)+
  theme_classic(base_size = 10) + ylim(0,35) +
  xlab("Mutation types") +
  ylab("Count") +
  theme(legend.position="bottom")

subset(DNMposition,!is.na(Phase) & PMNumber==1) %>% group_by(Type_5,CpG) %>% 
  count(Phase) %>% 
  mutate(prop = case_when(Phase == "M" ~ n/94,
                          Phase == "P" ~ n/120)) %>%
  ggplot(aes(Type_5, prop, color = Phase),alpha=0.5) +
  geom_bar(aes(fill=CpG),stat = 'identity', position = 'dodge', size = 0.5,width=0.8) +
  scale_color_manual(values=c("#3690c0", "#016450"),name=NULL) +
  scale_fill_manual(values=c("#81b1d3","#ffffff"),name=NULL)+
  theme_classic(base_size = 10) + 
  xlab("Mutation types") +
  ylab("Proportion") +
  theme(legend.position="bottom")

subset(DNMposition,PMNumber==1&!is.na(Phase)) %>% group_by(Type_5) %>% count(Phase)
chisq.test(c(6,21,8,18,7,21,13),p=c(10,16,16,19,10,18,31)/120)

chisq.test(c(81,13),p=c(89,31)/120) ## X-squared = 7.069, df = 1, p-value = 0.007843

### phase

phase <- subset(DNMposition,PMNumber==1) %>% group_by(Sample_ID) %>% count(Phase)
subset(DNMposition,PMNumber==1) %>% group_by(Parental_Mosaic) %>% count(Phase)

ggplot(data = subset(phase,!is.na(Phase)),aes(x=Phase,y=n)) +
  geom_boxplot(aes(fill=Phase),color="#226E9C",outlier.colour = NA) +
  geom_jitter(size=0.5,color="#226E9C") +
  theme_classic(base_size = 10) +
  xlab("Parental source") +
  ylab("Count") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#ffffff","#81b1d3"),name=NULL) +
  stat_compare_means(method = "t.test", label.y = 6,size=3) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.",size = 3,label.y = 5,hide.ns = T,show.legend = F)  

phase3 <- read_csv("phase.csv")
phase3.1 <- subset(phase3, (P+M) != 0)
mean(phase3.1$P/(phase3.1$M+phase3.1$P))
l3.model <- lm(P/(P+M) ~ 1, phase3.1)
confint(l3.model, level=0.95) 

phase3.lm <- lm(P ~ Generation + as.character(PedigreeType), phase3)
summary(phase3.lm)
anova(phase3.lm)
require(nlme)
glm.phase3.1 <- lme(P/(M+P)~Generation, random =~1|Pedigree,data=phase3.1)
anova(glm.phase3.1) ## p=0.5522

ggplot(data = phase3.1,aes(x=Generation,y=P/(P+M),color=Generation)) +
  geom_boxplot() +
  geom_jitter(size=0.5) +
  theme_classic(base_size = 10) +
  xlab("Generation") +
  ylab("Alpha") +
  stat_compare_means(method = "t.test", label.y = 1.25,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.1,size=3.5) +
  theme(legend.position="bottom")
t.test(P/(P+M)~Generation, data=phase3.1) ## p-value = 0.7612

ggplot(data = subset(phase3.1,(P/(P+M))!=0 & (P/(P+M))!=1),aes(x=Generation,y=P/(P+M),color=Generation)) +
  geom_boxplot() +
  geom_jitter(size=0.5) +
  theme_classic(base_size = 10) +
  xlab("Generation") +
  ylab("Alpha") +
  stat_compare_means(method = "t.test", label.y = 1.25,size =3.5) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 1.1,size=3.5) +
  theme(legend.position="bottom")

### mu comparison
prev.ns <- c(1.42E-08, 6.8E-08,3.7E-08,1E-08)
mu <- 0.456E-08
prev.ns/mu
subset(otherDNM, (CI_higherbound..E.08.>=0.401 & CI_higherbound..E.08.<=0.512) | (CI_lowerbound..E.08.>=0.401 & CI_lowerbound..E.08.<=0.512) | (CI_lowerbound..E.08.<=0.401 & CI_higherbound..E.08.>=0.512) | per.generation..E.08.>=0.401 &per.generation..E.08.<=0.512)
subset(otherDNM, per.generation..E.08.>=0.401 &per.generation..E.08.<=0.512)


