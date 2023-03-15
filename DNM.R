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
otherDNM <- read.csv("Otherpapers.csv",stringsAsFactors = FALSE)



###### filters #######################
######################################

Value <- c(mutations$Mendelian_0,mutations$indelremove_0,mutations$GQ80_0,mutations$DP_0,mutations$AD0_0,mutations$AB_0,mutations$DNM)
Filter <-  c(rep("1:Mendelian_Violation",each=106), rep("2:Indel_adjacent",each=106), rep("3:GQ",each=106),rep("4:DP",each=106),rep("5:AD",each=106),rep("6:AB",each=106),rep("8:IGV",each=106))
Filter <- data.frame(Value,Filter)
Filter$SampleID <- rep(mutations$Sample_label,7)
Filter$Family <- rep(mutations$Pedigree,7)
Filter$ID <- c(c(1:10),c(1:12),rep(c(1:10),4),c(1:24),c(1:10),c(1:10))
Filter$POP <- c(rep("TVA",each=52),rep("POR",each=54))


ggplot(Filter , aes(fill=Filter,x = factor(ID), y = Value)) + 
  geom_bar(position="stack",stat = "identity",width = 1) +
  scale_x_discrete(name = "Family ID") +
  scale_y_continuous(name = "Number of DNM candidates left") +
  theme_classic(base_size = 10) +
  facet_grid(. ~ Family, space = 'free_x', scales = 'free_x')+
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line=element_line(color="grey"),strip.background =element_rect(color="grey")) 


require(reshape2)
ggplot(melt(Filter), aes(x=Value, fill = Filter)) + 
  #call geom_histogram with position="dodge" to offset the bars and manual binwidth of 2
  geom_histogram(position = "dodge", binwidth=1)




###### mutation rates ################
######################################

mean(mutations$DNM_rate) ## 4.367264e-09
sd(mutations$DNM_rate) ##2.779147e-09
l2.model <- lm(DNM_rate ~ 1, mutations)
confint(l2.model, level=0.95)  ##3.832034e-09 4.902495e-09
confint(l2.model, level=0.95) /2 ##1.916017e-09 2.451247e-09
se <- sd(mutations$DNM_rate)/sqrt(106) ##2.699346e-10
qt(p=0.05/2, df=105,lower.tail=F) ##1.982815
4.367264e-09 - 1.982815*se ##3.884317e-09

mean(subset(mutations,Population=="TVA")$DNM_rate)  #4.082692e-09
l2.model <- lm(DNM_rate ~ 1, subset(mutations,Population=="TVA")) 
confint(l2.model, level=0.95) # 3.306172e-09 4.859213e-09

mean(subset(mutations,Population=="POR")$DNM_rate)  #4.641296e-09
l2.model <- lm(DNM_rate ~ 1, subset(mutations,Population=="POR")) 
confint(l2.model, level=0.95) # 3.885891e-09 5.396701e-09

mean(mutations$DNM)  #3.028302
t.test(DNM_rate~ParentGeneration, data=mutations)


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
  theme_classic(base_size = 15) +
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
count <- c(158,97,28,12)
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

ggplot(data = subset(DNMposition,PMNumber==1&!is.na(Exon_type))) +
  geom_bar(aes(x=Exon_type,fill=Type_5)) +
  theme_classic(base_size = 10) +
  ylab("Number of DNMs observed")+
  xlab("Subtitution types") 
subset(DNMposition,PMNumber==1&!is.na(Exon_type)) %>% group_by(Exon_type) %>% count(Type_5) 


### Ts vs Tv

ggplot(data = subset(DNMposition,PMNumber==1&!is.na(Type_4))) +
  geom_bar(aes(x=Type_4,fill=Type_4)) +
  theme_classic(base_size = 10) +
  ylab("Number of DNMs observed")+
  xlab("Subtitution types") +
  scale_fill_manual(values=c("#8b7765", "#000080"),name=NULL) +
  guides(fill="none") 


Tsvcount <- subset(DNMposition,PMNumber==1) %>% count(Type_4) 
Tsvcount <- c(167,128)
chisq.test(Tsvcount,p=c(1/3,2/3))


### CGI

t.test(CGIcount$DNMrate~CGIcount$CGI,alternative="two.sided")

ggplot(data = CGIcount,aes(x=CGI,y=DNMrate,fill=CGI,color=CGI)) +
  geom_boxplot(outlier.shape = NA)+
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank()) +
  ylim(0, 2.1e-08)+
  ylab("DNM rates") +
  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)+
  scale_color_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)+
  stat_compare_means(method = "t.test", label.y = 2e-08,size=3) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.", size = 3,label.y = 1.8e-08,hide.ns = T,show.legend = F) + 
  theme(legend.position="bottom")



### weak strong pairing

ggplot(data = WScount,aes(x=Type_3,y=DNMrate,fill=Type_3)) +
  geom_point(size=0.5)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  theme_classic(base_size = 10) +
  theme(axis.title.x = element_blank()) +
  ylim(0, 1.15e-08)+
  ylab("DNM rates (currated by FNR)") +
  #  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)+
  #  scale_color_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)+
  stat_compare_means(method = "kruskal.test", label.y = 1.1e-08,size=3) +
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.",size = 3,label.y = 1.0e-08,hide.ns = T,show.legend = F) + 
  theme(legend.position="bottom")+
  scale_fill_discrete(name="")

WSrateSW <- subset(WScount,Type_3 == "S>W" | Type_3 == "W>S")
t.test(WSrateSW$DNMrate~WSrateSW$Type_3,alternative="two.sided")
mean(subset(WScount,Type_3=="S>W")$DNMrate)
## 2.879577e-09 1.470535e-09 2.388338e-09 1.606703e-09
## S>W S>S W>S W>W
subset(DNMposition,PMNumber==1) %>% count(Type_3) 
t.test(WSrateSW$n~WSrateSW$Type_3,alternative="two.sided")
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

ggplot(data = subset(DNMposition,PMNumber==1), aes(x = Type_5)) +
  geom_bar(aes(y=(..count..),fill=CpG)) +
  geom_text(stat="count",size=3,aes(group=Type_2,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)

# not PM
ggplot(data = subset(DNMposition,is.na(Parental_Mosaic)), aes(x = Type_5)) +
  geom_bar(aes(y=(..count..),fill=CpG)) +
  geom_text(stat="count",size=3,aes(group=Type_5,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)

# PM
ggplot(data = subset(DNMposition,!is.na(Parental_Mosaic) & PMNumber==1), aes(x = Type_5)) +
  geom_bar(aes(y=(..count..),fill=CpG)) +
  geom_text(stat="count",size=3,aes(group=Type_5,label=..count..),vjust=-1)+
  theme_classic(base_size = 10) +
  xlab("Mutation types") +
  ylab("Count") +
  ylim(0,60) +
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#F39B7FB2", "#00A087B2"),name=NULL)


###### current available studies #####
######################################


otherDNM$Reference <- factor(otherDNM$Reference, levels = otherDNM$Reference)
ggplot(data = otherDNM, aes(x = Reference, y = per.generation..E.08., color=Type)) +
  geom_point() +
  geom_errorbar(aes(x = Reference, ymin = CI_lowerbound..E.08., ymax = CI_higherbound..E.08.),color="#999999",width=0.1,size=0.3) +
  geom_text(aes(label=Species),hjust=0, vjust=-0.5,size=2.5,angle=90)+
  ylab("DNM rate (E-08/bp/generation)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values=c("#669966", "wheat4","skyblue2","#FF9966","plum4"),name=NULL) +
  theme(legend.position="bottom") 



