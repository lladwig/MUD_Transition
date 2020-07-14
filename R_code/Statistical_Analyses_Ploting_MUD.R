#####Code for the analyses used in Ladwig et al 2020####

#####Begin analyses####
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")


library("ggplot2")
library("phyloseq")
library(dplyr)
library(vegan)
library(MASS) 
library(lme4)
library(Ternary) # making triangular plots
library(car) #for correct ANOVA table
library(ggrepel)
library(RVAideMemoire) # for tests after PERMANOVA
library(emmeans)
library("seqinr")
library(reshape2)
options(contrasts=c("contr.sum", "contr.poly"))
pa=function(x)(ifelse(x>0,1,0))
'%w/o%' <- function(x,y)!('%in%'(x,y))
#Load enzyme data
eea <-read.csv("R_files/EEA_MUD_transition.csv")

#Load soil characteristics
char <-read.csv("R_files/SoilCharacteristics_MUD_Transition.csv")

#Load in DESEQ2 normalized OTU matrix
MUD.Fung_vst=read.table("R_files/MUD_fungi_OTU_field_ITS_only_fung_VST.txt", header = T,row.names = 1)
fun_otu <- MUD.Fung_vst
fun_otu[1:10,1:10]
min(rowSums(fun_otu))

#Load in the mapping file
trans_map=read.csv("R_files/MUD_Transition_map.csv")
MUD.OTU = otu_table(fun_otu, taxa_are_rows = TRUE)
length(colnames(MUD.OTU))
#56
sort(colnames(MUD.OTU))
min(taxa_sums(MUD.OTU))
nrow(MUD.OTU)
ncol(MUD.OTU)
colnames(MUD.OTU)

#Load in the taxon table
MUD.fung.tax=as.matrix(read.table("R_files/MUD_fungi_fungi_taxon_ITS_trunc_phyl.txt",header=T))
nrow(MUD.fung.tax)
#3624
MUD.fung.TAX = tax_table(MUD.fung.tax)
head(eea)
length(row.names(eea))
#56
length(row.names(trans_map))
#90

#Merge the enzyme data with the soil charateristics
MUD.map_eea=merge(eea, trans_map,by="site")
head(MUD.map_eea)
length(row.names(MUD.map_eea))
#56
head(char)
length(row.names(char))
#56

#Remove colomuns that are repeats
char[,c("Location","Spp","Rep")]=NULL
head(char)
length(row.names(char))
#56

#Merge the mapping file with the enzymes ans soil charateristics
MUD.map_eea_char=merge(char, MUD.map_eea,by="site")
head(MUD.map_eea_char)
length(row.names(MUD.map_eea_char))
#56
rownames(MUD.map_eea_char)=MUD.map_eea_char$sampleID
MUD.map_eea_char$site_spp=with(MUD.map_eea_char, interaction(Site,Species))

#COmbine files into a phyloseq object
MUD.data=phyloseq(MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX)


ntaxa(MUD.data)
#2673
sum(otu_table(MUD.data))
#99908.89

min(sample_sums(MUD.data))
#1051.587

max(sample_sums(MUD.data))
#5041.467

sort(rank(sample_sums(MUD.data)))

min(taxa_sums(MUD.data))
#9.456247
sort(sample_sums(MUD.data))
MUD.data<-prune_taxa(taxa_sums(MUD.data) > 0, MUD.data)
ntaxa(MUD.data)
#2673
sum(otu_table(MUD.data))
#99908.89

min(sample_sums(MUD.data))
#1051.587

max(sample_sums(MUD.data))
#5041.467

length(sample_sums(MUD.data))
#56

#Save this phyloseq object so we can load it in each section
save(MUD.data, file = "R_files/MUD.data_dseq2_phyloseq_obj.RData")


####Diversity Analyses####
#Installing phyloseq. It was more than just "install.packages..." The below rdata won't load without it
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

load("R_files/MUD.data_dseq2_phyloseq_obj.RData")

alpha_meas = c("Shannon", "InvSimpson")
MUD.data.divfil=estimate_richness(MUD.data,measures=alpha_meas)
MUD.data_map=sample_data(MUD.data)

MUD.data.divfil=merge(MUD.data.divfil, MUD.data_map, by ="row.names")
nrow(MUD.data.divfil)
#56



#####Shannon Diversity Model####
hist(MUD.data.divfil$Shannon)
#Shannon
Shannon_mod <- lm((Shannon)^15 ~ Species + Site, data=MUD.data.divfil)
qqPlot(stdres(Shannon_mod))
hist(stdres(Shannon_mod))
shapiro.test(stdres(Shannon_mod))
#0.04533
summary(Shannon_mod)
Anova(Shannon_mod, type=3)
#nada sig


#####Graphing Shannon Diversity####
positions2=c("S","E","G")
spp_pos=c("Latr","Boer",  "Plja", "Bogr")
Shannon_g=ggplot(MUD.data.divfil, aes(x=Site, y=Shannon),
                    fill=Species)

(Shannon_box_p=Shannon_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Shannon_g$data, aes(x=Site, y=Shannon,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Shannon diversity")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                         values=c("turquoise3", "sienna", "bisque1","coral1" ),
                                                                                         labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=12,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

#####Inverse Simpson Diversity Model####
Simpson_mod <- lm((InvSimpson) ~ Species + Site, data=MUD.data.divfil)
qqPlot(stdres(Simpson_mod))
hist(stdres(Simpson_mod))
shapiro.test(stdres(Simpson_mod))
#0.2494
summary(Simpson_mod)
Anova(Simpson_mod, type=3)
#nada sig
#emmeans(Shannon_mod, pairwise~Site)

#####Graphing Inverse Simpson Diversity####
positions2=c("S","E","G")
spp_pos=c("Latr","Boer",  "Plja", "Bogr")
Simpson_g=ggplot(MUD.data.divfil, aes(x=Site, y=InvSimpson),
                 fill=Species)

(Simpson_box_p=Simpson_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Simpson_g$data, aes(x=Site, y=InvSimpson,
                                          fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Inverse Simpson diversity")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                values=c("turquoise3", "sienna", "bisque1", "coral1"),
                                                                                labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=12,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

### Graph for Chao is run around line 590  ######

# putting diversity graphs together
dev.off() #cleaning R, just in case
pdf("/Users/laura/Desktop/Writing Projects/MUD Transition/R 2020/MUD_Transition/results/Diversity_FigureS1_20200714.pdf", width = 5, height = 14) #This saves the pdf
ggarrange(Shannon_box_p, Simpson_box_p, Chao1_R_box_p,
          ncol = 1,
          nrow = 3,
          legend = "top", 
          common.legend = TRUE) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()


#####Phylum Level Stacked Bar Graph####
#How do phylum level community change over the transect

MUD.data_fact=merge_samples(MUD.data, "site_spp")
sample_names(MUD.data_fact)     

get_taxa_unique(MUD.data_fact, taxonomic.rank="Phylum")
#13
(MUD.data_fact.phylum<-tax_glom(MUD.data_fact, taxrank="Phylum"))




TopPHYL = names(sort(taxa_sums(MUD.data_fact.phylum), TRUE)[1:10])
mud_fung.T10 = prune_taxa(TopPHYL, MUD.data_fact.phylum)

plot_bar(mud_fung.T10, x= "site_spp", fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




MUD.data_fact.phylum.prop=transform_sample_counts(MUD.data_fact.phylum, function(x)x/sum(x))

taxa_phyla_names=sort(get_taxa_unique(mud_fung.T10, taxonomic.rank="Phylum"))
taxa_phyla_names_sep <- colsplit(taxa_phyla_names, ":", c("letter", "Phyl_name"))
taxa_phyla_names_sep[10,2]="Unclassified"
sample_names(MUD.data_fact.phylum.prop)
otu_table(MUD.data_fact.phylum.prop)
get_taxa_unique(MUD.data_fact.phylum.prop, taxonomic.rank="Phylum")

trans_positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")
fung.clayT10.prop = prune_taxa(TopPHYL, MUD.data_fact.phylum.prop)
(p_fung_T10=plot_bar(fung.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Spectral", label=c(taxa_phyla_names_sep$Phyl_name))+theme_bw()+theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                                                             axis.title=element_text(size=20),panel.grid.major=element_blank(),
                                                             panel.grid.minor=element_blank())+
    scale_x_discrete(limits = trans_positions))




#####Graphing the NMDS of DESEQ2 Normalized MAtrix####
load("R_files/MUD.data_dseq2_phyloseq_obj.RData")
#VST
MUD.data.ord <- ordinate(MUD.data, method="NMDS",distance = "bray")
#*** Solution reached
#Stress:     0.1480211 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'wisconsin(veganifyOTU(physeq))' 

plot_ordination(MUD.data, MUD.data.ord)
MUD.data_map=sample_data(MUD.data)
nrow(MUD.data_map)
#56
MUD.data_ord_points=merge(MUD.data.ord$points,MUD.data_map,by="row.names")
colnames(MUD.data_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(MUD.data_map$Site)
ggplot(data=MUD.data_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())


#####Community composition permANOVA####
MUD.data.dist<-phyloseq::distance(MUD.data,method = "bray")
MUD.data_map=sample_data(MUD.data)
test_perm <- adonis (MUD.data.dist ~ MUD.data_map$Location + MUD.data_map$Spp, permutations=10000)
print(test_perm)
#Pairwise analyses of communities between sites
pairwise.perm.manova(MUD.data.dist, MUD.data_map$Location, nperm=2000)
# RESULT: all locations different than one another

# airwise analyses of communities between species
pairwise.perm.manova(MUD.data.dist, MUD.data_map$Spp, nperm=2000)



#####SIMPER ANALYSES####
load("R_files/MUD.data_dseq2_phyloseq_obj.RData")
#need the OTU table
MUD.data_otu=t(otu_table(MUD.data))

colnames(MUD.data_otu)
row.names(MUD.data_otu)

#I also need the taxonomy table for extracting the classified species identity

MUD.data_taxa=as.data.frame(tax_table(MUD.data))

#let load in the OLD funguild results based on Sintax against UNITEv8
#####THIS FUNGUILD CLASSIFICATION WAS UPDATED BELOW####
funguild.80c=read.delim("R_files/MUD_fungi_fungi_OTU_ITS_trunc_phyl_for_FunGuildV2.guilds.txt",sep = "\t")
head(funguild.80c)
row.names(funguild.80c)=funguild.80c$OTU.ID
funguild.80c_names=funguild.80c[92:ncol(funguild.80c)]



#####LATR SIMPER####
#Now I need to extract out the samples that have LATR

MUD.data_map=sample_data(MUD.data)
MUD.data_map_LATR=subset(MUD.data_map, Spp=="LATR")
nrow(MUD.data_map_LATR)
MUD.data_LATR_otu=merge(MUD.data_map_LATR,MUD.data_otu, by="row.names")
nrow(MUD.data_LATR_otu)
colnames(MUD.data_LATR_otu)


MUD.data_LATR.simp <- with(MUD.data_LATR_otu, simper(MUD.data_LATR_otu[,29:ncol(MUD.data_LATR_otu)], Site,permutations=9999))
summary(MUD.data_LATR.simp,ordered = T)
MUD.data_LATR.simp_mat_num=as.data.frame(cbind(as.numeric(MUD.data_LATR.simp$E_S$average),as.numeric(MUD.data_LATR.simp$E_S$ava),
                                           as.numeric(MUD.data_LATR.simp$E_S$avb),as.numeric(MUD.data_LATR.simp$E_S$p)))
head(MUD.data_LATR.simp_mat_num)
row.names(MUD.data_LATR.simp_mat_num)=MUD.data_LATR.simp$E_S$species
summary(MUD.data_LATR.simp_mat_num)
colnames(MUD.data_LATR.simp_mat_num)[c(1:4)]=c("average","av_Ecotone","av_Shrub","pval")
(MUD.data_LATR.simp_mat_num[order(-MUD.data_LATR.simp_mat_num$average),])[1:10,]

#add in the taxonomy 

MUD.data_LATR.simp_taxa_mat=merge(MUD.data_LATR.simp_mat_num,MUD.data_taxa, by="row.names", all.x = T)
head(MUD.data_LATR.simp_taxa_mat)
#rename the species column
colnames(MUD.data_LATR.simp_taxa_mat)[1]="OTU"
rownames(MUD.data_LATR.simp_taxa_mat)=MUD.data_LATR.simp_taxa_mat$OTU
MUD.data_LATR.simp_taxa_mat_sig=subset(MUD.data_LATR.simp_taxa_mat, pval<0.05)
summary(MUD.data_LATR.simp_taxa_mat_sig)
nrow(MUD.data_LATR.simp_taxa_mat_sig)
MUD.data_LATR.simp_funguild_mat_sig=merge(MUD.data_LATR.simp_taxa_mat_sig,funguild.80c_names, by="row.names")
nrow(MUD.data_LATR.simp_funguild_mat_sig)
#56
write.csv(MUD.data_LATR.simp_funguild_mat_sig, "R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv")
MUD.data_LATR.simp_funguild_mat_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv", header = T)

#how many of the SIMPER id otus were in Funguild 
MUD.data_LATR.simp_funguild_mat_sig_class=subset(MUD.data_LATR.simp_funguild_mat_sig, Guild!="-")
nrow(MUD.data_LATR.simp_funguild_mat_sig_class)
#19


#I want to extract the rep sequences that did not classify well for LATR_Eco_V_Shrub
rep_set.fung<- read.fasta(file = "D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/USEARCHv11_files/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa", as.string = TRUE, set.attributes = FALSE)
colnames(MUD.data_LATR.simp_funguild_mat_sig)
nrow(MUD.data_LATR.simp_funguild_mat_sig)
MUD.data_LATR.simp_funguild_mat_sig_uncl=subset(MUD.data_LATR.simp_funguild_mat_sig,Taxon=="-")
nrow(MUD.data_LATR.simp_funguild_mat_sig_uncl)
#37

rep_set_LATR_Eco_V_Shrub_uncl_OTUs=rep_set.fung[names(rep_set.fung) %in% MUD.data_LATR.simp_funguild_mat_sig_uncl[,"OTU"]]
head(rep_set_LATR_Eco_V_Shrub_uncl_OTUs)
length(rep_set_LATR_Eco_V_Shrub_uncl_OTUs)
#37
write.fasta(sequences =rep_set_LATR_Eco_V_Shrub_uncl_OTUs, names = names(rep_set_LATR_Eco_V_Shrub_uncl_OTUs), file.out ="R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub_uncl_rep_set.fna")

#To start out, let's take the the top 10 that explain the most variation

MUD.data_LATR.simp_taxa_mat_sig_top10=(MUD.data_LATR.simp_taxa_mat_sig[order(-MUD.data_LATR.simp_taxa_mat_sig$average),])[1:10,]
rownames(MUD.data_LATR.simp_taxa_mat_sig_top10)



#####Ecotone SIMPER####

MUD.data_map=sample_data(MUD.data)
MUD.data_map_Ecotone=subset(MUD.data_map, Location=="Ecotone")
nrow(MUD.data_map_Ecotone)


MUD.data_map_Ecotone_otu=merge(MUD.data_map_Ecotone,MUD.data_otu, by="row.names")
nrow(MUD.data_map_Ecotone_otu)
colnames(MUD.data_map_Ecotone_otu)

MUD.data_map_Ecotone.simp <- with(MUD.data_map_Ecotone_otu, simper(MUD.data_map_Ecotone_otu[,29:ncol(MUD.data_map_Ecotone_otu)], Spp,permutations=9999))
summary(MUD.data_map_Ecotone.simp,ordered = T)

#In Ecotone BOER compared to LATR
MUD.data_map_Ecotone.simp_mat_BOER_LATR=as.data.frame(cbind(as.numeric(MUD.data_map_Ecotone.simp$BOER_LATR$average),as.numeric(MUD.data_map_Ecotone.simp$BOER_LATR$ava),
                                               as.numeric(MUD.data_map_Ecotone.simp$BOER_LATR$avb),as.numeric(MUD.data_map_Ecotone.simp$BOER_LATR$p)))
head(MUD.data_map_Ecotone.simp_mat_BOER_LATR)
row.names(MUD.data_map_Ecotone.simp_mat_BOER_LATR)=MUD.data_map_Ecotone.simp$BOER_LATR$species
summary(MUD.data_map_Ecotone.simp_mat_BOER_LATR)
colnames(MUD.data_map_Ecotone.simp_mat_BOER_LATR)[c(1:4)]=c("average","av_BEOR","av_LATR","pval")
(MUD.data_map_Ecotone.simp_mat_BOER_LATR[order(-MUD.data_map_Ecotone.simp_mat_BOER_LATR$average),])[1:10,]

#add in the taxonomy 

MUD.data_map_Ecotone.simp_taxa_BOER_LATR=merge(MUD.data_map_Ecotone.simp_mat_BOER_LATR,MUD.data_taxa, by="row.names", all.x = T)
head(MUD.data_map_Ecotone.simp_taxa_BOER_LATR)
#rename the species column
colnames(MUD.data_map_Ecotone.simp_taxa_BOER_LATR)[1]="OTU"
rownames(MUD.data_map_Ecotone.simp_taxa_BOER_LATR)=MUD.data_map_Ecotone.simp_taxa_BOER_LATR$OTU
MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig=subset(MUD.data_map_Ecotone.simp_taxa_BOER_LATR, pval<0.05)
summary(MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig)
nrow(MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig)
#66
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=merge(MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig, "R_files/MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR.csv")

#To start out, let's take the the top 10 that explain the most variation

MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top10=(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig[order(-MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$average),])[1:10,]
rownames(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top10)=MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top10$OTU
rownames(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top10)

#In Ecotone PLJA compared to LATR
MUD.data_map_Ecotone.simp_mat_LATR_PLJA=as.data.frame(cbind(as.numeric(MUD.data_map_Ecotone.simp$LATR_PLJA$average),as.numeric(MUD.data_map_Ecotone.simp$LATR_PLJA$ava),
                                                            as.numeric(MUD.data_map_Ecotone.simp$LATR_PLJA$avb),as.numeric(MUD.data_map_Ecotone.simp$LATR_PLJA$p)))
head(MUD.data_map_Ecotone.simp_mat_LATR_PLJA)
row.names(MUD.data_map_Ecotone.simp_mat_LATR_PLJA)=MUD.data_map_Ecotone.simp$LATR_PLJA$species
summary(MUD.data_map_Ecotone.simp_mat_LATR_PLJA)
colnames(MUD.data_map_Ecotone.simp_mat_LATR_PLJA)[c(1:4)]=c("average","av_LATR","av_PLJA","pval")
(MUD.data_map_Ecotone.simp_mat_LATR_PLJA[order(-MUD.data_map_Ecotone.simp_mat_LATR_PLJA$average),])[1:10,]

#add in the taxonomy 

MUD.data_map_Ecotone.simp_taxa_LATR_PLJA=merge(MUD.data_map_Ecotone.simp_mat_LATR_PLJA,MUD.data_taxa, by="row.names", all.x = T)
head(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA)
#rename the species column
colnames(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA)[1]="OTU"
rownames(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA)=MUD.data_map_Ecotone.simp_taxa_LATR_PLJA$OTU
MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig=subset(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA, pval<0.05)
summary(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig)
nrow(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig)
#337
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=merge(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig, "R_files/MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA.csv")


#####BEOR SIMPER####
#Now I need to extract out the samples that have BOER

MUD.data_map=sample_data(MUD.data)
MUD.data_map_BOER=subset(MUD.data_map, Spp=="BOER")
nrow(MUD.data_map_BOER)
MUD.data_map_BOER_otu=merge(MUD.data_map_BOER,MUD.data_otu, by="row.names")
nrow(MUD.data_map_BOER_otu)
colnames(MUD.data_map_BOER_otu)


MUD.data_BOER.simp <- with(MUD.data_map_BOER_otu, simper(MUD.data_map_BOER_otu[,29:ncol(MUD.data_map_BOER_otu)], Site,permutations=9999))
summary(MUD.data_BOER.simp,ordered = T)
MUD.data_BOER.simp_num=as.data.frame(cbind(as.numeric(MUD.data_BOER.simp$E_G$average),as.numeric(MUD.data_BOER.simp$E_G$ava),
                                               as.numeric(MUD.data_BOER.simp$E_G$avb),as.numeric(MUD.data_BOER.simp$E_G$p)))
head(MUD.data_BOER.simp_num)
row.names(MUD.data_BOER.simp_num)=MUD.data_BOER.simp$E_G$species
summary(MUD.data_BOER.simp_num)
colnames(MUD.data_BOER.simp_num)[c(1:4)]=c("average","av_Ecotone","av_Grassland","pval")
(MUD.data_BOER.simp_num[order(-MUD.data_BOER.simp_num$average),])[1:10,]

#add in the taxonomy 

MUD.data_BOER.simp_taxa=merge(MUD.data_BOER.simp_num,MUD.data_taxa, by="row.names", all.x = T)
head(MUD.data_BOER.simp_taxa)
#rename the species column
colnames(MUD.data_BOER.simp_taxa)[1]="OTU"
rownames(MUD.data_BOER.simp_taxa)=MUD.data_BOER.simp_taxa$OTU
MUD.data_BOER.simp_taxa_sig=subset(MUD.data_BOER.simp_taxa, pval<0.05)
summary(MUD.data_BOER.simp_taxa_sig)
nrow(MUD.data_BOER.simp_taxa_sig)
#119
MUD.data_BOER.simp_funguild_sig=merge(MUD.data_BOER.simp_taxa_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_BOER.simp_funguild_sig, "R_files/MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland.csv")


#####PLJA SIMPER####
#Now I need to extract out the samples that have PLJA

MUD.data_map=sample_data(MUD.data)
MUD.data_map_PLJA=subset(MUD.data_map, Spp=="PLJA")
nrow(MUD.data_map_PLJA)
#16
MUD.data_map_PLJA_otu=merge(MUD.data_map_PLJA,MUD.data_otu, by="row.names")
nrow(MUD.data_map_PLJA_otu)
colnames(MUD.data_map_PLJA_otu)


MUD.data_PLJA.simp <- with(MUD.data_map_PLJA_otu, simper(MUD.data_map_PLJA_otu[,29:ncol(MUD.data_map_PLJA_otu)], Site,permutations=9999))
summary(MUD.data_PLJA.simp,ordered = T)
MUD.data_PLJA.simp_num=as.data.frame(cbind(as.numeric(MUD.data_PLJA.simp$E_G$average),as.numeric(MUD.data_PLJA.simp$E_G$ava),
                                           as.numeric(MUD.data_PLJA.simp$E_G$avb),as.numeric(MUD.data_PLJA.simp$E_G$p)))
head(MUD.data_PLJA.simp_num)
row.names(MUD.data_PLJA.simp_num)=MUD.data_PLJA.simp$E_G$species
summary(MUD.data_PLJA.simp_num)
colnames(MUD.data_PLJA.simp_num)[c(1:4)]=c("average","av_Ecotone","av_Grassland","pval")
(MUD.data_PLJA.simp_num[order(-MUD.data_PLJA.simp_num$average),])[1:10,]

#add in the taxonomy 

MUD.data_PLJA.simp_taxa=merge(MUD.data_PLJA.simp_num,MUD.data_taxa, by="row.names", all.x = T)
head(MUD.data_PLJA.simp_taxa)
#rename the species column
colnames(MUD.data_PLJA.simp_taxa)[1]="OTU"
rownames(MUD.data_PLJA.simp_taxa)=MUD.data_PLJA.simp_taxa$OTU
MUD.data_PLJA.simp_taxa_sig=subset(MUD.data_PLJA.simp_taxa, pval<0.05)
summary(MUD.data_PLJA.simp_taxa_sig)
nrow(MUD.data_PLJA.simp_taxa_sig)
MUD.data_PLJA.simp_funguild_sig=merge(MUD.data_PLJA.simp_taxa_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_PLJA.simp_funguild_sig, "R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv")



#To start out, let's take the the top 10 that explain the most variation

MUD.data_LATR.simp_taxa_mat_sig_top10=(MUD.data_LATR.simp_taxa_mat_sig[order(-MUD.data_LATR.simp_taxa_mat_sig$average),])[1:10,]
rownames(MUD.data_LATR.simp_taxa_mat_sig_top10)

#We need a fasta file with the sequences from the SIMPER sig OTUs


#####Extract SIMPER identified OTUs for further Analyses####

#need to load in the csv files we created above

MUD.data_PLJA.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv")
MUD.data_BOER.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland.csv")
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA.csv")
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR.csv")
MUD.data_LATR.simp_funguild_mat_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv")

head(MUD.data_PLJA.simp_funguild_sig)
MUD.data_simper_OTU=c(as.character(MUD.data_PLJA.simp_funguild_sig$OTU),as.character(MUD.data_BOER.simp_funguild_sig$OTU),
                      as.character(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig$OTU),
                      as.character(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$OTU),as.character(MUD.data_LATR.simp_funguild_mat_sig$OTU))

length(MUD.data_simper_OTU)
#384

#remove duplicated OTUs

MUD.data_simper_OTU_unq=unique(MUD.data_simper_OTU)
length(MUD.data_simper_OTU_unq)
#290
anyDuplicated(MUD.data_simper_OTU_unq)

rep_set.fung<- read.fasta(file = "D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/USEARCHv11_files/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa", as.string = TRUE, set.attributes = FALSE)

head(rep_set.fung)
length(rep_set.fung)
#4833
rep_set_simper_OTUs=rep_set.fung[names(rep_set.fung) %in% MUD.data_simper_OTU_unq]
head(rep_set_simper_OTUs)
length(rep_set_simper_OTUs)
#290
write.fasta(sequences =rep_set_simper_OTUs, names = names(rep_set_simper_OTUs), file.out ="R_files/MUD.data_simper_OTU_rep_set.fna")


#####Community analyses on Untransformed data####
load("R_files/MUD.data_dseq2_phyloseq_obj.RData")
#write.table(otu_table(MUD.Fung_only_field), "D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/R_files/MUD_fungi_OTU_field_ITS_only_fung_untransformed.txt")
untrans.otu=read.table("R_files/MUD_fungi_OTU_field_ITS_only_fung_untransformed.txt")
colnames(untrans.otu)[colnames(untrans.otu)=="fungi.Shurb.LATR.2"]="fungi.Shrub.LATR.2"
untrans.MUD.OTU = otu_table(untrans.otu, taxa_are_rows = TRUE)
sum(taxa_sums(untrans.MUD.OTU))
#2288921
MUD.data_untrans=phyloseq(untrans.MUD.OTU,sample_data(MUD.data),tax_table(MUD.data))
MUD.data_untrans<-prune_taxa(taxa_sums(MUD.data_untrans) > 0, MUD.data_untrans)
sum(taxa_sums(MUD.data_untrans))
#2288921
ntaxa(MUD.data_untrans)
#2673

untrans.MUD.data.divfil=estimate_richness(MUD.data_untrans,measures="Chao1")
untrans.MUD.data_map=sample_data(MUD.data_untrans)

untrans.MUD.data.divfil=merge(untrans.MUD.data.divfil, untrans.MUD.data_map, by ="row.names")

#####TO BE DELETED####
MUD.data.divfil$Row.names=NULL

write.csv(merge(MUD.data.divfil,untrans.MUD.data.divfil[,c("sampleID","Chao1")], by="sampleID"),"R_files/Diversity_data_MUD.csv")
#####TO BE DELETED####


#####Chao1 Statistical Analyses#####
Chao1_mod <- lm((Chao1) ~ Species + Site, data=untrans.MUD.data.divfil)
qqPlot(stdres(Chao1_mod))
hist(stdres(Chao1_mod))
shapiro.test(stdres(Chao1_mod))
#0.9977
summary(Chao1_mod)
Anova(Chao1_mod, type=3)
#nada sig


#####Graphing Chao1####
positions2=c("S","E","G")
spp_pos=c("Latr","Boer",  "Plja", "Bogr")
Chao1_R_g=ggplot(untrans.MUD.data.divfil, aes(x=Site, y=Chao1),
                    fill=Species)

(Chao1_R_box_p=Chao1_R_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Chao1_R_g$data, aes(x=Site, y=Chao1,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Chao1 Richness")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                    values=c("turquoise3", "sienna", "bisque1", "coral1"),
                                                                                    labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=12,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))



#####Lee Taylor lead Protax and BLAST Taxonomy####
#Here Lee Tayler constructed a taxonmy with Protax and BLAST against NCBI to better characterize unclassified taxa.
#This anaylses was focused on the Top 500 OTUs by abundance
#Lee also used FunGuild to Classify the function of the OTUs 

#VST phyloseq obj beginning of analyses section 

load("R_files/MUD.data_dseq2_phyloseq_obj.RData")
#MUD.data
ntaxa(MUD.data)
#2673
sum(otu_table(MUD.data))
#99908.89

min(sample_sums(MUD.data))
#1051.587

max(sample_sums(MUD.data))
#5041.467

length(sample_sums(MUD.data))
#56

#extract OTU table
otu_table(MUD.data)[1:10,1:10]
otu_MUD.data=data.frame(otu_table(MUD.data))
otu_MUD.data[1:10,1:10]
otu_MUD.data$OTU.ID=row.names(otu_MUD.data)
head(otu_MUD.data$OTU.ID)
subset(otu_MUD.data, OTU.ID=="OTU10 ")


#load in the top 500 taxa classified with protax and funGuild

top500_funGuild_raw=read.csv("R_files/Top500Funguild_R.csv", header = T)
colnames(top500_funGuild_raw)
nrow(top500_funGuild_raw)
#[1] "Rank"                        "Final.Taxonomy.for.Funguild" "OTU.ID"                     
#[4] "taxonomy"                    "Taxon"                       "Taxon.Level"                
#[7] "Trophic.Mode"                "Guild"                       "Confidence.Ranking"         
#[10] "Growth.Morphology"           "Trait"                       "Notes"                      
#[13] "Citation.Source"      
head(top500_funGuild_raw$OTU.ID)
length(top500_funGuild_raw$OTU.ID)
#501
otu_MUD_top500_FG=merge(otu_MUD.data,top500_funGuild_raw, by="OTU.ID", all.y = T)
nrow(otu_MUD_top500_FG)
#501
subset(otu_MUD_top500_FG,is.na(fungi.Grass.BOGR.1))$OTU.ID
#"OTU001" "OTU10" are not in my original otu table 

#OTU10 is the Malassezia which I removed from previous analyses
#OTU001 is blank in the file from Lee? I think I need to drop both of them
otu_MUD_top500_FG_sub=subset(otu_MUD_top500_FG, OTU.ID!="OTU001"&OTU.ID!="OTU10")
nrow(otu_MUD_top500_FG_sub)
#499
colnames(otu_MUD_top500_FG_sub)

otu_MUD_top500_FG_sub$OTU_sum=rowSums(otu_MUD_top500_FG_sub[,2:57])
otu_MUD_top500_FG_sub_sort=otu_MUD_top500_FG_sub[order(-otu_MUD_top500_FG_sub$OTU_sum),]
head(otu_MUD_top500_FG_sub_sort$OTU.ID)
head(MUD.OTU_sum_taxa_funguild_sort$OTU)
unique(otu_MUD_top500_FG_sub$Trophic.Mode)
#[1] Saprotroph                        -                                 Pathotroph-Saprotroph-Symbiotroph
#[4] Saprotroph-Symbiotroph            Pathotroph                        Pathotroph-Saprotroph            
#[7] Symbiotroph                       Saprotroph-Pathotroph-Symbiotroph

unique(otu_MUD_top500_FG_sub$Guild)
#there are 53 seperate Guilds

#The prevalence and abundance of the Trophic modes
otu_MUD_top500_FG_sub %>% group_by(Trophic.Mode) %>% summarise_at(vars(OTU_sum),list(~sum(.),~n(),se=~sd(.)/sqrt(n()),~sd(.)))

#The prevalence and abundance of the Guilds
otu_MUD_top500_FG_sub_guild=otu_MUD_top500_FG_sub %>% group_by(Guild) %>% summarise_at(vars(OTU_sum),list(~sum(.),~n(),se=~sd(.)/sqrt(n()),~sd(.)))

#Let's extract the tava in order of abundance
otu_MUD_top500_FG_sub_guild[order(-otu_MUD_top500_FG_sub_guild$n),]


#####Trophic mode analyses####

#Let's look at main categories of
#Saprotroph Pathotroph Symbiotroph
colnames(otu_MUD_top500_FG_sub)

#Creation of the Trophic mode dataset by summing OTUs by Trophic mode
otu_MUD_top500_FG_sub_troph=otu_MUD_top500_FG_sub[,c(2:57,63)]%>%group_by(Trophic.Mode)%>%summarise_all(~sum(.))
otu_MUD_top500_FG_sub_main_trop=subset(otu_MUD_top500_FG_sub_troph, Trophic.Mode=="Saprotroph"|
                                          Trophic.Mode=="Pathotroph"|Trophic.Mode=="Symbiotroph")

otu_MUD_top500_FG_sub_main_trop_M=melt(otu_MUD_top500_FG_sub_main_trop)

#Merge with the metadata
otu_MUD_top500_FG_sub_main_trop_M_trt=merge(otu_MUD_top500_FG_sub_main_trop_M,sample_data(MUD.data),by.x = "variable",
                                            by.y = "row.names")


#####TO BE DELETED####
head(otu_MUD_top500_FG_sub_main_trop_M_trt)
write.csv(otu_MUD_top500_FG_sub_main_trop_M_trt, "R_files/Trophic_mode_FUNGuild_data.csv")
#####TO BE DELETED####

####Saprotroph Analyses####
otu_MUD_top500_FG_sub_main_trop_M_trt_sap=subset(otu_MUD_top500_FG_sub_main_trop_M_trt, Trophic.Mode=="Saprotroph")

Saprt_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_trop_M_trt_sap)
qqPlot(stdres(Saprt_abun_R_mod))
hist(stdres(Saprt_abun_R_mod))
shapiro.test(stdres(Saprt_abun_R_mod))
#0.0333
summary(Saprt_abun_R_mod)
Anova(Saprt_abun_R_mod, type=3)
#nada sig


####Graphing Saprotroph Abundance####
positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Saprotrophs_g=ggplot(otu_MUD_top500_FG_sub_main_trop_M_trt_sap, aes(x=Site, y=value),
                 fill=Species)

(Saprotrophs_T_p=Saprotrophs_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Saprotrophs_g$data, aes(x=Site, y=value,
                                          fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Saprotroph reads\n(VST normalized)")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                values=c("coral1", "sienna", "bisque1", "turquoise3"),
                                                                                labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))



#####Pathotroph Ananlyses####
otu_MUD_top500_FG_sub_main_trop_M_trt_path=subset(otu_MUD_top500_FG_sub_main_trop_M_trt, Trophic.Mode=="Pathotroph")

Patho_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_trop_M_trt_path)
qqPlot(stdres(Patho_abun_R_mod))
hist(stdres(Patho_abun_R_mod))
shapiro.test(stdres(Patho_abun_R_mod))
#0.02566
summary(Patho_abun_R_mod)
Anova(Patho_abun_R_mod, type=3)
#nada sig


####Graphing Pathotroph Abundance####
positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Pathotrophs_g=ggplot(otu_MUD_top500_FG_sub_main_trop_M_trt_path, aes(x=Site, y=value),
                     fill=Species)

(Pathotroph_T_p=Pathotrophs_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Pathotrophs_g$data, aes(x=Site, y=value,
                                              fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Pathotroph")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                                 values=c("coral1", "sienna", "bisque1", "turquoise3"),
                                                                                                 labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))




#####Symbiotroph Analyses####
otu_MUD_top500_FG_sub_main_trop_M_trt_symbio=subset(otu_MUD_top500_FG_sub_main_trop_M_trt, Trophic.Mode=="Symbiotroph")

Symbio_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_trop_M_trt_symbio)
qqPlot(stdres(Symbio_abun_R_mod))
hist(stdres(Symbio_abun_R_mod))
shapiro.test(stdres(Symbio_abun_R_mod))
#0.5152
summary(Symbio_abun_R_mod)
Anova(Symbio_abun_R_mod, type=3)
#Site          0.77  2    2.9378 0.06219 .
emmeans(Symbio_abun_R_mod, pairwise~Site)
#E - G      0.0833 0.128 50 0.650   0.7931 
#E - S      0.4229 0.181 50 2.335   0.0601 
#G - S      0.3396 0.222 50 1.531   0.2852 

#####Graphing Symbiotroph####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Symbiotroph_g=ggplot(otu_MUD_top500_FG_sub_main_trop_M_trt_symbio, aes(x=Site, y=value),
                     fill=Species)

(Symbiotroph_T_p=Symbiotroph_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Symbiotroph_g$data, aes(x=Site, y=value,
                                              fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Symbiotroph")+xlab(NULL)+
    annotate("text", x = 0.75, y = 10 , label = "p = 0.06", size = 5) +
    scale_fill_manual(limits = spp_pos, 
                                                                                                 values=c("turquoise3", "sienna", "bisque1", "coral1"),
                                                                                                 labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))





#####Guild Analyses####

#Need to first creat the data file by summing OTUs by Guild classification

#First let's combine similar guilds (i.e. Saprotroph groups)
otu_MUD_top500_FG_sub$Guild_comb=if_else(otu_MUD_top500_FG_sub$Guild=="Undefined Saprotroph","Saprotroph",
                                                       if_else(otu_MUD_top500_FG_sub$Guild=="Soil Saprotroph","Saprotroph",
                                                               if_else(otu_MUD_top500_FG_sub$Guild=="Dung Saprotroph-Plant Saprotroph-Soil Saprotroph",
                                                                       "Saprotroph",if_else(otu_MUD_top500_FG_sub$Guild=="Leaf Saprotroph","Saprotroph",
                                                                                            if_else(otu_MUD_top500_FG_sub$Guild=="Dung Saprotroph-Plant Saprotroph","Saprotroph",
                                                                                                      if_else(otu_MUD_top500_FG_sub$Guild=="Wood Saprotroph","Saprotroph",
                                                                                                              if_else(otu_MUD_top500_FG_sub$Guild=="Plant Saprotroph-Wood Saprotroph", "Saprotroph",
                                                                                                                      if_else(otu_MUD_top500_FG_sub$Guild=="Dung Saprotroph-Wood Saprotroph","Saprotroph",
                                                                                                                              if_else(otu_MUD_top500_FG_sub$Guild=="Dung Saprotroph-Soil Saprotroph", 
                                                                                                                                      "Saprotroph",as.character(otu_MUD_top500_FG_sub$Guild))))))))))

colnames(otu_MUD_top500_FG_sub)
otu_MUD_top500_FG_sub_GUILD=otu_MUD_top500_FG_sub[,c(2:57,71)]%>%group_by(Guild_comb)%>%summarise_all(~sum(.))
otu_MUD_top500_FG_sub_main_guild=subset(otu_MUD_top500_FG_sub_GUILD, Guild_comb=="Arbuscular Mycorrhizal"|
                                          Guild_comb=="Saprotroph"|Guild_comb=="Plant Pathogen"|Guild_comb=="Endophyte"
                                          )

otu_MUD_top500_FG_sub_main_guild_M=melt(otu_MUD_top500_FG_sub_main_guild)

otu_MUD_top500_FG_sub_main_guild_M_trt=merge(otu_MUD_top500_FG_sub_main_guild_M,sample_data(MUD.data),by.x = "variable",
                                            by.y = "row.names")


#####TO BE DELETED####
head(otu_MUD_top500_FG_sub_main_guild_M_trt)
write.csv(otu_MUD_top500_FG_sub_main_guild_M_trt, "R_files/Guil_comb_FUNGuild_data.csv")
#####TO BE DELETED####




#####Arbuscular Mycorrhizal Analyses####
otu_MUD_top500_FG_sub_main_guild_M_trt_arb=subset(otu_MUD_top500_FG_sub_main_guild_M_trt, Guild_comb=="Arbuscular Mycorrhizal")

Arbu_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_guild_M_trt_arb)
qqPlot(stdres(Arbu_abun_R_mod))
hist(stdres(Arbu_abun_R_mod))
shapiro.test(stdres(Arbu_abun_R_mod))
#0.9303
summary(Arbu_abun_R_mod)
Anova(Arbu_abun_R_mod, type=3)
#Site          1.57  2    3.8232 0.02851 * 
emmeans(Arbu_abun_R_mod, pairwise~Site)
#E - G       0.112 0.160 50 0.700   0.7648 
#E - S       0.605 0.226 50 2.675   0.0267 
#G - S       0.493 0.277 50 1.780   0.1866



#####Graphing Arbuscular Mycorrhizal#####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Arbuscular_g=ggplot(otu_MUD_top500_FG_sub_main_guild_M_trt_arb, aes(x=Site, y=value),
                     fill=Species)

(Arbuscular_T_p=Arbuscular_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Arbuscular_g$data, aes(x=Site, y=value, fill=factor(Species, levels=spp_pos)), outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+
    scale_colour_manual(values=c("black","black","black","black"), labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Arbuscular mycorrhizal")+
    xlab(NULL)+
    scale_fill_manual(limits = spp_pos, values=c("turquoise3", "sienna", "bisque1", "coral1"),                                          labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+
    annotate("text", x = 0.75, y = 10 , label = "p = 0.03", size = 5) +
    annotate("text", x = 1, y =100, label = "A", size = 6) +
    annotate("text", x = 2, y =250, label = "B", size = 6) +
    annotate("text", x = 3, y =200, label = "AB", size = 6) +
    theme(legend.title = element_blank(), legend.text=element_text(size=16,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))




#####Saprotroph Guild Analyses#####
otu_MUD_top500_FG_sub_main_guild_M_trt_Sapr=subset(otu_MUD_top500_FG_sub_main_guild_M_trt, Guild_comb=="Saprotroph")

Sapr_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_guild_M_trt_Sapr)
qqPlot(stdres(Sapr_abun_R_mod))
hist(stdres(Sapr_abun_R_mod))
shapiro.test(stdres(Sapr_abun_R_mod))
#0.01
summary(Sapr_abun_R_mod)
Anova(Sapr_abun_R_mod, type=3)
#nada sig
emmeans(Sapr_abun_R_mod, pairwise~Site)


#####Graphing Saprotroph Guild####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Sapro_g=ggplot(otu_MUD_top500_FG_sub_main_guild_M_trt_Sapr, aes(x=Site, y=value),
                    fill=Species)

(Saprotroph_G_p=Sapro_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Sapro_g$data, aes(x=Site, y=value,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Saprotroph")+xlab(NULL)+
    scale_fill_manual(limits = spp_pos, values=c("turquoise3", "sienna", "bisque1", "coral1"), labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+
    theme(legend.title = element_blank(), legend.text=element_text(size=16,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))



####Plant Pathogen Analyses####
otu_MUD_top500_FG_sub_main_guild_M_trt_Path=subset(otu_MUD_top500_FG_sub_main_guild_M_trt, Guild_comb=="Plant Pathogen")

Path_abun_R_mod <- lm(log(value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_guild_M_trt_Path)
qqPlot(stdres(Path_abun_R_mod))
hist(stdres(Path_abun_R_mod))
shapiro.test(stdres(Path_abun_R_mod))
#0.1994
summary(Path_abun_R_mod)
Anova(Path_abun_R_mod, type=3)
#nada sig
emmeans(Path_abun_R_mod, pairwise~Site)

####Graphing Plant Pathogen####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Plant_path_g=ggplot(otu_MUD_top500_FG_sub_main_guild_M_trt_Path, aes(x=Site, y=value),
               fill=Species)

(Plant_path_G_p=Plant_path_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Plant_path_g$data, aes(x=Site, y=value,
                                        fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Plant pathogen")+
    xlab(NULL)+
    scale_fill_manual(limits = spp_pos,  values=c("turquoise3", "sienna", "bisque1", "coral1"), labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))



#####Endophytes Analyses####
otu_MUD_top500_FG_sub_main_guild_M_trt_Endo=subset(otu_MUD_top500_FG_sub_main_guild_M_trt, Guild_comb=="Endophyte")

Endo_abun_R_mod <- lm((value) ~ Species + Site, data=otu_MUD_top500_FG_sub_main_guild_M_trt_Endo)
qqPlot(stdres(Endo_abun_R_mod))
hist(stdres(Endo_abun_R_mod))
shapiro.test(stdres(Endo_abun_R_mod))
#0.3664
summary(Endo_abun_R_mod)
Anova(Endo_abun_R_mod, type=3)
#nada sig
emmeans(Endo_abun_R_mod, pairwise~Site)

#####Graphing Endophytes####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Endop_g=ggplot(otu_MUD_top500_FG_sub_main_guild_M_trt_Endo, aes(x=Site, y=value),
                    fill=Species)

(Endophyte_G_p=Endop_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Endop_g$data, aes(x=Site, y=value,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+
    scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Endophyte")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                                     values=c("turquoise3", "sienna", "bisque1", "coral1"),
                                                                                                     labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16,face = "italic"), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))


## Trying to organize the graphs into a panel
library(ggpubr)

## TO DO:  ajust font size to look good, 
dev.off() #cleaning R, just in case
pdf("/Users/laura/Desktop/Writing Projects/MUD Transition/R 2020/MUD_Transition/results/FUNGuild_Figure6_20200714.pdf", width = 10, height = 14) #This saves the pdf

ggarrange(Endophyte_G_p, Symbiotroph_T_p, Plant_path_G_p, Saprotroph_G_p, Arbuscular_T_p,
          ncol = 2,
          nrow = 3,
          legend = "top", 
          common.legend = TRUE)
dev.off()




#####SIMPER analyses with the updated funguild results#####
top500_funGuild_raw=read.csv("R_files/Top500Funguild_R.csv", header = T)
colnames(top500_funGuild_raw)
nrow(top500_funGuild_raw)
head(top500_funGuild_raw)

#####PLJA SIMPER FUNGuild####
#PLJA Between Grassland and Ecotone
MUD.data_PLJA.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv", row.names = "Row.names")
head(MUD.data_PLJA.simp_funguild_sig)
nrow(MUD.data_PLJA.simp_funguild_sig)
#70
#need to remove the previous classifications 
MUD.data_PLJA.simp_funguild_sig[,c("X","Domain","Phylum","Class","Order","Family","Genus","Species","taxonomy","Taxon","Taxon.Level","Trophic.Mode","Guild",
                                   "Growth.Morphology","Trait","Confidence.Ranking","Notes","Citation.Source")]=NULL

#Let's create a column that categories the simper direction 
MUD.data_PLJA.simp_funguild_sig$SIM_dir=if_else(MUD.data_PLJA.simp_funguild_sig$av_Ecotone<MUD.data_PLJA.simp_funguild_sig$av_Grassland, "Grassland","Ecotone")


#Merge with the new classification
MUD.data_PLJA.simp_funguild_sig_top500=merge(MUD.data_PLJA.simp_funguild_sig,top500_funGuild_raw,by.x="OTU", by.y="OTU.ID")
head(MUD.data_PLJA.simp_funguild_sig_top500)
nrow(MUD.data_PLJA.simp_funguild_sig_top500)
#45

MUD.data_PLJA.simp_funguild_sig_top500 %>% group_by(SIM_dir,Trophic.Mode,Guild)%>%summarise_at(vars(av_Ecotone,av_Grassland),
                                                                                               list(~sum(.),~n()))


write.csv(MUD.data_PLJA.simp_funguild_sig_top500,"R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland_T500.csv", row.names = T)

#####BOER SIMPER FUNGuild####
#BOER in Grassland v Ecotone
MUD.data_BOER.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland.csv")
head(MUD.data_BOER.simp_funguild_sig)
nrow(MUD.data_BOER.simp_funguild_sig)
#119
#need to remove the previous classifications 
MUD.data_BOER.simp_funguild_sig[,c("X","Domain","Phylum","Class","Order","Family","Genus","Species","taxonomy","Taxon","Taxon.Level","Trophic.Mode","Guild",
                                   "Growth.Morphology","Trait","Confidence.Ranking","Notes","Citation.Source")]=NULL

#Let's create a column that categories the simper direction 
MUD.data_BOER.simp_funguild_sig$SIM_dir=if_else(MUD.data_BOER.simp_funguild_sig$av_Ecotone<MUD.data_BOER.simp_funguild_sig$av_Grassland, "Grassland","Ecotone")


#Merge with the new classification
MUD.data_BOER.simp_funguild_sig_top500=merge(MUD.data_BOER.simp_funguild_sig,top500_funGuild_raw,by.x="OTU", by.y="OTU.ID")
head(MUD.data_BOER.simp_funguild_sig_top500)
nrow(MUD.data_BOER.simp_funguild_sig_top500)
#70

MUD.data_BOER.simp_funguild_sig_top500 %>% group_by(SIM_dir,Trophic.Mode,Guild)%>%summarise_at(vars(av_Ecotone,av_Grassland),
                                                                                               list(~sum(.),~n()))


write.csv(MUD.data_BOER.simp_funguild_sig_top500,"R_files/MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland_T500.csv", row.names = T)

#####Ecotone SIMPER FUNGuild####
#In Ecotone LATR versus PLJA
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA.csv")
head(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig)
nrow(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig)
#73
#need to remove the previous classifications 
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig[,c("X","Domain","Phylum","Class","Order","Family","Genus","Species","taxonomy","Taxon","Taxon.Level","Trophic.Mode","Guild",
                                   "Growth.Morphology","Trait","Confidence.Ranking","Notes","Citation.Source")]=NULL

#Let's create a column that categories the simper direction 
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig$SIM_dir=if_else(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig$av_LATR<MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig$av_PLJA, "PLJA","LATR")


#Merge with the new classification
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig_top500=merge(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig,top500_funGuild_raw,by.x="OTU", by.y="OTU.ID")
head(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig_top500)
nrow(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig_top500)
#36

MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig_top500 %>% group_by(SIM_dir,Trophic.Mode,Guild)%>%summarise_at(vars(av_LATR,av_PLJA),
                                                                                               list(~sum(.),~n()))


write.csv(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig_top500,"R_files/MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA_T500.csv", row.names = T)



#In Ecotone LATR versus BOER
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR.csv")
head(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig)
nrow(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig)
#66
#need to remove the previous classifications 
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig[,c("X","Domain","Phylum","Class","Order","Family","Genus","Species","taxonomy","Taxon","Taxon.Level","Trophic.Mode","Guild",
                                                    "Growth.Morphology","Trait","Confidence.Ranking","Notes","Citation.Source")]=NULL

#Let's create a column that categories the simper direction 
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$SIM_dir=if_else(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$av_LATR<MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$av_BEOR, "BEOR","LATR")


#Merge with the new classification
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top500=merge(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig,top500_funGuild_raw,by.x="OTU", by.y="OTU.ID")
head(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top500)
nrow(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top500)
#40

MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top500 %>% group_by(SIM_dir,Trophic.Mode,Guild)%>%summarise_at(vars(av_LATR,av_BEOR),
                                                                                                                list(~sum(.),~n()))


write.csv(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top500,"R_files/MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR_T500.csv", row.names = T)

#####BOER SIMPER FUNGuild####
#BOER in Shrubland v Ecotone
MUD.data_LATR.simp_funguild_mat_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv")
head(MUD.data_LATR.simp_funguild_mat_sig)
nrow(MUD.data_LATR.simp_funguild_mat_sig)
#56
#need to remove the previous classifications 
MUD.data_LATR.simp_funguild_mat_sig[,c("X","Domain","Phylum","Class","Order","Family","Genus","Species","taxonomy","Taxon","Taxon.Level","Trophic.Mode","Guild",
                                                    "Growth.Morphology","Trait","Confidence.Ranking","Notes","Citation.Source")]=NULL

#Let's create a column that categories the simper direction 
MUD.data_LATR.simp_funguild_mat_sig$SIM_dir=if_else(MUD.data_LATR.simp_funguild_mat_sig$av_Ecotone<MUD.data_LATR.simp_funguild_mat_sig$av_Shrub, "Shrubland","Ecotone")


#Merge with the new classification
MUD.data_LATR.simp_funguild_mat_sig_top500=merge(MUD.data_LATR.simp_funguild_mat_sig,top500_funGuild_raw,by.x="OTU", by.y="OTU.ID")
head(MUD.data_LATR.simp_funguild_mat_sig_top500)
nrow(MUD.data_LATR.simp_funguild_mat_sig_top500)
#34

MUD.data_LATR.simp_funguild_mat_sig_top500 %>% group_by(SIM_dir,Trophic.Mode,Guild)%>%summarise_at(vars(av_Ecotone,av_Shrub),
                                                                                                                list(~sum(.),~n()))


write.csv(MUD.data_LATR.simp_funguild_mat_sig_top500,"R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub_T500.csv", row.names = T)
