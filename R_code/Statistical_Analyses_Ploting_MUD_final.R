#####Code for the analyses used in Ladwig et al 2020####

#####Begin analyses####
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")


library("ggplot2")
library(phyloseq)
library(dplyr)
library(vegan)
library(MASS) 
# I don;t know what package I just accidentally deleted here, but "undo" is not working... Maybe be lme4 and lmerTest ?
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
MUD.fung.tax=as.matrix(read.table("R_files/MUD_fungi_taxon_PROTAX_ITS_trunc_phyl.txt",header=T))
nrow(MUD.fung.tax)
#2673
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

#Make a grouping factor for mixed effects models


sample_data(MUD.data)$transect_grp=with(sample_data(MUD.data),interaction(Location,Rep))

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
Shannon_mod <- lmer((Shannon)^15 ~ Species + Site+(1|transect_grp), data=MUD.data.divfil)
plot(Shannon_mod)
qqPlot(resid(Shannon_mod))
hist(resid(Shannon_mod))
shapiro.test(resid(Shannon_mod))
#0.03886

anova(Shannon_mod, type=3)
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
Simpson_mod <- lmer((InvSimpson) ~ Species + Site+(1|transect_grp), data=MUD.data.divfil)
plot(Simpson_mod)
qqPlot(resid(Simpson_mod))
hist(resid(Simpson_mod))
shapiro.test(resid(Simpson_mod))
#0.2116

anova(Simpson_mod, type=3)
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
pdf("/Users/laura/Desktop/Desktop/Writing Projects/MUD Transition/R 2020/MUD_Transition/results/Diversity_FigureS1_20210304.pdf", width = 5, height = 14) #This saves the pdf
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

## Copied over graphing code from other script below 
#unique(MUD.data_map$Site)
#ggplot(data=MUD.data_ord_points, aes(x=MDS1,y=MDS2))+
#  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
#  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
#                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
#  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
#  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
#                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
#                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())




# By site; Colored based on location
dev.off() #cleaning R, just in case
pdf("/Users/laura/Desktop/Desktop/Writing Projects/MUD Transition/R 2020/MUD_Transition/results/NMDS_Figure3_20210304.pdf", width = 7, height = 5.5) #This saves the pdf
ord_fig <- ggplot((data = MUD.data_ord_points), mapping = aes(x = MDS2, y = MDS1)) +
  geom_point(aes(color = Site, fill = Site, shape = Species), size = 4) + 
  scale_shape_manual(values = c(25, 22, 21, 24)) +
  #geom_segment(data=vec_sort_c,aes(x=0,xend=MDS1/2,y=0,yend=MDS2/2),arrow = arrow(length = unit(0.5, "cm")),colour="grey",inherit_aes=FALSE) + 
  #geom_text(data=vec_sort_d ,aes(x=MDS1/2,y=MDS2/2,label=Lees_id),size=4)+
  theme_classic() +
  scale_fill_manual(values = c("saddlebrown","brown1", "darkturquoise")) + #CHANGE!
  scale_color_manual(values = c("black", "black", "black")) +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),#black line all around graph
        text = element_text(size = 12),
        #legend.key.size = unit(4, "lines"),
        #legend.spacing = unit(0.5, "lines"),
        legend.title = element_blank(),
        legend.position = "top"
  )
ord_fig
dev.off()


#####Community composition permANOVA####
MUD.data.dist<-phyloseq::distance(MUD.data,method = "bray")
MUD.data_map=sample_data(MUD.data)
test_perm <- adonis (MUD.data.dist ~ MUD.data_map$Location + MUD.data_map$Spp, strata = MUD.data_map$transect_grp,permutations=10000)
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

#let load in the OLD funguild results based on Protax against UNITEv8

funguild.PROTAX=read.delim("R_files/MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.guilds.txt",sep = "\t")
head(funguild.PROTAX)
row.names(funguild.PROTAX)=funguild.PROTAX$OTU.ID
colnames(funguild.PROTAX)
funguild.PROTAX_names=funguild.PROTAX[58:ncol(funguild.PROTAX)]



#####LATR SIMPER####
#Now I need to extract out the samples that have LATR

MUD.data_map=sample_data(MUD.data)
MUD.data_map_LATR=subset(MUD.data_map, Spp=="LATR")
nrow(MUD.data_map_LATR)
MUD.data_LATR_otu=merge(MUD.data_map_LATR,MUD.data_otu, by="row.names")
nrow(MUD.data_LATR_otu)
colnames(MUD.data_LATR_otu)

#march 4, 2021 changed 29 to 30, but I'm not sure why an extra column would have been added with this update - LML
MUD.data_LATR.simp <- with(MUD.data_LATR_otu, simper(MUD.data_LATR_otu[,30:ncol(MUD.data_LATR_otu)], Site,permutations=9999))
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
MUD.data_LATR.simp_funguild_mat_sig=merge(MUD.data_LATR.simp_taxa_mat_sig,funguild.PROTAX_names, by="row.names")
nrow(MUD.data_LATR.simp_funguild_mat_sig)
#56
write.csv(MUD.data_LATR.simp_funguild_mat_sig, "R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv")
MUD.data_LATR.simp_funguild_mat_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv", header = T)

#how many of the SIMPER id otus were in Funguild 
MUD.data_LATR.simp_funguild_mat_sig_class=subset(MUD.data_LATR.simp_funguild_mat_sig, Guild!="-")
nrow(MUD.data_LATR.simp_funguild_mat_sig_class)
#19
#updated to 34; march 2021


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

MUD.data_map_Ecotone.simp <- with(MUD.data_map_Ecotone_otu, simper(MUD.data_map_Ecotone_otu[,30:ncol(MUD.data_map_Ecotone_otu)], Spp,permutations=9999))
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
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=merge(MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig,funguild.PROTAX_names, by="row.names")
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
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=merge(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig,funguild.PROTAX_names, by="row.names")
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
MUD.data_BOER.simp_funguild_sig=merge(MUD.data_BOER.simp_taxa_sig,funguild.PROTAX_names, by="row.names")
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
MUD.data_PLJA.simp_funguild_sig=merge(MUD.data_PLJA.simp_taxa_sig,funguild.PROTAX_names, by="row.names")
write.csv(MUD.data_PLJA.simp_funguild_sig, "R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv")



#To start out, let's take the the top 10 that explain the most variation

MUD.data_LATR.simp_taxa_mat_sig_top10=(MUD.data_LATR.simp_taxa_mat_sig[order(-MUD.data_LATR.simp_taxa_mat_sig$average),])[1:10,]
rownames(MUD.data_LATR.simp_taxa_mat_sig_top10)

#We need a fasta file with the sequences from the SIMPER sig OTUs


#####Extract SIMPER identified OTUs for further Analyses####

#need to load in the csv files we created above

MUD.data_PLJA.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv",row.names = 1)
MUD.data_BOER.simp_funguild_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland.csv",row.names = 1)
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA.csv",row.names = 1)
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR.csv",row.names = 1)
MUD.data_LATR.simp_funguild_mat_sig=read.csv("R_files/MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv",row.names = 1)



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
Chao1_mod <- lmer((Chao1) ~ Species + Site+ (1|transect_grp), data=untrans.MUD.data.divfil)
plot(Chao1_mod)
qqPlot(resid(Chao1_mod))
hist(resid(Chao1_mod))
shapiro.test(resid(Chao1_mod))
#0.9987

anova(Chao1_mod, type=3)
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



#####FunGuild to Classify the function of the OTUs####



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



#load in the top 500 taxa classified with protax and funGuild

funguild.PROTAX=read.delim("R_files/MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.guilds.txt",sep = "\t")
colnames(funguild.PROTAX)
nrow(funguild.PROTAX)
funguild.PROTAX_dt=funguild.PROTAX[,c(1,58:ncol(funguild.PROTAX))]
head(funguild.PROTAX_dt)
otu_MUD_PrTax_FG=merge(otu_MUD.data,funguild.PROTAX_dt, by="OTU.ID", all.y = T)
nrow(otu_MUD_PrTax_FG)
#2673
head(otu_MUD_PrTax_FG)



colnames(otu_MUD_PrTax_FG)
otu_MUD_PrTax_FG$OTU_sum=rowSums(otu_MUD_PrTax_FG[,2:57])
otu_MUD_PrTax_FG_sort=otu_MUD_PrTax_FG[order(-otu_MUD_PrTax_FG$OTU_sum),]
head(otu_MUD_PrTax_FG_sort$OTU.ID)

unique(otu_MUD_PrTax_FG$Trophic.Mode)
#[1] "Saprotroph"                        "-"                                 "Pathotroph"                        "Pathotroph-Saprotroph-Symbiotroph"
#[5] "Pathotroph-Saprotroph"             "Saprotroph-Symbiotroph"            "Symbiotroph"                       "Pathotroph-Symbiotroph"           
#[9] "Saprotroph; Symbiotroph"                       Saprotroph-Pathotroph-Symbiotroph

unique(otu_MUD_PrTax_FG$Guild)
#there are 73 seperate Guilds

#The prevalence and abundance of the Trophic modes
otu_MUD_PrTax_FG %>% group_by(Trophic.Mode) %>% summarise_at(vars(OTU_sum),list(~sum(.),~n(),se=~sd(.)/sqrt(n()),~sd(.)))

#The prevalence and abundance of the Guilds
otu_MUD_PrTax_FG_guild=otu_MUD_PrTax_FG %>% group_by(Guild) %>% summarise_at(vars(OTU_sum),list(~sum(.),~n(),se=~sd(.)/sqrt(n()),~sd(.)))

#Let's extract the tava in order of abundance
otu_MUD_PrTax_FG_guild[order(-otu_MUD_PrTax_FG_guild$n),]


unique(otu_MUD_PrTax_FG$Confidence.Ranking)
#Limit the classification to only the "Probable" and "Highly Probable" Confidence

otu_MUD_PrTax_FG_sub=subset(otu_MUD_PrTax_FG,Confidence.Ranking=="Probable"|Confidence.Ranking=="Highly Probable")
#Confidence.Ranking=="Probable"|Confidence.Ranking=="Highly Probable"
#How many of the OTUs does this represent?

nrow(otu_MUD_PrTax_FG_sub)
#1215
nrow(otu_MUD_PrTax_FG)
#2673

nrow(otu_MUD_PrTax_FG_sub)/nrow(otu_MUD_PrTax_FG)
#0.4545455

#How many of the VST reads does this represent?
sum(rowSums(otu_MUD_PrTax_FG_sub[,2:57]))
#44129.13
sum(rowSums(otu_MUD_PrTax_FG[,2:57]))
#99908.89

sum(rowSums(otu_MUD_PrTax_FG_sub[,2:57]))/sum(rowSums(otu_MUD_PrTax_FG[,2:57]))
#0.4416938



#####Trophic mode analyses####

#Let's look at main categories of
#Saprotroph Pathotroph Symbiotroph
colnames(otu_MUD_PrTax_FG_sub)

#Creation of the Trophic mode dataset by summing OTUs by Trophic mode
otu_MUD_PrTax_FG_sub_troph=otu_MUD_PrTax_FG_sub[,c(2:57,61)]%>%group_by(Trophic.Mode)%>%summarise_all(~sum(.))
otu_MUD_PrTax_FG_sub_main_trop=subset(otu_MUD_PrTax_FG_sub_troph, Trophic.Mode=="Saprotroph"|
                                          Trophic.Mode=="Pathotroph"|Trophic.Mode=="Symbiotroph")

otu_MUD_PrTax_FG_sub_main_trop_M=melt(otu_MUD_PrTax_FG_sub_main_trop)

#Merge with the metadata
otu_MUD_PrTax_FG_sub_main_trop_M_trt=merge(otu_MUD_PrTax_FG_sub_main_trop_M,sample_data(MUD.data),by.x = "variable",
                                            by.y = "row.names")


#####TO BE DELETED####
head(otu_MUD_PrTax_FG_sub_main_trop_M_trt)
write.csv(otu_MUD_PrTax_FG_sub_main_trop_M_trt, "R_files/Trophic_mode_FUNGuild_data.csv")
#####TO BE DELETED####

####Saprotroph Analyses####
otu_MUD_PrTax_FG_sub_main_trop_M_trt_sap=subset(otu_MUD_PrTax_FG_sub_main_trop_M_trt, Trophic.Mode=="Saprotroph")

Saprt_abun_R_mod <- lmer((value)^-1 ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_trop_M_trt_sap)
plot(Saprt_abun_R_mod)
qqPlot(resid(Saprt_abun_R_mod))
hist(resid(Saprt_abun_R_mod))
shapiro.test(resid(Saprt_abun_R_mod))
#0.9915

anova(Saprt_abun_R_mod, type=3)
#nada sig


####Graphing Saprotroph Abundance####
positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Saprotrophs_g=ggplot(otu_MUD_PrTax_FG_sub_main_trop_M_trt_sap, aes(x=Site, y=value),
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
otu_MUD_PrTax_FG_sub_main_trop_M_trt_path=subset(otu_MUD_PrTax_FG_sub_main_trop_M_trt, Trophic.Mode=="Pathotroph")

Patho_abun_R_mod <- lmer(log(value) ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_trop_M_trt_path)
plot(Patho_abun_R_mod)
qqPlot(resid(Patho_abun_R_mod))
hist(resid(Patho_abun_R_mod))
shapiro.test(resid(Patho_abun_R_mod))
#0.1749

anova(Patho_abun_R_mod, type=3)
#nada sig


####Graphing Pathotroph Abundance####
positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Pathotrophs_g=ggplot(otu_MUD_PrTax_FG_sub_main_trop_M_trt_path, aes(x=Site, y=value),
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
otu_MUD_PrTax_FG_sub_main_trop_M_trt_symbio=subset(otu_MUD_PrTax_FG_sub_main_trop_M_trt, Trophic.Mode=="Symbiotroph")

Symbio_abun_R_mod <- lmer(log(value) ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_trop_M_trt_symbio)
plot(Symbio_abun_R_mod)
qqPlot(resid(Symbio_abun_R_mod))
hist(resid(Symbio_abun_R_mod))
shapiro.test(resid(Symbio_abun_R_mod))
#0.2965

anova(Symbio_abun_R_mod, type=3)
#nada sig


#####Graphing Symbiotroph####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Symbiotroph_g=ggplot(otu_MUD_PrTax_FG_sub_main_trop_M_trt_symbio, aes(x=Site, y=value),
                     fill=Species)

(Symbiotroph_T_p=Symbiotroph_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Symbiotroph_g$data, aes(x=Site, y=value,
                                              fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Shrubland","Ecotone","Grassland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    scale_y_continuous(name = "Symbiotroph")+xlab(NULL)+
    annotate("text", x = 0.75, y = 10 , label = "p = 0.25", size = 5) +
    scale_fill_manual(limits = spp_pos, 
                                                                                                 values=c("turquoise3", "sienna", "bisque1", "coral1"),
                                                                                                 labels= c("L. tridentata","B. eriopoda","P. jamesii", "B. gracilis"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_blank(), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))





#####Guild Analyses####

#Need to first creat the data file by summing OTUs by Guild classification

#First let's combine similar guilds (i.e. Saprotroph groups)
otu_MUD_PrTax_FG_sub$Guild_comb=if_else(otu_MUD_PrTax_FG_sub$Guild=="Undefined Saprotroph","Saprotroph",
                                                       if_else(otu_MUD_PrTax_FG_sub$Guild=="Soil Saprotroph","Saprotroph",
                                                               if_else(otu_MUD_PrTax_FG_sub$Guild=="Dung Saprotroph-Plant Saprotroph-Soil Saprotroph",
                                                                       "Saprotroph",if_else(otu_MUD_PrTax_FG_sub$Guild=="Leaf Saprotroph","Saprotroph",
                                                                                            if_else(otu_MUD_PrTax_FG_sub$Guild=="Dung Saprotroph-Plant Saprotroph","Saprotroph",
                                                                                                      if_else(otu_MUD_PrTax_FG_sub$Guild=="Wood Saprotroph","Saprotroph",
                                                                                                              if_else(otu_MUD_PrTax_FG_sub$Guild=="Plant Saprotroph-Wood Saprotroph", "Saprotroph",
                                                                                                                      if_else(otu_MUD_PrTax_FG_sub$Guild=="Dung Saprotroph-Wood Saprotroph","Saprotroph",
                                                                                                                              if_else(otu_MUD_PrTax_FG_sub$Guild=="Dung Saprotroph-Soil Saprotroph", 
                                                                                                                                      "Saprotroph",as.character(otu_MUD_PrTax_FG_sub$Guild))))))))))

colnames(otu_MUD_PrTax_FG_sub)
ncol(otu_MUD_PrTax_FG_sub)
otu_MUD_PrTax_FG_sub_GUILD=otu_MUD_PrTax_FG_sub[,c(2:57,69)]%>%group_by(Guild_comb)%>%summarise_all(~sum(.))
otu_MUD_PrTax_FG_sub_main_guild=subset(otu_MUD_PrTax_FG_sub_GUILD, Guild_comb=="Arbuscular Mycorrhizal"|
                                          Guild_comb=="Saprotroph"|Guild_comb=="Plant Pathogen"|Guild_comb=="Endophyte"
                                          )

otu_MUD_PrTax_FG_sub_main_guild_M=melt(otu_MUD_PrTax_FG_sub_main_guild)

otu_MUD_PrTax_FG_sub_main_guild_M_trt=merge(otu_MUD_PrTax_FG_sub_main_guild_M,sample_data(MUD.data),by.x = "variable",
                                            by.y = "row.names")


#####TO BE DELETED####
head(otu_MUD_PrTax_FG_sub_main_guild_M_trt)
write.csv(otu_MUD_PrTax_FG_sub_main_guild_M_trt, "R_files/Guil_comb_FUNGuild_data.csv")
#####TO BE DELETED####




#####Arbuscular Mycorrhizal Analyses####
otu_MUD_PrTax_FG_sub_main_guild_M_trt_arb=subset(otu_MUD_PrTax_FG_sub_main_guild_M_trt, Guild_comb=="Arbuscular Mycorrhizal")

Arbu_abun_R_mod <- lmer(log(value) ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_guild_M_trt_arb)
plot(Arbu_abun_R_mod)
qqPlot(resid(Arbu_abun_R_mod))
hist(resid(Arbu_abun_R_mod))
shapiro.test(resid(Arbu_abun_R_mod))
#0.2492

anova(Arbu_abun_R_mod, type=3)
#Site    0.80881  0.4044     2    50  2.9878 0.05947 .
emmeans(Arbu_abun_R_mod, pairwise~Site)
#E - G      0.0575 0.130 31.0 0.442   0.8981 
#E - S      0.4422 0.184 50.0 2.404   0.0513 
#G - S      0.3847 0.225 47.3 1.708   0.2130 



#####Graphing Arbuscular Mycorrhizal#####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Arbuscular_g=ggplot(otu_MUD_PrTax_FG_sub_main_guild_M_trt_arb, aes(x=Site, y=value),
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
otu_MUD_PrTax_FG_sub_sub_main_guild_M_trt_Sapr=subset(otu_MUD_PrTax_FG_sub_main_guild_M_trt, Guild_comb=="Saprotroph")

Sapr_abun_R_mod <- lmer((value)^-1 ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_sub_main_guild_M_trt_Sapr)
plot(Sapr_abun_R_mod)
qqPlot(resid(Sapr_abun_R_mod))
hist(resid(Sapr_abun_R_mod))
shapiro.test(resid(Sapr_abun_R_mod))
#0.571

anova(Sapr_abun_R_mod, type=3)
#nada sig
emmeans(Sapr_abun_R_mod, pairwise~Site)


#####Graphing Saprotroph Guild####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Sapro_g=ggplot(otu_MUD_PrTax_FG_sub_sub_main_guild_M_trt_Sapr, aes(x=Site, y=value),
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
otu_MUD_PrTax_FG_sub_main_guild_M_trt_Path=subset(otu_MUD_PrTax_FG_sub_main_guild_M_trt, Guild_comb=="Plant Pathogen")

Path_abun_R_mod <- lmer((value)^-1 ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_guild_M_trt_Path)
plot(Path_abun_R_mod)
qqPlot(resid(Path_abun_R_mod))
hist(resid(Path_abun_R_mod))
shapiro.test(resid(Path_abun_R_mod))
#0.1218

anova(Path_abun_R_mod, type=3)
#nada sig
emmeans(Path_abun_R_mod, pairwise~Site)

####Graphing Plant Pathogen####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Plant_path_g=ggplot(otu_MUD_PrTax_FG_sub_main_guild_M_trt_Path, aes(x=Site, y=value),
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
otu_MUD_PrTax_FG_sub_main_guild_M_trt_Endo=subset(otu_MUD_PrTax_FG_sub_main_guild_M_trt, Guild_comb=="Endophyte")

Endo_abun_R_mod <- lmer(log(value) ~ Species + Site + (1|transect_grp), data=otu_MUD_PrTax_FG_sub_main_guild_M_trt_Endo)
plot(Endo_abun_R_mod)
qqPlot(resid(Endo_abun_R_mod))
hist(resid(Endo_abun_R_mod))
shapiro.test(resid(Endo_abun_R_mod))
#0.1391

anova(Endo_abun_R_mod, type=3)
#nada sig
emmeans(Endo_abun_R_mod, pairwise~Site)

#####Graphing Endophytes####
positions2=c("S","E","G")
spp_pos=c("Latr" ,"Boer",  "Plja", "Bogr")
Endop_g=ggplot(otu_MUD_PrTax_FG_sub_main_guild_M_trt_Endo, aes(x=Site, y=value),
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




#####Soil Charateristics ######
load("R_files/MUD.data_dseq2_phyloseq_obj.RData")
MUD.data_map=data.frame(sample_data(MUD.data))

#SAND          CLAY  pH P_ppm K_ppm OM_percent NO3_N_ppm
#% sand, 

Sand_R_mod <- lmer((SAND) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(Sand_R_mod)
qqPlot(resid(Sand_R_mod))
hist(resid(Sand_R_mod))
shapiro.test(resid(Sand_R_mod))
#0.8724

anova(Sand_R_mod, type=3)
#Site    581.01 290.507     2    49 11.1113 0.0001049 ***

emmeans(Sand_R_mod, pairwise~Site)


#% clay, 

CLAY_R_mod <- lmer((CLAY) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(CLAY_R_mod)
qqPlot(resid(CLAY_R_mod))
hist(resid(CLAY_R_mod))
shapiro.test(resid(CLAY_R_mod))
#0.6392

anova(CLAY_R_mod, type=3)
#Site    183.791  91.896     2    49  8.8651 0.0005175 ***

emmeans(CLAY_R_mod, pairwise~Site)

#% silt, 

SILT_R_mod <- lmer((SILT) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(SILT_R_mod)
qqPlot(resid(SILT_R_mod))
hist(resid(SILT_R_mod))
shapiro.test(resid(SILT_R_mod))
# 0.7855

anova(SILT_R_mod, type=3)
#Site    89.217  44.609     2 36.186  5.5168 0.008105 **

emmeans(SILT_R_mod, pairwise~Site)


#% organic matter (OM) 

OM_percent_R_mod <- lmer((OM_percent)^-1 ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(OM_percent_R_mod)
qqPlot(resid(OM_percent_R_mod))
hist(resid(OM_percent_R_mod))
shapiro.test(resid(OM_percent_R_mod))
#0.2242

anova(OM_percent_R_mod, type=3)
#Site    1.33705 0.66852     2    50 28.9058 4.546e-09 ***

emmeans(OM_percent_R_mod, pairwise~Site)

#ppm of P

P_ppm_R_mod <- lmer((P_ppm) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(P_ppm_R_mod)
qqPlot(resid(P_ppm_R_mod))
hist(resid(P_ppm_R_mod))
shapiro.test(resid(P_ppm_R_mod))
#0.6392

anova(P_ppm_R_mod, type=3)
#Species  140.89   46.96     3 33.647  3.0576   0.04152 *  
#Site    2685.54 1342.77     2 37.458 87.4251 7.747e-15 ***

emmeans(P_ppm_R_mod, pairwise~Site)
emmeans(P_ppm_R_mod, pairwise~Species)


#ppm of K

K_ppm_R_mod <- lmer((K_ppm) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(K_ppm_R_mod)
qqPlot(resid(K_ppm_R_mod))
hist(resid(K_ppm_R_mod))
shapiro.test(resid(K_ppm_R_mod))
#0.3247

anova(K_ppm_R_mod, type=3)
#Species  22097    7366     3 31.690  7.1289 0.0008592 ***
#Site     68551   34275     2 38.522 33.1737 4.194e-09 ***

emmeans(K_ppm_R_mod, pairwise~Site)
emmeans(K_ppm_R_mod, pairwise~Species)

#ppm of NO3-N

NO3_N_ppm_R_mod <- lmer((NO3_N_ppm)^-1 ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(NO3_N_ppm_R_mod)
qqPlot(resid(NO3_N_ppm_R_mod))
hist(resid(NO3_N_ppm_R_mod))
shapiro.test(resid(NO3_N_ppm_R_mod))
#0.09364

anova(NO3_N_ppm_R_mod, type=3)
#Site    0.50710 0.253552     2 38.869  4.9591 0.01207 *

emmeans(NO3_N_ppm_R_mod, pairwise~Site)



#####OM and Soil enzymes####

#OM NAG

NAG_OM_R_mod <- lmer(log(NAG_OM) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(NAG_OM_R_mod)
qqPlot(resid(NAG_OM_R_mod))
hist(resid(NAG_OM_R_mod))
shapiro.test(resid(NAG_OM_R_mod))
#0.4409


anova(NAG_OM_R_mod, type=3)
#nada sig

emmeans(NAG_OM_R_mod, pairwise~Site)

#OM BG 

BG_OM_R_mod <- lmer(log(BG_OM) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(BG_OM_R_mod)
qqPlot(resid(BG_OM_R_mod))
hist(resid(BG_OM_R_mod))
shapiro.test(resid(BG_OM_R_mod))
#0.9771

anova(BG_OM_R_mod, type=3)
#nada sig

emmeans(BG_OM_R_mod, pairwise~Site)

#OM AAP
MUD.data_map$AAP_OM=as.numeric(MUD.data_map$AAP_OM)
AAP_OM_R_mod <- lmer((AAP_OM) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(AAP_OM_R_mod)
qqPlot(resid(AAP_OM_R_mod))
hist(resid(AAP_OM_R_mod))
shapiro.test(resid(AAP_OM_R_mod))
#0.1833

anova(AAP_OM_R_mod, type=3)
#nada sig 

emmeans(AAP_OM_R_mod, pairwise~Site)


#OM AlkP

AlkP_OM_R_mod <- lmer(log(AlkP_OM) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(AlkP_OM_R_mod)
qqPlot(resid(AlkP_OM_R_mod))
hist(resid(AlkP_OM_R_mod))
shapiro.test(resid(AlkP_OM_R_mod))
#0.2341

anova(AlkP_OM_R_mod, type=3)
#nada sig 

emmeans(AlkP_OM_R_mod, pairwise~Site)

#soil NAG

NAG_soil_R_mod <- lmer(log(NAG_soil) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(NAG_soil_R_mod)
qqPlot(resid(NAG_soil_R_mod))
hist(resid(NAG_soil_R_mod))
shapiro.test(resid(NAG_soil_R_mod))
# 0.6768


anova(NAG_soil_R_mod, type=3)
#Site    2.2498 1.12491     2 33.773  4.8641 0.01392 *

emmeans(NAG_soil_R_mod, pairwise~Site)

#soil BG 

BG_soil_R_mod <- lmer(log(BG_soil) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(BG_soil_R_mod)
qqPlot(resid(BG_soil_R_mod))
hist(resid(BG_soil_R_mod))
shapiro.test(resid(BG_soil_R_mod))
#0.847

anova(BG_soil_R_mod, type=3)
#nada sig

emmeans(BG_soil_R_mod, pairwise~Site)

#soil AAP
MUD.data_map$AAP_soil=as.numeric(MUD.data_map$AAP_soil)
AAP_soil_R_mod <- lmer((AAP_soil) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(AAP_soil_R_mod)
qqPlot(resid(AAP_soil_R_mod))
hist(resid(AAP_soil_R_mod))
shapiro.test(resid(AAP_soil_R_mod))
#0.8596

anova(AAP_soil_R_mod, type=3)
#Site    13231.4  6615.7     2 31.802  3.1625 0.05589 .

emmeans(AAP_soil_R_mod, pairwise~Site)


#soil AlkP

AlkP_soil_R_mod <- lmer(log(AlkP_soil) ~ Species + Site + (1|transect_grp), data=MUD.data_map)
plot(AlkP_soil_R_mod)
qqPlot(resid(AlkP_soil_R_mod))
hist(resid(AlkP_soil_R_mod))
shapiro.test(resid(AlkP_soil_R_mod))
#0.5347

anova(AlkP_soil_R_mod, type=3)
#Site    0.80754 0.40377     2 28.172  3.0833 0.06153 .

emmeans(AlkP_soil_R_mod, pairwise~Site)
