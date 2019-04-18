rm(list=ls(all=TRUE))
library(phyloseq)
#####Processing conducted before by hand#####
#fungi

MUD.fung=subset_samples(MUD.data, funct_grp=="fungi")

MUD.FUNG<-prune_taxa(taxa_sums(MUD.fung) > 0, MUD.fung)
ntaxa(MUD.FUNG)
#4960

sum(otu_table(MUD.FUNG))
#4238112

MUD.Fung=subset_taxa(MUD.FUNG, Domain=="d:Fungi")
ntaxa(MUD.Fung)
#4951
sum(otu_table(MUD.Fung))
#4235030

sort(sample_sums(MUD.Fung))
#Min=7815
#Max=79429

write.table(tax_table(MUD.Fung), "MUD_fungi_taxon_ITS_trunc_phyl.txt") 
write.table(otu_table(MUD.Fung), "MUD_fungi_OTU_ITS_trunc_phyl.txt") 





#nomalizing using DESEQ2
MUD.Fung.deseq=phyloseq_to_deseq2(MUD.Fung, ~X.SampleID)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
MUD.Fung.vsd<-varianceStabilizingTransformation(MUD.Fung.deseq)#may result in the below error
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

head(MUD.Fung.vsd)
plotSparsity(assay(MUD.Fung.vsd))
meanSdPlot(assay(MUD.Fung.vsd))
#Fix for above common error from the support website https://support.bioconductor.org/p/63229/


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x>0]),na.rm = na.rm)/length(x))
}
geoMeans1=apply(counts(MUD.Fung.deseq), 1, gm_mean)
dseq_f1=estimateSizeFactors(MUD.Fung.deseq, geoMeans=geoMeans1)
head(dseq_f1)
plotSparsity(dseq_f1)
dseq_f2=varianceStabilizingTransformation(dseq_f1, fitType = "local")
head(assay(dseq_f2))
dseq_f3=varianceStabilizingTransformation(dseq_f1, blind = TRUE)
head(assay(dseq_f3))
dseq_f4=rlog(dseq_f1, blind = TRUE)#this is much slower than VST 
#rlog requires fitting a shrinkage term for each sample and each gene
#if your dataset has a lot of zeros, you will get the below message

#Warning message:
#In sparseTest(counts(object, normalized = TRUE), 0.9, 100, 0.1) :
#  the rlog assumes that data is close to a negative binomial distribution, an assumption
#which is sometimes not compatible with datasets where many genes have many zero counts
#despite a few very large counts.
#In this data, for 14.1% of genes with a sum of normalized counts above 100, it was the case 
#that a single sample's normalized count made up more than 90% of the sum over all samples.
#the threshold for this warning is 10% of genes. See plotSparsity(dds) for a visualization of this.
#We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).
head(assay(dseq_f4))

#below are visualizations of the sparsity of the data across the normalization strategies
plotSparsity(assay(dseq_f4))
plotSparsity(assay(dseq_f3))
plotSparsity(assay(dseq_f2))

#below are the visualizations of the standard deviation as read counts increase 
#https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
meanSdPlot(log2(counts(dseq_f1,normalized=TRUE) + 1))
meanSdPlot(assay(dseq_f4))
meanSdPlot(assay(dseq_f3))
meanSdPlot(assay(dseq_f2))
#"Note that the vertical axis in such plots is the square root of the variance over all samples, 
#so including the variance due to the experimental conditions. While a flat curve of the square root 
#of variance over the mean may seem like the goal of such transformations, this may be unreasonable 
#in the case of datasets with many true differences due to the experimental conditions."


MUD.Fung_vst=assay(MUD.Fung.vsd)#variance stabilized transformation OTU table in matrix form
MUD.Fung_vst=MUD.Fung_vst+abs(min(MUD.Fung_vst))#add a constant to remove any zeros
#the normalization results in negatives
#I have not found a consensus on how to handle them 

#now we need to turn it back into a phyloseq object
OTU.table.vst = otu_table(MUD.Fung_vst, taxa_are_rows = TRUE)
write.table(otu_table(OTU.table.vst), "MUD_fungi_OTU_ITS_VST.txt") 
phyloseq_obj_vst=phyloseq(OTU.table.vst,TAX.file, sample_data(map.file))
######End of Processing#####

setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")

#library(reshape)
#library(reshape2)
#library(stringr)
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
options(contrasts=c("contr.sum", "contr.poly"))
pa=function(x)(ifelse(x>0,1,0))
eea <-read.csv("EEA_MUD_transition.csv")
char <-read.csv("SoilCharacteristics_MUD_Transition.csv")
fun_otu <- read.csv("MUD_fungi_OTU_ITS_VST.csv",header=T, row.names = 1)
#untrans.otu=read.table("MUD_fungi_OTU_ITS_trunc_phyl.txt",header=T, row.names = 1)
#MUD.map=import_qiime_sample_data("MUD_Transition_map_qiime.txt")
trans_map=read.csv("MUD_Transition_map.csv")
MUD.OTU = otu_table(fun_otu, taxa_are_rows = TRUE)
MUD.fung.tax=as.matrix(read.table("MUD_fungi_taxon_ITS_trunc_phyl.txt",header=T))
MUD.fung.TAX = tax_table(MUD.fung.tax)
head(eea)
length(row.names(eea))
#56
length(row.names(trans_map))
#90
MUD.map_eea=merge(eea, trans_map,by="site")
head(MUD.map_eea)
length(row.names(MUD.map_eea))
head(char)
length(row.names(char))
char[,c("Location","Spp","Rep")]=NULL
head(char)
length(row.names(char))
MUD.map_eea_char=merge(char, MUD.map_eea,by="site")
head(MUD.map_eea_char)
length(row.names(MUD.map_eea_char))
rownames(MUD.map_eea_char)=MUD.map_eea_char$sampleID
MUD.map_eea_char$site_spp=with(MUD.map_eea_char, interaction(Site,Species))
MUD.data=phyloseq(MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX)


ntaxa(MUD.data)
#4951
sum(otu_table(MUD.data))
#120193.4

min(sample_sums(MUD.data))
#1114.705

max(sample_sums(MUD.data))
#3088.89

sort(rank(sample_sums(MUD.data)))

sort(sample_sums(MUD.data))
MUD.data<-prune_taxa(taxa_sums(MUD.data) > 0, MUD.data)
ntaxa(MUD.data)
#3640
sum(otu_table(MUD.data))
#120193.4

min(sample_sums(MUD.data))
#1114.705

max(sample_sums(MUD.data))
#3088.89


#Diversity 

alpha_meas = c("Shannon", "InvSimpson")
MUD.data.divfil=estimate_richness(MUD.data,measures=alpha_meas)
MUD.data_map=sample_data(MUD.data)

MUD.data.divfil=merge(MUD.data.divfil, MUD.data_map, by ="row.names")


#Shannon
Shannon_mod <- lm((Shannon) ~ Species + Site, data=MUD.data.divfil)
qqPlot(stdres(Shannon_mod))
hist(stdres(Shannon_mod))
shapiro.test(stdres(Shannon_mod))
#0.3025
summary(Shannon_mod)
Anova(Shannon_mod, type=3)
#nada sig
#emmeans(Shannon_mod, pairwise~Site)


positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Shannon_g=ggplot(MUD.data.divfil, aes(x=Site, y=Shannon),
                    fill=Species)

(Shannon_box_p=Shannon_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Shannon_g$data, aes(x=Site, y=Shannon,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Shannon diversity")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                         values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                                                                                         labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

#Simpson
Simpson_mod <- lm((InvSimpson) ~ Species + Site, data=MUD.data.divfil)
qqPlot(stdres(Simpson_mod))
hist(stdres(Simpson_mod))
shapiro.test(stdres(Simpson_mod))
#0.7644
summary(Simpson_mod)
Anova(Simpson_mod, type=3)
#nada sig
#emmeans(Shannon_mod, pairwise~Site)


positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Simpson_g=ggplot(MUD.data.divfil, aes(x=Site, y=InvSimpson),
                 fill=Species)

(Simpson_box_p=Simpson_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Simpson_g$data, aes(x=Site, y=InvSimpson,
                                          fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Inverse Simpson diversity")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                                                                                labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))



#How do phylum level community change over the transect

MUD.data_fact=merge_samples(MUD.data, "site_spp")
sample_names(MUD.data_fact)     

get_taxa_unique(MUD.data_fact, taxonomic.rank="Phylum")
#12
(MUD.data_fact.phylum<-tax_glom(MUD.data_fact, taxrank="Phylum"))




TopPHYL = names(sort(taxa_sums(MUD.data_fact.phylum), TRUE)[1:10])
mud_fung.T10 = prune_taxa(TopPHYL, MUD.data_fact.phylum)

plot_bar(mud_fung.T10, x= "site_spp", fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




MUD.data_fact.phylum.prop=transform_sample_counts(MUD.data_fact.phylum, function(x)x/sum(x))
sort(get_taxa_unique(MUD.data_fact, taxonomic.rank="Phylum"))

sample_names(MUD.data_fact.phylum.prop)

trans_positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")
fung.clayT10.prop = prune_taxa(TopPHYL, MUD.data_fact.phylum.prop)
(p_fung_T10=plot_bar(fung.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Spectral")+theme_bw()+theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                                                             axis.title=element_text(size=20),panel.grid.major=element_blank(),
                                                             panel.grid.minor=element_blank())+
    scale_x_discrete(limits = trans_positions))


#VST
MUD.data.ord <- ordinate(MUD.data, method="NMDS",distance = "bray")
#*** Solution reached
#Stress:     0.1418 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'wisconsin(veganifyOTU(physeq))' 

plot_ordination(MUD.data, MUD.data.ord)
MUD.data_map=sample_data(MUD.data)
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
MUD.data.dist<-distance(MUD.data,method = "bray")
MUD.data_map=sample_data(MUD.data)
test_perm <- adonis (MUD.data.dist ~ MUD.data_map$Location + MUD.data_map$Spp, permutations=10000)
print(test_perm)
# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(MUD.data.dist, MUD.data_map$Location, nperm=2000)
# RESULT: all locations different than one another

# This code compares species
pairwise.perm.manova(MUD.data.dist, MUD.data_map$Spp, nperm=2000)


mmm1 <- metaMDS(MUD.data.dist, k = 3, trymax = 500,
              autotransform =FALSE,  
              wascores = TRUE, expand = TRUE,
              trace = 1, plot = TRUE)
#0.133854

#let run the simper

#need the OTU table
MUD.data_otu=t(otu_table(MUD.data))

colnames(MUD.data_otu)
row.names(MUD.data_otu)

#I also need the taxonomy table for extracting the classified species identity

MUD.data_taxa=as.data.frame(tax_table(MUD.data))

#let load in the funguild results
funguild.80c=read.delim("MUD_fungi_OTU_ITS_VST_80c.guilds.txt",sep = "\t")
head(funguild.80c)
row.names(funguild.80c)=funguild.80c$X
funguild.80c_names=funguild.80c[92:ncol(funguild.80c)]


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
write.csv(MUD.data_LATR.simp_funguild_mat_sig, "MUD_sig_simper_fungi_80c_LATR_Eco_V_Shrub.csv")

#To start out, let's take the the top 10 that explain the most variation

MUD.data_LATR.simp_taxa_mat_sig_top10=(MUD.data_LATR.simp_taxa_mat_sig[order(-MUD.data_LATR.simp_taxa_mat_sig$average),])[1:10,]
rownames(MUD.data_LATR.simp_taxa_mat_sig_top10)





#Now I need to extract out the samples that have Ecotone

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
#504
MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig=merge(MUD.data_map_Ecotone.simp_taxa_BOER_LATR_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig, "MUD_sig_simper_fungi_80c_Ecotone_BEOR_LATR.csv")

#To start out, let's take the the top 10 that explain the most variation

MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig_top10=(MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig[order(-MUD.data_map_Ecotone.simp_funguild_BOER_LATR_sig$average),])[1:10,]
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
#449
MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig=merge(MUD.data_map_Ecotone.simp_taxa_LATR_PLJA_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_map_Ecotone.simp_funguild_LATR_PLJA_sig, "MUD_sig_simper_fungi_80c_Ecotone_LATR_PLJA.csv")



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
MUD.data_BOER.simp_funguild_sig=merge(MUD.data_BOER.simp_taxa_sig,funguild.80c_names, by="row.names")
write.csv(MUD.data_BOER.simp_funguild_sig, "MUD_sig_simper_fungi_80c_BOER_Eco_V_Grassland.csv")

#Now I need to extract out the samples that have PLJA

MUD.data_map=sample_data(MUD.data)
MUD.data_map_PLJA=subset(MUD.data_map, Spp=="PLJA")
nrow(MUD.data_map_PLJA)
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
write.csv(MUD.data_PLJA.simp_funguild_sig, "MUD_sig_simper_fungi_80c_PLJA_Eco_V_Grassland.csv")



#To start out, let's take the the top 10 that explain the most variation

MUD.data_LATR.simp_taxa_mat_sig_top10=(MUD.data_LATR.simp_taxa_mat_sig[order(-MUD.data_LATR.simp_taxa_mat_sig$average),])[1:10,]
rownames(MUD.data_LATR.simp_taxa_mat_sig_top10)


############### MORE SOIL CHARACTERISTICS ##################################
MUD.data_map=sample_data(MUD.data)
### Want to compare the species matrix to an environmental matrix
soil_list=c("SAND","SILT","CLAY","pH","P_ppm","K_ppm","OM_percent","NO3_N_ppm")

soil <- MUD.data_map[,soil_list]
str(soil)
soil=as.matrix(soil)
soil=as.data.frame(soil)



(vec <-envfit(MUD.data.ord, soil, perm=9999, na.rm=TRUE)) #I think i added the environmental data
vec.df <-as.data.frame(vec$vectors$arrows*sqrt(vec$vectors$r))
vec.df$huh <-rownames(vec.df)
scores(vec, "vectors")
summary(vec)

(vec2 <-envfit(mmm1, soil.1, perm=9999, na.rm=TRUE)) #I think i added the environmental data
##Followed this info below for graphing arrows: https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
## PROBLEM: this drops the catagorical traits from the analysis/graph and I don't know why
plot_ordination(MUD.data, MUD.data.ord)+
  geom_point(size=3, aes(color=Site,shape=Species))+
  scale_color_brewer(palette = "Dark2",labels=c("Ecotone",
                                                "Grassland",
                                                "Shrub"))+
  scale_shape_manual(values=c(16,17,18,15),name="",labels=c("BEOR","BOGR","LATR","PLJA"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_segment(data=vec.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length=unit(0.5, "cm")), colour="grey")+
  geom_text_repel (data=vec.df, aes(x=NMDS1, y=NMDS2, label=huh), size=5)


(fit <- envfit(MUD.data.ord, soil, perm = 999, na.rm=TRUE))
scores(fit, "vectors")
plot(m1)
plot(fit)
plot(fit, p.max = 0.05, col = "red")
summary(fit)




eea_list <- c("NAG_soil", "NAG_OM", "AlkP_soil","AlkP_OM","AAP_soil","AAP_OM", "BG_soil", "BG_OM")
eea <- MUD.data_map[,eea_list]
str(eea)
eea$AAP_soil <- as.numeric(levels(eea$AAP_soil)) [eea$AAP_soil] # changing factor to numeric
eea$AAP_OM <- as.numeric(levels(eea$AAP_OM)) [eea$AAP_OM] # changing factor to numeric
str(eea)
eea=as.matrix(eea)
eea=as.data.frame(eea)
str(eea)





(eea_vec <-envfit(MUD.data.ord, eea, perm=9999, na.rm=TRUE)) #I think i added the environmental data
eea_vec.df <-as.data.frame(eea_vec$vectors$arrows*sqrt(eea_vec$vectors$r))
eea_vec.df$huh <-rownames(eea_vec.df)
scores(eea_vec, "vectors")
summary(eea_vec)

plot_ordination(MUD.data, MUD.data.ord)+
  geom_point(size=3, aes(color=Site,shape=Species))+
  scale_color_brewer(palette = "Dark2",labels=c("Ecotone",
                                                "Grassland",
                                                "Shrub"))+
  scale_shape_manual(values=c(16,17,18,15),name="",labels=c("BEOR","BOGR","LATR","PLJA"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_segment(data=eea_vec.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
               arrow = arrow(length=unit(0.5, "cm")), colour="grey")+
  geom_text_repel (data=eea_vec.df, aes(x=NMDS1, y=NMDS2, label=huh), size=5)


#let's run a mantel test between the soil EEA and community dist
#first let's normalize so values are between 1 and 0


eea.nom=decostand(eea, "range",na.rm=T)


eea.dist=vegdist(eea.nom,method="euclidean")

mantel(MUD.data.dist,eea.dist)
#Mantel statistic r: 0.2045 
#Significance: 0.001 

#EEA SOils
eea_soil <- c("NAG_soil", "AlkP_soil","AAP_soil","BG_soil")
eea.nom_soil=eea.nom[,eea_soil]
trt_map=MUD.data_map[,c("Site","Species")]
eea.dist_soil=vegdist(eea.nom_soil,method="euclidean")

eea.nom_soil_trt=merge(eea.nom_soil,trt_map,by="row.names")
eea.nom_soil.ord=metaMDS(eea.dist_soil)

eea.nom_soi_ord_points=merge(eea.nom_soil.ord$points,MUD.data_map,by="row.names")
colnames(eea.nom_soi_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(MUD.data_map$Site)
ggplot(data=eea.nom_soi_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())


eea_soil=eea[,eea_soil]
trt_map=MUD.data_map[,c("Site","Species")]
eea.dist_soil=vegdist(eea_soil,method="euclidean")

eea.nom_soil_trt=merge(eea.nom_soil,trt_map,by="row.names")
eea.nom_soil.ord=metaMDS(eea.dist_soil)

eea.nom_soi_ord_points=merge(eea.nom_soil.ord$points,MUD.data_map,by="row.names")
colnames(eea.nom_soi_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(MUD.data_map$Site)
ggplot(data=eea.nom_soi_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())

test_perm <- adonis (eea.dist_soil ~ MUD.data_map$Location + MUD.data_map$Spp, permutations=10000)
print(test_perm)
# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(eea.dist_soil, MUD.data_map$Location, nperm=2000)



#OM 
eea_OM_names <- c("NAG_OM", "AlkP_OM","AAP_OM","BG_OM")
eea_OM=eea[,eea_OM_names]
trt_map=MUD.data_map[,c("Site","Species")]
eea.dist_OM=vegdist(eea_OM,method="euclidean")

eea.nom_OM_trt=merge(eea_OM,trt_map,by="row.names")
eea.nom_OM.ord=metaMDS(eea.dist_OM)

eea.nom_OM_ord_points=merge(eea.nom_OM.ord$points,MUD.data_map,by="row.names")
colnames(eea.nom_OM_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(MUD.data_map$Site)
ggplot(data=eea.nom_OM_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())


test_perm <- adonis (eea.dist_OM ~ MUD.data_map$Location + MUD.data_map$Spp, permutations=10000)
print(test_perm)
# I *think* this code is a quasi post-hoc test that will compare the locations
#pairwise.perm.manova(MUD.data.dist, MUD.data_map$Location, nperm=2000)

#let's run a mantel test between the soil nutrients and community dist
#first let's normalize so values are between 1 and 0



char.nom=decostand(soil, "range",na.rm=T)
char.dist=vegdist(char.nom,method="euclidean")
mantel(MUD.data.dist,char.dist)
#Mantel statistic r: 0.3958 
#Significance: 0.001 

#finally the soil nut versus soil eea 

mantel(eea.dist,char.dist)
#Mantel statistic r: 0.138 
#Significance: 0.029 



char.dist=vegdist(char.dist,method="euclidean")

char.nom_trt=merge(eea_OM,trt_map,by="row.names")
char.nom.ord=metaMDS(char.dist)

char.nom_ord_points=merge(char.nom.ord$points,MUD.data_map,by="row.names")
colnames(char.nom_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(MUD.data_map$Site)
ggplot(data=char.nom_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())


test_perm <- adonis (char.dist ~ MUD.data_map$Location + MUD.data_map$Spp, permutations=10000)
print(test_perm)
# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(char.dist, MUD.data_map$Location, nperm=2000)

pairwise.perm.manova(char.dist, MUD.data_map$Species, nperm=2000)

#let try a partial mantel 

mantel.partial(MUD.data.dist,eea.dist,char.dist)
#Mantel statistic r: -0.05462 
#Significance: 0.805 

#all rank Taxon
MUD.fung.tax_allRank=as.matrix(read.table("its_MUD_sintax_all_rank_simple.txt",header=T))
MUD.fung.TAX_allRank = tax_table(MUD.fung.tax_allRank)
MUD.data_allRank=phyloseq(otu_table(MUD.data),sample_data(MUD.data),MUD.fung.TAX_allRank)

MUD.data_fact_allRank=merge_samples(MUD.data_allRank, "site_spp")
sample_names(MUD.data_fact_allRank)     

get_taxa_unique(MUD.data_fact_allRank, taxonomic.rank="Phylum")
#12
(MUD.data_fact_allRank.phylum<-tax_glom(MUD.data_fact_allRank, taxrank="Phylum"))




TopPHYL_allRank = names(sort(taxa_sums(MUD.data_fact_allRank.phylum), TRUE)[1:10])
mud_fung_allRank.T10 = prune_taxa(TopPHYL_allRank, MUD.data_fact_allRank.phylum)

plot_bar(mud_fung_allRank.T10, x= "site_spp", fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




MUD.data_fact_allRank.phylum.prop=transform_sample_counts(MUD.data_fact_allRank.phylum, function(x)x/sum(x))
sort(get_taxa_unique(MUD.data_fact_allRank, taxonomic.rank="Phylum"))

sample_names(MUD.data_fact_allRank.phylum.prop)

trans_positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")
fung_allRank.clayT10.prop = prune_taxa(TopPHYL_allRank, MUD.data_fact_allRank.phylum.prop)
(p_fung_T10_allRank=plot_bar(fung_allRank.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Spectral")+theme_bw()+theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                                                             axis.title=element_text(size=20),panel.grid.major=element_blank(),
                                                             panel.grid.minor=element_blank())+
    scale_x_discrete(limits = trans_positions))


#untranformed for the analyses
untrans.otu=read.table("MUD_fungi_OTU_ITS_trunc_phyl.txt")
untrans.MUD.OTU = otu_table(untrans.otu, taxa_are_rows = TRUE)
MUD.data_untrans=phyloseq(untrans.MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX)


Rich_meas = c("Observed", "Chao1")
untrans.MUD.data.divfil=estimate_richness(MUD.data_untrans,measures=Rich_meas)
untrans.MUD.data_map=sample_data(MUD.data_untrans)

untrans.MUD.data.divfil=merge(untrans.MUD.data.divfil, untrans.MUD.data_map, by ="row.names")


#Observed Richness
Observed_R_mod <- lm((Observed) ~ Species + Site, data=untrans.MUD.data.divfil)
qqPlot(stdres(Observed_R_mod))
hist(stdres(Observed_R_mod))
shapiro.test(stdres(Observed_R_mod))
#0.8072
summary(Observed_R_mod)
Anova(Observed_R_mod, type=3)
#nada sig
#emmeans(Shannon_mod, pairwise~Site)


positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Observed_R_g=ggplot(untrans.MUD.data.divfil, aes(x=Site, y=Observed),
                 fill=Species)

(Observed_R_box_p=Observed_R_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Observed_R_g$data, aes(x=Site, y=Observed,
                                          fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Observed raw richness")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                                                                                labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))


#Chao1
Chao1_mod <- lm((Chao1) ~ Species + Site, data=untrans.MUD.data.divfil)
qqPlot(stdres(Chao1_mod))
hist(stdres(Chao1_mod))
shapiro.test(stdres(Chao1_mod))
#0.9112
summary(Chao1_mod)
Anova(Chao1_mod, type=3)
#nada sig
#emmeans(Shannon_mod, pairwise~Site)


positions2=c("G","E","S")
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
Chao1_R_g=ggplot(untrans.MUD.data.divfil, aes(x=Site, y=Chao1),
                    fill=Species)

(Chao1_R_box_p=Chao1_R_g+stat_boxplot(geom = "errorbar", aes(color=factor(Species, levels=spp_pos)),position = position_dodge2(width = 0.5, padding = 0.5,preserve = "single"))+
    geom_boxplot(data=Chao1_R_g$data, aes(x=Site, y=Chao1,
                                             fill=factor(Species, levels=spp_pos)),
                 outlier.shape = 19, outlier.size = 2.5,position = position_dodge2(preserve = "single") )+
    scale_x_discrete(limits = positions2,labels= c("Grassland","Ecotone","Shrubland"))+scale_colour_manual(values=c("black","black","black","black"),
                                                                                                           labels= c("BOGR","BOER","PLJA","LATR"))+
    scale_y_continuous(name = "Chao1 Richness")+xlab(NULL)+scale_fill_manual(limits = spp_pos, 
                                                                                    values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                                                                                    labels= c("BOGR","BOER","PLJA","LATR"))+
    theme_bw()+theme(legend.title = element_blank(), legend.text=element_text(size=16), axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title=element_text(size=20),panel.grid.major=element_blank(),panel.grid.minor=element_blank()))

#presence absence
MUD.data_pa=transform_sample_counts(MUD.data_untrans,pa)
max(otu_table(MUD.data_pa))

MUD.data_pa.ord <- ordinate(MUD.data_pa, method="NMDS",distance = "jaccard")
#*** Solution reached
#Dimensions: 2 
#Stress:     0.1793214 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'veganifyOTU(physeq)' 

plot_ordination(MUD.data_pa, MUD.data_pa.ord)
untrans.MUD.data_map=sample_data(MUD.data_untrans)
MUD.data_pa_ord_points=merge(MUD.data_pa.ord$points,untrans.MUD.data_map,by="row.names")
colnames(MUD.data_pa_ord_points)
spp_pos=c("Bogr","Boer",  "Plja", "Latr")
positions2=c("G","E","S")

unique(untrans.MUD.data_map$Site)
ggplot(data=MUD.data_pa_ord_points, aes(x=MDS1,y=MDS2))+
  geom_point(size=3, aes(color=factor(Species, levels=spp_pos),shape=factor(Site, levels=positions2)))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     name="",labels= c("BOGR","BOER","PLJA","LATR"))+
  scale_shape_manual(values=c(16,17,15),name="",labels=c("Grassland","Ecotone","Shrubland"))+
  theme_bw()+theme(axis.text=element_text(size=10),legend.text = element_text(size = 14),
                   axis.title.x=element_text(size=12),axis.title.y=element_blank(),
                   panel.grid.major=element_blank(),panel.grid.minor=element_blank())
MUD.data_pa.dist<-distance(MUD.data_pa,method = "jaccard")
MUD.data_pa_map=sample_data(MUD.data_pa)
test_perm.pa <- adonis (MUD.data_pa.dist ~ MUD.data_pa_map$Location + MUD.data_pa_map$Spp, permutations=10000)
print(test_perm.pa)
# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(MUD.data_pa.dist, MUD.data_pa_map$Location, nperm=2000)
# RESULT: all locations different than one another

# This code compares species
pairwise.perm.manova(MUD.data_pa.dist, MUD.data_pa_map$Spp, nperm=2000)



#Funguild seems like an interesting analysis 
head(MUD.map_eea_char$sampleID)
length(row.names(MUD.map_eea_char))
funguild.allRank=read.delim("MUD_fungi_OTU_ITS_VST_allRank.guilds.txt",sep = "\t")
head(colnames(funguild.allRank))
head(row.names(funguild.allRank))
names(funguild.allRank) <- gsub(" ", "_", names(funguild.allRank))
funguild.allRank=funguild.allRank %>% group_by(Guild)
funguild.allRank_sum=summarise_at(funguild.allRank,vars(row.names(MUD.map_eea_char)),sum)
row.names(funguild.allRank_sum)=funguild.allRank_sum$Guild
head(colnames(funguild.allRank_sum))
head(row.names(funguild.allRank_sum))
length(row.names(funguild.allRank_sum))
colnames(funguild.allRank_sum[,2:57] )
row.names(sort(rowSums(funguild.allRank_sum[,2:57] )))
row.names(funguild.allRank_sum)
funguild.allRank_sum[,1]

length(row.names(funguild.allRank_sum))
funguild.allRank_sum.1=mutate(funguild.allRank_sum, total_reads=rowSums(funguild.allRank_sum[,2:57] ))
row.names(funguild.allRank_sum.1)=funguild.allRank_sum.1$Guild
row.names(funguild.allRank_sum.1)
funguild.allRank_sum.1[,"total_reads"]

newdata <- funguild.allRank_sum.1[order(funguild.allRank_sum.1$total_reads),]
tail(newdata[,c("Guild","total_reads")])


"# A tibble: 6 x 2
  Guild                  total_reads
  <fct>                        <dbl>
1 Arbuscular Mycorrhizal       6455.
2 Lichenized                   9051.
3 Ectomycorrhizal              9732.
4 Plant Pathogen              11622.
5 Undefined Saprotroph        20659.
6 -                           29076."

newdata2=as.data.frame(newdata[,c("Guild","total_reads")])
"79                                                                                       Endophyte  2946.438332
80                                                                                 Wood Saprotroph  4229.818863
81                                                                          Arbuscular Mycorrhizal  6455.031877
82                                                                                      Lichenized  9051.089935
83                                                                                 Ectomycorrhizal  9731.956105
84                                                                                  Plant Pathogen 11622.078225
85                                                                            Undefined Saprotroph 20659.044194
86                                                                                               - 29075.937469"

funguild.allRank_sum[,"Guild"]=NULL
funguild.allRank_sum_t=t(funguild.allRank_sum)
row.names(funguild.allRank_sum_t)
colnames(funguild.allRank_sum_t)

funguild.allRank_sum_T=merge(MUD.map_eea_char,funguild.allRank_sum_t, by="row.names")

funguild.allRank_sum_T=funguild.allRank_sum_T %>% group_by(Site,Species)
funguild.allRank_sum_T_sum=summarise_at(funguild.allRank_sum_T, vars("Undefined Saprotroph","Plant Pathogen","Ectomycorrhizal","Lichenized","Arbuscular Mycorrhizal","Wood Saprotroph",Endophyte),funs(mean,se=sd(.)/sqrt(n()),sd))
funguild.allRank_sum_T_sum$site_spp=with(funguild.allRank_sum_T_sum, interaction(Site,Species))
names(funguild.allRank_sum_T_sum) <- gsub(" ", "_", names(funguild.allRank_sum_T_sum))
positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")

Eco_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Ectomycorrhizal_mean, ymin = Ectomycorrhizal_mean-Ectomycorrhizal_se, ymax = Ectomycorrhizal_mean+Ectomycorrhizal_se))

(Eco_fun_p=Eco_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Ectomycorrhizae")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


AMF_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Arbuscular_Mycorrhizal_mean, 
                                                 ymin = Arbuscular_Mycorrhizal_mean-Arbuscular_Mycorrhizal_se,
                                                 ymax = Arbuscular_Mycorrhizal_mean+Arbuscular_Mycorrhizal_se))

(AMF_fun_p=AMF_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Arbuscular mycorrhizae")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))

UnSap_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Undefined_Saprotroph_mean, 
                                                 ymin = Undefined_Saprotroph_mean-Undefined_Saprotroph_se,
                                                 ymax = Undefined_Saprotroph_mean+Undefined_Saprotroph_se))

(UnSap_fun_p=UnSap_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Unidentified Saprotroph")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


PlantPath_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Plant_Pathogen_mean, 
                                                   ymin = Plant_Pathogen_mean-Plant_Pathogen_se,
                                                   ymax = Plant_Pathogen_mean+Plant_Pathogen_se))

(PlantPath_fun_p=PlantPath_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Plant Pathogen")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))

#Lichenized

Lichen_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Lichenized_mean, 
                                                       ymin = Lichenized_mean-Lichenized_se,
                                                       ymax = Lichenized_mean+Lichenized_se))

(Lichen_fun_p=Lichen_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Lichenized")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))

#Wood_Saprotroph

WSap_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Wood_Saprotroph_mean, 
                                                    ymin = Wood_Saprotroph_mean-Wood_Saprotroph_se,
                                                    ymax = Wood_Saprotroph_mean+Wood_Saprotroph_se))

(WSap_fun_p=WSap_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Wood Saprotroph")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


#Endophyte
Endo_fun_g=ggplot(funguild.allRank_sum_T_sum, aes(site_spp, Endophyte_mean, 
                                                  ymin = Endophyte_mean-Endophyte_se,
                                                  ymax = Endophyte_mean+Endophyte_se))

(Endo_fun_p=Endo_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Endophyte")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


colnames(funguild.allRank)
funguild.allRank=funguild.allRank %>% group_by(Growth.Morphology)
funguild.allRank_sumG=summarise_at(funguild.allRank,vars(row.names(MUD.map_eea_char)),sum)
row.names(funguild.allRank_sumG)=funguild.allRank_sumG$Growth.Morphology
head(colnames(funguild.allRank_sumG))
head(row.names(funguild.allRank_sumG))
length(row.names(funguild.allRank_sumG))
colnames(funguild.allRank_sumG[,2:57] )
row.names(funguild.allRank_sumG)
funguild.allRank_sumG[,1]

length(row.names(funguild.allRank_sumG))
funguild.allRank_sumG.1=mutate(funguild.allRank_sumG, total_reads=rowSums(funguild.allRank_sumG[,2:57] ))
row.names(funguild.allRank_sumG.1)=funguild.allRank_sumG.1$Growth.Morphology
row.names(funguild.allRank_sumG.1)
funguild.allRank_sumG.1[,"total_reads"]

newdataG <- funguild.allRank_sumG.1[order(funguild.allRank_sumG.1$total_reads),]
tail(newdataG[,c("Growth.Morphology","total_reads")])


"# A tibble: 6 x 2
  Growth.Morphology total_reads
<fct>                   <dbl>
1 Corticioid              4746.
2 Thallus                 9051.
3 Agaricoid               9074.
4 Microfungus            12763.
5 -                      29076.
6 NULL                   38313."

newdataG=as.data.frame(newdataG[,c("Growth.Morphology","total_reads")])
"21                   Dark Septate Endophyte  2787.863352
22                        Facultative Yeast  4574.220823
23                               Corticioid  4746.403880
24                                  Thallus  9051.089935
25                                Agaricoid  9074.490509
26                              Microfungus 12762.999394
27                                        - 29075.937469
28                                     NULL 38312.572352"

funguild.allRank_sumG[,"Growth.Morphology"]=NULL
funguild.allRank_sumG_t=t(funguild.allRank_sumG)
row.names(funguild.allRank_sumG_t)
colnames(funguild.allRank_sumG_t)

funguild.allRank_sumG_T=merge(MUD.map_eea_char,funguild.allRank_sumG_t, by="row.names")

funguild.allRank_sumG_T=funguild.allRank_sumG_T %>% group_by(Site,Species)
funguild.allRank_sumG_T_sum=summarise_at(funguild.allRank_sumG_T, vars("Microfungus","Agaricoid","Thallus",
                                                                       "Corticioid","Facultative Yeast","Dark Septate Endophyte"),
                                         funs(mean,se=sd(.)/sqrt(n()),sd))
funguild.allRank_sumG_T_sum$site_spp=with(funguild.allRank_sumG_T_sum, interaction(Site,Species))
names(funguild.allRank_sumG_T_sum) <- gsub(" ", "_", names(funguild.allRank_sumG_T_sum))
positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")

Dark_sept_fun_g=ggplot(funguild.allRank_sumG_T_sum, aes(site_spp, Dark_Septate_Endophyte_mean, 
                                                        ymin = Dark_Septate_Endophyte_mean-Dark_Septate_Endophyte_se, 
                                                        ymax = Dark_Septate_Endophyte_mean+Dark_Septate_Endophyte_se))

(Dark_sept_fun_p=Dark_sept_fun_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Dark septate endophytes")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))



#Most recent UNITEv7.2 1201017
MUD.fung.tax.1201017=as.matrix(read.table("its_MUD_sintax_simple_80_v72.1201017.txt",header=T))
MUD.fung.TAX.1201017 = tax_table(MUD.fung.tax)
MUD.data.1201017=phyloseq(MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX.1201017)


sum(otu_table(MUD.data.1201017))
#4238112
#120193.4
MUD.Fung=subset_taxa(MUD.data.1201017, Domain=="d:Fungi")
ntaxa(MUD.Fung)
#4951
sum(otu_table(MUD.Fung))
#4235030
#120193.4

sort(sample_sums(MUD.Fung))


MUD.Fung_fact=merge_samples(MUD.Fung, "site_spp")
sample_names(MUD.Fung_fact)     

get_taxa_unique(MUD.Fung_fact, taxonomic.rank="Phylum")
#12
(MUD.Fung_fact.phylum<-tax_glom(MUD.Fung_fact, taxrank="Phylum"))




TopPHYL.f = names(sort(taxa_sums(MUD.Fung_fact.phylum), TRUE)[1:10])
mud_fung.T10.f = prune_taxa(TopPHYL.f, MUD.data_fact.phylum)

plot_bar(mud_fung.T10.f, x= "site_spp", fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




MUD.data_fact.phylum.prop=transform_sample_counts(MUD.data_fact.phylum, function(x)x/sum(x))
sort(get_taxa_unique(MUD.data_fact, taxonomic.rank="Phylum"))

sample_names(MUD.data_fact.phylum.prop)

trans_positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")
fung.clayT10.prop = prune_taxa(TopPHYL, MUD.data_fact.phylum.prop)
(p_fung_T10=plot_bar(fung.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Spectral")+theme_bw()+theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                                                             axis.title=element_text(size=20),panel.grid.major=element_blank(),
                                                             panel.grid.minor=element_blank())+
    scale_x_discrete(limits = trans_positions))




#Simple 60c
MUD.fung.tax_60=as.matrix(read.table("ITS_MUD_sintax_simple_60.txt",header=T))
MUD.fung.TAX_60 = tax_table(MUD.fung.tax_60)
MUD.data_60=phyloseq(MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX_60)


sum(otu_table(MUD.data_60))
#4238112
#120193.4
MUD.Fung=subset_taxa(MUD.data_60, Domain=="d:Fungi")
ntaxa(MUD.Fung)
#4951
sum(otu_table(MUD.Fung))
#4235030
#120193.4

sort(sample_sums(MUD.Fung))


MUD.Fung_fact=merge_samples(MUD.Fung, "site_spp")
sample_names(MUD.Fung_fact)     

get_taxa_unique(MUD.Fung_fact, taxonomic.rank="Phylum")
#12
(MUD.Fung_fact.phylum<-tax_glom(MUD.Fung_fact, taxrank="Phylum"))




TopPHYL.f = names(sort(taxa_sums(MUD.Fung_fact.phylum), TRUE)[1:10])
mud_fung.T10.f = prune_taxa(TopPHYL.f, MUD.Fung_fact.phylum)

plot_bar(mud_fung.T10.f, x= "site_spp", fill="Phylum")+ 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")




MUD.data_fact.phylum.prop=transform_sample_counts(MUD.Fung_fact.phylum, function(x)x/sum(x))
sort(get_taxa_unique(MUD.data_fact.phylum.prop, taxonomic.rank="Phylum"))

sample_names(MUD.data_fact.phylum.prop)

trans_positions=c("G.Bogr","G.Boer", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")
fung.clayT10.prop = prune_taxa(TopPHYL.f, MUD.data_fact.phylum.prop)
(p_fung_T10=plot_bar(fung.clayT10.prop, fill="Phylum")+ylab("Proportion")+ 
    geom_bar(aes( fill=factor(Phylum)), stat="identity", position="stack",color="black")+xlab(NULL)+
    scale_fill_brewer(palette = "Spectral")+theme_bw()+theme(axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                                                             axis.title=element_text(size=20),panel.grid.major=element_blank(),
                                                             panel.grid.minor=element_blank())+
    scale_x_discrete(limits = trans_positions,labels= c("BOGR","BOER","PLJA","BOER","PLJA","LATR","LATR")))



# Online script to generate cohesion metrics for a set of samples 
# CMH 26Apr17; cherren@wisc.edu

# User instructions: read in a sample table (in absolute or relative abundance) as object "b".
# If using a custom correlation matrix, read in that matrix at the designated line.
# Run the entire script, and the 4 vectors (2 of connectedness and 2 of cohesion) are generated for each sample at the end.
# Parameters that can be adjusted include pers.cutoff (persistence cutoff for retaining taxa in analysis), iter (number of iterations for the null model), tax.shuffle (whether to use taxon shuffle or row shuffle randomization), and use.custom.cors (whether to use a pre-determined correlation matrix)

####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################

## Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis
pers.cutoff <- 0.10
## Decide the number of iterations to run for each taxon. (>= 200 is recommended)
# Larger values of iter mean the script takes longer to run
iter <- 100
## Decide whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
tax.shuffle <- T
## Option to input your own correlation table
# Note that your correlation table MUST have the same number of taxa as the abundance table. There should be no empty (all zero) taxon vectors in the abundance table. 
# Even if you input your own correlation table, the persistence cutoff will be applied
use.custom.cors <- F

###################################################################
###################################################################

# Read in dataset
## Data should be in a matrix where each row is a sample.
#MUD.data.otu<-otu_table(MUD.data)
#MUD.data.Otu<-as.matrix(MUD.data.otu)
#b<-t(as.data.frame(MUD.data.otu))

untrans.otu=t(read.table("MUD_fungi_OTU_ITS_trunc_phyl.txt"))
attributes(untrans.otu)
b<-(as.data.frame(untrans.otu))
row.names(b)
#b <- read.csv("your_path_here.csv", header = T, row.names = 1)

# Read in custom correlation matrix, if desired. Must set "use.custom.cors" to TRUE
if(use.custom.cors == T) {
  custom.cor.mat <- read.csv("your_path_here.csv", header = T, row.names = 1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #Check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2] == dim(custom.cor.mat)[2])
}


# Suggested steps to re-format data. At the end of these steps, the data should be in a matrix "c" where there are no empty samples or blank taxon columns. 
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]

# Optionally re-order dataset to be in chronological order. Change date format for your data. 
#c <- c[order(as.Date(rownames(c), format = "%m/%d/%Y")), ]

# Save total number of individuals in each sample in the original matrix. This will be 1 if data are in relative abundance, but not if matrix c is count data
rowsums.orig <- rowSums(c)

# Based on persistence cutoff, define a cutoff for the number of zeroes allowed in a taxon's distribution
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])

# Remove taxa that are below the persistence cutoff
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
# Remove any samples that no longer have any individuals, due to removing taxa
d <- d[rowSums(d) > 0, ]

#If using custom correlation matrix, need to remove rows/columns corresponding to the taxa below persistence cutoff
if(use.custom.cors == T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c, 2, zero) < (dim(c)[1]-zero.cutoff), apply(c, 2, zero) < (dim(c)[1]-zero.cutoff)]
}

# Create relative abundance matrix.  
rel.d <- d / rowsums.orig
# Optionally, check to see what proportion of the community is retained after cutting out taxa
hist(rowSums(rel.d))

# Create observed correlation matrix
cor.mat.true <- cor(rel.d)

# Create vector to hold median otu-otu correlations for initial otu
med.tax.cors <- vector()

# Run this loop for the null model to get expected pairwise correlations
# Bypass null model if the option to input custom correlation matrix is TRUE
if(use.custom.cors == F) {
  ifelse(tax.shuffle, {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create empty matrix of same dimension as rel.d
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #For each otu
        for(j in 1:dim(rel.d)[2]){ 
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        
        # Do not randomize focal column 
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } , {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          
          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
  )
}

"Error in ans[test & ok] <- rep(yes, length.out = length(ans))[test & ok] : 
  replacement has length zero
In addition: Warning message:
In rep(yes, length.out = length(ans)) :
  'x' is NULL so the result will be NULL"

# Save observed minus expected correlations. Use custom correlations if use.custom.cors = TRUE
ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)

diag(obs.exp.cors.mat) <- 0

#### 
#### Produce desired vectors of connectedness and cohesion 

# Calculate connectedness by averaging positive and negative observed - expected correlations
connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

# Calculate cohesion by multiplying the relative abundance dataset by associated connectedness
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

####
#### Combine vectors into one list and print 
output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")

print(output)

mud_untrans_cohension=cbind(output$`Negative Cohesion`, output$`Positive Cohesion`)
colnames(mud_untrans_cohension)=c("Neg_Cohesion","Pos_Cohesion")

write.csv(mud_untrans_cohension, "mud_untrans_cohension_pos_neg_200int.csv")

length(row.names(mud_untrans_cohension))
#need to fix the one misspelled row name
row.names(mud_untrans_cohension)[match("fungi.Shurb.LATR.2",row.names(mud_untrans_cohension))] ="fungi.Shrub.LATR.2"

MUD_untrans_data=merge(MUD.map_eea_char,mud_untrans_cohension, by= "row.names")
length(row.names(MUD.map_eea_char))
length(row.names(MUD_untrans_data))


MUD_untrans_data=MUD_untrans_data %>% group_by(Site,Species)
mud_coh_sum=summarise_at(MUD_untrans_data, vars("Neg_Cohesion","Pos_Cohesion"),funs(mean,se=sd(.)/sqrt(n()),sd))
mud_coh_sum$site_spp=with(mud_coh_sum, interaction(Site,Species))
positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")

Neg_coh_g=ggplot(mud_coh_sum, aes(site_spp, Neg_Cohesion_mean, ymin = Neg_Cohesion_mean-Neg_Cohesion_se, ymax = Neg_Cohesion_mean+Neg_Cohesion_se))

(Neg_coh_p=Neg_coh_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Negative Cohension")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


Neg_coh_mod <- lm(Neg_Cohesion ~ Species + Site, data=MUD_untrans_data)
qqPlot(stdres(Neg_coh_mod))
hist(stdres(Neg_coh_mod))
shapiro.test(stdres(Neg_coh_mod))
summary(Neg_coh_mod)
Anova(Neg_coh_mod, type=3)
#nada sig


Pos_coh_g=ggplot(mud_coh_sum, aes(site_spp, Pos_Cohesion_mean, ymin = Pos_Cohesion_mean-Pos_Cohesion_se, ymax = Pos_Cohesion_mean+Pos_Cohesion_se))

(Pos_coh_p=Pos_coh_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "Positive Cohension")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))


Pos_coh_mod <- lm((Pos_Cohesion)^3 ~ Species + Site, data=MUD_untrans_data)
qqPlot(stdres(Pos_coh_mod))
hist(stdres(Pos_coh_mod))
shapiro.test(stdres(Pos_coh_mod))
#0.05848
summary(Pos_coh_mod)
Anova(Pos_coh_mod, type=3)
#nada sig
emmeans(Pos_coh_mod, pairwise~Site)
