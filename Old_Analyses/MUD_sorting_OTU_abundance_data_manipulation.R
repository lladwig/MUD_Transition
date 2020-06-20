library("phyloseq")
library(seqinr)
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")

#meta-data
eea <-read.csv("EEA_MUD_transition.csv")
char <-read.csv("SoilCharacteristics_MUD_Transition.csv")
#OTU table that was VST normalized


MUD.Fung_vst=read.table("MUD_fungi_fungi_OTU_ITS_VST.txt", header = T,row.names = 1)
#there is a miss-spelled sample name
colnames(MUD.Fung_vst)[colnames(MUD.Fung_vst)=="fungi.Shurb.LATR.2"]="fungi.Shrub.LATR.2"
MUD.Fung_vst[1:10,1:10]
fun_otu <- MUD.Fung_vst
fun_otu[1:10,1:10]
min(rowSums(fun_otu))
#Mapping file
trans_map=read.csv("MUD_Transition_map.csv")
MUD.OTU = otu_table(fun_otu, taxa_are_rows = TRUE)
length(colnames(MUD.OTU))
#90
sort(colnames(MUD.OTU))
min(taxa_sums(MUD.OTU))
nrow(MUD.OTU)
ncol(MUD.OTU)

colnames(MUD.OTU)

#trimmed phylogeny 
MUD.fung.tax=as.matrix(read.table("MUD_fungi_fungi_taxon_ITS_trunc_phyl.txt",header=T))
nrow(MUD.fung.tax)
#3750
MUD.fung.TAX = tax_table(MUD.fung.tax)
head(eea)
length(row.names(eea))
#56
length(row.names(trans_map))
#90
MUD.map_eea=merge(eea, trans_map,by="site")
head(MUD.map_eea)
length(row.names(MUD.map_eea))
#56
head(char)
length(row.names(char))
#56
char[,c("Location","Spp","Rep")]=NULL
head(char)
length(row.names(char))
#56
MUD.map_eea_char=merge(char, MUD.map_eea,by="site")
head(MUD.map_eea_char)
length(row.names(MUD.map_eea_char))
#56
rownames(MUD.map_eea_char)=MUD.map_eea_char$sampleID
MUD.map_eea_char$site_spp=with(MUD.map_eea_char, interaction(Site,Species))
MUD.data=phyloseq(MUD.OTU,sample_data(MUD.map_eea_char),MUD.fung.TAX)


ntaxa(MUD.data)
#3750
sum(otu_table(MUD.data))
#170867.8

min(sample_sums(MUD.data))
#1012.763

max(sample_sums(MUD.data))
#5212.293

sort(rank(sample_sums(MUD.data)))

min(taxa_sums(MUD.data))

sort(sample_sums(MUD.data))
MUD.data<-prune_taxa(taxa_sums(MUD.data) > 0, MUD.data)
ntaxa(MUD.data)
#2758
sum(otu_table(MUD.data))
#170867.8

min(sample_sums(MUD.data))
#1012.763

max(sample_sums(MUD.data))
#5212.293

length(sample_sums(MUD.data))

MUD.OTU_sum<-data.frame(taxa_sums(MUD.data))
colnames(MUD.OTU_sum)="OTU_sum"
MUD.OTU_sum_taxa=data.frame(merge(MUD.OTU_sum,data.frame(tax_table(MUD.data)), by="row.names", all.x = T))
head(MUD.OTU_sum_taxa)
colnames(MUD.OTU_sum_taxa)[1]="OTU"


#Let's add in the funguild classifications
funguild.80c=read.delim("MUD_fungi_OTU_ITS_VST_80c.guilds.txt",sep = "\t")
head(funguild.80c)
#one sample name is misspelled
colnames(funguild.80c)[colnames(funguild.80c)=="fungi.Shurb.LATR.2"]="fungi.Shrub.LATR.2"
colnames(funguild.80c)

MUD.OTU_sum_taxa_funguild=data.frame(merge(MUD.OTU_sum_taxa,funguild.80c[,c(1,93:101)], by.x="OTU", by.y = "X",all.x = T))

MUD.OTU_sum_taxa_funguild_sort=MUD.OTU_sum_taxa[order(-MUD.OTU_sum_taxa_funguild$OTU_sum),]
head(MUD.OTU_sum_taxa_funguild_sort)

#Ouput this file for sharing with Lee

write.csv(MUD.OTU_sum_taxa_funguild_sort,"MUD_taxa_funguild_OTU_sum_desending_sort.csv")


#Subsetting the representative sequences 
unfilter_rep_set_seq<- read.fasta("uniques_fwd_reads_demux_nophix_cut_fil_otus.fa", as.string = TRUE, set.attributes = FALSE)
MUD.sort_rep_set_seq=unfilter_rep_set_seq[names(unfilter_rep_set_seq) %in% MUD.OTU_sum_taxa_funguild_sort$OTU]
MUD.sort_rep_set_seq[names(MUD.sort_rep_set_seq)]
length(MUD.sort_rep_set_seq)
#2758
MUD.sort_rep_set_seq=MUD.sort_rep_set_seq[order(names(MUD.sort_rep_set_seq))]
head(MUD.sort_rep_set_seq)
length(MUD.sort_rep_set_seq)
MUD.sort_rep_set_seq_sort=MUD.sort_rep_set_seq[order(-MUD.OTU_sum_taxa_funguild$OTU_sum)]
head(MUD.sort_rep_set_seq_sort)
write.fasta(sequences =MUD.sort_rep_set_seq_sort, names = names(MUD.sort_rep_set_seq_sort), file.out ="MUD_rep_set_desending_sort.fna")
