library("phyloseq")
library(seqinr)
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")

#OTU table that was VST normalized
MUD.Fung_vst=read.table("MUD_fungi_fungi_OTU_ITS_VST.txt", header = T,row.names = 1)

#there is a miss-spelled sample name in the OTU Table 
colnames(MUD.Fung_vst)[colnames(MUD.Fung_vst)=="fungi.Shurb.LATR.2"]="fungi.Shrub.LATR.2"

MUD.OTU = otu_table(MUD.Fung_vst, taxa_are_rows = TRUE)
nrow(MUD.OTU)
#3750
#taxonomy based off of UNITE
MUD.fung.tax=as.matrix(read.table("MUD_fungi_fungi_taxon_ITS_trunc_phyl.txt",header=T))
nrow(MUD.fung.tax)
#3750

MUD.phylo=phyloseq(MUD.OTU,MUD.fung.tax)
ntaxa(MUD.phylo)
#3750

MUD.OTU_sum<-data.frame(taxa_sums(MUD.phylo))
colnames(MUD.OTU_sum)="OTU_sum"
MUD.OTU_sum_taxa=data.frame(merge(MUD.OTU_sum,data.frame(MUD.fung.tax), by="row.names", all.x = T))
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
#3750
MUD.sort_rep_set_seq=MUD.sort_rep_set_seq[order(names(MUD.sort_rep_set_seq))]
head(MUD.sort_rep_set_seq)
length(MUD.sort_rep_set_seq)
MUD.sort_rep_set_seq_sort=MUD.sort_rep_set_seq[order(-MUD.OTU_sum_taxa_funguild$OTU_sum)]
head(MUD.sort_rep_set_seq_sort)
write.fasta(sequences =MUD.sort_rep_set_seq_sort, names = names(MUD.sort_rep_set_seq_sort), file.out ="MUD_rep_set_desending_sort.fna")
