#####Code for the processing of Fungal Community before analyses####

#Load Packages
library(phyloseq)
library(seqinr)
library(taxize)
library(tidyr)
library(stringr)
#Set this for querying data from the NCBI
#Sys.setenv(ENTREZ_KEY = )


#Set your working directory to where the github was cloned
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")

#####Processing conducted before by hand#####


#Load in OTU table


MUD.otu=read.table("USEARCHv11_files/rep_set_fwd_reads_demux_phix_filtered_fil_OTU_table.txt", header=T, row.names = 1)
colnames(MUD.otu)

#there is a sample that is mislabeled
colnames(MUD.otu)[colnames(MUD.otu)=="fungi.Shurb.LATR.2"]="fungi.Shrub.LATR.2"
head(row.names(MUD.otu))
#Convert to a phyloseq formate
MUD.OTU = otu_table(MUD.otu, taxa_are_rows = TRUE)

#Load in mapping file
#this file only has fungal samples in it
MUD.map=sample_data(read.csv("R_files/MUD_Transition_map.csv", row.names = "sampleID"))
head(MUD.map)
nrow(MUD.map)

#Load in taxon table
MUD.fung.taxa_raw_UNITE8= 
  read.table("USEARCHv11_files/ITS_fwd_reads_demux_phix_filtered_fil_OTU_tax_v8.sintax",sep = "\t")
nrow(MUD.fung.taxa_raw_UNITE8)
#4833
head(MUD.fung.taxa_raw_UNITE8)

#Format the taxon table for phyloseq
MUD.fung.taxa_80C_UNITE8=MUD.fung.taxa_raw_UNITE8[,c(1,4)]
head(MUD.fung.taxa_80C_UNITE8)

MUD.fung.taxa_80C_UNITE8_sep=MUD.fung.taxa_80C_UNITE8 %>% separate(V4, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ",")
row.names(MUD.fung.taxa_80C_UNITE8_sep)=MUD.fung.taxa_80C_UNITE8_sep$V1
MUD.fung.taxa_80C_UNITE8_sep$V1=NULL
MUD.fung.taxa_80C_UNITE8_sep[is.na(MUD.fung.taxa_80C_UNITE8_sep)] <- "UNKNOWN"
MUD.fung.taxa_80C_UNITE8_sep_mat=as.matrix(MUD.fung.taxa_80C_UNITE8_sep)
head(MUD.fung.taxa_80C_UNITE8_sep_mat)
TAXA_80C_UNITE8_MUD_all=tax_table(MUD.fung.taxa_80C_UNITE8_sep_mat)





#Combine the files to make the phyloseq object
MUD.fung=phyloseq(MUD.OTU,sample_data(MUD.map),TAXA_80C_UNITE8_MUD_all)
nsamples(MUD.fung)


MUD.FUNG<-prune_taxa(taxa_sums(MUD.fung) > 0, MUD.fung)
ntaxa(MUD.FUNG)
#4807

sum(otu_table(MUD.FUNG))
#4134944

#Filter out OTUs that are not classified to Fungi at >=80% confidence
MUD.Fung=subset_taxa(MUD.FUNG, Domain=="d:Fungi")
ntaxa(MUD.Fung)
#3762
sum(otu_table(MUD.Fung))
#3949775

sort(sample_sums(MUD.Fung))
#Min=7241
#Max=78152


#Let's make a rep set fasta for fungal dataset so we search can for non-fungal reads outside of R
rep_set.fung<- read.fasta(file = "USEARCHv11_files/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa", as.string = TRUE, set.attributes = FALSE)
head(rep_set.fung)

MUD_rep_set.fung_names=rep_set.fung[names(rep_set.fung) %in% row.names(otu_table(MUD.Fung))]
write.fasta(sequences =MUD_rep_set.fung_names, names = names(MUD_rep_set.fung_names), 
            file.out ="R_files/rep_seq_MUD.Fung_OTU.fna")


#####KRACKEN Bact Search####
#Using Galaxy online bioinformatics and Kraken to find and filter out Bacterial sequences
#Kracken was used with default settings 
#https://toolshed.g2.bx.psu.edu/view/devteam/kraken/aec58624706f

MUD_poss_bact=read.delim("Galaxy_generated_taxonomy/Galaxy9-[Kraken-filter_on_data_5].tabular",
                         sep = "\t", header = F)
head(MUD_poss_bact)
summary(MUD_poss_bact)
#Filter out OTUs that did not match to bacteria
MUD_poss_bact_hit_B=subset(MUD_poss_bact, V5!="P=0.000")
nrow(MUD_poss_bact_hit_B)



#I am going to export the rep set from this and blast them to double check that there are not fungi
rep_set_Kracken_bacteria=rep_set.fung[names(rep_set.fung) %in% MUD_poss_bact_hit_B[,"V2"]]
head(rep_set_Kracken_bacteria)
length(rep_set_Kracken_bacteria)
#2
write.fasta(sequences =rep_set_Kracken_bacteria, names = names(rep_set_Kracken_bacteria), file.out ="Galaxy_generated_taxonomy/MUD_Kranken_bacteria.fna")

#Top hits for these taxa were 

#OTU4472  Darksidea sp. isolate DS913 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and large subunit ribosomal RNA gene, partial sequence	174	174	62%	7e-40	100.00%	MK808996.1
#OTU4819  Pezizales sp. isolate CL40PH internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence	167	167	62%	1e-37	98.92%	KU612291.1


#now I want to look at the taxa classification that were id as bacteria using KRACKEN
head(MUD_poss_bact_hit_B)
MUD.Fung_bact=prune_taxa(as.character(MUD_poss_bact_hit_B$V2), MUD.Fung)
tax_table(MUD.Fung_bact)
#TOTU4819 "d:Fungi" "p:Ascomycota" "c:Pezizomycetes" "o:Pezizales" "f:Ascobolaceae" "UNKNOWN" "UNKNOWN"
#OTU4472 "d:Fungi" "UNKNOWN"      "UNKNOWN"         "UNKNOWN"     "UNKNOWN"        "UNKNOWN" "UNKNOWN"

#Now I want to export unknown sequences and blast them to see if there are any crappy sequences or lingering bacteria
MUD.Fung_unk=subset_taxa(MUD.Fung, Phylum=="UNKNOWN")
ntaxa(MUD.Fung_unk)
#938
sum(otu_table(MUD.Fung_unk))
#716934


#I am going to export the rep set from this and blast them 
rep_set_unknown=rep_set.fung[names(rep_set.fung) %in% taxa_names(MUD.Fung_unk)]
head(rep_set_unknown)
length(rep_set_unknown)
#938
write.fasta(sequences =rep_set_unknown, names = names(rep_set_unknown), file.out ="UNKNOWN_taxa/MUD_unknown_phyla_seq.fna")


#####NCBI UNKNOWN OTUs Search####
#Representative sequences were BLAST against NCBI nucleotide DB with default settings
#Now load in the NCBI classifications to look for any OTUs that did not blast to anything or were more likley bacteria
#https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome



#Read in hit table 


ncbi_UNKNOWN=read.csv("UNKNOWN_taxa/U4MJ00GY014-Alignment-HitTable.csv",
                      header = F)
head(ncbi_UNKNOWN)
nrow(ncbi_UNKNOWN)
#32212
ncbi_UNKNOWN[,c(1,2)]

#We only care about the top value so I need to remove the duplicate hits
ncbi_UNKNOWN_uq <- ncbi_UNKNOWN[!duplicated(ncbi_UNKNOWN[,"V1"]),]
nrow(ncbi_UNKNOWN_uq)
#800


#Now we need to download the classifications from NCBI usingthe IDs
for (i in seq(nrow(ncbi_UNKNOWN_uq))) {
  ncbi_UNKNOWN_uq$uid[i]=genbank2uid(id = (as.character(ncbi_UNKNOWN_uq[i,2])))
  
  
}
classification((ncbi_UNKNOWN_uq$uid[[1]]))
attributes(ncbi_UNKNOWN_uq$uid[[1]])



ncbi_UNKNOWN_uq_taxID=genbank2uid(id = c(as.character(ncbi_UNKNOWN_uq[,2])))
length(ncbi_UNKNOWN_uq_taxID)
ncbi_UNKNOWN_uq_taxID[[1]]
attributes(ncbi_UNKNOWN_uq_taxID[[1]])
classification(ncbi_UNKNOWN_uq_taxID[[1]])
length(ncbi_UNKNOWN_uq_taxID)



list_new=c()
for (i in seq(length((ncbi_UNKNOWN_uq$uid)))) {
  list_new[i]=classification((ncbi_UNKNOWN_uq$uid[[i]]))
  for (j in seq(length(list_new))) {
    ncbi_UNKNOWN_uq$King[i] =list_new[[j]][2,1]
    
  }
  
}


ncbi_UNKNOWN_uq$King


#subset to only Eukaryota
ncbi_UNKNOWN_uq_Euk=subset(ncbi_UNKNOWN_uq,King=="Eukaryota")
nrow(ncbi_UNKNOWN_uq_Euk)
#800

#Now I want to remove any reads that blasted as non fungi
head(ncbi_UNKNOWN_uq_Euk)

list_new=c()
for (i in seq(length((ncbi_UNKNOWN_uq_Euk$uid)))) {
  list_new[i]=classification((ncbi_UNKNOWN_uq_Euk$uid[[i]]))
  for (j in seq(length(list_new))) {
    ncbi_UNKNOWN_uq_Euk$King_2[i] =list_new[[j]][4,1]
    
  }
  
}

unique(ncbi_UNKNOWN_uq_Euk$King_2)
#Subset taxa to only include fungi
ncbi_UNKNOWN_uq_Euk_fungi=subset(ncbi_UNKNOWN_uq_Euk,King_2=="Fungi")
nrow(ncbi_UNKNOWN_uq_Euk_fungi)
#800


#Now create a list of OTUs that need to be removed from the dataset
unk_but_fungi=ncbi_UNKNOWN_uq_Euk_fungi$V1

MUD.Fung_unk_OTU_names=taxa_names(MUD.Fung_unk)
length(MUD.Fung_unk_OTU_names)
#938
MUD.Fung_tax_names_unk_not_fungi<-MUD.Fung_unk_OTU_names[MUD.Fung_unk_OTU_names %w/o% unk_but_fungi]
length(MUD.Fung_tax_names_unk_not_fungi)
#138

MUD.Fung_tax_names=taxa_names(MUD.Fung)
length(MUD.Fung_tax_names)
#3762

MUD.Fung_no_bact_tax_names_fungi<-MUD.Fung_tax_names[MUD.Fung_tax_names %w/o% MUD.Fung_tax_names_unk_not_fungi]
length(MUD.Fung_no_bact_tax_names_fungi)
#3624
save(MUD.Fung_no_bact_tax_names_fungi, file = "MUD_Transition/UNKNOWN_taxa/List_NCBI_fungi_OTUs.Rdata")
#####NCBI UNKNOWN OTUs Search####

#####Filtering out non-fungi from the Phyloseq Object####
load(file = "UNKNOWN_taxa/List_NCBI_fungi_OTUs.Rdata")


#Now filter out the taxa from the phyloseq object
MUD.Fung_fungi=prune_taxa(MUD.Fung_no_bact_tax_names_fungi, MUD.Fung)


ntaxa(MUD.Fung_fungi)
#3624
sum(otu_table(MUD.Fung_fungi))
#3931394

sort(sample_sums(MUD.Fung_fungi))
#Min=7229 
#Max=78152

sum(otu_table(MUD.Fung_fungi))
#3931394



#####Filtering out non-fungi from the Phyloseq Object####


#Output the newly formatted OTU table and taxon table
write.table(tax_table(MUD.Fung_fungi), "R_files/MUD_fungi_fungi_taxon_ITS_trunc_phyl.txt") 
write.table(otu_table(MUD.Fung_fungi), "R_files/MUD_fungi_fungi_OTU_ITS_trunc_phyl.txt") 

#####Creation of OTU sums####
#Let's first remove the mock community and SPEGAR (an unrelated project)

MUD.Fung_fungi_map=sample_data(MUD.Fung_fungi)
head(MUD.Fung_fungi_map)
unique(MUD.Fung_fungi_map$site)
nrow(MUD.Fung_fungi_map)
MUD.Fung_fungi_field=subset_samples(MUD.Fung_fungi, project!="mock_com"&project!="spegar")
unique(sample_data(MUD.Fung_fungi_field)$project)

nrow(sample_data(MUD.Fung_fungi_field))
#56
MUD.Fung_fungi_field<-prune_taxa(taxa_sums(MUD.Fung_fungi_field) > 0, MUD.Fung_fungi_field)
ntaxa(MUD.Fung_fungi_field)
#2674
sum(otu_table(MUD.Fung_fungi_field))
#2321517
#let's look at the top taxa

MUD.OTU_sum<-data.frame(taxa_sums(MUD.Fung_fungi_field))
colnames(MUD.OTU_sum)="total_OTU_sum"
sum(MUD.OTU_sum)

#####OTU sums for each ecotone species combination####

head(sample_data(MUD.Fung_fungi_field))
unique(sample_data(MUD.Fung_fungi_field)$Spp)
unique(sample_data(MUD.Fung_fungi_field)$Location)
MUD.OTU_spp_loc_sum<-data.frame(Grass_BOGR_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="BOGR"&Location=="Grass")),
                                Grass_BOER_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="BOER"&Location=="Grass")),
                                Grass_PLJA_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="PLJA"&Location=="Grass")),
                                Ecotone_BOER_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="BOER"&Location=="Ecotone")),
                                Ecotone_PLJA_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="PLJA"&Location=="Ecotone")),
                                Ecotone_LATR_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="LATR"&Location=="Ecotone")),
                                Shrub_LATR_OTU_sum=taxa_sums(subset_samples(MUD.Fung_fungi_field, Spp=="LATR"&Location=="Shrub")))
sum(MUD.OTU_spp_loc_sum)


MUD.OTU_sum_comb=data.frame(merge(MUD.OTU_sum,MUD.OTU_spp_loc_sum, by="row.names", all.x = T))
head(MUD.OTU_sum_comb)

MUD.OTU_sum_comb_taxa=data.frame(merge(MUD.OTU_sum_comb,data.frame(tax_table(MUD.Fung_fungi_field)),by.x = "Row.names" ,by.y="row.names", all.x = T))
head(MUD.OTU_sum_comb_taxa)
colnames(MUD.OTU_sum_comb_taxa)[1]="OTU"
nrow(MUD.OTU_sum_comb_taxa)


#Subsetting the representative sequences to inlcude only fungal sequences
unfilter_rep_set_seq<- read.fasta(file = "USEARCHv11_files/rep_set_fwd_reads_demux_phix_filtered_fil_OTU.fa", as.string = TRUE, set.attributes = FALSE)

MUD.sort_rep_set_seq=unfilter_rep_set_seq[names(unfilter_rep_set_seq) %in% MUD.OTU_sum_comb_taxa$OTU]
MUD.sort_rep_set_seq[names(MUD.sort_rep_set_seq)]
length(MUD.sort_rep_set_seq)
#2674
MUD.sort_rep_set_seq=MUD.sort_rep_set_seq[order(names(MUD.sort_rep_set_seq))]
head(MUD.sort_rep_set_seq)
length(MUD.sort_rep_set_seq)
MUD.sort_rep_set_seq_sort=MUD.sort_rep_set_seq[order(-MUD.OTU_sum_comb_taxa$total_OTU_sum)]
head(MUD.sort_rep_set_seq_sort)
write.fasta(sequences =MUD.sort_rep_set_seq_sort, names = names(MUD.sort_rep_set_seq_sort), file.out ="R_files/MUD_rep_set_desending_sort.fna")

#Spot check of the top 20 most abundant

#OTU10 Eukaryota; Fungi; Dikarya; Basidiomycota; Ustilaginomycotina; Malasseziomycetes; Malasseziales; Malasseziaceae; Malassezia.
#Lee looked through the genus-level IDs on the revised top 500 OTUs, only OTU10 seemed to be a contaminant





#####Final filtered OTU table#####
MUD.Fung_fungi_field=subset_samples(MUD.Fung_fungi, project!="mock_com"&project!="spegar")
MUD.Fung_fungi_field<-prune_taxa(taxa_sums(MUD.Fung_fungi_field) > 0, MUD.Fung_fungi_field)
tax_table(MUD.Fung_fungi_field)
ntaxa(MUD.Fung_fungi_field)
#2674
sum(taxa_sums(MUD.Fung_fungi_field))
#2321517
#####Protax Classified Taxa to replace syntax####

MUD.fung.taxa_PROTAX_UNITE8=read.csv("R_files/MUD_taxa_PROTAX.csv", header = T)

head(MUD.fung.taxa_PROTAX_UNITE8)

#Clean up the names by remove unnecessary characters

MUD.fung.taxa_PROTAX_UNITE8$PROTAX=str_replace_all(MUD.fung.taxa_PROTAX_UNITE8$PROTAX,"_[1234567890]+","")
head(MUD.fung.taxa_PROTAX_UNITE8)


MUD.fung.taxa_PROTAX_UNITE8_sep=MUD.fung.taxa_PROTAX_UNITE8 %>% separate(PROTAX, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep = ",")
row.names(MUD.fung.taxa_PROTAX_UNITE8_sep)=MUD.fung.taxa_PROTAX_UNITE8_sep$OTU
MUD.fung.taxa_PROTAX_UNITE8_sep$OTU=NULL
MUD.fung.taxa_PROTAX_UNITE8_sep[is.na(MUD.fung.taxa_PROTAX_UNITE8_sep)] <- "unk"
unique(MUD.fung.taxa_PROTAX_UNITE8_sep$Domain)
MUD.fung.taxa_PROTAX_UNITE8_sep$Domain=replace(MUD.fung.taxa_PROTAX_UNITE8_sep$Domain,MUD.fung.taxa_PROTAX_UNITE8_sep$Domain=="ungi","Fungi")
MUD.fung.taxa_PROTAX_UNITE8_sep_mat=as.matrix(MUD.fung.taxa_PROTAX_UNITE8_sep)
head(MUD.fung.taxa_PROTAX_UNITE8_sep_mat)

TAXA_PROTAX_UNITE8_MUD_all=tax_table(MUD.fung.taxa_PROTAX_UNITE8_sep_mat)


#Combine with the phyloseq object
MUD.Fung_fungi_field_P=phyloseq(otu_table(MUD.Fung_fungi_field),sample_data(MUD.Fung_fungi_field),TAXA_PROTAX_UNITE8_MUD_all)
ntaxa(MUD.Fung_fungi_field_P)
#2674
get_taxa_unique(MUD.Fung_fungi_field_P, "Domain")
#"Fungi" "" 
#Two OTUs are unknown at Domain level

taxa_sums(subset_taxa(MUD.Fung_fungi_field_P,Domain==""))
#OTU4136 OTU2554 
#6639      24

#They both BLAST to fungi on NCBI

#Let's remove the common contaminant Malassezia 

MUD.Fung_fungi_field_P=subset_taxa(MUD.Fung_fungi_field_P, Genus!="Malassezia")
ntaxa(MUD.Fung_fungi_field_P)
#2654
sum(taxa_sums(MUD.Fung_fungi_field_P))
#2274136

save(MUD.Fung_fungi_field_P, file = "R_files/MUD.Fung_only_field_untransformed_phyloseq.RData")
write.table(otu_table(MUD.Fung_fungi_field_P), "R_files/MUD_fungi_OTU_field_ITS_only_fung_untransformed.txt")
write.table(tax_table(MUD.Fung_fungi_field_P), "R_files/MUD_fungi_taxon_PROTAX_ITS_trunc_phyl.txt") 

#####FunGuild File Creation#####
#Example formating

#OTU ID	sample1	sample2	sample3	sample4	sample5	taxonomy
#OTU_100	0	1	0	0	0	93.6%|Laetisaria_fuciformis|EU118639|SH012042.06FU|reps_singleton|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Corticiales;f__Corticiaceae;g__Laetisaria;s__Laetisaria_fuciformis
load(file = "R_files/MUD.Fung_only_field_untransformed_phyloseq.RData")
MUD.Fung_only_field_PROTAX = read.table("R_files/MUD_fungi_taxon_PROTAX_ITS_trunc_phyl.txt",header = T)
head(MUD.Fung_only_field_PROTAX)
nrow(MUD.Fung_only_field_PROTAX)
#2654


#Make a Taxonmy column 
MUD.Fung_only_field_PROTAX$taxonomy=paste("k__",MUD.Fung_only_field_PROTAX$Domain,";p__",MUD.Fung_only_field_PROTAX$Phylum,";c__",MUD.Fung_only_field_PROTAX$Class,
                                          ";o__",MUD.Fung_only_field_PROTAX$Order,";f__",MUD.Fung_only_field_PROTAX$Family,";g__",MUD.Fung_only_field_PROTAX$Genus,
                                          ";s__",MUD.Fung_only_field_PROTAX$Species,sep = "")
MUD.Fung_only_field_PROTAX$OTUs=row.names(MUD.Fung_only_field_PROTAX)
head(MUD.Fung_only_field_PROTAX)
MUD.Fung_only_field_OTU_PROTAX=merge(otu_table(MUD.Fung_fungi_field_P),MUD.Fung_only_field_PROTAX[,c("OTUs","taxonomy")], 
                                    by.x="row.names",by.y="OTUs")
head(MUD.Fung_only_field_OTU_PROTAX)


colnames(MUD.Fung_only_field_OTU_PROTAX)[1]="OTUID"


head(MUD.Fung_only_field_OTU_PROTAX$taxonomy)

write.table(MUD.Fung_only_field_OTU_PROTAX, 
            "R_files/MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.txt", row.names = F, sep = "\t") 
write.csv(MUD.Fung_only_field_OTU_PROTAX, "R_files/MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.csv", row.names = F) 
#I had to modify this in excel (i.e. turn "OTU.ID "to "OTU ID")

#RUN in shell 

#cd HardDrive/MUD_SequenceData/Analyses_collaboration/MUD_Transition/R_files/
#python Guilds_v1.1.py -otu MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.csv -db fungi -m -u

#Found 3422 matching taxonomy records in the database.
#Dereplicating and sorting the result...
#FunGuild tried to assign function to 2654 OTUs in 'MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.csv'.
#FUNGuild made assignments on 1628 OTUs.
#Result saved to 'MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.guilds.txt'

#Additional output:
#  FUNGuild made assignments on 1628 OTUs, these have been saved to MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.guilds_matched.txt.
#1026 OTUs were unassigned, these are saved to MUD.Fung_only_field_OTU_PROTAX_for_FunGuild.guilds_unmatched.txt.

#Total calculating time: 15.91 seconds.





#####Top 250 Taxa####

load("R_files/MUD.Fung_only_field_untransformed_phyloseq.RData")

sum(sort(taxa_sums(MUD.Fung_fungi_field_P),decreasing = T)[1:250])
#1769527
sum(taxa_sums(MUD.Fung_fungi_field_P))
#2274136
1769527/2274136
#0.7781096



#####Normalize the OTU matrix with DESEQ2####
#nomalizing using DESEQ2
load("R_files/MUD.Fung_only_field_untransformed_phyloseq.RData")
library(DESeq2)
library(vsn)
head(sample_data(MUD.Fung_fungi_field_P))
MUD.Fung.deseq=phyloseq_to_deseq2(MUD.Fung_fungi_field_P, ~site)#need to convert the phyloseq object to deseq2 dataset
#your grouping factor should be you samples

#First try the direct transformations
MUD.Fung.vsd<-varianceStabilizingTransformation(MUD.Fung.deseq)
#-- note: fitType='parametric', but the dispersion trend was not well captured by the
#function: y = a/x + b, and a local regression fit was automatically substituted.
#specify fitType='local' or 'mean' to avoid this message next time.

head(MUD.Fung.vsd)
plotSparsity(assay(MUD.Fung.vsd))
meanSdPlot(assay(MUD.Fung.vsd))

#Fix for common error from the support website https://support.bioconductor.org/p/63229/


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
#Warning message:
#  In sparseTest(counts(object, normalized = TRUE), 0.9, 100, 0.1) :
#  the rlog assumes that data is close to a negative binomial distribution, an assumption
#which is sometimes not compatible with datasets where many genes have many zero counts
#despite a few very large counts.
#In this data, for 20.2% of genes with a sum of normalized counts above 100, it was the case 
#that a single sample's normalized count made up more than 90% of the sum over all samples.
#the threshold for this warning is 10% of genes. See plotSparsity(dds) for a visualization of this.
#We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).

#Below are visualizations of the sparsity of the data across the normalization strategies
#We can use this information to choose the best transformation

plotSparsity(assay(dseq_f3))
plotSparsity(assay(dseq_f2))
plotSparsity(assay(dseq_f4))
plotSparsity(assay(MUD.Fung.vsd))
#below are the visualizations of the standard deviation as read counts increase 
#https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
meanSdPlot(log2(counts(dseq_f1,normalized=TRUE) + 1))

meanSdPlot(assay(dseq_f3))
meanSdPlot(assay(dseq_f2))
meanSdPlot(assay(dseq_f4))
meanSdPlot(assay(MUD.Fung.vsd))
#"Note that the vertical axis in such plots is the square root of the variance over all samples, 
#so including the variance due to the experimental conditions. While a flat curve of the square root 
#of variance over the mean may seem like the goal of such transformations, this may be unreasonable 
#in the case of datasets with many true differences due to the experimental conditions."

min(assay(MUD.Fung.vsd))

MUD.Fung_vst=assay(MUD.Fung.vsd)#variance stabilized transformation OTU table in matrix form
min(MUD.Fung_vst)
#0.1839722

MUD.Fung_vst[1:10,1:10]
min(rowSums(MUD.Fung_vst))


#Wow we need to turn it back into a phyloseq object
OTU.table.vst = otu_table(MUD.Fung_vst, taxa_are_rows = TRUE)
write.table(otu_table(OTU.table.vst), "R_files/MUD_fungi_OTU_field_ITS_only_fung_VST.txt") 

######End of Processing#####
