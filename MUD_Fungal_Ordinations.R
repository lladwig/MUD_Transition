## MUD: Microbes Under Dominants
## LML
## march 5, 2018

# Reset Rs brain
rm(list=ls())

# getwd tells you where R is currently looking
getwd()

# # setwd tells R where to look
setwd("/Users/laura/Desktop/Writing Projects/MUD transition/R/Data/")

# use getwd to confirm that R is now looking here
getwd()

#install.packages("gsubfn")
library(ggplot2)
library(gsubfn)
library(reshape)
library(reshape2)
library(stringr)
library(cowplot)
library(dplyr)
library(vegan)
library(MASS) 
library(RVAideMemoire) # for tests after PERMANOVA

# This is the fungal data that Lukas pipelined for me and I'm getting to look at for the first time! The bacterial sequesces were a dud, so this project will mainly be on fungi. For methods on how data were cleaned and pipelined, see Lukas

#fun_otu2 <- read.table("MUD_fungi_OTU_ITS_trunc_phyl.txt", header = TRUE, sep = "\t")
#fun_otu <- read.delim("MUD_fungi_OTU_ITS_trunc_phyl.txt") #another import method
#str(fun_otu)
#head(fun_otu2)
#g <- as.data.frame(fun_otu)
## Found this online, but it seems like too much...
## https://stackoverflow.com/questions/20297564/converting-text-file-into-data-frame-in-r

## Lukas sent me the data in a .txt but for the life of me I can't get it to import correctly, (it was separted with spaces and I didn't now how to code that) so I pasted the data into excel and made it into a CSV because I know how to work with that

#fun_otu <- read.csv("MUD_fungi_OTU_ITS_trunc_phyl.csv")# Initial data that Lukas sent
fun_otu <- read.csv("MUD_fungi_OTU_ITS_VST.csv")# This has been transformed/normalized
#str(fun_otu)
#head(fun_otu)

#First need to separate data based on project, so I need to rearrange things to pull out SPEGAR

# A simple transformation with t fucks up the data format
#fun_otu_t <-t(fun_otu) #transposing data so I can separate sample names by .
#fun_otu_t <- as.data.frame(fun_otu.t)

fuckit <- cast(melt(fun_otu, id.var = "OTU"), variable ~ OTU)
#head(fuckit)
#str(fuckit)
fuckit <- as.data.frame(fuckit)

#Now separate out "variable" based on .

#First it would probably be best to do this separating in a separte dataset
fuckit_first <-fuckit$variable
fuckit_first <- as.data.frame(fuckit_first)
fuckit_first$split_me <- fuckit$variable
#split_me <- as.character(fuckit_first$split_me)

fuckit_first <- str_split_fixed (fuckit_first$split_me, "[.]", 4) #splits the one column into 4 based on where the periods are. preiod has to be in hard brackets otherwise it doesn't know it's a dot and just thinks it's any letter. the result is not a list, so it's easier to deal with

colnames(fuckit_first) <-c("taxon", "site", "spp", "rep")
fuckit_first <- as.data.frame(fuckit_first)

fuckit_first$variable <- paste(fuckit_first$taxon, fuckit_first$site, fuckit_first$spp, fuckit_first$rep, sep=".") #making the variable again so it can be merged back with the real data

fuckit_first$grass_rep <-paste(fuckit_first$spp, fuckit_first$rep, sep="_")#making a unique variable for the ordination

## Ok, now I have to merge fuckit_first back with the actual data in fuckit
fun_otu.b <- merge(fuckit_first, fuckit, all=TRUE)

fuckit_first$site_spp_rep <-paste(fuckit_first$site, fuckit_first$grass_rep, sep="_")#making a unique variable for the ordination
fun_otu.c <- merge(fuckit_first, fuckit, all=TRUE)
transect_fun_otu <- fun_otu.c %>%
	filter(site != "SPEGAR") #Filtering out SPEGAR data

# Need to get rid of the "mock" data. Not exactly sure how I will do this the best way, so for now I guess I'll filter it out
transect_fun_otu <- transect_fun_otu %>%
	filter(site !="mock1") #filtering out mock1 data
transect_fun_otu <- transect_fun_otu %>%
	filter(site !="mock2") #filtering out mock2 data

## Formatting dataset so it is ready for the nmds
transect_fun_otu_k <-transect_fun_otu [,-c(1,2,3,4,5,6)] #taking out extra columns
long.t <- melt(transect_fun_otu_k) # making wide to long
colnames(long.t) <-c("site", "spp", "value")
long.t_here <- long.t %>%
	filter(value > 0)

short=dcast(long.t_here, site~spp) #there shouldn't be duplicates here, so I don't know why sum is needed...
short[is.na(short)] = 0 #replaces NAs with 0s
rownames(short)<-short$site #set rownames to site names
short=short[,-1] 		#drop site names column

m1 <- metaMDS(short, distance = "bray", k = 3, trymax = 500,
        autotransform =FALSE,  
        wascores = TRUE, expand = TRUE,
        trace = 1, plot = TRUE,)

#plot ordination
plot(m1, display = c("sites"), choices = c(1,2),
     type = "t", shrink = TRUE)
plot(m1, display = c("sites"), choices = c(1,3),
     type = "t", shrink = TRUE)
plot(m1, display = c("sites"), choices = c(2,3),
     type = "t", shrink = TRUE)
scores(m1, display = c("sites", "species"), shrink = FALSE)
summary(m1)
 
#write.table(short,"/Users/laura/Desktop/Writing Projects/Savannas - WI/R/Results/OrdinationScoresAllCommunities_20170509.txt", sep="\t")       
print (m1)
str(m1)     
m1$species
m1$points
m1$ndim    

## Want to make a dataset of NMDS points that can be used for graphing
nmds.abun <-m1$points #pulling out just the scores
nmds.abun <- as.data.frame(nmds.abun) # making it a dataframe
nmds.abun <- add_rownames(nmds.abun, "sample") #making rownames into a column
nmds.abun <- as.data.frame(nmds.abun) #making it a dataframe again
nmds.abun<- nmds.abun %>%
	separate(sample, c("habitat", "spp", "rep"),"_") #split grouping column

write.csv(nmds.abun, "/Users/laura/Desktop/Writing Projects/MUD transition/R/Results/NMDS_Abundance_Scores_MUDtrans_20180621.csv")

##---------------------------------------------------
## Graphing the NMDS of fungal soil abundance
##---------------------------------------------------

#### ~*~*~ Figure Edit Wish List ~*~*~
## I'd like to control (manually pick) the colors of the habitats
## control (manually pick) the shapes of species
## Make a single legend that represents both species and habitat
## Make a similar graph with MDS1 and MDS2 and put the two graphs in a panel together


abun <-ggplot(nmds.abun, aes(x=MDS1, y=MDS2, fill = habitat, shape = spp)) +
	geom_point(aes(colour=habitat)) +
	theme_classic () +
	ylab("NMDS 2") +
	xlab("NMDS 1") +
	scale_fill_manual(values = c("red", "black", "orange")) + #this doesn't work
	theme(axis.text.y = element_text(size=14),
		axis.text.x = element_text(size = 14),
		axis.title.y = element_text(size = 16, face = "plain"),
		axis.title.x = element_text(size = 16, face = "plain"),
		panel.border = element_rect(colour ="black", size = 1, fill=NA)
		)
abun


#--------------------ADONIS--------------------------
# Testing if fungal communities vary by location or species
# First, need to make an environmental matrix to compare sequneces to. Basically just has to have the envornmental variables that I want to use for comparisons in the ANODIS, so spp and location

# Making the environmental matrix
data.enviro <- 0
data.enviro$matching <- transect_fun_otu$site_spp_rep #this is the column that will match up with the species dataset. I'll have to change the name of it to match what it is called in the species matrix, but for now I will call it "matching"
data.enviro$split <- transect_fun_otu$site_spp_rep #This needs to be split into site, spp and rep columns

data.enviro.b <- str_split_fixed (data.enviro$split, "_", 3) #splits the one column into 3 based on where the underscores are. the result is not a list, so it's easier to deal with
colnames(data.enviro.b) <- c("location", "spp", "rep") # naming the columns that were just created
data.enviro.b=as.data.frame(data.enviro.b) #changing it back to a dataframe
data.enviro.b$match = data.enviro$matching #adding the matching codes back in
enviro <- unique (data.enviro.b) # this shouldn't change it, since everything should be unique anyways, but it's a good check anyway. If numbers change, soemthing was wrong

#Now making species matrix
data.spp <- 0
data.spp$match <- long.t$site
data.spp$spp <- long.t$spp
data.spp$abun <-long.t$value
data.spp=as.data.frame(data.spp)
data.spp <-data.spp [,-c(1)] # dropping first column of zeros
data.sppwide <- dcast(data.spp, match~spp, sum) # turn to wide format, the "sum" statement says to return the value in the cell, otherwise it defaults to length, which just turns everything into "1", not quite true, it only turns things to a count if there are duplicates, since there are no duplicates in my data it should just enter the values fine. But that said, I don't think it hurts to have "sum" in the code
## Not sure if the matix has to be all abundance values or not, so for now, taking out the sample id column in the lines of code below
rownames(data.sppwide)<-data.sppwide$match #sets rownames to sample names
data.sppwide=data.sppwide[,-1] 		#drop site names column
#data.sppwide <- as.data.frame (data.sppwide) #changes back to a dataframe


## The statistical test, a 2 way permanova without an interaction
test <- adonis (data.sppwide ~ enviro$location + enviro$spp, permutations=10000)
print(test)
# RESULTS: fungal communities varied with species and location. 

# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(dist(data.sppwide, "euclidian"), enviro$location, nperm=2000)
# RESULT: all locations different than one another

# This code compares species
pairwise.perm.manova(dist(data.sppwide, "euclidian"), enviro$spp, nperm=2000)
# RESULTS: Nearly all species combinations are different, with the exception of BOER/PLJA which p=0.061

#####~*~*~*~*~*~*~*~*~ Presense/absense ~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
short.pa <- short
short.pa[short.pa>0] <-1 #This code makes it presense/absence data. 

m2 <- metaMDS(short.pa, distance = "jaccard", k = 3, trymax = 500,
        autotransform =FALSE,  
        wascores = TRUE, expand = TRUE,
        trace = 1, plot = TRUE,)

plot(m2, display = c("sites"), choices = c(1,2),
     type = "t", shrink = TRUE)
plot(m2, display = c("sites"), choices = c(1,3),
     type = "t", shrink = TRUE)
plot(m2, display = c("sites"), choices = c(2,3),
     type = "t", shrink = TRUE)
scores(m2, display = c("sites", "species"), shrink = FALSE)
summary(m2)

print (m2)
str(m2)     
m2$species
m2$points
m2$ndim

## Making a presense/absense dataframe to use for the ADONIS
data.sppwide.pa <- data.sppwide #First making a new dataframe so I don't overwrite the old one
data.sppwide.pa[data.sppwide.pa>0] <-1 #Then changing any value over 0 to a 1

## The statistical test, a 2 way permanova without an interaction
test.pa <- adonis (data.sppwide.pa ~ enviro$location + enviro$spp, permutations=10000)
print(test.pa)
# RESULTS: fungal OTU presense/absence statistically varies with both location (<0.0001) and species (p=0.0015)

# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(dist(data.sppwide.pa, "euclidian"), enviro$location, nperm=2000)
# RESULT: all locations different than one another

# This code compares species
pairwise.perm.manova(dist(data.sppwide.pa, "euclidian"), enviro$spp, nperm=2000)
# RESULTS: Nearly all species combinations are different, with the exception of BOER/PLJA which p=0.176


############################# more visuals #####################################

## It would be nice if I could look at the number of shared taxa











############### MORE SOIL CHARACTERISTICS ##################################

### Want to compare the species matrix to an environmental matrix
soil <- read.csv("SoilCharacteristics_MUD_Transition.csv")
str(soil)

#m1 <-as.data.frame(m1)
#soil <-as.data.frame(soil)

soil <- soil[,-c(4,5)] #getting rid of non-soil data columns

vec <-envfit(m1, soil, perm=999, na.rm=TRUE) #I think i added the environmental data
vec.df <-as.data.frame(vec$vectors$arrows*sqrt(vec$vectors$r))
vec.df$huh <-rownames(vec.df)

spp.matrix =data.frame(MDS1 = m1$points[,1], MDS2 = m1$points[,2])#making something to graph with

##Followed this info below for graphing arrows: https://stackoverflow.com/questions/14711470/plotting-envfit-vectors-vegan-package-in-ggplot2
## PROBLEM: this drops the catagorical traits from the analysis/graph and I don't know why
ggplot(data = spp.matrix, aes(MDS1, MDS2)) +
	#geom_point(aes(data=))
	geom_segment(data=vec.df, aes(x=0,xend=NMDS1,y=0,yend=NMDS2),
		arrow = arrow(length=unit(0.5, "cm")), colour="grey",inherit_aes=FALSE) +
	geom_text (data=vec.df, aes(x=NMDS1, y=NMDS2, label=huh), size=5)+
	coord_fixed()










####### CODE EXAMPLES######################

# I don't need this now, but helpful code to have for merging variables, basically creating the variables that I just split up in the code above.
#data.n$g <-paste (data.n$new.site, data.n$timing, sep= "_")

