# MUD Shrubland to Grassland Transition
# Soil Enzyme Activity

# May 21, 2018
# LML

# Reset Rs brain
rm(list=ls())

# getwd tells you where R is currently looking
getwd()

# # setwd tells R where to look
setwd("D:/MUD_SequenceData/Analyses_collaboration/MUD_Transition/")

# use getwd to confirm that R is now looking here
getwd()

# Install and load packages
#install.packages("gsubfn")
library(ggplot2)
#library(gsubfn)
library(reshape)
library(reshape2)
library(stringr)
#library(cowplot)
library(dplyr)
library(vegan)
library(MASS) 
library(lme4)
#library(vcd)
library(Ternary) # making triangular plots
library(car) #for correct ANOVA table
library(RVAideMemoire) # for tests after PERMANOVA
options(contrasts=c("contr.sum", "contr.poly"))
#Importing the data
eea <-read.csv("EEA_MUD_transition.csv")
str(eea)
eea$AAP_soil <- as.numeric(levels(eea$AAP_soil)) [eea$AAP_soil] # changing factor to numeric
eea$AAP_OM <- as.numeric(levels(eea$AAP_OM)) [eea$AAP_OM] # changing factor to numeric
colnames(eea)[1]="Site"

# ~*~*~*~*~*~Statistical Tests to see if actvity varies by species or location~*~*~*~*~*~
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

# Measured enzyme activity levels per gram soil and also per gram organic matter (OM). Not sure which I want to present, but it should be consistent across enzymes 
#RESULTS: There are significant differences in enzyme activities for nearly all enzymes whne looking at activity per gram soil. Any patterns seem to crap out when looking at activity per OM. So it looks like when there is OM in the soil, plants and microbes deploy the same suite of enzymes in the same concentrations to harvest it. But, since OM varies between locations (as seen from the soil OM analysis in the other file) we see higher enzyme activies in soils that have higher amounts of OM

# Does EEA vary between species or location?
## RESULTS: Location seems much more important than species. So not nesessarily who you are, but where you are matters. But that's a bit tricky, since we are only compairing two points for this, both speices with the ecotone, since grasses aren't in shrublands and shrubs aren't in grasslands

######## EEA Per Gram Soil Tests ##############
###############################################
mod_nag_soil <- lm(log(NAG_soil) ~ Species + Site, data=eea)
qqPlot(stdres(mod_nag_soil))
hist(stdres(mod_nag_soil))
shapiro.test(stdres(mod_nag_soil))
summary(mod_nag_soil)
Anova(mod_nag_soil, type=3)
# RESULTS: varies by location
a1 <- aov(log(eea$NAG_soil) ~ eea$Species + eea$Site)
summary(a1)
posthoc1 <-TukeyHSD(x=a1, 'eea$Site', conf.level=0.95)
print(posthoc1)
# RESULTS: grassland different than ecotone

#mod_alkp_soil <- lm(AlkP_soil ~ Species + Site, data=eea)
#summary(mod_alkp_soil)
#Anova(mod_alkp_soil, type=3)
## RESULT: ns, p=0.062 (but if sig, would vary by site)

mod_alkp_soil_log <- lm(log(AlkP_soil) ~ Species + Site, data=eea)
qqPlot(stdres(mod_alkp_soil_log))
hist(stdres(mod_alkp_soil_log))
shapiro.test(stdres(mod_alkp_soil_log))
summary(mod_alkp_soil_log)
Anova(mod_alkp_soil_log, type=3)
# RESULT: p=0.04, varied by location
a2 <- aov(log(eea$AlkP_soil) ~ eea$Species + eea$Site)
summary(a2)
posthoc2 <-TukeyHSD(x=a2, 'eea$Site', conf.level=0.95)
print(posthoc2)
# RESULTS: grassland and ecotone different

mod_aap_soil <- lm(AAP_soil ~ Species + Site, data=eea)
qqPlot(stdres(mod_aap_soil))
hist(stdres(mod_aap_soil))
shapiro.test(stdres(mod_aap_soil))
summary(mod_aap_soil)
Anova(mod_aap_soil, type=3)
# RESULTS: ns, p=0.0987

a3 <- aov((eea$AAP_soil) ~ eea$Species + eea$Site)
summary(a3 )
posthoc3 <-TukeyHSD(x=a3, 'eea$Site', conf.level=0.95)
print(posthoc3)

#mod_bg_soil <- lm(BG_soil ~ Species + Site, data=eea)
#summary(mod_bg_soil)
#Anova(mod_bg_soil, type=3)
## RESULTS: ns, p=0.68

mod_bg_soil_log <- lm(log(BG_soil) ~ Species + Site, data=eea)
qqPlot(stdres(mod_bg_soil_log))
hist(stdres(mod_bg_soil_log))
shapiro.test(stdres(mod_bg_soil_log))
summary(mod_bg_soil_log)
Anova(mod_bg_soil_log, type=3)
# RESULTS: ns, p=0.58



# Checking model assumptions
par(mfrow=c(2,2))

plot(mod_nag_soil) #looks okay
plot(mod_alkp_soil) # looks okayish
plot(mod_alkp_soil_log) #yes, I think it fits the line better
plot(mod_aap_soil) #looks great
#plot(mod_bg_soil) # ends are starting to come off the line
plot(mod_bg_soil_log) #log does make the data fit better, should use it


######## EEA Per Gram Organic Matter Tests ##########
#####################################################

#mod_nag_om <- lm(NAG_OM ~ Species + Site, data=eea)
#summary(mod_nag_om)
#Anova(mod_nag_om, type=3)
## RESULTS: ns, p=0.63

mod_nag_om_log <- lm(log(NAG_OM) ~ Species + Site, data=eea)
summary(mod_nag_om_log)
Anova(mod_nag_om_log, type=3)
# RESULTS: ns, p=0.57

mod_alkp_om <- lm(AlkP_OM ~ Species + Site, data=eea)
summary(mod_alkp_om)
Anova(mod_alkp_om, type=3)
# RESULTS: ns, p=0.68

mod_aap_om <- lm(AAP_OM ~ Species + Site, data=eea)
summary(mod_aap_om)
Anova(mod_aap_om, type=3)
# RESULTS: ns, p=0.515

#mod_bg_om <- lm(BG_OM ~ Species + Site, data=eea)
#summary(mod_bg_om)
#Anova(mod_bg_om, type=3)
## RESULTS: ns, p=0.81

mod_bg_om_log <- lm(log(BG_OM) ~ Species + Site, data=eea)
summary(mod_bg_om_log)
Anova(mod_bg_om_log, type=3)
# RESULTS: ns, p=0.93

# Checking model assumptions
par(mfrow=c(2,2))

#plot(mod_nag_om) #not the greatest, consider log transformation
plot(mod_nag_om_log) #log fits much better
plot(mod_alkp_om) #looks okay
plot(mod_aap_om) #looks great
#plot(mod_bg_om) #okay, but consider checking log transformation
plot(mod_bg_om_log) #yes, log fits it better


# ~*~*~*~*~*~                 Graphing enzyme activities             ~*~*~*~*~*~
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

# Ideal graphs: Separate graph for each enzyme. Scatter plot. Activity level on the y-axis, location on x-axis. 
#For grassland and ecotone where several species are present, make the points much closer together so it's 
#clear they are all in that location. Different shape point for each species, different color based on growth 
#habit (grass = green, shrub = brown). Error bars

eea=eea %>% group_by(Site,Species)
eea_sum=summarise_at(eea, vars(NAG_soil,AlkP_soil,AAP_soil,BG_soil),funs(mean,se=sd(.)/sqrt(n()),sd))
eea_sum$site_spp=with(eea_sum, interaction(Site,Species))
positions=c("G.Boer", "G.Bogr", "G.Plja", "E.Boer", "E.Plja", "E.Latr", "S.Latr")

NAG_g=ggplot(eea_sum, aes(site_spp, NAG_soil_mean, ymin = NAG_soil_mean-NAG_soil_se, ymax = NAG_soil_mean+NAG_soil_se))

(NAG_scat_p=NAG_g+geom_errorbar(width=0.25)+geom_line(aes(group = Site))+
    geom_point(size=7,aes(color=Species))+scale_color_brewer(palette="Dark2")+
    scale_y_continuous(name = "NAG activity (??mol,h,g soil)")+scale_x_discrete(limits = positions)+
    theme_bw()+theme(legend.position = "none",axis.text.y=element_text(size=18),axis.text.x=element_text(size=18), 
                     axis.title.y=element_text(size=20),axis.title.x=element_blank(),panel.grid.major=element_blank(),
                     panel.grid.minor=element_blank()))



## ~*~*~*~*~*~*~*~ GRAPHING SOIL TEXTURE~*~*~*~*~*~*~*~ 
# This works, just have to get different colors
## Attempting to make a triangular (ternary) plot
#First creating a unique naming column, then I might need to sub-setting data so it is just the soil texture variables
#char$unique <- paste(char$Location, char$Spp, char$Rep, sep="_")
char <-read.csv("SoilCharacteristics_MUD_Transition.csv")
char$LocSpp <-paste(char$Location, char$Spp, sep="_")
print(char)

#subsetting soil texture. have to figure out how to keep "unique" in there
tex <- char[c(5:7)]
str(tex)

#get rid of NAs
tex<-tex[complete.cases(tex), ]
str(tex)

# Ternaryplot with vcd takes matrices but not data frames, but with Ternary data frames can be used, too
# Future possibilities: Only graph averages of spp by location
#xM <-as.matrix()
TernaryPlot(alab='sand', blab ='silt', clab='clay')
TernaryPoints(tex, col="red")



#~*~*~*~*~*~*~STATISTICAL TESTS FOR ENVIROMENTAL VARIABLES*~*~*~*~*~*~*~
#lmer mixed model with random effects

#Fixed effect models
# separate model for each environmental factor
# The Anova statement gives me an anova table, just like I'd get in SAS. I took out the species by location interaction since not all species were at all locations. 

# Running data test
mod_clay <- lm(CLAY ~ Spp + Location, data=char)
summary(mod_clay)
#anova(mod_clay) #Just anova gives the wrong numbers, but I think the code below is correct
Anova(mod_clay, type=3)
# RESULT: varies with location
a3 <- aov((char$CLAY) ~ char$Spp + char$Location)
posthoc3 <-TukeyHSD(x=a3, 'char$Location', conf.level=0.95)
print(posthoc3)
# RESULTS: Ecotone different than shrub and grassland; grassland and shrubland not different.

mod_silt <- lm(SILT ~ Spp + Location, data=char)
summary(mod_silt)
Anova(mod_silt, type=3)
#RESULT: differs by location
a4 <- aov((char$SILT) ~ char$Spp + char$Location)
posthoc4 <-TukeyHSD(x=a4, 'char$Location', conf.level=0.95)
print(posthoc4)
# RESULTS: shrub and ecotone differ, all else similar

mod_sand <- lm(SAND ~ Spp + Location, data=char)
summary(mod_sand)
Anova(mod_sand, type=3)
#RESULT: differ by location
a5 <- aov((char$SAND) ~ char$Spp + char$Location)
posthoc5 <-TukeyHSD(x=a5, 'char$Location', conf.level=0.95)
print(posthoc5)
# RESULTS: Ecotone different than shrub and grassland; grassland and shrubland not different.

mod_pH <- lm(pH ~ Spp + Location, data=char)
summary(mod_pH)
Anova(mod_pH, type=3)
#RESULT: Differ by location and species
a6 <- aov((char$pH) ~ char$Spp + char$Location)
posthoc6 <-TukeyHSD(x=a6, 'char$Location', conf.level=0.95)
print(posthoc6)
#RESULT: No difference between sites in post hoc test, location just looks like it's getting taken along for the significant ride
a7 <- aov(char$pH ~ char$Spp + char$Location)
posthoc7 <-TukeyHSD(x=a7, 'char$Spp', conf.level=0.95)
print(posthoc7)
#RESULTS: LATR and BOGR different, all else not




## Note: this data distribution doens't really match model assuptions very well. But there are a decent amount of zeros in the data so I also cannot log transform the data and run the test without having to throw out those data. I'd like to keep the zeros in, so for now I am sticking with the less than ideal model until I can find something better
mod_p <- lm(P_ppm ~ Spp + Location, data=char)
summary(mod_p)
Anova(mod_p, type=3)
#RESULT: differe by location
a8 <- aov(char$P_ppm ~ char$Spp + char$Location)
posthoc8 <-TukeyHSD(x=a8, 'char$Location', conf.level=0.95)
print(posthoc8)
# RESULTS: all locations different

#char$P_ppm[which(!is.finite(char$P_ppm))] = NA
#mod_p_log <- lm(log(P_ppm) ~ Spp + Location, data=char)
#summary(mod_p_log)
#Anova(mod_p_log, type=3)


mod_k <- lm(K_ppm ~ Spp + Location, data=char)
summary(mod_k)
Anova(mod_k, type=3)
#RESULT: Differ by species and location
a9 <- aov(char$K_ppm ~ char$Spp + char$Location)
posthoc9 <-TukeyHSD(x=a9, 'char$Location', conf.level=0.95)
print(posthoc9)
#RESULT: Grass different than ecotone or shrub; ecotone and shurb similar
a10 <- aov(char$K_ppm ~ char$Spp + char$Location)
posthoc10 <-TukeyHSD(x=a10, 'char$Spp', conf.level=0.95)
print(posthoc10)
#RESULT: different BOGR-BOER LATR-BOGR PLJA-BOGR; similar LATR-BOER PLJA-BOER PLJA-LATR

mod_om <- lm(OM_percent ~ Spp + Location, data=char)
summary(mod_om)
#anova(mod_om)
Anova(mod_om, type=3)
#RESULTS: differ my location
a11 <- aov(char$OM_percent ~ char$Spp + char$Location)
posthoc11 <-TukeyHSD(x=a11, 'char$Location', conf.level=0.95)
print(posthoc11)
#RESULTS: Econtone differ from grass and srhurb; grass and srhrub similar

mod_n <- lm(NO3_N_ppm ~ Spp + Location, data=char)
summary(mod_n)
Anova(mod_n, type=3)
#RESULT: differ by location
a12 <- aov(char$NO3_N_ppm ~ char$Spp + char$Location)
posthoc12 <-TukeyHSD(x=a12, 'char$Location', conf.level=0.95)
print(posthoc12)
#RESULTS: grass and shrub similar; ecotone different than everyone else


## This was code I found online to make a good anova table, but I don't understand it. Instead, I'm using Anova with type =3 as Jon Henn suggested
#results <- anova(mod_om)
#Df <- results$Df
#SumSq <-results$"Sum Sq"
#MeanSq <- results$"Mean Sq"
#Fvalue <- results$"F value"
#Pvalue <- results$"Pr(>F)"
#Error.Term <- MeanSq[3]
#df.error <- Df[3]

#Fvalue[1] <- MeanSq[1]/Error.Term
#Pvalue[1] <- 1 - pf(Fvalue[1], Df[1], df.error)

#Fvalue[2] <- MeanSq[2]/Error.Term
#Pvalue[2] <- 1 - pf(Fvalue[2], Df[2], df.error)

#Ftable <- cbind(Df, SumSq, MeanSq, Fvalue, Pvalue)
#rownames(Ftable) <- c("Spp", "Location", "Spp:Location", "Residuals")
#print(Ftable)


# Checking model assumptions
# In general, they look okay
par(mfrow=c(2,2))

plot(mod_clay)
plot (mod_silt)
plot (mod_sand)
plot (mod_pH)
plot (mod_p) #maybe wonky? Don't remember how to log transform
plot (mod_k)
plot (mod_om)
plot (mod_n)
par(mfrow=c(1,1))
#plot(fitted(mod_clay), resid(mod_clay))
#qqnorm(resid(mod_clay))
#qqnorm(ranef(mod_clay)[,1])
#plot (x, resid(mod_clay))

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

#First need to separate data based on project, so pulling out all the SPEGAR first and working with them might be best.

# A simple transformation with t fucks up the data format
#fun_otu_t <-t(fun_otu) #transposing data so I can separate sample names by .
#fun_otu_t <- as.data.frame(fun_otu.t)

fuckit <- dcast(melt(fun_otu, id.var = "OTU"), variable ~ OTU)
#head(fuckit)
#str(fuckit)


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

spegar_fun_otu <- fun_otu.b %>%
  filter(site == "SPEGAR") #Filtering out just SPEGAR data

spegar_fun_otu <-spegar_fun_otu [,-c(1,2,3,4,5)] # Taking out excess variables for ordination

long <- melt(spegar_fun_otu) # making wide to long
colnames(long) <-c("site", "spp", "value")

# I need to get rid of OTUs that were present in the other project but not SPEGAR. Not sure what the best way do to it is, but my plan is to take out all the zeros in the long format, then turn it wide and put the zeros back in. Note: it might be worth it to run an ordination with presense/absence data too.

long_here <- long %>%
  filter(value > 0)

#wide_here <- dcast(long_here, site~spp) #turned data wide 

##This is code for the package "reshape2". 
#Cast (wtf) abundance data
short=dcast(long_here, site~spp) #there shouldn't be duplicates here, so I don't know why sum is needed...
short[is.na(short)] = 0 #replaces NAs with 0s
rownames(short)<-short$site #set rownames to site names
short=short[,-1] 		#drop site names column

#library(vegan)
#library(MASS)     
m1 <- metaMDS(short, distance = "bray", k = 3, trymax = 500,
              autotransform =FALSE,  
              wascores = TRUE, expand = TRUE,
              trace = 1, plot = TRUE,)

"Data:     short 
Distance: bray 

Dimensions: 3 
Stress:     0.1231396 
Stress type 1, weak ties
Two convergent solutions found after 20 tries
Scaling: centring, PC rotation, halfchange scaling 
Species: expanded scores based on 'short' "


#This is the code I should need for running an ordination
#plot ordination
plot(m1, display = c("sites"), choices = c(1,2),
     type = "t", shrink = TRUE)
plot(m1, display = c("sites"), choices = c(1,3),
     type = "t", shrink = TRUE)
plot(m1, display = c("sites"), choices = c(2,3),
     type = "t", shrink = TRUE)
scores(m1, display = c("sites", "species"), shrink = FALSE)

#write.table(short,"/Users/laura/Desktop/Writing Projects/Savannas - WI/R/Results/OrdinationScoresAllCommunities_20170509.txt", sep="\t")       
print (m1)
str(m1)     
m1$species
m1$points
m1$ndim    


##### Now lets look at the transect data (eg, non-SPEGAR)
# --------------------------------------------------------------------------------
fuckit_first$site_spp_rep <-paste(fuckit_first$site, fuckit_first$grass_rep, sep="_")#making a unique variable for the ordination
fun_otu.c <- merge(fuckit_first, fuckit, all=TRUE)
transect_fun_otu <- fun_otu.c %>%
  filter(site != "SPEGAR") #Filtering out SPEGAR data

# Need to get rid of the "mock" data. Not exactly sure how I will do this the best way, so for now I guess I'll filter it out
transect_fun_otu <- transect_fun_otu %>%
  filter(site !="mock1") #filtering out mock1 data
transect_fun_otu <- transect_fun_otu %>%
  filter(site !="mock2") #filtering out mock2 data

transect_fun_otu_k <-transect_fun_otu [,-c(1,2,3,4,5,6)] 
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
par(mfrow=c(1,1))
#This is the code I should need for running an ordination
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
untrans.otu=t(read.table("MUD_fungi_OTU_ITS_trunc_phyl.txt"))
row.names(untrans.OTU)
colnames(untrans.OTU)
tans_map=read.csv("MUD_Transition_map.csv",row.names = 1)
untrans.OTU=merge(tans_map,untrans.otu, by="row.names")

row.names(untrans.OTU)=untrans.OTU[,"Row.names"]
untrans.OTU_trans=subset(untrans.OTU, project=="transitions")
untrans.OTU_trans_o=untrans.OTU_trans[,6:length(colnames(untrans.OTU))]
untrans.OTU_map=untrans.OTU_trans[,1:6]
untrans.OTU_trans_o[,1]
short.pa <- untrans.OTU_trans_o
short.pa[short.pa>0] <-1 #This code makes it presense/absence data. 
max(short.pa)
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
test.pa <- adonis (short.pa ~ untrans.OTU_map$Location + untrans.OTU_map$Spp, permutations=10000)
print(test.pa)
# RESULTS: fungal OTU presense/absence statistically varies with both location (<0.0001) and species (p=0.0015)

# I *think* this code is a quasi post-hoc test that will compare the locations
pairwise.perm.manova(dist(short.pa, "euclidian"), untrans.OTU_map$Location, nperm=2000)
# RESULT: all locations different than one another

# This code compares species
pairwise.perm.manova(dist(short.pa, "euclidian"), untrans.OTU_map$Spp, nperm=2000)
# RESULTS: Nearly all species combinations are different, with the exception of BOER/PLJA which p=0.176


############################# more visuals #####################################

## It would be nice if I could look at the number of shared taxa











############### MORE SOIL CHARACTERISTICS ##################################

### Want to compare the species matrix to an environmental matrix
soil <- read.csv("SoilCharacteristics_MUD_Transition.csv")
str(soil)

#m1 <-as.data.frame(m1)
#soil <-as.data.frame(soil)
soil$site
rownames(short)

soil <- soil[,-c(4,5)] #getting rid of non-soil data columns

(vec <-envfit(m1, soil, perm=999, na.rm=TRUE)) #I think i added the environmental data
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



(fit <- envfit(m1, soil, perm = 999, na.rm=TRUE))
scores(fit, "vectors")
plot(m1)
plot(fit)
plot(fit, p.max = 0.05, col = "red")
summary(fit)



#let's run a mantel test between the soil EEA and community dist
#first let's normalize so values are between 1 and 0
eea <-read.csv("EEA_MUD_transition.csv")
str(eea)
eea$AAP_soil <- as.numeric(levels(eea$AAP_soil)) [eea$AAP_soil] # changing factor to numeric
eea$AAP_OM <- as.numeric(levels(eea$AAP_OM)) [eea$AAP_OM] # changing factor to numeric
row.names(eea)=eea$site
eea_o=eea[,c("NAG_soil", "NAG_OM","AlkP_soil","AlkP_OM","AAP_soil","AAP_OM","BG_soil","BG_OM")]
str(eea_o)
summary(eea_o)
eea.nom=decostand(eea_o, "range",na.rm=T)


eea.dist=vegdist(eea.nom,method="euclidean")
short_dist=vegdist(short,method="bray")
mantel(short_dist,eea.dist)
#Mantel statistic r: 0.2045 
#Significance: 0.003 

#let's run a mantel test between the soil nutrients and community dist
#first let's normalize so values are between 1 and 0

char <-read.csv("SoilCharacteristics_MUD_Transition.csv")
char$LocSpp <-paste(char$Location, char$Spp, sep="_")
print(char)


soil_nut=c("SAND","SILT","CLAY","pH","P_ppm","K_ppm","OM_percent","NO3_N_ppm")

char_o=char[,soil_nut]

char.nom=decostand(char_o, "range",na.rm=T)
char.dist=vegdist(char.nom,method="euclidean")
short_dist=vegdist(short,method="bray")
mantel(short_dist,char.dist)
#Mantel statistic r: 0.3313 
#Significance: 0.001 

#finally the soil nut versus soil eea 

mantel(eea.dist,char.dist)
#Mantel statistic r: 0.03963 
#Significance: 0.248 


#let try a partial mantel 

mantel.partial(short_dist,eea.dist,char.dist)
#Mantel statistic r: 0.203 
#Significance: 0.002 

####### CODE EXAMPLES######################

# I don't need this now, but helpful code to have for merging variables, basically creating the variables that I just split up in the code above.
#data.n$g <-paste (data.n$new.site, data.n$timing, sep= "_")

