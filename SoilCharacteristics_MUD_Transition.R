#MUD Shrubland to Grassland transition
#General Soil Characteristics

#March 30, 2018
#LML

# Reset Rs brain
rm(list=ls())

# getwd tells you where R is currently looking
getwd()

# # setwd tells R where to look
setwd("/Users/laura/Desktop/Writing Projects/MUD transition/R/Data/")

# use getwd to confirm that R is now looking here
getwd()

#install.packages("gsubfn")
#install.packages("vcd")
#install.packages("Ternary")
#install.packages("car")
library(ggplot2)
#library(gsubfn)
#library(reshape)
#library(reshape2)
#library(stringr)
#library(cowplot)
#library(dplyr)
library(vegan)
#library(MASS) 
library(lme4)
#library(vcd)
library(Ternary)
library(car)

char<-read.csv("SoilCharacteristics_MUD_Transition.csv")

## ~*~*~*~*~*~*~*~ GRAPHING SOIL TEXTURE~*~*~*~*~*~*~*~ 
# This works, just have to get different colors
## Attempting to make a triangular (ternary) plot
#First creating a unique naming column, then I might need to sub-setting data so it is just the soil texture variables
#char$unique <- paste(char$Location, char$Spp, char$Rep, sep="_")
char$LocSpp <-paste(char$Location, char$Spp, sep="_")
print(char)

#subsetting soil texture. have to figure out how to keep "unique" in there
tex <- char[c(6:8)]
str(tex)

#get rid of NAs
tex<-tex[complete.cases(tex), ]
str(tex)

# Ternaryplot with vcd takes matrices but not data frames, but with Ternary data frames can be used, too
# Future possibilities to make graph prettier: 
# Only graph averages of spp by location
# color-code points to match other graphs
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

#plot(fitted(mod_clay), resid(mod_clay))
#qqnorm(resid(mod_clay))
#qqnorm(ranef(mod_clay)[,1])
#plot (x, resid(mod_clay))



