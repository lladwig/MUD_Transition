# MUD Shrubland to Grassland Transition
# Soil Enzyme Activity

# May 21, 2018
# LML

# Reset Rs brain
rm(list=ls())

# getwd tells you where R is currently looking
getwd()

# # setwd tells R where to look
setwd("/Users/laura/Desktop/Writing Projects/MUD transition/R/Data/")

# use getwd to confirm that R is now looking here
getwd()

# Install and load packages
#install.packages("gsubfn")
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
#library(Ternary) # making triangular plots
library(car) #for correct ANOVA table

#Importing the data
eea <-read.csv("EEA_MUD_transition.csv")
str(eea)
eea$AAP_soil <- as.numeric(levels(eea$AAP_soil)) [eea$AAP_soil] # changing factor to numeric
eea$AAP_OM <- as.numeric(levels(eea$AAP_OM)) [eea$AAP_OM] # changing factor to numeric


# ~*~*~*~*~*~Statistical Tests to see if actvity varies by species or location~*~*~*~*~*~
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

# Measured enzyme activity levels per gram soil and also per gram organic matter (OM). Not sure which I want to present, but it should be consistent across enzymes 
#RESULTS: There are significant differences in enzyme activities for nearly all enzymes whne looking at activity per gram soil. Any patterns seem to crap out when looking at activity per OM. So it looks like when there is OM in the soil, plants and microbes deploy the same suite of enzymes in the same concentrations to harvest it. But, since OM varies between locations (as seen from the soil OM analysis in the other file) we see higher enzyme activies in soils that have higher amounts of OM

# Does EEA vary between species or location?
## RESULTS: Location seems much more important than species. So not nesessarily who you are, but where you are matters. But that's a bit tricky, since we are only compairing two points for this, both speices with the ecotone, since grasses aren't in shrublands and shrubs aren't in grasslands

######## EEA Per Gram Soil Tests ##############
###############################################
mod_nag_soil <- lm(NAG_soil ~ Species + Site, data=eea)
summary(mod_nag_soil)
Anova(mod_nag_soil, type=3)
# RESULTS: varies by location
a1 <- aov(log(eea$NAG_soil) ~ eea$Species + eea$Site)
posthoc1 <-TukeyHSD(x=a1, 'eea$Site', conf.level=0.95)
print(posthoc1)
# RESULTS: grassland different than ecotone

#mod_alkp_soil <- lm(AlkP_soil ~ Species + Site, data=eea)
#summary(mod_alkp_soil)
#Anova(mod_alkp_soil, type=3)
## RESULT: ns, p=0.062 (but if sig, would vary by site)

mod_alkp_soil_log <- lm(log(AlkP_soil) ~ Species + Site, data=eea)
summary(mod_alkp_soil_log)
Anova(mod_alkp_soil_log, type=3)
# RESULT: p=0.04, varied by location
a2 <- aov(log(eea$AlkP_soil) ~ eea$Species + eea$Site)
posthoc2 <-TukeyHSD(x=a2, 'eea$Site', conf.level=0.95)
print(posthoc2)
# RESULTS: grassland and ecotone different

mod_aap_soil <- lm(AAP_soil ~ Species + Site, data=eea)
summary(mod_aap_soil)
Anova(mod_aap_soil, type=3)
# RESULTS: ns, p=0.0987

#mod_bg_soil <- lm(BG_soil ~ Species + Site, data=eea)
#summary(mod_bg_soil)
#Anova(mod_bg_soil, type=3)
## RESULTS: ns, p=0.68

mod_bg_soil_log <- lm(log(BG_soil) ~ Species + Site, data=eea)
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

# Ideal graphs: Separate graph for each enzyme. Scatter plot. Activity level on the y-axis, location on x-axis. For grassland and ecotone where several species are present, make the points much closer together so it's clear they are all in that location. Different shape point for each species, different color based on growth habit (grass = green, shrub = brown). Error bars


