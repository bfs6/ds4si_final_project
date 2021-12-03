###  
### A Meta-Analysis of Country Level Studies on Environmental Change and Migration
###
### AUTHORS: Hoffmann, R., Dimitrova, A., Muttarak, R., Crespo-Cuaresma, J., Peisker, J. 
###

# >>> Set your working directory here
#setwd(" ")
# >>> Save the main RData files in the working directory to access them with this code

### **************************************************************
### OUTLINE ----------------------------------------------------
### **************************************************************

# 1. PREPARATION
# 2. OUTPUT TABLES AND FIGURES FOR MAIN PAPER
# 3. OUTPUT TABLES AND FIGURES FOR EXTENDED DATA
# 4. OUTPUT TABLES AND FIGURES FOR SUPPLEMENTARY MATERIAL


#
##
### **************************************************************
### 1. PREPARATIONS ---------------------------------------------- 
### **************************************************************
##
#


##
## DATA AND VARIABLE DESCRIPTIONS  ******************************** 
## 

## SL = study line / individual model

# DATA SETS       | DESCRIPTION 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# meta            | main meta data file, contains information on  1803 study line cases
# meta.plus       | main meta file plus coefficients estimated for environmental changes in *destination* countries
# country         | country data set with information on country characteristics and environmental change observed in countries for predictions
# country_rep     | NUmber of times countries are included in the country samples

# VARIABLE        | DESCRIPTION                       | SCALE   | CATEGORIES 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - # "paper_short"   | paper identifier                  | factor
# "paper"         | paper identifier                  | factor
# "author_id      | author pair identifier            | factor
# "Table"         | table identifier                  | factor
# "year"          | year of publication               | continuous

# stancoef        | standardized coefficient          | continuous
# stancoef.abs    | absolute standardized coefficient | continuous
# stancoef.w      | weighted standardized coefficient | continuous
# stanse          | standardized standard errors      | continuous
# coef            | orig. coefficient from study line | continuous
# serror          | orig. stan. error from study line | continuous
# mig_sd          | stan. deviation of migration var. | continuous
# env_sd          | stan. deviation of env. var.      | continuous

# "fe_time"       | study uses time fixed effects     | dummy   | 0=no fixed effects, 1=fixed effects
# "fe_spatial"    | study uses spatial fixed effects  | dummy   | 0=no fixed effects, 1=fixed effects
# interaction"    | effect size based on interation   | factor  | 0=no interaction, 1=interaction with dummy, effect is main coefficient, 2=interaction with dummy, combined effect, 3=interaction with continuous variable, combined effect
# "tab"           | table count within paper          | continuous
# "countrysample" | number of countries in sample     | continuous 
# "yearscovered"  | number of years in panel          | continuous
# "lin"           | all variables in linear form      | dummy
# "robustse"      | estimates (cluster) robust SE     | dummy  

# "rapidonset"    | SL considers rapid-onset event... | dummy   | 0=other env. factor, 1=rapid-onset event
# "pre_lev"       | ... precipitation level change    | dummy   | 0=other env. factor, 1=prec. level change
# "pre_var"       | ... precipitation variability     | dummy   | 0=other env. factor, 1=prec. variability change/anomaly
# "tem_lev"       | ... temperature level change      | dummy   | 0=other env. factor, 1=temp level change
# "tem_var",      | ... temperature variability       | dummy   | 0=other env. factor, 1=temp. variability change/anomaly
# "precipitation" | SL considers prec. level/variab   | dummy   | 0=other env. factor, 1=precipitation 
# "temperature"   |  SL considers temp. level/variab  | dummy   | 0=other env. factor, 1=temperature 
# "env_lag"       | lag between ind. and dep. var.    | continuous
# "env_lag_dum    | lag between ind. and dep. var.    | dummy   | 0=no lag, 1 = lag > 0
# "env_other_sum" | other env. factor controlled for  | count
# "env_other_dum" | other env. factor controlled for  | dummy
# "env_timespan"  | time interval to measure env. var | factor  | 1=1 year, 2=5 years, 3=10years
# "env_timespan_dum" | interval to measure env. var   | dummy   | 0=1 year, 1=>1 year 
# "env_others"    | other env. var. controlled for    | string
# "o_tem"         | other env. factor: temperature    | dummy   | 0=no other temp. measure, 1=other temp. measure controlled for
# "o_pre"         | other env. factor: precipitation  | dummy   | 0=no other prec. measure, 1=other prec. measure controlled for
# "o_sho"         | other env. factor: rapidonset     | dummy   | 0=no other rapid onset measure, 1=other rapid onset measure controlled for

# "internal",     |  SL considers internal migr. only | dummy   | 0=other migration, 1=internal migration 
# "dest_high"     | SL considers migr. to high-inc. c.| dummy   | 0=other migration, 1=international migration to high-income countries
# "dest_low",     | SL considers migr. to low-inc. c. | dummy   | 0=other migration, 1=international migration to low-income countries
# "dest_world"    | worldwide migr.                   | dummy   | 0=other migration, 1=international migration, worldwide
# "dest_ambi"     | destination cannot be categorized  | dummy   | 0=other migration, 1=international migration, destination ambiguous

# "period start"  | Start year of panel               | continuous
# "period end"    | End year of panel                 | continuous
# "income_channel"    | controls for income channel   | dummy   | 0=no control, 1= control
# "conflict_channel"  | controls for conflict channel | dummy   | 0=no control, 1= control 
# "weights"           | SL uses weights               | dummy   | 0=no weights, 1= weights 
# "robustness_check"  | SL is a robustness check      | dummy   | 0=no robustness check, 1= robustness check 
# "lin"               | SL uses a lin-lin estimation  | dummy   | 0=no lin-lin, 1= lin-lin 
# "gravity"           | SL uses a gravity model       | dummy   | 0=no gravity, 1= gravity 

# "published1"        | study was published           | dummy   | 0=unpublished (working paper, mimeo), 1= published 
# "published2"        | study was published by IF     | factor   | 0=unpublished, 1= published with IF <=2,2 = published with IF >2
# "pubimp"            | impact factor of publication  | continuous with 0 if not published by journal

# * Compositional shares *
# "nonoecd", "H","L" ,"LM" ,"UM", "M" ,"agr" "conflict_mepv_5" , "europena" ,"asia","lac" , "ssa" ,"mena"   
# > non oecd , high-income, low-income, lower-middle-income, upper-middle-income, agriculturally dependent, conflict, Europe/North America, Asia, Latin America/Caribbean, Sub Saharan Africa, Middle East and North Africa

# * Alternative conflict measures * 
# "conflict_mepv_5civ", "conflict_mepv_10", "conflict_mepv_10civ","conflict_ucdp_25_8" , "conflict_ucdp_0_8","conflict_ucdp_100_8", "conflict_ucdp_50_8"

# * Further controls included in specification *
# "controls_political" "controls_population" "controls_pastmigr"   "controls_econlevel"  "controls_culture" "controls_geo"   
#  "other_controls"  -> list of all control variables in specification 

# * Further study line model related information *
# "R2", "t.stats" ,  "se_detail" , "sample.size"  , "m_type", "m_detail"

# * Further information on migration outcome *
# "mig_destination", "mig_measure" , "mig_detail"  , "region"  

# * Further information on environmental factor/key independent variable *
# "env_type",  "env_location"  



##
## LOADING DATA AND PACKAGES *************************************
## 

rm(list=ls())

# >>> The following packages have to be installed for the code to run

require(meta)
require(metafor)
require(lfe)
require(lme4)
require(tidyverse)
require(stargazer)
require(ggpubr)
require(jtools)
require(RColorBrewer) 
require(rgeos)
require(rnaturalearth)
require(svglite)

# get world map
world <- ne_countries(returnclass = "sf")
# set ggplot theme
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

citation("meta")
citation("metafor")
citation("lfe")
citation("lme4")
citation("tidyverse")
citation("stargazer")
citation("ggpubr")
citation("jtools")
citation("RColorBrewer") 
citation("rgeos")
citation("rnaturalearth")

# load data
load("data-raw/Environmental Change and Migration_Meta and Country Data.rdata")

##  
## DESCRIPTIVE STATISTICS  ***************************************
##  

prop.table(table(meta$temperature));  table(meta$temperature)
prop.table(table(meta$precipitation));table(meta$precipitation)
prop.table(table(meta$rapidonset)); table(meta$rapidonset)
prop.table(table(meta$internal));table(meta$internal)


## ESTIMATING AVERAGE EFFECTS ACROSS ALL COEFFICIENTS Using Meta and Metafor packages
## see also: https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/

m.pm <- metagen(stancoeff,
                stanse,
                data = meta,
                byvar=paper,
                studlab = paste(paper),
                comb.fixed = FALSE,
                comb.random = TRUE,
                method.tau = "PM",
                hakn = T,
                prediction = TRUE,
                sm = "COR")
m.pm


##
##  DEFINING VARIABLES FOR MODELING ******************************
##

## MODELS IN MAIN TEXT

# Baseline controls
con_base <- c("interaction", "fe_time", "fe_spatial", "control_sum") 

# Environmental variables 1
con_clim1 <- c("pre_var", "rapidonset", "tem_lev", "tem_var")

# Environmental variables 2
con_clim2 <- c("env_lag_dum", "env_timespan_dum", "env_other_dum")

# Migration variables
con_mig <- c("internal", "dest_low", "dest_high", "dest_ambi") # , 

# Compositional measures 1: Regional composition
con_context1 <- "nonoecd"

# Compositional measures 2: Composition by country characteristics
con_context2 <-c("L", "LM", "UM", "agr", "conflict_mepv_5")


## MODELS IN SUPPLEMENTARY MATERIALS

# control variables used in models
con_control <- c("income_channel", "conflict_channel", "controls_political", "controls_population", "controls_pastmigr", "controls_econlevel", "controls_culture", "controls_geo")

# Sample size of models
con_samplesize <- c("countrysample", "yearscovered")

# Time and spatial fixed effects
con_fixedeffects <- c("fe_time", "fe_spatial")

# Weighting and further estimation features
con_estimation <- c("weights", "lin", "robustness_check")

# Publication characteristics of studies
con_pub1 <- c( "published1")

# Additional controls for focus on differential time periods
con_time2 <- c("period_start" , "period_end")

# Defining data and weights
dat <- meta
wt <- 1/meta$stanse^2


## 
## FUNCTIONS TO IMPLEMENT MODELS *****************************
##  

func.felm <- function(yname,xnames, fixedeffect.cluster) {
  controls <- paste(xnames, collapse = "+")
  fe.c <- paste(fixedeffect.cluster,  collapse = "|" )
  formula <- as.formula(paste(yname,"~",controls, "|", fe.c, sep = ""))
  return(formula)
}

func.lmer <- function(yname,xnames) {
  controls <- paste(xnames, collapse = "+")
  formula <- as.formula(paste(yname,"~",controls,"+(1|paper)", sep = ""))
  return(formula)
}



#
##
### **************************************************************
### 2. OUTPUT TABLES AND FIGURES FOR MAIN PAPER ------------------
### **************************************************************
##
#


##
## FIGURE 1 - CODE AVAILABLE UPON REQUEST
##


## 
## TABLE 1 WEIGHTED FELM MODEL: STANDARDIZED COEFF AS OUTCOME -----------------
## 

## COLUMN 1
m1a <- func.felm("stancoeff", c(con_base, con_clim1), c("paper", "0", "paper"))
m1b <- felm(m1a , data = dat, weights=wt)
summary(m1b)

## COLUMN 2
m2a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2), c("paper", "0", "paper"))
m2b <- felm(m2a , data = dat, weights=wt)
summary(m2b)

## COLUMN 3
m3a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig), c("paper", "0", "paper"))
m3b <- felm(m3a , data = dat, weights=wt)
summary(m3b)

## COLUMN 4
m4a <- func.felm("stancoeff", c(con_base, con_clim1, con_clim2, con_mig, con_context1), c("paper", "0" , "paper"))
m4b <- felm(m4a , data = dat, weights=wt) 
summary(m4b) 

## COLUMN 5
m5a <- func.felm("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2), c("paper", "0" , "paper"))
m5b <- felm(m5a ,  data = dat, weights=wt)
summary(m5b)


## TABLE 1
stargazer(m1b, m2b, m3b, m4b, m5b,
          type="html",
          out="output/table 1_baseline_m1-m5.doc",
          ci=F, 
          notes="Test",
          model.names = T,
          single.row = T,
          omit = c("interaction", "fe_time", "fe_spatial", "control_sum", "yearscovered", "countrysample"))

rm(m1a, m2a, m3a, m4a, m5a, m1b, m2b, m3b, m4b, m5b)

