###  
### A Meta-Analysis of Country Level Studies on Environmental Change and Migration
###
### AUTHORS: Hoffmann, R., Dimitrova, A., Muttarak, R., Crespo-Cuaresma, J., Peisker, J. 
###


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
library(sjPlot)

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
citation("svglite")
citation("sjPlot")

# load data
load("data-raw/Environmental Change and Migration_Meta and Country Data.rdata")


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
con_control <- c("income_channel", "conflict_channel", "controls_political", 
                 "controls_population", "controls_pastmigr", "controls_econlevel",
                 "controls_culture", "controls_geo")

# Sample size of models
con_samplesize <- c("countrysample", "yearscovered")

# Time and spatial fixed effects
con_fixedeffects <- c("fe_time", "fe_spatial")

# Weighting and further estimation features
con_estimation <- c("weights", "lin", "robustness_check")

# Publication characteristics of studies
con_pub1 <- c("published1")

# Additional controls for focus on differential time periods
con_time2 <- c("period_start" , "period_end")

# Defining data and weights
# dat <- meta
wt <- 1/meta$stanse^2

## 
## FUNCTIONS TO IMPLEMENT MODELS *****************************
##  

func.felm <- function(yname,xnames, fixedeffect.cluster) {
  controls <- paste(xnames, collapse = " + ")
  fe.c <- paste(fixedeffect.cluster,  collapse = "|" )
  formula <- as.formula(paste(yname,"~",controls, "|", fe.c, sep = ""))
  return(formula)
}

func.lmer <- function(yname,xnames) {
  controls <- paste(xnames, collapse = "+")
  formula <- as.formula(paste(yname,"~",controls,"+(1|paper)", sep = ""))
  return(formula)
}

#################################################
####Reproducing Model and Figure for Figure 3####
#################################################

####Clean the Data####
dat <- 
  meta %>% 
  mutate(region = 
           region %>%
           str_to_title(),
         region = 
           case_when(region == "East Asia And The Pacific" ~ "East Asia Pacific",
                     region %in% c("Latin America Caribbean", "Latin America", "Latin America And Caribbean") ~ "LAC",
                     region == "Themiddle East And Northafrica" ~ "MENA",
                     region == "Nonmiddle East And Northafrica" ~ "non-MENA",
                     region == "Sub-Saharan Africa" ~ "SSA",
                     region == "Nonlatin America And Caribbean" ~ "non-LAC",
                     region == "Nonsouth Asia" ~ "non-South Asia",
                     region == "Noneurope And Central Asia" ~ "non-Europe and Central Asia",
                     region == "Noneast Asia And The Pacific" ~ "non-East Asia Pacific",
                     TRUE ~ as.character(region)),
         region = 
           region %>% 
           str_replace_all("East Asia Pacific", "East Asia/Pacific") %>% 
           str_replace_all("Europe And Central Asia", "Europe/Central Asia") %>% 
           str_replace_all("Lac", "LAC") %>% 
           str_replace_all("Mena", "MENA") %>% 
           str_replace_all("Oecd", "OECD") %>% 
           str_replace_all("Non", "non") %>% 
           str_replace_all("Ssa", "SSA") %>% 
           str_replace_all("Countries In The South", "Global South") %>%
           str_replace_all("World Without Central And Eastern Europe As Well As The Balkan States And Cyprus", 
                           "non-Central and Eastern Europe"))

######################
####Start Analyses####
######################

####Analysis 1: Add Predictors to Mixed-Effects Model####
##Create Model
pred1a <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2))
  ## The new mixed effects model includes all of the Environmental drivers, Further environmental controls, 
  ## Migration destination, and Sample composition. The variables used are all of the variables used in the m5b model 
  ## from the original reproducibility code. 
pred1a_og <- func.lmer("stancoeff", c(con_base, con_context2))
pred1b_new <- lmer(pred1a, data = dat , weights = wt) ## This is the new mixed effects model 
pred1b_og <- lmer(pred1a_og, data = dat, weights = wt) ## This is the old mixed effects model
summary(pred1b_new)

##Plots and Predictions
dat.predict <- 
  dat %>% 
  mutate(predict = predict(pred1b_new, dat))

  # Same Plots from Figure 3 Rerun for New Mixed Effects Model
  # plot1, plot2, and plot3 in the below code are the same as the plots from the original paper, but 
  # they are all based on the new predictions from the pred1b_new mixed effects model we changed
plot1 <- 
  dat.predict %>% 
  filter(agr<0.5 & conflict_mepv_5<0.5) %>% 
  ggplot(aes(y=predict, x=M))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.05, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color = "black") +
  labs(
    x = "% middle-income countries in sample", 
    y = "Predicted environmental effect"
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))
plot1

plot2 <- 
  dat.predict %>%
  filter(L>0.8) %>% 
  ggplot(aes(y=predict, x=agr))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm",se=T, color="Black")+
  labs(
    x = "% agriculturally dependent countries in sample",
    y = ""
  )+
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.3),
        axis.text=element_text(size=11.3))
plot2

plot3 <- dat.predict %>%
  filter(agr<0.5) %>% 
  ggplot(aes(y=predict, x=conflict_mepv_5))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color="black")+
  labs(
    x = "% conflict countries in sample",
    y = ""
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  ) +
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))

plot3

plot123 <- ggarrange(plot1, plot2, plot3, align="h", nrow = 1, labels = "auto")

plot123    
ggsave("plots/new-analyses/plots_mixed_effects_model_1.jpg", width=15, height=5, dpi=600)


####Analysis 2: Add Region Variables to Mixed Effects Model####
##Create Model
con_location <- c("asia", "lac", "mena", "ssa", "europena")
me_pred_formula_region <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_location, con_context2))
me_pred_formula_region <- 
  paste0(
    as.character(me_pred_formula_region)[2], 
    as.character(me_pred_formula_region)[1], 
    as.character(me_pred_formula_region)[3], sep = " ") %>% 
  paste0(" + (1 | region)") %>% 
  as.formula()
  #The model includes location variables but also adds region as a random effect
me_lmer_region <- lmer(me_pred_formula_region, data = dat, weights = wt)

##Plots and Predictions  
dat.predict2 <- 
  dat %>% 
  mutate(predict = predict(me_lmer_region, dat))

# Same Plots from Figure 3 Rerun for New Mixed Effects Model
# plot1, plot2, and plot3 in the below code are the same as the plots from the original paper, but 
# they are all based on the new predictions from the pred1b_new mixed effects model we changed
plot1 <- 
  dat.predict2 %>% 
  filter(agr<0.5 & conflict_mepv_5<0.5) %>% 
  ggplot(aes(y=predict, x=M))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.05, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color = "black") +
  labs(
    x = "% middle-income countries in sample", 
    y = "Predicted environmental effect"
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))
plot1

plot2 <- 
  dat.predict2 %>%
  filter(L>0.8) %>% 
  ggplot(aes(y=predict, x=agr))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm",se=T, color="Black")+
  labs(
    x = "% agriculturally dependent countries in sample",
    y = ""
  )+
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.3),
        axis.text=element_text(size=11.3))
plot2

plot3 <- 
  dat.predict2 %>%
  filter(agr<0.5) %>% 
  ggplot(aes(y=predict, x=conflict_mepv_5))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color="black")+
  labs(
    x = "% conflict countries in sample",
    y = ""
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  ) +
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))

plot3

plot123 <- ggarrange(plot1, plot2, plot3, align="h", nrow = 1, labels = "auto")

plot123    
ggsave("plots/new-analyses/plots_mixed_effects_model_2.jpg", width=15, height=5, dpi=600)

##Visualize Region and Location Predictors 
effects_plots2 <- 
  plot_model(me_lmer_region, terms = con_location, 
             axis.labels=c("Asia", "Latin America/Caribbean", "Middle East/North Africa", "Sub-Saharan Africa"),
             show.values=TRUE, show.p=TRUE,
             title="Effects of Regional Proportions in Meta-Analysis")
ggsave(plot = effects_plots2, "plots/new-analyses/effects_plots2.jpg", width=7, height=7, dpi=600)


####Analysis 3: Filter for Conflict > 0.5####
##Filter and Create Model
conflict_data <- 
  dat %>% 
  filter(conflict_mepv_5 > 0.5)
wt_conflict <- 1/conflict_data$stanse^2
conflict_formula <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_location, con_context2))
conflict_model <- lmer(conflict_formula, data = conflict_data, weights = wt_conflict)
summary(conflict_model) 

##Plots and Predictions 
dat.predict3 <- 
  conflict_data %>% 
  mutate(predict = predict(conflict_model, conflict_data))

# Same Plots from Figure 3 Rerun for New Mixed Effects Model
# plot1, plot2, and plot3 in the below code are the same as the plots from the original paper, but 
# they are all based on the new predictions from the pred1b_new mixed effects model we changed
plot1 <- 
  dat.predict3 %>% 
  filter(agr<0.5) %>% 
  ggplot(aes(y=predict, x=M))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.05, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color = "black") +
  labs(
    x = "% middle-income countries in sample", 
    y = "Predicted environmental effect"
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))
plot1

plot2 <- 
  dat.predict3 %>%
  filter(L>0.8) %>% 
  ggplot(aes(y=predict, x=agr))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm",se=T, color="Black")+
  labs(
    x = "% agriculturally dependent countries in sample",
    y = ""
  )+
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  )+
  scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.3),
        axis.text=element_text(size=11.3))
plot2

plot3 <- 
  dat.predict3 %>%
  filter(agr<0.5) %>% 
  ggplot(aes(y=predict, x=conflict_mepv_5))+
  geom_hline(yintercept = 0) +
  geom_jitter(width=0.03, height=0.01, alpha=0.25)+
  geom_smooth(method="lm", color="black")+
  labs(
    x = "% conflict countries in sample",
    y = ""
  ) +
  coord_cartesian(
    ylim = c(-0.05, 0.05), 
    xlim = c(0, 1)
  ) +
  scale_y_continuous(breaks=seq(-0.05,0.05,0.025))+
  scale_x_continuous(labels=scales::percent)+
  theme(axis.title=element_text(size=13.8),
        axis.text=element_text(size=11.3))

plot3

plot123 <- ggarrange(plot1, plot2, plot3, align="h", nrow = 1, labels = "auto")

plot123    
ggsave("plots/new-analyses/plots_mixed_effects_model_3.jpg", width=15, height=5, dpi=600)

##Visualize Regional Effects
effects_plots3 <- 
  plot_model(conflict_model, terms = con_location, 
             axis.labels=c("Asia", "Latin America/Caribbean", "Middle East/North Africa", "Sub-Saharan Africa"),
             show.values=TRUE, show.p=TRUE,
             title="Effects of Regional Proportions in Meta-Analysis, Conflict Model")
ggsave(plot = effects_plots3, "plots/new-analyses/effects_plots3.jpg", width=7, height=7, dpi=600)


####Analysis 4: Add Interaction Terms to Mixed Effects Model#### 
##Create Model
interaction_model_formula <- func.lmer("stancoeff", c(con_base,  con_clim1, con_clim2, con_mig, con_context2, con_location, "agr:conflict_mepv_5", "agr:L", "conflict_mepv_5:L", "agr:conflict_mepv_5:L"))
interaction_model_formula <- 
  paste0(
    as.character(interaction_model_formula)[2], 
    as.character(interaction_model_formula)[1], 
    as.character(interaction_model_formula)[3], sep = " ") %>% 
  paste0(" + (1 | region)") %>% 
  as.formula()
interaction_model <- lmer(interaction_model_formula, data = dat, weights = wt)
summary(interaction_model)

##Visualize Interaction Terms
full_interaction_term_effects_plot4 <- 
  plot_model(interaction_model, 
             title = "Effects of Model w/ Interaction Terms")
ggsave(plot = full_interaction_term_effects_plot4, "plots/new-analyses/full_interaction_term_effects_plot4.jpg", width=5, height=5, dpi=600)

interaction_effects_model_plot4 <- 
  plot_model(interaction_model, 
             terms = c("agr:conflict_mepv_5", "L:agr", "L:conflict_mepv_5", "L:agr:conflict_mepv_5"),
             show.values = TRUE, show.p = TRUE,
             axis.labels = c("Agricultural Dependence x Conflict", 
                             "Low-income Countries x Agricultural Dependence",
                             "Low-income Countries x Conflict", 
                             "Low-income Countries x Agricultural Dependence x Conflict"),
             title = "Interaction Term Effects")
ggsave(plot = interaction_effects_model_plot4, "plots/new-analyses/interaction_effects_model_plot4.jpg", width=7, height=7, dpi=600)


####Analysis 5: Build Stepwise LM Model####
##Create Model
lm_pred_formula <- paste0("stancoeff ~ ", paste0(c(con_base,  con_clim1, con_clim2, con_mig, con_context2, con_location, "agr:conflict_mepv_5", "agr:L", "conflict_mepv_5:L", "agr:conflict_mepv_5:L"), collapse = " + "), " + paper + region")
lm_version <- lm(formula = lm_pred_formula, data = dat, weights = wt)
step_lm_version <- step(lm_version)

##Visualize Model
effects_plots5 <- 
  plot_model(step_lm_version, 
             title = "Effects of Stepwise Linear Model", 
             rm.terms = c("paperGhimire et al. 2015"))
effects_plots5_v2 <- 
  plot_model(step_lm_version, 
             title = "Effects of Stepwise Linear Model w/o Region and Paper",
             terms = str_subset(names(step_lm_version$coefficients), 
                                "^region|^paper|^fe|^interaction|^control", negate = TRUE),
             axis.labels = c("Precipitation variability/anomalies", "Rapid-Onset event", "Temperature level change",
                             "Temperature variability/anomalies", "Measurement timeframe > 1 year", 
                             "Internal migration", "International, destination only low/middle-income countries",
                             "% low-income countries in sample", "% low/middle-income countries in sample", 
                             "% agriculturally dependent countries in sample", "% conflict countries in sample",
                             "% Asian countries in sample", "% Latin American/Caribbean countries in sample", 
                             "% Middle Eastern/North African countries in sample", 
                             "Agricultural Dependence x Conflict", "Low-income Countries x Agricultural Dependence",
                             "Low-income Countries x Conflict", "Low-income Countries x Agricultural Dependence x Conflict"))
ggsave(plot = effects_plots5, "plots/new-analyses/effects_plots5.jpg", width = 8, height = 7, dpi = 600)
ggsave(plot = effects_plots5_v2, "plots/new-analyses/effects_plots5_v2.jpg", width = 7, height = 7, dpi = 600)


####Visualize All Models in Stargazer Table####
stargazer(pred1b_og, pred1b_new, me_lmer_region, conflict_model, interaction_model,
          type="html",
          out="output/new-analyses/stargazer_table.doc",
          ci=FALSE,
          title = "Comparison of Mixed Effects Models",
          model.names = TRUE,
          single.row = TRUE, 
          omit = c("interaction", "fe_time", "fe_spatial", "control_sum", "yearscovered", "countrysample"), 
          covariate.labels=c("Precipitation variability/anomalies", "Rapid-Onset event", "Temperature level change",
                             "Temperature variability/anomalies", "Environment-migration lag in years",
                             "Measurement timeframe > 1 year", "Other environmental factors controlled for in original model",
                             "Internal migration", "International, destination only low/middle-income countries",
                             "International, destination only high-income countries", "International, destination ambiguous",
                             "% Asian countries in sample", "% Latin American/Caribbean countries in sample", 
                             "% Middle Eastern/North African countries in sample", "% Sub-Saharan African countries in sample",
                             "Agricultural Dependence x Conflict", "Low-income Countries x Agricultural Dependence",
                             "Low-income Countries x Conflict", "Low-income Countries x Agricultural Dependence x Conflict",
                             "% low-income countries in sample", "% low/middle-income countries in sample", 
                             "% upper/middle-income countries in sample", "% agriculturally dependent countries in sample", 
                             "% conflict countries in sample", "Intercept"),
          dep.var.labels.include = TRUE
)




# ####Plots with Updated Filters and Outcomes - Analysis 3 (3 analyses in this one)####
# # On top of changes in filters and variables visualized, these predictions are also based on the new 
# # mixed-effects model
# 
# # plot_1_new keeps the same filters but looks at the percentage of high-income and 
# # upper middle-income countries in the sample, instead of low-income countries. 
# # Even in richer countries, enivronmental effects of migration can still have an increasing effect. 
# plot1_new <- 
#   dat.predict %>% 
#   mutate(MnL = L + M,
#          HnUM = H + UM) %>% 
#   filter(agr<0.5 & conflict_mepv_5<0.5) %>% 
#   ggplot(aes(y=predict, x=HnUM))+
#   geom_hline(yintercept = 0) +
#   geom_jitter(width=0.05, height=0.01, alpha=0.25)+
#   geom_smooth(method="lm", color = "black") +
#   labs(
#     x = "% High-income and upper middle-income countries in sample", 
#     y = "Predicted environmental effect"
#   ) +
#   coord_cartesian(
#     ylim = c(-0.05, 0.05), 
#     xlim = c(0, 1)
#   )+
#   scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
#   scale_x_continuous(labels=scales::percent)+
#   theme(axis.title=element_text(size=13.8),
#         axis.text=element_text(size=11.3))
# plot1_new
# 
# # plot_2_new looks at the same variables, but the filter changes from the low-income country sample % being > 0.8 to < 0.8. 
# # The obvious takeaway is that the relationships are completely inverse when the sample changes to be mostly non-low-income countries.
# # It seems that the environmental effects decrease. Higher-income countries that are dependent on agriculture, see less of an 
# # environmental effect on migration. There are more opportunities for economic opportunity, separate from environmental factors. 
# plot2_new <- 
#   dat.predict %>%
#   filter(L < 0.8) %>% 
#   ggplot(aes(y=predict, x=agr))+
#   geom_hline(yintercept = 0) +
#   geom_jitter(width=0.03, height=0.01, alpha=0.25)+
#   geom_smooth(method="lm",se=T, color="Black")+
#   labs(
#     x = "% agriculturally dependent countries in sample",
#     y = ""
#   )+
#   coord_cartesian(
#     ylim = c(-0.05, 0.05), 
#     xlim = c(0, 1)
#   )+
#   scale_y_continuous(breaks = seq(-0.05,0.05,0.025))+
#   scale_x_continuous(labels=scales::percent)+
#   theme(axis.title=element_text(size=13.3),
#         axis.text=element_text(size=11.3))
# plot2_new
# 
# 
# # plot3_new is the same as plot3, except that the filter is now agr > 0.5, meaning that we are only focused on the samples where 
# # agriculturally dependent countries is greater than 50%. This emphasizes the idea that agriculturally dependent countries that are 
# # also conflict countries, tend to be very strongly affected by environmental effects on mirgrations 
# plot3_new <- 
#   dat.predict %>%
#   filter(agr > 0.5) %>% 
#   ggplot(aes(y=predict, x=conflict_mepv_5))+
#   geom_hline(yintercept = 0) +
#   geom_jitter(width=0.03, height=0.01, alpha=0.25)+
#   geom_smooth(method="lm", color="black")+
#   labs(
#     x = "% conflict countries in sample",
#     y = ""
#   ) +
#   coord_cartesian(
#     ylim = c(-0.05, 0.05), 
#     xlim = c(0, 1)
#   ) +
#   scale_y_continuous(breaks=seq(-0.05,0.05,0.025))+
#   scale_x_continuous(labels=scales::percent)+
#   theme(axis.title=element_text(size=13.8),
#         axis.text=element_text(size=11.3))
# 
# plot3_new
# 
# plot123_new <- ggarrange(plot1_new, plot2_new, plot3_new, align="h", nrow = 1, labels = "auto")
# 
# plot123_new
# 
# ggsave("plots/new-analyses/figure-3_line_plots_sample_composition_effects_NEW_PLOTS.jpg", width=15, height=5, dpi=600)
# 
