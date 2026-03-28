## assessing the shapes of relations between temperature and trait and estimating the par-rs for these relations

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)
library(minpack.lm)
library(lme4)
library(glmmTMB)
library(stringr)
library(patchwork)

## source functions
source('./R/fit_shape.R')

## read in the studies using air temperature
fil <- list.files(path = './output/output_temp/', full.names = TRUE)

airt_df <- bind_rows(lapply(fil, FUN = function(x){readRDS(file = x)}))
airt_df <- airt_df %>%
  mutate(ID_fac = as.factor(ID), temper = 'air')

length(unique(airt_df$ID_fac))  #167

# read in the studies with sea surface temperature
fil_sst <- list.files(path = './output/output_temp_SeaB/', full.names = TRUE)

sst_df <- bind_rows(lapply(fil_sst, FUN = function(x){readRDS(file = x)}))
sst_df <- sst_df %>%
  mutate(ID_fac = as.factor(ID), temper = 'sst')

length(unique(sst_df$ID_fac))  # 45

# Merge both datasets
all <- bind_rows(airt_df, sst_df)
length(unique(all$ID_fac))  # 212

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                analyses of temperature over years across studies    ####

# 1. all dataset: both air temperature and SST
## Check how temperature across studies across years look like

ggplot(all, aes(x = Year, y = mean_Temp)) +
  geom_point()

ggplot(all, aes(x = Year, y = mean_Temp, col = ID_fac)) +
  geom_point() + theme(legend.position = 'none')

hist(all$mean_Temp)
quantile(all$mean_Temp, na.rm = T, probs = c(0.25, 0.5, 0.75))
mean(all$mean_Temp, na.rm = T)

## centering temperature across studies across years
all <- all %>%
  group_by(ID) %>%
  mutate(temp_center = scale(mean_Temp, center = TRUE, scale = FALSE)) %>%
  ungroup() %>%
  mutate(temp_center = as.numeric(temp_center))

ggplot(all, aes(x = Year, y = temp_center, col = ID_fac)) +
  geom_point() + theme(legend.position = 'none')
hist(all$temp_center)

## checking what ID has the temp centered beow and above -10 and 10
all$ID_fac[all$temp_center < -10]
all$ID_fac[all$temp_center > 10] ## 471 and 472

check <- subset(all, ID_fac %in% c('471', '472'))
head(check)
unique(check$Species)  # it is fish. I will just exclude these studies...

all <- all %>%
  filter(temp_center > -4 & temp_center < 5) %>%  ## as this leaves 1 or 2 data points in 3 studies, remove
  filter (! ID_fac %in% c("447", "471", "472"))# those studies altogether
length(unique(all$ID_fac)) ## 209 studies

sd(all$temp_center, na.rm = TRUE)  ##  SD: 0.631277

# number of phenological vs morphological
all_unique <- all %>%
  distinct(pick(ID), .keep_all = TRUE)

table(all_unique$Trait_Categ)

# Phenological Morphological
# 94           115

## duration per trait categ
all_dur_perTrait <- all %>%
  dplyr::summarize(.by = c(ID, Trait_Categ), Dur = dplyr::n()) %>%
  group_by(Trait_Categ) %>%
  dplyr::summarise(meanDur = mean(Dur),
                   medianDur = median(Dur),
                   minDur = min(Dur),
                   maxDur = max(Dur))

# Trait_Categ   meanDur medianDur minDur maxDur
# <fct>           <dbl>     <dbl>  <int>  <int>
#   1 Phenological     27.2        24      9     63
# 2 Morphological    19.8        15      9     60


all_Var <- all %>%
  dplyr:: group_by(ID) %>%
  summarize(SD = sd(temp_center, na.rm = TRUE))

hist(all_Var$SD, breaks = 20)

# consider including this figure in SI
mean(all_Var$SD, na.rm = TRUE)  ## 0.5991938  -> the value to be used as SD in baseline scenario
median(all_Var$SD, na.rm = TRUE)  ## 0.5988186
quantile(all_Var$SD, probs = c(0.25, 0.75), na.rm = TRUE)
#       25%       75%   ## updated values: 0.5197058 0.7351683
# 0.4777274 0.7021184
## with Julie we agreed to use 0.5 as baseline / burn-in SD for temp; and argue that this is the temperature for
# which the intercepts in the pop model for all LHS were par-rised CHECK in the par-risation file if that is the case


## but maybe 0.6 makes more sense now as median????
## Eva used 0.58 and it seemed to be fine. SO perhaps we
# can also go with 0.6 ---> but then also the slope of
## SD change have to be updated accordingly
min(all_Var$SD, na.rm = TRUE)  ##  0.07761635

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(all_Var$SD) ## mode 0.6028734

# explore relation between temperature and year
ggplot(all, aes(x = Year, y = temp_center, col = ID_fac)) +
  geom_point() + theme(legend.position = 'none')

## fit mixed-effects model with the study as random intercept, year as explanatory adn centered tmep as response
## to see what is the slope of temperature increase over time in this dataset
# singular fit...
mod_T_randi <- lmer(temp_center ~ Year + (1|ID_fac), data= all, REML = FALSE)
summary(mod_T_randi)  ## slope of temp on year: 0.006603


## include autocor
all$Year_fac <- as.factor(all$Year)
mod_T_randi_aut <- glmmTMB(temp_center ~ Year + (1 |ID_fac) +  ar1(Year_fac - 1 | ID_fac),
                   data = all)
summary(mod_T_randi_aut)  ## slope of temp on year:  0.007658, autocor = 0.42

# fitting no autocor with GLMMTmb - also issues with convergence
mod_T_randi_tmb <- glmmTMB(temp_center ~ Year + (1 |ID_fac),
                           data = all)
summary(mod_T_randi_tmb)   ## slope: 0.0066
### here because temp is centered the intercept is forced to 0 for all of them
# so there is no variation in random_intercept, and that is why this model
# has difficulties with singularity

# try with raw data, without centering
mod_T_randi_tmb_raw <- glmmTMB(mean_Temp ~ Year + (1 |ID_fac),
                           data = all)
summary(mod_T_randi_tmb_raw)   ## slope: 0.014
# AIC: 10885.1

# and the same with raw data but also autocor
mod_T_randi_aut_raw <- glmmTMB(mean_Temp ~ Year + (1 |ID_fac) +  ar1(Year_fac - 1 | ID_fac),
                           data = all)
summary(mod_T_randi_aut_raw)  ## slope of temp on year:  0.013932, autocor = 0.61
# AIC: 10423.2


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                prepare the data for fitting diff. shapes    ####


## add scaled values of trait  per study (i.e. z-score)
all_trans <- all %>%
  group_by(ID) %>%
  mutate(Trait_z = scale(Trait_mean)) %>%
  ungroup() %>%
  mutate(Trait_z= as.numeric(Trait_z))

hist(all_trans$Trait_z)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                plots of z-scored traits vs temperature    ####

## plot relations between trait and raw temperature and also of trait and temperature residuals per study. I actually think for traits I would first have to standardise them, and only then detrend
## (or maybe even no detrending at all is needed)
dat_abbrev <- all_trans %>%
  dplyr::group_by(., ID) %>%
  dplyr::distinct(., ID, Study_Authors, Species, .keep_all = TRUE)


abbrev_d <- dat_abbrev %>%
  dplyr::group_by(., ID) %>%
  tidyr::nest(.data =., cols = ID) %>%
  dplyr::mutate(
    fLet = purrr::map_chr(dat_abbrev$Species, ~substr(., start = 1, stop = 1)),
    sName = purrr::map_chr(dat_abbrev$Species, ~substr(strsplit(as.character(.), split = ' ')[[1]][2], 1, 4)),
    label = paste(Study_Authors, paste(fLet, sName, sep = '.'), Demog_rate, sep = ': ')
  ) %>%
  tidyr::unnest(cols = c(cols))


## plot of raw data
Temp_trait_raw <- ggplot(all_trans, aes(x = mean_Temp, y = Trait_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, raw') +
  ylab('Trait, raw')
print(Temp_trait_raw)

ggforce::n_pages(Temp_trait_raw)
pdf('./output/output_nonL/TraitvsTemperature_Raw.pdf')
for(i in 1:11){
  Temp_trait_raw <- ggplot(all_trans, aes(x = mean_Temp, y = Trait_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Temperature, raw') +
    ylab('Trait, raw')
  print(Temp_trait_raw)
}
dev.off()


## plot of raw traits vs centered temp
TempCenter_traitRaw <- ggplot(all_trans, aes(x = temp_center, y = Trait_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, centered') +
  ylab('Trait, raw')
print(TempCenter_traitRaw)

ggforce::n_pages(TempCenter_traitRaw)
pdf('./output/output_nonL/TraitRaw_vs_TemperatureCentered.pdf')
for(i in 1:11){
  TempCenter_traitRaw <- ggplot(all_trans, aes(x = temp_center, y = Trait_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')  + xlab('Temperature, centered') +
    ylab('Trait, raw')
  print(TempCenter_traitRaw)
}
dev.off()



## plot of scaled traits vs centered temp
TempCenter_traitZ <- ggplot(all_trans, aes(x = temp_center, y = Trait_z)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, centered') +
  ylab('Trait, z-score')
print(TempCenter_traitZ)

ggforce::n_pages(TempCenter_traitZ)
pdf('./output/output_nonL/TraitZ_vs_TemperatureCentered.pdf')
for(i in 1:11){
  TempCenter_traitZ <- ggplot(all_trans, aes(x = temp_center, y = Trait_z)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')  + xlab('Temperature, centered') +
    ylab('Trait, z-score')
  print(TempCenter_traitZ)
}
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####              compare diff shapes               ####

## test fit_shape for trait-temperature relations
check <- fit_shape(data = all_trans,
                   x = 'temp_center',
                   y = 'Trait_z',
                   ID = 1,
                   Thresh = 5,
                   out_folder = './output/output_nonL/shapes_temtrait/')

# suppress warnings? It is about the fit of the null model
# run function with detrended temperaturr eand standardized trait across studies
## trying with a rather low threshold of 4
shapes_fit <- do.call('rbind',
                      lapply(unique(all_trans$ID),
                             FUN = function(x){fit_shape(data = all_trans, ID = x,
                                                         Thresh = 5,  ## a rather low thresh
                                                         x = 'temp_center',
                                                         y = 'Trait_z',
                                                         out_folder = './output/output_nonL/shapes_temtrait/')}))
nrow(shapes_fit)  # 204, so out of 209 studies for which the data were available for several of them
# convergence with a particular (or several shapes) was not achieved
5/209 # 0.024 - so 2% of the studies
table(shapes_fit$Selected)
## linear  linear/sigmoid      quadratic        sigmoid
#     1            201              2

# check if any relations at all are better than the null model
hist(shapes_fit$Delta_Best_null)

nrow(shapes_fit[shapes_fit$Delta_Best_null > 0, ]) ## 147
nrow(shapes_fit)  ## 204

# So, not so many studies for which we see any shape fit better than the intercept only model (57 out of 204)
# 28%; however, the intercept only model can actually be considered as one type of a linear model (where slope is 0; a
# as slopes in this climate on trait model can be rather shallow). So, actually, I can consider dropping that
# bit with the null model from the function. see later

table(shapes_fit$mod_minAIC)
##  linear quadratic   sigmoid
##   96        22       86

##we will  use the minimum AIC to select the type of relations

# a figure like this will have to go to SI to back up the claim that it is very difficutl to differentiate
# linear from sigmoid relations in the data
hist(shapes_fit$Delta1_AIC[shapes_fit$Selected == 'linear/sigmoid'])  ## those usually are very very close...

## save the output
saveRDS(shapes_fit, file = './output/output_nonL/shapes_temtrait/Shapes_climTrait.RDS')

# read the saved file
#shapes_fit <- readRDS(file =  './output/output_nonL/shapes_temtrait/Shapes_climTrait.RDS')

## So: We will draw direct combis of slopes and and intercepts for the simulation model, using mod_minAIC to separate the three groups of relations
## relations selected based on the MinAIC
shapes_fitMin <- shapes_fit %>%
  mutate(sel_minAIC = case_when(mod_minAIC == 'linear' ~ 'linear',
                                mod_minAIC == 'sigmoid' ~  'sigmoid',
                                mod_minAIC == 'quadratic' ~ 'quadratic'))

table(shapes_fitMin$sel_minAIC)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####              explore fitted relationships                               ####

## 1. linear relations
Lin_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'linear') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'trait')

# explore the intercepts and slopes
hist(Lin_minAIC$Int_Lin)
min(Lin_minAIC$Int_Lin)  # -0.03894848
max(Lin_minAIC$Int_Lin)  # 0.03723515
# the intercept does not really differ from 0 as we have centered the tempeature
# prior to analyses; and we can fix it to 0 for linear relations

hist(Lin_minAIC$Beta_Lin)
median(Lin_minAIC$Beta_Lin)  # -0.2732398
mean(Lin_minAIC$Beta_Lin)  # -0.2726592
## climate affects the trait negatively because majority of the traits are phenologicla ones

# 2. quadratic relations
Quad_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'quadratic') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'trait')

# explore the intercepts and slopes
hist(Quad_minAIC$Int_Quad)
min(Quad_minAIC$Int_Quad)  # -0.4503806
max(Quad_minAIC$Int_Quad)  # 0.5538111
# given the sample is small, it i sokay to assume that the intercpet is 0

hist(Quad_minAIC$Beta_Quad)
median(Quad_minAIC$Beta_Quad)  # 0.04721768
mean(Quad_minAIC$Beta_Quad)  # -0.2625724

# 3. sigmoid relations
Sigm_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'sigmoid') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'trait')

# explore the intercepts and slopes
hist(Sigm_minAIC$Int_Sigm)
min(Sigm_minAIC$Int_Sigm)  # -0.007025333
max(Sigm_minAIC$Int_Sigm)  # 0.01258764
# Okay, again: we centered the predictor to 0; fine to assume that these are hardly different from 0

hist(Sigm_minAIC$Beta_Sigm)
median(Sigm_minAIC$Beta_Sigm)  #  -0.04716691
mean(Sigm_minAIC$Beta_Sigm)  # -0.05037761
## these slopes for sigmoid are very shallow, making it difficult to distinguish from a linear relation


## plots of slopes vs intercepts
ggplot(Lin_minAIC, aes(x = Int_Lin, y = Beta_Lin)) + geom_point()
ggplot(Quad_minAIC, aes(x = Int_Quad, y = Beta_Quad)) + geom_point()  ## remove the study with intercept < -15 !!!

ggplot(Quad_minAIC, aes(x = Int_Quad, y = Beta2_Quad)) + geom_point()
ggplot(Sigm_minAIC, aes(x = Int_Sigm, y = Beta_Sigm)) + geom_point()


# look ath the histogram of intercepts across all relations

shapes <-
  Lin_minAIC %>%
  ## import + join data
  bind_rows(Sigm_minAIC)  %>%
    bind_rows(Quad_minAIC) %>%

  ## add beta2 for linear and sigmoid
  mutate(Beta2_Lin = 0, Beta2_Sigm = 0) %>%

  ## only keep inputs for shape of interest, all in column with same name:
  ## final data will have 6 columns: ID, shape, relationship, int, beta, beta2
  pivot_longer(cols = -c(ID, shape, relationship), names_to = "var", values_to = "val") %>%
  separate(var, into = c("var", "type"), sep = "_") %>%
  mutate(
    type = case_when(
      type == "Lin"  ~ "linear",
      type == "Quad" ~ "quadratic",
      type == "Sigm" ~ "sigmoid"
    ),
    var = str_to_lower(var)
  ) %>%
  filter(shape == type) %>%
  pivot_wider(id_cols = c(ID, shape, relationship), names_from = var, values_from = val) %>%

  ## sort data
  arrange(ID)

shapes

hist(shapes$int)
quantile(shapes$int, proba = c(0.05, 0.1, 0.09, 0.95))
median(shapes$int)  # 2.556024e-10
mean(shapes$int)   #  -0.009556259

# here still add the mean / median - a plot like this, but cosmetics improved - in the main text
ggplot(shapes, aes(beta)) +
  geom_histogram(aes(y = after_stat(density))) +
  facet_wrap(vars(shape)) + theme_bw()

# actually better yet would be to directly plot the fits as lines;
# I think with abline or so this should be possible, check !!! TO DO
data <- data.frame(Temperature = seq(-3, 3, by = 0.1), Trait = seq(-3, 3, by = 0.1))

subs_linear <- subset(shapes, shape == 'linear')
lin_rel_temp_trait <- ggplot(data, aes(x = Temperature, y = Trait)) +
  lims(x = c(-3, 3), y =  c(-3, 3)) +
  geom_blank() +
  geom_abline(data = subs_linear, aes(intercept = 0, slope = beta),
                       alpha = 0.7, col = 'grey', lwd = 2) +
  geom_abline(aes(intercept = 0, slope = median(subs_linear$beta)
              ), col = 'black', lwd = 4) +
  scale_colour_brewer(palette = 'Dark2') +
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
lin_rel_temp_trait

pdf('./output/output_nonL/shapes_temtrait/plot_lin_relations_temp_trait.pdf')
lin_rel_temp_trait
dev.off()

## and I will need plots like this also for sigm and for the quadratic rel
# for them I would have to project over a rnage of the predictor
# using the betas

subs_quad <- subset(shapes, shape == 'quadratic')
subs_quad <- subs_quad[subs_quad$beta > -5, ]  ## remove the outlier
subs_quad <- subs_quad[subs_quad$ID != 304, ]
data_quad <- data.frame(Temperature = rep(seq(-2, 2, by = 0.1), nrow(subs_quad)),
                        ID =  rep(unique(subs_quad$ID), each = length(seq(-2, 2, by = 0.1))),
                        Trait = 0)

for(i in unique(subs_quad$ID)){
  data_quad$Trait[data_quad$ID == i] <- data_quad$Temperature[data_quad$ID == i]*subs_quad$beta[subs_quad$ID == i] +
    data_quad$Temperature[data_quad$ID == i]^2*subs_quad$beta2[subs_quad$ID == i]
}

# data_quad$ID[data_quad$Trait < -100]  ### 557
# subs_quad[subs_quad$ID == 557, ]  ## hthis study has very outlier-like looking estimates
# exclude
# data_quad$ID[data_quad$Trait < -15]  # 304
# subs_quad[subs_quad$ID == 304, ]  # also remove

data_quad$ID[data_quad$Trait < -15]

quad_rel_temp_trait <- ggplot(data_quad, aes(x = Temperature, y = Trait, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd =2) + ylim(-3, 3) + #xlim(-2, 2) +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))

quad_rel_temp_trait

pdf('./output/output_nonL/shapes_temtrait/plot_quad_relations_temp_trait.pdf')
quad_rel_temp_trait
dev.off()

# and, finally, the plot for the sigmoid relation
subs_sigm <- subset(shapes, shape == 'sigmoid')
data_sigm <- data.frame(Temperature = rep(seq(-2, 2, by = 0.1), nrow(subs_sigm)),
                        ID =  rep(unique(subs_sigm$ID), each = length(seq(-2, 2, by = 0.1))),
                        Trait = 0)

for(i in unique(subs_sigm$ID)){
  data_sigm$Trait[data_sigm$ID == i] <-
  (1/(1+exp(-5*(0 + subs_sigm$beta[subs_sigm$ID == i] * data_sigm$Temperature[data_sigm$ID == i]))) - 0.5)*4

}

# diff colours
sigm_rel_temp_trait <- ggplot(data_sigm, aes(x = Temperature, y = Trait, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd = 2) +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
sigm_rel_temp_trait

pdf('./output/output_nonL/shapes_temtrait/plot_sigm_relations_temp_trait.pdf')
sigm_rel_temp_trait
dev.off()


# grey shade
sigm_rel_temp_trait_grey <- ggplot(data_sigm, aes(x = Temperature, y = Trait, group = as.factor(ID))) +
  geom_line(col = 'darkgrey') +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 10))
sigm_rel_temp_trait_grey


# full plot with all shapes together
(lin_rel_temp_trait + sigm_rel_temp_trait + quad_rel_temp_trait) +
  plot_layout(nrow = 1)
