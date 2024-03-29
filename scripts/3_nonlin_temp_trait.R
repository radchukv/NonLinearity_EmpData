## assessing the shapes of relations between temperature and trait and estimating the par-rs for these relations

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)
library(minpack.lm)
library(lme4)
library(glmmTMB)

## source functions
source('./R/fit_shape.R')

## read in the studies
fil <- list.files(path = './output/output_temp/', full.names = TRUE)

all<- bind_rows(lapply(fil, FUN = function(x){readRDS(file = x)}))
all <- all %>%
  mutate(ID_fac = as.factor(ID))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                analyses of temperature over years across studies    ####

## Check how temperature across studies across years look
ggplot(all, aes(x = Year, y = mean_Temp)) +
  geom_point()

ggplot(all, aes(x = Year, y = mean_Temp, col = ID_fac)) +
  geom_point() + theme(legend.position = 'none')

hist(all$mean_Temp)
quantile(all$mean_Temp, na.rm = T, probs = c(0.25, 0.5, 0.75))
mean(all$mean_Temp, na.rm = T)

## centering temperature across studies across years
all <- all %>%
  mutate(temp_center = as.numeric(scale(mean_Temp, center = TRUE, scale = FALSE)))

ggplot(all, aes(x = Year, y = temp_center, col = ID_fac)) +
  geom_point() + theme(legend.position = 'none')
hist(all$temp_center)

sd(all$temp_center, na.rm = TRUE)  ##  SD: 5.814646  -> the value to be used as SD in baseline scenario


## fit mixed-effects model with the study as random intercept, year as explanatory adn centered tmep as response
## to see what is the slope of temperature increase over time in this dataset
mod_T_randi <- lmer(temp_center ~ Year + (1|ID_fac), data= all, REML = FALSE)
summary(mod_T_randi)  ## slope of temp on year: 0.021691


## include autocor
all$Year_fac <- as.factor(all$Year)
mod_T_randi_aut <- glmmTMB(temp_center ~ Year + (1 |ID_fac) +  ar1(Year_fac - 1 | ID_fac),
                   data = all)
summary(mod_T_randi_aut)  ## slope of temp on year:  0.040716, autocor = 0.94



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                prepare the data for fitting diff. shapes    ####


## add scaled values of trait  per study (i.e. z-score)
all_trans <- all %>%
  group_by(ID) %>%
  mutate(Trait_z = scale(Trait_mean))

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
for(i in 1:8){
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
for(i in 1:13){
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
for(i in 1:13){
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
                   Thresh = 2,
                   out_folder = './output/output_nonL/shapes_temtrait/')

# run function with detrended temperaturr eand standardized trait across studies
## trying with a rather low threshold of 4
shapes_fit <- do.call('rbind', lapply(unique(all_trans$ID), FUN = function(x){fit_shape(data = all_trans, ID = x,
                                                                                        Thresh = 4,  ## a rather low thresh
                                                                                        x = 'temp_center',
                                                                                        y = 'Trait_z',
                                                                                        out_folder = './output/output_nonL/shapes_temtrait/')}))

table(shapes_fit$Selected)
## linear  linear/sigmoid      quadratic        sigmoid
# 3            236              1              1

table(shapes_fit$mod_minAIC)
##  linear quadratic   sigmoid
##   97        16       128

##we will  use the minimum AIC to select the type of relations

hist(shapes_fit$Delta1_AIC[shapes_fit$Selected == 'linear/sigmoid'])  ## those usually are very very close...

## if using a threshol of 2:
## linear linear/sigmoid      quadratic        sigmoid
## 47            188              5              1

## save the output with the deltaAIC of 4
#saveRDS(shapes_fit, file = './output/output_nonL/shapes_temtrait/Shapes_4DeltaAIC.RDS')

# read the saved file
shapes_fit <- readRDS(file =  './output/output_nonL/shapes_temtrait/Shapes_4DeltaAIC.RDS')
table(shapes_fitMin$Selected)

# some exploration of the results
Lin <- shapes_fit %>%
  filter(Selected %in% c('linear/sigmoid', 'linear'))

Sigm <- shapes_fit %>%
  filter(Selected %in% c('linear/sigmoid', 'sigmoid'))

hist(Lin$Beta_Lin)
hist(Lin$Int_Lin)

# check of the correlation
ggplot(Lin, aes(x = Int_Lin, y = Beta_Lin)) + geom_point()
cor.test(Lin$Int_Lin, Lin$Beta_Lin)  ##        cor  -0.4361498

## So: We will draw direct combis of slopes and and intercepts for the simulation model, using mod_minAIC to separate the three groups of relations
## relations selected based on the MinAIC
shapes_fitMin <- shapes_fit %>%
  mutate(sel_minAIC = case_when(mod_minAIC == 'linear' ~ 'linear',
                                mod_minAIC == 'sigmoid' ~  'sigmoid',
                                mod_minAIC == 'quadratic' ~ 'quadratic'))

table(shapes_fitMin$sel_minAIC)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####              extract var-cov matrices               ####

## fit the mixed-effect model for linear relation
Lin_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'linear')

Quad_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'quadratic')

Sigm_minAIC <- shapes_fitMin %>%
  filter(sel_minAIC == 'sigmoid')

# first create a subset of the studies in the raw data set for linear relations
lin_raw <- all_trans %>%
  filter(ID %in% Lin_minAIC$ID) %>%
  mutate(Trait_z = as.numeric(Trait_z))

lin_mix <- glmmTMB(Trait_z ~ temp_center + (1|ID_fac), data = lin_raw, REML = FALSE)
summary(lin_mix)

# these will be parameters for the rmvnorm
vcov(lin_mix)
fixef(lin_mix)

lin_mix_rSl <- glmmTMB(Trait_z ~ temp_center + (temp_center|ID_fac), data = lin_raw, REML = FALSE)
## In fitTMB(TMBStruc) :
## Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')


# a subset of the raw data for studied for which quadratic relations seem to be most supported
quad_raw <- all_trans %>%
  filter(ID %in% Quad_minAIC$ID) %>%
  mutate(Trait_z = as.numeric(Trait_z)) %>%
  mutate(temp_center2 = temp_center ^2)

lin_quad <- glmmTMB(Trait_z ~ temp_center + temp_center2 + (1|ID_fac), data = quad_raw, REML = FALSE)
summary(lin_quad)

# par-ms for the quadratic model
vcov(lin_quad)
fixef(lin_quad)


# a subset of the raw data for studied for which quadratic relations seem to be most supported
sigm_raw <- all_trans %>%
  filter(ID %in% Sigm_minAIC$ID) %>%
  mutate(Trait_z = as.numeric(Trait_z))

# Not sure how to fit sigmoid now within the mixed-effect framework...

