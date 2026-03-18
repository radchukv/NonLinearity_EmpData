library(ggplot2)
library(magrittr)
library(minpack.lm)
library(patchwork)
library(tidyverse)


# source functions
source('./R/fit_shape.R')

# 1. Read in the data -----------------------------------------------------


## reading in the data per study returned from merging all the climwin outputs prepared for SEM
temp_SEM <- readRDS(file = './data-raw/all_SEM.RDS')

# check that we always have one unique combi per trait, location, species and dem rate STILL TO DO
temp_SEM_uniCombi <- temp_SEM %>%
  distinct(pick(Species, Location, Trait, Demog_rate), .keep_all = TRUE)

nrow(temp_SEM_uniCombi)      ## 303
length(unique(temp_SEM$ID))  ## 309
## 198, 202; here 201 & 202 should be kept (of  201, 198, 202, 199, 203, 200, 204,  205, 206, 207, 208) ; 207 & 201
check_162_163 <- subset(temp_SEM, ID %in% c(201, 198, 202, 199)) %>%
  distinct(pick(ID), .keep_all = TRUE)

# 203, 200, 204, 205, 206, 207, 208
check_162_163 <- subset(temp_SEM, ID %in% c(201, 198, 202, 199, 203, 200, 204)) %>%
  distinct(pick(ID), .keep_all = TRUE)


 # cleaning: some replicate studies on Hydrobates pelagicus, remove
cl_tr <- temp_SEM %>%
  filter(! ID %in% c(198, 202, 199, 203, 200, 204, 205, 206, 208)) %>%
  mutate(ID = as_factor(ID))
## some cleaning: updating the values for percentage to be shown in proportions
cl_tr$Demog_rate_mean[cl_tr$ID == 187] <- cl_tr$Demog_rate_mean[cl_tr$ID == 187] / 100

## check for the median study duration across both trait categories
temp_SEM_dur <- cl_tr %>%
  dplyr::group_by(ID) %>%
  dplyr::summarize(Dur = dplyr::n())
median(temp_SEM_dur$Dur)  ## 20
min(temp_SEM_dur$Dur); max(temp_SEM_dur$Dur) # 7; 63
hist(temp_SEM_dur$Dur)

# number of phenological vs morphological
temp_SEM_uniqID <- cl_tr %>%
  distinct(pick(ID), .keep_all = TRUE)

table(temp_SEM_uniqID$Trait_Categ)
#
# Phenological Morphological
# 145           155

## duration per trait categ
temp_SEM_dur_perTrait <- cl_tr %>%
  dplyr::summarize(.by = c(ID, Trait_Categ), Dur = dplyr::n()) %>%
  group_by(Trait_Categ) %>%
  dplyr::summarise(meanDur = mean(Dur), medianDur = median(Dur))

# Trait_Categ   meanDur medianDur
# 1 Phenological     25.1        24
# 2 Morphological    19.6        16

# 2. Data prep ----------------------------------------------------
rep <- subset(cl_tr, Demog_rate_Categ == 'Reproduction')  # do not drop the leevls here, because otherwise will be difficult to merge the datasets into one again later
surv <- subset(cl_tr, Demog_rate_Categ == 'Survival')
rec <- subset(cl_tr, Demog_rate_Categ == 'Recruitment')

hist(surv$Demog_rate_mean)  ## ok
hist(rep$Demog_rate_mean)  # we would have to round or so to be able to assume a poisson... Not sure.
table(rep$Demog_rate)  ## litter success may be 0/1... Same for HatchingSuccess, Nest Success and Fledgling success. I think it may be best to remove those

# check the potentially binary (reproductive success variables)
rep_success <- subset(rep, Demog_rate %in% c('FledgingSuccess', 'HatchingSuccess', 'LitterSuccess', 'NestSuccess'))
hist(rep_success$Demog_rate_mean)  ## hmm, not only
## check one by one to see which ones are 'true' success to remove those from the data for fitting
rep_Fledlings <- subset(rep, Demog_rate == 'FledgingSuccess') #ok
rep_LitSuc <- subset(rep, Demog_rate == 'LitterSuccess')
hist(rep_LitSuc$Demog_rate_mean)  # this is a true binary, to be removed from the dataset used for fitting !!!
# + in documenting of analyses in the Methods will have to mention that the studies on litter success were removed, thus also
# a lower number of studies in the end

rep_HatchSuc <- subset(rep, Demog_rate == 'HatchingSuccess')
hist(rep_HatchSuc$Demog_rate_mean)  ## also a binary, so remove - as for the litter success
rep_NestSuc <- subset(rep, Demog_rate == 'NestSuccess')
hist(rep_NestSuc$Demog_rate_mean)  ## also a binary, so remove

rep_noSuc <- rep %>%
  filter(! Demog_rate %in% c('HatchingSuccess', 'NestSuccess', 'LitterSuccess'))

hist(rec$Demog_rate_mean)  ## recruitment goes from 0 to ~ 2. See whether I can include these studies in the
## analyses

unique(rec$Demog_rate) ## those that are on juv. survival, and FirstYear survival
# can definitely be included under surv - combine

subs_rec_surv <- subset(rec, Demog_rate %in% c('JuvenileSurvival', 'JuvenileSurvival_Early',
                                               'JuvenileSurvival_Male', 'JuvenileSurvival_Female',
                                               'ChickSurvival', 'FirstYearSurvival',
                                               'OffspringSurvival', 'SurvivalTo6Mnths',
                                               'PupSurvival', 'SurvivalToRecruitment'))
hist(subs_rec_surv$Demog_rate_mean)  ## okay, merge with the surv data

surv_comb <- rbind(surv, subs_rec_surv)
nrow(surv_comb)
surv_comb$Demog_rate_Categ <- 'Survival'
hist(surv_comb$Demog_rate_mean)


tot <- rbind(rep_noSuc, surv_comb)
length(unique(tot$ID))  ## 266


# reproduction and survival with scaled traits separately
# by visual inspection conducted afterwards I saw some of the studies actually have fecundities between 0 and 1,
# remove those: 86, 121, 122, 123, 140, 143, 159, 190, 201, 209, 210, 211, 212, 213, 308, 309, 310, 311, 430, 433,
## 434, 435, 439, 440, 447, 460, 487, 556, 557, 558, 559, 560, 561, 562, 563
rep_noSuc_scaled <- droplevels(rep_noSuc %>%
                                 filter(! is.na(Demog_rate_mean) & ! is.nan(Demog_rate_mean)) %>%
                                 filter(! ID %in% c(86, 121, 122, 123, 140, 143, 159, 190, 201, 209, 210, 211, 212,
                                                    213, 308, 309, 310, 311, 430, 433, 434, 435, 439, 440, 447, 460,
                                                    487, 556, 557, 558, 559, 560, 561, 562, 563)) %>%
                                 mutate(ID = as_factor(ID)) %>%
                                 dplyr::group_by(ID) %>%
                                 dplyr::mutate(Trait_z = scale(Trait_mean)) %>%
                                 dplyr::ungroup())
length(unique(rep_noSuc_scaled$ID)) ## 99

surv_scaled <- droplevels(surv_comb %>%
                            dplyr::filter(! is.na(Demog_rate_mean) & ! is.nan(Demog_rate_mean)) %>%
                                 mutate(ID = as_factor(ID)) %>%
                                 dplyr::group_by(ID) %>%
                                 dplyr::mutate(Trait_z = scale(Trait_mean)) %>%
                                 dplyr::ungroup())
length(unique(surv_scaled$ID))  # 132

# 3.  Exploratory plots --------------------------------------------------

# plotting separately for reproduction
## and now plotting the scaled dem. rate values against the trait
dat_abbrev_rep <- rep_noSuc_scaled %>%
  dplyr::group_by(., ID) %>%
  dplyr::distinct(., ID, Study_Authors, Species, .keep_all = TRUE)


abbrev_d_rep <- dat_abbrev_rep %>%
  dplyr::group_by(., ID) %>%
  tidyr::nest(.data =., cols = ID) %>%
  dplyr::mutate(
    fLet = purrr::map_chr(dat_abbrev_rep$Species, ~substr(., start = 1, stop = 1)),
    sName = purrr::map_chr(dat_abbrev_rep$Species, ~substr(strsplit(as.character(.), split = ' ')[[1]][2], 1, 4)),
    label = paste(Study_Authors, paste(fLet, sName, sep = '.'), Demog_rate, sep = ': ')
  ) %>%
  tidyr::unnest(cols = c(cols))


## plot of non-scaled data
trait_rep <- ggplot(rep_noSuc_scaled, aes(x = Trait_mean, y = Demog_rate_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d_rep,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Fecundity')
print(trait_rep)

ggforce::n_pages(trait_rep)
pdf('./output/output_nonL/Repro_vsTraitRaw.pdf')
for(i in 1:5){
  trait_rep_scale <- ggplot(rep_noSuc_scaled, aes(x = Trait_mean, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_rep,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Fecundity')
  print(trait_rep_scale)
}
dev.off()

# and plot of the scaled data
trait_rep_scale <- ggplot(rep_noSuc_scaled, aes(x = Trait_z, y = Demog_rate_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d_rep,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Fecundity')
print(trait_rep_scale)

ggforce::n_pages(trait_rep_scale)
pdf('./output/output_nonL/Repro_vsTrait_zScaled.pdf')
for(i in 1:5){
  trait_rep_scale <- ggplot(rep_noSuc_scaled, aes(x = Trait_z, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_rep,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Fecundity')
  print(trait_rep_scale)
}
dev.off()



# plotting separately for survival
## and now plotting the scaled dem. rate values against the trait
dat_abbrev_s <- surv_scaled %>%
  dplyr::group_by(., ID) %>%
  dplyr::distinct(., ID, Study_Authors, Species, .keep_all = TRUE)


abbrev_d_surv <- dat_abbrev_s %>%
  dplyr::group_by(., ID) %>%
  tidyr::nest(.data =., cols = ID) %>%
  dplyr::mutate(
    fLet = purrr::map_chr(dat_abbrev_s$Species, ~substr(., start = 1, stop = 1)),
    sName = purrr::map_chr(dat_abbrev_s$Species, ~substr(strsplit(as.character(.), split = ' ')[[1]][2], 1, 4)),
    label = paste(Study_Authors, paste(fLet, sName, sep = '.'), Demog_rate, sep = ': ')
  ) %>%
  tidyr::unnest(cols = c(cols))


## plot of non-scaled data
trait_surv <- ggplot(surv_scaled, aes(x = Trait_mean, y = Demog_rate_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d_surv,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Survival rate')
print(trait_surv)

ggforce::n_pages(trait_surv)
pdf('./output/output_nonL/Surv_vsTraitRaw.pdf')
for(i in 1:7){
  trait_surv_scale <- ggplot(surv_scaled, aes(x = Trait_mean, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_surv,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Survival rate')
  print(trait_surv_scale)
}
dev.off()


# also the same plot with 0-1 for y axis
pdf('./output/output_nonL/Surv_vsTraitRaw_01_limits.pdf')
for(i in 1:7){
  trait_surv_scale <- ggplot(surv_scaled, aes(x = Trait_mean, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_surv,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Survival rate') + ylim(0,1)
  print(trait_surv_scale)
}
dev.off()

# and plot of the scaled data
trait_surv_scale <- ggplot(surv_scaled, aes(x = Trait_z, y = Demog_rate_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d_surv,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Survival rate')
print(trait_surv_scale)

ggforce::n_pages(trait_surv_scale)
pdf('./output/output_nonL/Surv_vsTrait_zScaled.pdf')
for(i in 1:7){
  trait_surv_scale <- ggplot(surv_scaled, aes(x = Trait_z, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_surv,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Survival rate')
  print(trait_surv_scale)
}
dev.off()


# also the same plot with 0-1 for y axis
pdf('./output/output_nonL/Surv_vsTrait_zScaled__01_limits.pdf')
for(i in 1:7){
  trait_surv_scale <- ggplot(surv_scaled, aes(x = Trait_z, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d_surv,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Survival rate') + ylim(0,1)
  print(trait_surv_scale)
}
dev.off()

# 4. Fitting the three shapes --------------------------------------------

## check the function
test <- fit_shape(data = surv_scaled,
          ID = 2, Thresh = 4,
          x = 'Trait_z',
          y = 'Demog_rate_mean',
          surv_rate = TRUE,
          out_folder = './output/output_nonL/shapes_traitdem/')


## applying the function to all repro studies
shapes_fit_rep <- do.call('rbind', lapply(unique(rep_noSuc_scaled$ID), FUN = function(x){fit_shape(data = rep_noSuc_scaled,
                                                                                     ID = x, Thresh = 4,
                                                                                     x = 'Trait_z',
                                                                                     y = 'Demog_rate_mean',
                                                                                     surv_rate  = FALSE,
                                                                                     out_folder = './output/output_nonL/shapes_traitdem/reprod/')}))

table(shapes_fit_rep$mod_minAIC)
## linear quadratic   sigmoid
## 75         8        14

# another way to look at it, though this one more free in taht it may also include cases where one relation clearly
# outperformed the other
pdf('./output/output_nonL/shapes_traitdem/Plot_DeltaAIC_Repro_Linear&Sigmoid.pdf')
hist(shapes_fit_rep$Delta1_AIC[shapes_fit_rep$mod_minAIC %in%
                                 c('linear', 'sigmoid')],
     xlab = 'Delta AICc',
     main = '', col = 'darkgrey')  ## ok - as expected, more spread
dev.off()


## save the output with the deltaAIC of 4, for drawing par-rs in simu model
saveRDS(shapes_fit_rep, file = './output/output_nonL/shapes_traitdem/reprod/Shapes_traitRepro_4DeltaAIC.RDS')

## applying the function to all studies on survival
shapes_fit_surv <- do.call('rbind', lapply(unique(surv_scaled$ID), FUN = function(x){fit_shape(data = surv_scaled,
                                                                                         ID = x, Thresh = 4,
                                                                                         x = 'Trait_z',
                                                                                         y = 'Demog_rate_mean',
                                                                                         surv_rate  = TRUE,
                                                                                         out_folder = './output/output_nonL/shapes_traitdem/survival/')}))



table(shapes_fit_surv$mod_minAIC)
## linear quadratic   sigmoid
##    56        10        66

# it may also include cases where one relation celarly outperformed the other
pdf('./output/output_nonL/shapes_traitdem/Plot_DeltaAIC_Survival_Linear&Sigmoid.pdf')
hist(shapes_fit_surv$Delta1_AIC[shapes_fit_surv$mod_minAIC %in%
                                  c('linear', 'sigmoid')],
     col = 'darkgrey',
     main = '',
     xlab = 'Delta AICc')
dev.off()

## save the output with the deltaAIC of 4, for drawing par-rs in simu model
saveRDS(shapes_fit_surv, file = './output/output_nonL/shapes_traitdem/survival/Shapes_traitSurv_4DeltaAIC.RDS')


## So, we generally have difficulties discriminating sigmoid from linear given shallow slopes for the latter



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####              explore fitted relationships                               ####

# I. for survival
## 1. linear relations
Lin_minAIC_s <- shapes_fit_surv %>%
  filter(mod_minAIC == 'linear') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'survival')

# explore the intercepts and slopes
hist(Lin_minAIC_s$Int_Lin)
min(Lin_minAIC_s$Int_Lin)  # 0.1159365
max(Lin_minAIC_s$Int_Lin)  # 0.9413462
mean(Lin_minAIC_s$Int_Lin)  # 0.6085007

# These intercepts cannot be 0 -the survival rates are not scaled...
# So the intercept here reflects the mean survival at '0' temperature
# So, it may be better to center the survival rates prior to fitting
# these relations? To have the slopes while intercept is aroiund 0????

# Same for reproduction: if I center the response variables before fitting the
#relaitons, then the intercept is at 0, but the slopes would be preserved

hist(Lin_minAIC_s$Beta_Lin)
median(Lin_minAIC_s$Beta_Lin)  # 0.000769066 - very small as it is ont he scale of survival rate!
mean(Lin_minAIC_s$Beta_Lin)  # 0.003906411
## POSITIVE mainly???

# 2. quadratic relations
Quad_minAIC_s <- shapes_fit_surv %>%
  filter(mod_minAIC == 'quadratic') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'survival')

# explore the intercepts and slopes
hist(Quad_minAIC_s$Int_Quad)
min(Quad_minAIC_s$Int_Quad)  # -0.4781995
max(Quad_minAIC_s$Int_Quad)  # 2.433744

hist(Quad_minAIC_s$Beta_Quad)
median(Quad_minAIC_s$Beta_Quad)  # -0.04615922
mean(Quad_minAIC_s$Beta_Quad)  # -0.02990958

# 3. sigmoid relations
Sigm_minAIC_s <- shapes_fit_surv %>%
  filter(mod_minAIC == 'sigmoid') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'survival')

# explore the intercepts and slopes
hist(Sigm_minAIC_s$Int_Sigm)
mean(Sigm_minAIC_s$Int_Sigm)  ## 0.04397168 = okay, that is because here we have the logit transform....
min(Sigm_minAIC_s$Int_Sigm)  # -2.565511
max(Sigm_minAIC_s$Int_Sigm)  # 2.985382
# Okay, fine to assume that these are hardly different from 0

hist(Sigm_minAIC_s$Beta_Sigm)
median(Sigm_minAIC_s$Beta_Sigm)  #  -0.005771454
mean(Sigm_minAIC_s$Beta_Sigm)  # -0.04991454



# II. for reproduction
## 1. linear relations
Lin_minAIC_r <- shapes_fit_rep %>%
  filter(mod_minAIC == 'linear') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'fecundity')

# explore the intercepts and slopes
hist(Lin_minAIC_r$Int_Lin)
min(Lin_minAIC_r$Int_Lin)  # 0.1159365
max(Lin_minAIC_r$Int_Lin)  # 0.9413462
mean(Lin_minAIC_r$Int_Lin)
# of course, the intercept dpeends on the study - the mean repro when temp = 0...

hist(Lin_minAIC_r$Beta_Lin)
median(Lin_minAIC_r$Beta_Lin)  # 0.02141382
mean(Lin_minAIC_r$Beta_Lin)  # 0.04087355
## POSITIVE in half of the cases - depends on the trait, likely - may be worth checking for the interpretation
# of the overall results later on

# 2. quadratic relations
Quad_minAIC_r <- shapes_fit_rep %>%
  filter(mod_minAIC == 'quadratic') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'fecundity')

# explore the intercepts and slopes
hist(Quad_minAIC_r$Int_Quad)
min(Quad_minAIC_r$Int_Quad)  # 0.3810891
max(Quad_minAIC_r$Int_Quad)  # 10
mean(Quad_minAIC_r$Int_Quad)  # 5.11


hist(Quad_minAIC_r$Beta_Quad)
median(Quad_minAIC_r$Beta_Quad)  # 0.04068422
mean(Quad_minAIC_r$Beta_Quad)  # 0.05084021

# 3. sigmoid relations
Sigm_minAIC_r <- shapes_fit_rep %>%
  filter(mod_minAIC == 'sigmoid') %>%
  dplyr::select(ID, shape = mod_minAIC, Int_Lin:Beta_Sigm) %>%
  mutate(relationship = 'fecundity')

# explore the intercepts and slopes
hist(Sigm_minAIC_r$Int_Sigm)
mean(Sigm_minAIC_r$Int_Sigm)  ## 0.314497
min(Sigm_minAIC_r$Int_Sigm)  # 0.133682
max(Sigm_minAIC_r$Int_Sigm)  # 0.6208231


hist(Sigm_minAIC_r$Beta_Sigm)
median(Sigm_minAIC_r$Beta_Sigm)  #  0.001588862
mean(Sigm_minAIC_r$Beta_Sigm)  # 0.007202677
nrow(Sigm_minAIC_r)


# look ath the histogram of slopes across all relations
shapes <-
  Lin_minAIC_r %>%
  ## import + join data
  bind_rows(Sigm_minAIC_r)  %>%
  bind_rows(Quad_minAIC_r) %>%
  bind_rows(Lin_minAIC_s) %>%
  bind_rows(Sigm_minAIC_s) %>%
  bind_rows(Quad_minAIC_s) %>%

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

summaries_slopes <- shapes %>%
  group_by(relationship, shape) %>%
  summarize(mean = mean(beta, na.rm = TRUE),
            median = median(beta, na.rm = TRUE))

# here still add the mean / median - a plot like this, but cosmetics improved - in the main text
pdf('./output/output_nonL/shapes_traitdem/Histogramme_Slopes_Surv&Repro.pdf',
    width = 10)
ggplot(shapes, aes(beta)) +
  geom_histogram(aes(y = after_stat(density))) +
  facet_grid(rows = vars(relationship), cols = vars(shape)) + theme_bw() +
  geom_vline(data = summaries_slopes, aes(xintercept = mean),
             col = 'grey', lwd = 1.5) +
  geom_vline(data = summaries_slopes, aes(xintercept = median),
             col = 'red', lwd = 1.5) +
  theme(strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()


summaries_interc <- shapes %>%
  group_by(relationship, shape) %>%
  summarize(mean = mean(int, na.rm = TRUE),
            median = median(int, na.rm = TRUE))


shapes_surv <- subset(shapes, relationship == 'survival')
summaries_surv_int <- subset(summaries_interc, relationship == 'survival')
# for the plot it makes more sense to actually have the survival and
# repro plotted separately otherwise the free_x does not help much...
# here still add the mean / median - a plot like this, but cosmetics improved - in the main text
pdf('./output/output_nonL/shapes_traitdem/Histogramme_Intercepts_Surv.pdf',
    width = 10)
ggplot(shapes_surv, aes(int)) +
  geom_histogram(aes(y = after_stat(density))) +
  facet_wrap(vars(shape), scales = 'free_x') + theme_bw() +
  geom_vline(data = summaries_surv_int, aes(xintercept = mean),
             col = 'grey', lwd = 1.5) +
  geom_vline(data = summaries_surv_int, aes(xintercept = median),
             col = 'red', lwd = 1.5) +
  theme(strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()


# and reproduction
shapes_rep <- subset(shapes, relationship == 'fecundity')
summaries_rep_int <- subset(summaries_interc, relationship == 'fecundity')

pdf('./output/output_nonL/shapes_traitdem/Histogramme_Intercepts_Rep.pdf',
    width = 10)
ggplot(shapes_rep, aes(int)) +
  geom_histogram(aes(y = after_stat(density)), bins = 7) +
  facet_wrap(vars(shape), scales = 'free') + theme_bw() +
  geom_vline(data = summaries_rep_int, aes(xintercept = mean),
             col = 'grey', lwd = 1.5) +
  geom_vline(data = summaries_rep_int, aes(xintercept = median),
             col = 'red', lwd = 1.5) +
  theme(strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
dev.off()


# 5. Correlations between intercepts and slopes ---------------------------

##  check correlations among the slopes and intercepts
Lin_traitSurv <- shapes_fit_surv %>%
  filter(mod_minAIC == 'linear')
ggplot(Lin_traitSurv, aes(x = Int_Lin, y = Beta_Lin)) + geom_point()
cor.test(Lin_traitSurv$Int_Lin, Lin_traitSurv$Beta_Lin)  ## uncorrelated

Quad_traitSurv <- shapes_fit_surv %>%
  filter(mod_minAIC == 'quadratic')
ggplot(Quad_traitSurv, aes(x = Int_Quad, y = Beta_Quad)) + geom_point()
cor.test(Quad_traitSurv$Int_Quad, Quad_traitSurv$Beta_Quad)  ## uncorrelated but power is low

ggplot(Quad_traitSurv, aes(x = Int_Quad, y = Beta2_Quad)) + geom_point()
cor.test(Quad_traitSurv$Int_Quad, Quad_traitSurv$Beta2_Quad)  ## correlated


Sigm_traitSurv <- shapes_fit_surv %>%
  filter(mod_minAIC == 'sigmoid')
ggplot(Sigm_traitSurv, aes(x = Int_Sigm, y = Beta_Sigm)) + geom_point()
cor.test(Sigm_traitSurv$Int_Sigm, Sigm_traitSurv$Beta_Sigm)  ## Strongly correlated

# remove the very different value from the set
Sigm_traitSurv <- Sigm_traitSurv[Sigm_traitSurv$Int_Sigm > -100, ]
ggplot(Sigm_traitSurv, aes(x = Int_Sigm, y = Beta_Sigm)) + geom_point()
cor.test(Sigm_traitSurv$Int_Sigm, Sigm_traitSurv$Beta_Sigm)  ## Still strongly correlated but te direction now flipped


## and now correlations for the set with reproduction
Lin_traitRep <- shapes_fit_rep %>%
  filter(mod_minAIC == 'linear')
ggplot(Lin_traitRep, aes(x = Int_Lin, y = Beta_Lin)) + geom_point()
cor.test(Lin_traitRep$Int_Lin, Lin_traitRep$Beta_Lin)  ## uncorrelated

Quad_traitRep <- shapes_fit_rep %>%
  filter(mod_minAIC == 'quadratic')
ggplot(Quad_traitRep, aes(x = Int_Quad, y = Beta_Quad)) + geom_point()
cor.test(Quad_traitRep$Int_Quad, Quad_traitRep$Beta_Quad)  ## uncorrelated but power is low

ggplot(Quad_traitRep, aes(x = Int_Quad, y = Beta2_Quad)) + geom_point()
cor.test(Quad_traitRep$Int_Quad, Quad_traitRep$Beta2_Quad)  ## correlated


Sigm_traitRep <- shapes_fit_rep %>%
  filter(mod_minAIC == 'sigmoid')
ggplot(Sigm_traitRep, aes(x = Int_Sigm, y = Beta_Sigm)) + geom_point()
cor.test(Sigm_traitRep$Int_Sigm, Sigm_traitRep$Beta_Sigm)

# remove the very different value from the dataset
Sigm_traitRep <- Sigm_traitRep[Sigm_traitRep$Int_Sigm < 2, ]
ggplot(Sigm_traitSurv, aes(x = Int_Sigm, y = Beta_Sigm)) + geom_point()
cor.test(Sigm_traitSurv$Int_Sigm, Sigm_traitSurv$Beta_Sigm)  ## Still strongly correlated



# 5. Plots of the fitted relations  ------------------------------------

# 1. survival
# linear relation
data <- data.frame(Trait = seq(-3, 3, by = 0.1),
                   Surv = seq(0, 1, length.out = length(seq(-3, 3, by = 0.1))))
subs_surv_linear <- subset(shapes, shape == 'linear' & relationship == 'survival')
lin_rel_trait_surv <- ggplot(data, aes(x = Trait, y = Surv)) +
  lims(x = c(-3, 3), y =  c(0, 1)) + ylab('Survival') +
  geom_blank() +
  geom_abline(data = subs_surv_linear, aes(intercept = int, slope = beta),
              alpha = 0.7, col = 'grey', lwd = 2) +
  geom_abline(aes(intercept = median(subs_surv_linear$int),
                  slope = median(subs_surv_linear$beta)
  ), col = 'black', lwd = 4) +
  scale_colour_brewer(palette = 'Dark2') +
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
lin_rel_trait_surv

pdf('./output/output_nonL/shapes_traitdem/plot_lin_relations_trait_surv.pdf')
lin_rel_trait_surv
dev.off()

# quadratic relation

subs_surv_quad <- subset(shapes, shape == 'quadratic' & relationship == 'survival')
data_quad_s <- data.frame(Trait = rep(seq(-2, 2, by = 0.1), nrow(subs_surv_quad)),
                        ID =  rep(unique(subs_surv_quad$ID),
                                  each = length(seq(-2, 2, by = 0.1))),
                        Surv = 0)

for(i in unique(subs_surv_quad$ID)){
  data_quad_s$Surv[data_quad_s$ID == i] <- 1 / (1 + exp(subs_surv_quad$int[subs_surv_quad$ID == i] +
                                                          data_quad_s$Trait[data_quad_s$ID == i]*subs_surv_quad$beta[subs_surv_quad$ID == i] +
    data_quad_s$Trait[data_quad_s$ID == i]^2*subs_surv_quad$beta2[subs_surv_quad$ID == i]))
}

quad_rel_trait_surv <- ggplot(data_quad_s, aes(x = Trait, y = Surv, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd =2) + ylab("Survival") +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))

quad_rel_trait_surv

pdf('./output/output_nonL/shapes_traitdem/plot_quad_relations_trait_surv.pdf')
quad_rel_trait_surv
dev.off()

# and sigmoidal
# and, finally, the plot for the sigmoid relation
subs_sigm_s <- subset(shapes, shape == 'sigmoid' & relationship == 'survival')
data_sigm_s <- data.frame(Trait = rep(seq(-2, 2, by = 0.1), nrow(subs_sigm)),
                        ID =  rep(unique(subs_sigm$ID), each = length(seq(-2, 2, by = 0.1))),
                        Surv = 0)

for(i in unique(subs_sigm_s$ID)){
  data_sigm_s$Surv[data_sigm_s$ID == i] <-
    (1/(1+exp(subs_sigm_s$beta[subs_sigm_s$ID == i] * data_sigm_s$Trait[data_sigm_s$ID == i])))

}

# diff colours
sigm_rel_trait_surv <- ggplot(data_sigm_s, aes(x = Trait, y = Surv, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd = 2) + ylab("Survival") +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
sigm_rel_trait_surv

pdf('./output/output_nonL/shapes_traitdem/plot_sigm_relations_trait_surv.pdf')
sigm_rel_trait_surv
dev.off()


# grey shade
sigm_rel_trait_surv_grey <- ggplot(data_sigm_s, aes(x = Trait, y = Surv, group = as.factor(ID))) +
  geom_line(col = 'darkgrey') + ylab('Survival') +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 12),
                     axis.text = element_text(size = 10))
sigm_rel_trait_surv_grey


pdf('./output/output_nonL/shapes_traitdem/plot_sigm_relations_trait_surv_GREY.pdf')
sigm_rel_trait_surv_grey
dev.off()



# 2. fecundity
# linear relation
data <- data.frame(Trait = seq(-3, 3, by = 0.1),
                   Fecundity = seq(0, 1, length.out = length(seq(-3, 3, by = 0.1))))
subs_rep_linear <- subset(shapes, shape == 'linear' & relationship == 'fecundity')
lin_rel_trait_fec <- ggplot(data, aes(x = Trait, y = Fecundity)) +
  lims(x = c(-3, 3), y =  c(0, 14)) + ylab('Fecundity') +
  geom_blank() +
  geom_abline(data = subs_rep_linear, aes(intercept = int, slope = beta),
              alpha = 0.7, col = 'darkgrey', lwd = 2) +
  geom_abline(aes(intercept = median(subs_rep_linear$int),
                  slope = median(subs_rep_linear$beta)
  ), col = 'black', lwd = 4) +
#  scale_colour_brewer(palette = 'Dark2') +
  theme_bw() + theme(axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
lin_rel_trait_fec

pdf('./output/output_nonL/shapes_traitdem/plot_lin_relations_trait_fec.pdf')
lin_rel_trait_fec
dev.off()

# quadratic relation

subs_rep_quad <- subset(shapes, shape == 'quadratic' & relationship == 'fecundity')
data_quad_fec <- data.frame(Trait = rep(seq(-2, 2, by = 0.1), nrow(subs_rep_quad)),
                          ID =  rep(unique(subs_rep_quad$ID),
                                    each = length(seq(-2, 2, by = 0.1))),
                          Fecundity = 0)

for(i in unique(subs_rep_quad$ID)){
  data_quad_fec$Fecundity[data_quad_fec$ID == i] <- subs_rep_quad$int[subs_rep_quad$ID == i] +
                                                          data_quad_fec$Trait[data_quad_fec$ID == i]*subs_rep_quad$beta[subs_rep_quad$ID == i] +
                                                          data_quad_fec$Trait[data_quad_fec$ID == i]^2*subs_rep_quad$beta2[subs_rep_quad$ID == i]
}

quad_rel_trait_fec <- ggplot(data_quad_fec, aes(x = Trait, y = Fecundity, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd =2) + ylab("Fecundity") +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))

quad_rel_trait_fec

pdf('./output/output_nonL/shapes_traitdem/plot_quad_relations_trait_fec.pdf')
quad_rel_trait_fec
dev.off()

# and sigmoidal
# and, finally, the plot for the sigmoid relation
subs_sigm_f <- subset(shapes, shape == 'sigmoid' & relationship == 'fecundity')
data_sigm_f <- data.frame(Trait = rep(seq(-2, 2, by = 0.1), nrow(subs_sigm_f)),
                          ID =  rep(unique(subs_sigm_f$ID), each = length(seq(-2, 2, by = 0.1))),
                          Fecundity = 0)

for(i in unique(subs_sigm_f$ID)){
  data_sigm_f$Fecundity[data_sigm_f$ID == i] <-
    (1/(1+exp(-5*(subs_sigm_f$int[subs_sigm_f$ID == i] +
                    subs_sigm_f$beta[subs_sigm_f$ID == i] * data_sigm_f$Trait[data_sigm_f$ID == i]))) - 0.5)*4

}

# diff colours
sigm_rel_trait_fec <- ggplot(data_sigm_f, aes(x = Trait, y = Fecundity, group = as.factor(ID))) +
  geom_line(aes(col = as.factor(ID)), lwd = 2) + ylab("Fecundity") +
  theme_bw() + theme(legend.position = 'none',
                     axis.title = element_text(size = 18),
                     axis.text = element_text(size = 12))
sigm_rel_trait_fec

pdf('./output/output_nonL/shapes_traitdem/plot_sigm_relations_trait_fec.pdf')
sigm_rel_trait_fec
dev.off()

pdf('./output/output_nonL/shapes_traitdem/Plot_allFitted_relations_Surv&Rep.pdf', width = 10)
(lin_rel_trait_surv + sigm_rel_trait_surv + quad_rel_trait_surv +
    lin_rel_trait_fec + sigm_rel_trait_fec + quad_rel_trait_fec) +
  plot_layout(nrow = 2) +  plot_annotation(tag_levels = "a")
dev.off()

