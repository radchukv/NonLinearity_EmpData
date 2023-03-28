library(ggplot2)
library(magrittr)
library(minpack.lm)
library(patchwork)


# source functions
source('./R/fit_shape.R')

# 1. Read in the data -----------------------------------------------------


## reading in the data per study returned from merging all the climwin outputs prepared for SEM
temp_SEM <- readRDS(file = './data-raw/all_SEM.RDS')

## some cleaning: updating the values for percentage to be shown in proportions
temp_SEM$Demog_rate_mean[temp_SEM$ID == 187] <- temp_SEM$Demog_rate_mean[temp_SEM$ID == 187] / 100

## check for the median study duration across both trait categories
temp_SEM_dur <- temp_SEM %>%
  dplyr::group_by(ID) %>%
  dplyr::summarize(Dur = dplyr::n())
median(temp_SEM_dur$Dur)  ## 20
min(temp_SEM_dur$Dur); max(temp_SEM_dur$Dur)
hist(temp_SEM_dur$Dur)

# 2. Data prep ----------------------------------------------------
## maybe split data prep from exploratory plots

rep <- subset(temp_SEM, Demog_rate_Categ == 'Reproduction')
surv <- subset(temp_SEM, Demog_rate_Categ == 'Survival')
rec <- subset(temp_SEM, Demog_rate_Categ == 'Recruitment')

hist(surv$Demog_rate_mean)  ## ok
hist(rep$Demog_rate_mean)
hist(rec$Demog_rate_mean)  ## recruitment goes from 0 to ~ 2. See whether I can include these studies in the
## analyses

unique(rec$Demog_rate) ## those thta are on juv. survival, and FirstYear survival
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

tot <- rbind(rep, surv_comb)
length(unique(tot$ID))  ## 288


## transform the demog_rate value to have it on the linear scale (for probabilities apply logit
## and for counts apply log

## I still have to subset the dataset, so as to have unique values per\
## combi trait-dem rate-location-species

nrow(tot)  ## 288
check_unique <- tot %>%
  dplyr::distinct(Species, Location, Trait, Demog_rate, .keep_all = TRUE)

nrow(check_unique)  ## 282 - so these should be fine to use

tot_unique <- tot %>%
  dplyr::filter(ID %in% check_unique$ID)

nrow(tot_unique)
tot_scale <- tot_unique %>%
  dplyr::filter(! is.na(Demog_rate_mean) & ! is.nan(Demog_rate_mean)) %>%
  dplyr::mutate(DemR_value = dplyr::case_when(
    Demog_rate_Categ == 'Survival' ~ Demog_rate_mean/ (1- Demog_rate_mean),
    Demog_rate_Categ == 'Reproduction' ~ Demog_rate_mean),
    DemR_value = log(DemR_value)) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Trait_z = scale(Trait_mean)) %>%
  dplyr::ungroup()## otherwise in one go applying log directly on odds and on fecundities
## leads to errors when log of 0 is asked for...


# 3.  Exploratory plots --------------------------------------------------

## and now plotting the scaled dem. rate values against the trait
dat_abbrev <- tot_scale %>%
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
  tidyr::unnest()


## plot of scaled data
trait_dem_scale <- ggplot(tot_scale, aes(x = Trait_z, y = DemR_value)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Demographic rate')
print(trait_dem_scale)

ggforce::n_pages(trait_dem_scale)
pdf('./output/output_nonL/DemRateScaled_vsTraitRaw.pdf')
for(i in 1:15){
  trait_dem_scale <- ggplot(tot_scale, aes(x = Trait_mean, y = DemR_value)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Demographic value, scaled')
  print(trait_dem_scale)
}
dev.off()


## using both scaled trait and scaled dem. rate
## plot of scaled data
trait_dem_Bothscale <- ggplot(tot_scale, aes(x = Trait_z, y = DemR_value)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Demographic rate')
print(trait_dem_Bothscale)

ggforce::n_pages(trait_dem_Bothscale)
pdf('./output/output_nonL/DemRatevsTrait_BothScaled.pdf')
for(i in 1:15){
  trait_dem_scale <- ggplot(tot_scale, aes(x = Trait_z, y = DemR_value)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, scaled') +
    ylab('Demographic value, scaled')
  print(trait_dem_scale)
}
dev.off()


## plot of untransformed, raw data
trait_dem <- ggplot(tot_scale, aes(x = Trait_mean, y = Demog_rate_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1.5, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Trait') +
  ylab('Demographic rate')
print(trait_dem)

pdf('./output/output_nonL/DemRatevsTrait_BothRawUnscaled.pdf')
for(i in 1:15){
  trait_dem <- ggplot(tot_scale, aes(x = Trait_mean, y = Demog_rate_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none') + xlab('Trait, raw') +
    ylab('Demographic rate, raw')
  print(trait_dem)
}
dev.off()



# 4. Fitting the three shapes --------------------------------------------

## check the function
test <- fit_shape(data = tot_scale,
          ID = 86, Thresh = 4,
          x = 'Trait_z',
          y = 'DemR_value',
          classic_sigm = FALSE,
          out_folder = './output/output_nonL/shapes_traitdem/')


## applying the function to all repro studies
reprod <- tot_scale %>%
 filter(Demog_rate_Categ =='Reproduction')
shapes_fit_rep <- do.call('rbind', lapply(unique(reprod$ID), FUN = function(x){fit_shape(data = reprod,
                                                                                     ID = x, Thresh = 4,
                                                                                     x = 'Trait_z',
                                                                                     y = 'DemR_value',
                                                                                     classic_sigm = FALSE,
                                                                                     out_folder = './output/output_nonL/shapes_traitdem/reprod/')}))

table(shapes_fit_rep$Selected)
## linear  linear/sigmoid      quadratic        sigmoid
# 10             120              5              2

table(shapes_fit_rep$mod_minAIC)
## linear quadratic   sigmoid
## 71        14        52

## save the output with the deltaAIC of 4, for drawing par-rs in simu model
saveRDS(shapes_fit_rep, file = './output/output_nonL/shapes_traitdem/reprod/Shapes_traitRepro_4DeltaAIC.RDS')

## applying the function to all studies on survival
surv <- tot_scale %>%
  filter(Demog_rate_Categ =='Survival')

shapes_fit_surv <- do.call('rbind', lapply(unique(surv$ID), FUN = function(x){fit_shape(data = surv,
                                                                                         ID = x, Thresh = 4,
                                                                                         x = 'Trait_z',
                                                                                         y = 'DemR_value',
                                                                                         classic_sigm = TRUE,
                                                                                         out_folder = './output/output_nonL/shapes_traitdem/survival/')}))


table(shapes_fit_surv$Selected)  # not a single time for sigmoid. Quadratic rather rare (11 out of 272)
##  linear linear/sigmoid      quadratic
##   54             71              2

table(shapes_fit_surv$mod_minAIC)
## linear quadratic   sigmoid
##   99        10        18

## save the output with the deltaAIC of 4, for drawing par-rs in simu model
saveRDS(shapes_fit_surv, file = './output/output_nonL/shapes_traitdem/survival/Shapes_traitSurv_4DeltaAIC.RDS')


## So, we generally have difficulties discriminating sigmoid from linear given shallow slopes for the latter


# 5. some histograms of obtained par-s ------------------------------------

## this section has to be revised --> drop altogether?  Or do the same but for the models with min AIC
# still visualize these par-rs to grasp how they look like


## a plot for the MS with the estimated betas and intercepts
## sigm
pl_Sigm_b <- ggplot(subset(shapes_fit_rep, Selected == 'linear/sigmoid'),
       aes(x = Beta_Sigm)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.y = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Beta_Sigm[shapes_fit_rep$Selected
                                                        == 'linear/sigmoid'],
                                probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Slope')

quantile(shapes_fit_rep$Beta_Sigm[shapes_fit_rep$Selected == 'linear/sigmoid'],
           probs = c(0.05, 0.95))  ## -0.09577593  0.13938209

pl_Sigm_int <- ggplot(subset(shapes_fit_rep, Selected == 'linear/sigmoid'),
       aes(x = Int_Sigm)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.y = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Int_Sigm[shapes_fit_rep$Selected
                                                        == 'linear/sigmoid'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Intercept') + labs(title = 'Sigmoid')


## linear
pl_Lin_b <- ggplot(subset(shapes_fit_rep, Selected %in% c('linear/sigmoid', 'linear')),
                    aes(x = Beta_Lin)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.x = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Beta_Lin[shapes_fit_rep$Selected
                                                        %in% c('linear/sigmoid', 'linear')],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2)

quantile(shapes_fit_rep$Beta_Lin[shapes_fit_rep$Selected %in% c('linear/sigmoid', 'linear')],
         probs = c(0.05, 0.95))  ## -0.2043896  0.2638683

pl_Lin_int <- ggplot(subset(shapes_fit_rep, Selected %in% c('linear/sigmoid', 'linear')),
                      aes(x = Int_Lin)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.x = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Int_Lin[shapes_fit_rep$Selected
                                                      %in% c('linear/sigmoid', 'linear')],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  labs(title = 'Linear')


## quadr
pl_Quad_b <- ggplot(subset(shapes_fit_rep, Selected == 'quadratic'),
                   aes(x = Beta_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Beta_Quad[shapes_fit_rep$Selected
                                                       == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2)

quantile(shapes_fit_rep$Beta_Quad[shapes_fit_rep$Selected == 'quadratic'],
         probs = c(0.05, 0.95))  ## -0.3502028  0.1509256

pl_Quad_int <- ggplot(subset(shapes_fit_rep, Selected == 'quadratic'),
                     aes(x = Int_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Int_Quad[shapes_fit_rep$Selected == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  labs(title = 'Quadratic')



pl_Quad_b2 <- ggplot(subset(shapes_fit_rep, Selected == 'quadratic'),
                    aes(x = Beta2_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2))) +
  geom_vline(xintercept = quantile(shapes_fit_rep$Beta2_Quad[shapes_fit_rep$Selected
                                                        == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Quadratic slope')



## plot all plots together

des_lay <- "
123#
1237
4567
456#"

pdf('./output_nonL/EstimatedParS_DemTraitRel_sTraitChange_DeltaAIC4.pdf',
    width = 10)
pl_Lin_int + pl_Sigm_int + pl_Quad_int +
  pl_Lin_b + pl_Sigm_b + pl_Quad_b + pl_Quad_b2 +
  plot_layout(design = des_lay) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = '(',
                  tag_suffix = ')') &
  theme(plot.tag = element_text(face = 'bold'),
        plot.margin = margin(c(0.2, 0.4, 0.2, 0.2)))
dev.off()

## subset those that have linear/sigmoid as identified relation, to see what deltaAIC is
lin_sigm <- droplevels(subset(shapes_fit_rep, Selected == 'linear/sigmoid'))

pdf('./output_nonL/SupplFig_DistrAIC_Sigmoid_Linear.pdf')
hist(lin_sigm$Delta1_AIC, main ='', xlab = expression(paste(Delta, 'AIC')),
     col = 'grey')
abline(v = mean(lin_sigm$Delta1_AIC), col = 'blue', lwd = 2)
abline(v = median(lin_sigm$Delta1_AIC), col = 'blue', lwd = 3, lty = 2)
dev.off()
hist(lin_sigm$Delta2_AIC[lin_sigm$mod_minAIC == 'quadratic'])  ## so, function works fine
table(lin_sigm$mod_minAIC)

