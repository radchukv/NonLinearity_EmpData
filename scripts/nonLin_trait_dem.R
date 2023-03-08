library(ggplot2)
library(magrittr)
library(minpack.lm)
library(patchwork)

## setting the ggplot theme
# ggplot2::theme_set(theme_classic())

## test fit_shape for trait-temperature relations
check <- fit_shape(data = all_trans,
          x = 'Temp_z',
          y = 'Trait_z',
          ID = 1,
          Thresh = 2,
          out_folder = './output/output_nonL/shapes/')


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

# I. Data prep ----------------------------------------------------
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

## not sure if I still have to pre-sieve the dataset, so as to have unique values per
## study_authors of species, location, trait, dem. rate - can check... I think this won't be necessary
## because at the moment the uniqueness of a study is defined by a combi trait-dem rate-location-species
## but can just check

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
  dplyr::mutate(Trait_value = scale(Trait_mean)) %>%
  dplyr::ungroup()## otherwise in one go applying log directly on odds and on fecundities
## leads to errors when log of 0 is asked for...


# II.  Exploratory plots --------------------------------------------------

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
trait_dem_scale <- ggplot(tot_scale, aes(x = Trait_mean, y = DemR_value)) +
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
pdf('./output_nonL/DemRatevsTrait_scaled.pdf')
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
    theme(legend.position = 'none')
  print(trait_dem_scale)
}
dev.off()


## suing both scaled trait and scaled dem. rate
## plot of scaled data
trait_dem_Bothscale <- ggplot(tot_scale, aes(x = Trait_value, y = DemR_value)) +
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
pdf('./output/output_nonL/DemRatevsTrait_Bothscaled.pdf')
for(i in 1:15){
  trait_dem_scale <- ggplot(tot_scale, aes(x = Trait_value, y = DemR_value)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')
  print(trait_dem_scale)
}
dev.off()


## plot of unscaled data
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

pdf('./output_nonL/DemRatevsTrait_rawUnscaled.pdf')
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
    theme(legend.position = 'none')
  print(trait_dem)
}
dev.off()




# III. Fitting the three shapes --------------------------------------------


## I should scale the traits to mean  = 0 and sd = 2,  ## so that it is closer to how things are simulated...
##  a function to assess the fit of diff. shapes to the data
fit_shape <- function(data = tot_scale,
                      ID = 1,
                      Thresh = 2,
                      out_folder = './output_nonL/'){
  sub_data <- droplevels(data[data$ID ==  ID, ])

  ## have to catch the errors while fitting
  tt.error.linRel <- tryCatch(linRel <- nlsLM(DemR_value ~ interc +  beta * Trait_value,
                                              start = list(interc = 0, beta = 1), upper = c(10, 10),
                                              data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                              error=function(e) e)
  if(is(tt.error.linRel,"error")){
    warning(cat('error when fitting linear relation \n',
                tt.error.linRel[1]$message, '\n'))
  }

  if(! is(tt.error.linRel,"error")){
    if(is.infinite(MuMIn::AICc(linRel))){
    warning(cat('infinite AIC for linear relation \n'))
    }
  }
  tt.error.quadRel <- tryCatch(quadRel <- nlsLM(DemR_value ~ interc + beta * Trait_value + beta2 * Trait_value^2,
                                                start = list(interc = 0, beta = 1, beta2 = 1), upper = c(10, 10, 10),
                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                               error=function(e) e)
  if(is(tt.error.quadRel,"error")){
    warning(cat('error when fitting quadratic relation \n',
                tt.error.quadRel[1]$message, '\n'))
  }
  if(! is(tt.error.quadRel,"error")) {
    if(is.infinite(MuMIn::AICc(quadRel))){
    warning(cat('infinite AIC for quadratic relation \n'))
    }
  }

  tt.error.sigmRel <- tryCatch(sigmRel <- nlsLM(DemR_value ~ (1/(1+exp(-5*beta*Trait_value)) - 0.5)*4 + interc,
                                                start = list(interc = 0, beta = 0), upper = c(10, 10),
                                                data = sub_data, algorithm  = "LM",control = list(maxiter = 200)),
                               error=function(e) e)
  if(is(tt.error.sigmRel,"error")){
    warning(cat('error when fitting sigmoid relation \n',
                tt.error.sigmRel[1]$message, '\n'))
  }
  if(! is(tt.error.sigmRel,"error")) {
    if(is.infinite(MuMIn::AICc(sigmRel))){
    warning(cat('infinite AIC for sigmoid relation \n'))
    }
  }

  ## saving a plot as a pdf
  if(! (is(tt.error.linRel,"error") | is(tt.error.quadRel,"error") |
        is(tt.error.sigmRel,"error") )){
    pdf(paste0(out_folder, ID, '_', unique(sub_data$Study_Authors),
               '_', unique(sub_data$Species), '.pdf'))
    plot(sub_data$Trait_value, sub_data$DemR_value,
         ylab = paste0(unique(sub_data$Demog_rate)),
         xlab = paste0(unique(sub_data$Trait)))
    points(sub_data$Trait_value, fitted(linRel), pch = 19, col = 'red')
    points(sub_data$Trait_value, fitted(sigmRel), pch = 21, col = 'grey')
    points(sub_data$Trait_value, fitted(quadRel), pch = 19, col = 'blue')
    legend('bottomright', c('linear', 'sigmoid', 'quadratic'),
           col = c('red', 'grey', 'blue'), pch = c(19, 21, 19))
    dev.off()

    surv_sub <- subset(sub_data, Demog_rate_Categ == 'Survival')
    if(nrow(surv_sub) > 1){
    pdf(paste0(out_folder, ID, '_', unique(sub_data$Study_Authors),
               '_', unique(sub_data$Species), 'Surv_probab.pdf'))
    plot(surv_sub$Trait_value, surv_sub$Demog_rate_mean,
         ylab = paste0(unique(surv_sub$Demog_rate)),
         xlab = paste0(unique(surv_sub$Trait)))
    points(surv_sub$Trait_value, exp(fitted(linRel)) / (1 + exp(fitted(linRel))), pch = 19, col = 'red')
    points(surv_sub$Trait_value, exp(fitted(sigmRel)) / (1 + exp(fitted(sigmRel))), pch = 21, col = 'grey')
    points(surv_sub$Trait_value, exp(fitted(quadRel)) / (1 + exp(fitted(quadRel))), pch = 19, col = 'blue')
    legend('bottomright', c('linear', 'sigmoid', 'quadratic'),
           col = c('red', 'grey', 'blue'), pch = c(19, 21, 19))
    dev.off()
    }
  }

  if(! (is(tt.error.linRel,"error") | is(tt.error.quadRel,"error") |
        is(tt.error.sigmRel,"error") )){
    AIC_lin <- MuMIn::AICc(linRel)
    AIC_quad <- MuMIn::AICc(quadRel)
    AIC_sigm <- MuMIn::AICc(sigmRel)
    if(! (is.infinite(AIC_lin) | is.infinite(AIC_quad) | is.infinite(AIC_sigm))){

      res <- data.frame('ID' = unique(sub_data$ID),
                        'AIC_Lin' = AIC_lin,
                        'AIC_Quad' = AIC_quad,
                        'AIC_Sigm' = AIC_sigm,
                        'Int_Lin' = coef(linRel)[names(coef(linRel)) == 'interc'],
                        'Beta_Lin' = coef(linRel)[names(coef(linRel)) == 'beta'],
                        'Int_Quad' = coef(quadRel)[names(coef(quadRel)) == 'interc'],
                        'Beta_Quad' = coef(quadRel)[names(coef(quadRel)) == 'beta'],
                        'Beta2_Quad' = coef(quadRel)[names(coef(quadRel)) == 'beta2'],
                        'Int_Sigm' = coef(sigmRel)[names(coef(sigmRel)) == 'interc'],
                        'Beta_Sigm' = coef(sigmRel)[names(coef(sigmRel)) == 'beta'])
      vect_AIC <- c(AIC_lin, AIC_quad, AIC_sigm)
      names(vect_AIC) <- c('linear', 'quadratic', 'sigmoid')
      ord <- order(c(AIC_lin, AIC_quad, AIC_sigm))
      min_AIC <- min(vect_AIC)
      Delta1 <- vect_AIC[ord][2] - min_AIC
      Delta2 <- vect_AIC[ord][3] - min_AIC


      res$minAIC <- min_AIC
      res$mod_minAIC <- names(which(vect_AIC == min_AIC))
      res$Delta1_AIC <- Delta1
      res$Delta2_AIC <- Delta2
      if(Delta1 > Thresh){res$Selected <- names(vect_AIC[ord][1])
      } else {
        if(Delta1 <= Thresh & Delta2 <= Thresh){
          # res$Selected <- 'linear'
          ## this still has to be adjusted accordingly
            if(length(intersect(names(vect_AIC[ord][c(1,2,3)]),
                                c('linear', 'sigmoid', 'quadratic'))) == 3){
              res$Selected <- 'linear/sigmoid'
            }
            # if(length(intersect(names(vect_AIC[ord][c(1,2)]),
            #                     c('linear', 'quadratic'))) == 2){
            #   res$Selected <- 'linear'
            # }
            # if(length(intersect(names(vect_AIC[ord][c(1,2)]),
            #                     c('sigmoid', 'quadratic'))) == 2){
            #   res$Selected <- 'sigmoid'
            # }
            #
                  } else {
          if(Delta1 <= Thresh & Delta2 > Thresh){
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('linear', 'sigmoid'))) == 2){
              res$Selected <- 'linear/sigmoid'
            }
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('linear', 'quadratic'))) == 2){
              res$Selected <- 'linear'
            }
            if(length(intersect(names(vect_AIC[ord][c(1,2)]),
                                c('sigmoid', 'quadratic'))) == 2){
              res$Selected <- 'sigmoid'
            }
          }
        }
      }
      return(res)
    } else {
      message(paste('infinite AIC for at least one of the models'))
    }
  } else {
    message(paste('error in convergence of at least one of the models'))
  }
}

## okay, so if AIC is infinite, it is infinite for all three relations, no point to fiddle further there

## check the function
fit_shape()
fit_shape(data = tot_scale,
          ID = 86, Thresh = 2)


## applying the function to all studies
shapes_fit <- do.call('rbind', lapply(unique(tot_scale$ID), FUN = function(x){fit_shape(data = tot_scale, ID = x, Thresh = 5)}))

table(shapes_fit$Selected)  # not a single time for sigmoid. Quadratic rather rare (11 out of 272)
##  linear/sigmoid      quadratic
##            260             7

table(shapes_fit$mod_minAIC)
## linear quadratic   sigmoid
## 116        25       126
## So, we have difficulties discriminating sigmoid from linear given shallow slopes
hist(shapes_fit$Beta_Sigm[shapes_fit$Selected == 'linear/sigmoid']) ## yes, so they mainly are very small values, therefore approximates linear
hist(shapes_fit$Int_Sigm[shapes_fit$Selected == 'linear/sigmoid'])
## well, these shallow slopes also mean that we simulated more extreme shapes and perhaps those
## are not common in nature??? + unclear whether under such more shallow slopes the findings
## on stability will hold...

hist(shapes_fit$Int_Lin[shapes_fit$Selected %in% c('linear/sigmoid', 'linear')])
hist(shapes_fit$Beta_Lin[shapes_fit$Selected %in% c('linear/sigmoid', 'linear')])  ## interesting, these slopes are more shallow than what we have tested for

hist(shapes_fit$Int_Quad[shapes_fit$Selected %in% c('quadratic')])
hist(shapes_fit$Beta_Quad[shapes_fit$Selected %in% c('quadratic')])
hist(shapes_fit$Beta2_Quad[shapes_fit$Selected %in% c('quadratic')]) ## these are also somewhat lower than the assumed abs value of 0.5 that we use in simus


## a plot for the MS with the estimated betas and intercepts
## sigm
pl_Sigm_b <- ggplot(subset(shapes_fit, Selected == 'linear/sigmoid'),
       aes(x = Beta_Sigm)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.y = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit$Beta_Sigm[shapes_fit$Selected
                                                        == 'linear/sigmoid'],
                                probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Slope')

quantile(shapes_fit$Beta_Sigm[shapes_fit$Selected == 'linear/sigmoid'],
           probs = c(0.05, 0.95))  ## -0.04087630  0.05204491

pl_Sigm_int <- ggplot(subset(shapes_fit, Selected == 'linear/sigmoid'),
       aes(x = Int_Sigm)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.y = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit$Int_Sigm[shapes_fit$Selected
                                                        == 'linear/sigmoid'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Intercept') + labs(title = 'Sigmoid')


## linear
pl_Lin_b <- ggplot(subset(shapes_fit, Selected %in% c('linear/sigmoid', 'linear')),
                    aes(x = Beta_Lin)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.x = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit$Beta_Lin[shapes_fit$Selected
                                                        %in% c('linear/sigmoid', 'linear')],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2)

quantile(shapes_fit$Beta_Lin[shapes_fit$Selected %in% c('linear/sigmoid', 'linear')],
         probs = c(0.05, 0.95))  ## -0.2043896  0.2638683

pl_Lin_int <- ggplot(subset(shapes_fit, Selected %in% c('linear/sigmoid', 'linear')),
                      aes(x = Int_Lin)) + geom_histogram() +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title.x = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit$Int_Lin[shapes_fit$Selected
                                                      %in% c('linear/sigmoid', 'linear')],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  labs(title = 'Linear')


## quadr
pl_Quad_b <- ggplot(subset(shapes_fit, Selected == 'quadratic'),
                   aes(x = Beta_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title = element_blank()) +
  geom_vline(xintercept = quantile(shapes_fit$Beta_Quad[shapes_fit$Selected
                                                       == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2)

quantile(shapes_fit$Beta_Quad[shapes_fit$Selected == 'quadratic'],
         probs = c(0.05, 0.95))  ## -0.3348143  0.1584980

pl_Quad_int <- ggplot(subset(shapes_fit, Selected == 'quadratic'),
                     aes(x = Int_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2)),
                     axis.title = element_blank(),
                     plot.title = element_text(hjust=0.5)) +
  geom_vline(xintercept = quantile(shapes_fit$Int_Quad[shapes_fit$Selected == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  labs(title = 'Quadratic')



pl_Quad_b2 <- ggplot(subset(shapes_fit, Selected == 'quadratic'),
                    aes(x = Beta2_Quad)) + geom_histogram(bins = 10) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.text = element_text(color = 'black'),
                     plot.margin = margin(c(0, 0.2, 0, 0.2))) +
  geom_vline(xintercept = quantile(shapes_fit$Beta2_Quad[shapes_fit$Selected
                                                        == 'quadratic'],
                                   probs = c(0.05, 0.95)), col = 'blue', lty = 2) +
  xlab('Quadratic slope')



## plot all plots together

des_lay <- "
123#
1237
4567
456#"

pdf('./output_nonL/EstimatedParS_DemTraitRel_sTraitChange_DeltaAIC5.pdf',
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
lin_sigm <- droplevels(subset(shapes_fit, Selected == 'linear/sigmoid'))

pdf('./output_nonL/SupplFig_DistrAIC_Sigmoid_Linear.pdf')
hist(lin_sigm$Delta1_AIC, main ='', xlab = expression(paste(Delta, 'AIC')),
     col = 'grey')
abline(v = mean(lin_sigm$Delta1_AIC), col = 'blue', lwd = 2)
abline(v = median(lin_sigm$Delta1_AIC), col = 'blue', lwd = 3, lty = 2)
dev.off()
hist(lin_sigm$Delta2_AIC[lin_sigm$mod_minAIC == 'quadratic'])  ## so, function works fine
table(lin_sigm$mod_minAIC)
## linear quadratic   sigmoid
##  116        18       126
## so, even if sometimes quadratic has lower AIC, overall the deltaAIC to either linear or sigmoid is low...

# III.  Quantify pop stability and link to shapes -------------------------

## calculate CV of popsize (and maybe also Growth rate?) I think these two explorations
## are enough (CV and the trend in popsize), no need to go for GR. Otherwise the reviewers
## might ask why we have not used GR as the measure of pop stability in the paper (which we
## did but do not report)
CV_pop <- tot_scale %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(CV = sd(Pop_mean, na.rm = TRUE) /
                  mean(Pop_mean, na.rm = TRUE)) %>%
  dplyr::distinct(ID, .keep_all = TRUE)


## fitting simple linear model to get the tredn in pop size over time
fit_lm <- function(df){lm(Pop_mean ~ Year, data = df)}

trendsPop <- tot_scale %>%
  dplyr::group_by(., ID) %>%
  dplyr::filter(., ! is.na(Pop_mean)) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    LM = purrr::map(purrr::map(data, fit_lm),
                    broom::tidy)) %>%
  tidyr::unnest(cols= c(LM), names_sep = '_') %>%
  dplyr::filter(., LM_term == 'Year') %>%
  dplyr::select(., -c(LM_term, LM_statistic)) %>%
  dplyr::mutate(., Trend = dplyr::case_when(LM_p.value >= 0.05 ~ 'stable',
                                            LM_p.value < 0.05 & LM_estimate > 0 ~ 'positive',
                                            LM_p.value < 0.05 & LM_estimate < 0 ~ 'negative'))

hist(trendsPop$LM_estimate, breaks = 20, col = 'lightgrey',
     xlab = 'Slope of population size over years', cex.lab = 1.4,
     main ='')
abline(v = median(trendsPop$LM_estimate), lwd = 2, col = 'blue')
table(trendsPop$Trend)

stab_dat <- merge(CV_pop, subset(trendsPop, select = -data),
                  by = 'ID')
stab_shape <- merge(stab_dat, shapes_fit, by = 'ID')
nrow(stab_shape)
nrow(shapes_fit)


## a simple linear model (hmmm, probably beta-reg would be better)
## ofr analyses remove the signle linear one

stab_shape_noL <- droplevels(subset(stab_shape, Selected != 'linear'))

## save this plot for the SI
pdf('./output_nonL/fig_SI_CV_vsDifferentShapes_sTraitChange.pdf')
ggplot(stab_shape_noL, aes(x = Selected, y = CV)) +
  geom_boxplot() + theme_bw() + ylab('CV in population size') +
  xlab('Supported shape of trait-demography relation') +
  theme(axis.text = element_text(color = 'black'),
        panel.grid.minor = element_blank())
dev.off()


lm_CV <- lm(CV ~ Selected, data = stab_shape_noL)
summary(lm_CV)  ## slope is not signif.
lm_int <- lm(CV ~ 1, data = stab_shape_noL)
anova(lm_int, lm_CV, test = 'F')  ## nothing

## have to remove the outliers otherwise nothing is really visible
stab_shape_filt <- stab_shape_noL[stab_shape_noL$LM_estimate < 2000 &
                                stab_shape_noL$LM_estimate > -1000, ]
ggplot(stab_shape_filt, aes(x = Selected, y = LM_estimate)) +
  geom_boxplot()

tab_Shape_trend <-table(stab_shape_filt$Selected, stab_shape_filt$Trend)
prop.table(tab_Shape_trend, margin = 1)  ## does not look like there are diff. BUT: numbers are too small
