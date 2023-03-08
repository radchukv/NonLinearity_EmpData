## detrending the temperature and trait values per study ID (still unsure what to do with the dem. rate)

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(broom)

## source functions
source('./R/fit_shape.R')

## read in the studies
fil <- list.files(path = './output/output_temp/', full.names = TRUE)

all<- bind_rows(lapply(fil, FUN = function(x){readRDS(file = x)}))


## detrending temperature oevr years
data_Tempdetr <- all %>%
  group_by(., ID) %>%
  nest() %>%
  mutate(fitm = map(data, ~lm(mean_Temp ~ Year,
                                            data = .)),
                res = map(fitm, augment)) %>%
  tidyr::unnest(cols = c(res)) %>%
  dplyr::rename(Temp_resid = .std.resid)
nrow(data_Tempdetr) ## 2853
nrow(all) ## 2865   -> hmm some data rows are lost...???

length(unique(data_Tempdetr$ID))


all_detrT <- merge(all, subset(data_Tempdetr, select = c(ID, Year, Temp_resid)), by = c('ID', 'Year'))
head(all_detrT)
hist(all_detrT$mean_Temp)
hist(all_detrT$Temp_resid)


## detrending trait over years
data_TTraitdetr <- all_detrT %>%
  group_by(., ID) %>%
  nest() %>%
  mutate(fitm = map(data, ~lm(Trait_mean ~ Year,
                              data = .)),
         res = map(fitm, augment)) %>%
  tidyr::unnest(cols = c(res)) %>%
  dplyr::rename(Trait_resid = .std.resid)
nrow(data_TTraitdetr) ## 2809
nrow(all_detrT) ## 2853   -> so again some rows are excluded

length(unique(data_TTraitdetr$ID))


all_detr <- merge(all_detrT, subset(data_TTraitdetr, select = c(ID, Year, Trait_resid)), by = c('ID', 'Year'))
head(all_detr)
hist(all_detr$Trait_mean)
hist(all_detr$Trait_resid)


## also add scaled values of trait and temperature per study (i.e. z-score)
all_trans <- all_detr %>%
  group_by(ID) %>%
  mutate(Temp_z = scale(mean_Temp),
         Trait_z = scale(Trait_mean))

hist(all_trans$Trait_z)
hist(all_trans$Temp_z)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                plots of residuals, raw and z-scores    ####

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


## plot of raw traits vs residual temp
TempResid_traitRaw <- ggplot(all_trans, aes(x = Temp_resid, y = Trait_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, residuals') +
  ylab('Trait, raw')
print(TempResid_traitRaw)

ggforce::n_pages(TempResid_traitRaw)
pdf('./output/output_nonL/TraitRaw_vs_TemperatureResid.pdf')
for(i in 1:8){
  TempResid_traitRaw <- ggplot(all_trans, aes(x = Temp_resid, y = Trait_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')  + xlab('Temperature, residuals') +
    ylab('Trait, raw')
  print(TempResid_traitRaw)
}
dev.off()
### Hmm, with residuals of temperature, the picture is not the same as with the raw temperature... Now I am not sure anymore whether it is
## more meaningful to just scale the variabels or indeed to detrend them


## plot of raw traits vs std temp
TempSTD_traitRaw <- ggplot(all_trans, aes(x = Temp_z, y = Trait_mean)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, z-score') +
  ylab('Trait, raw')
print(TempSTD_traitRaw)

ggforce::n_pages(TempSTD_traitRaw)
pdf('./output/output_nonL/TraitRaw_vs_TemperatureSTD.pdf')
for(i in 1:8){
  TempSTD_traitRaw <- ggplot(all_trans, aes(x = Temp_z, y = Trait_mean)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')  + xlab('Temperature, z-score') +
    ylab('Trait, raw')
  print(TempSTD_traitRaw)
}
dev.off()


## plot of residual traits vs residual temp
TempResid_traitResid <- ggplot(all_trans, aes(x = Temp_resid, y = Trait_resid)) +
  geom_point() +
  ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                               scales = 'free', nrow = 5, page = 1) +
  geom_smooth(method = 'lm', se = FALSE, col = 'black') +
  geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
  geom_text(size = 1, data = abbrev_d,
            aes(y = Inf, x = -Inf, label = label, vjust = 1,
                hjust = 0, colour  = 3)) +
  theme(legend.position = 'none') + xlab('Temperature, residuals') +
  ylab('Trait, residulas')
print(TempResid_traitResid)

ggforce::n_pages(TempResid_traitResid)
pdf('./output/output_nonL/TraitResid_vs_TemperatureResid.pdf')
for(i in 1:8){
  TempResid_traitResid <- ggplot(all_trans, aes(x = Temp_resid, y = Trait_resid)) +
    geom_point() +
    ggforce::facet_wrap_paginate(~ ID, ncol = 4,
                                 scales = 'free', nrow = 5, page = i) +
    geom_smooth(method = 'lm', se = FALSE, col = 'black') +
    geom_smooth(col = 'chocolate', fill = 'darkorange2', linetype = 'dashed') +
    geom_text(size = 1.5, data = abbrev_d,
              aes(y = Inf, x = -Inf, label = label, vjust = 1,
                  hjust = 0, colour  = 3)) +
    theme(legend.position = 'none')  + xlab('Temperature, residuals') +
    ylab('Trait, residuals')
  print(TempResid_traitResid)
}
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####              compare diff shapes               ####

## test fit_shape for trait-temperature relations
check <- fit_shape(data = all_trans,
                   x = 'Temp_z',
                   y = 'Trait_z',
                   ID = 1,
                   Thresh = 2,
                   out_folder = './output/output_nonL/shapes/')

# run function with detrended temperaturr eand standardized trait across studies
## trying with a rather low threshold of 4
shapes_fit <- do.call('rbind', lapply(unique(all_trans$ID), FUN = function(x){fit_shape(data = all_trans, ID = x,
                                                                                        Thresh = 4,  ## a rather low thresh
                                                                                        x = 'Temp_resid',
                                                                                        y = 'Trait_z',
                                                                                        out_folder = './output/output_nonL/shapes/')}))
