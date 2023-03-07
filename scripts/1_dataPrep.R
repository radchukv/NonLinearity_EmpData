#### Climwin analyses of sTraitChange data
#### Authors: V Radchuk & C Jones

devtools::install_github('radchukv/sTraitChange', auth_token = '4a4aba816767714248758ae04b068ff7d2d92097')
devtools::install_github("LiamDBailey/climwin", ref = "devel")
library(sTraitChange)
library(dplyr)
library(ggplot2)
library(magrittr)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather EUROPE                        ####
meanT <- raster::stack('./data-raw/tg_ens_mean_0.1deg_reg_v18.0e.nc', varname = 'tg')
meanP <- raster::stack('./data-raw/rr_ens_mean_0.1deg_reg_v18.0e.nc', varname = 'rr')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather USA                           ####
# precip
precip_CPC_US <- list.files('./data-raw/weather/precip_USA/', full.names = T)
totP_US <- raster::stack(lapply(1:length(precip_CPC_US), FUN = function(x){
  raster::stack(precip_CPC_US[x])}))
## need to bring all the coordinate systems to the common denominator...
## climatic model use longitude coordinates that run from 0 to 360, rather than -180 to 180
## see here: http://r-sig-geo.2731867.n2.nabble.com/Raster-projection-alignment-issues-td7587564.html
## needs rotation
totP_US <- raster::rotate(totP_US)

temp_CPC <- list.files('./data-raw/weather/temp_world/', full.names = T)


## only read min_CPC, to assign the names of the layers for the mean stack
min_CPC <- grep('tmin', temp_CPC, value = TRUE)
minT <- raster::stack(lapply(1:length(min_CPC), FUN = function(x){
  raster::stack(min_CPC[x])}))


## read the prepared mean T stack
meanT_world <- raster::stack('./data-raw/weather/temp_world/meanT_world_CPC.nc', varname = 'variable') ## does not preserve the names of the layers
names(meanT_world) <- names(minT)
# pdf('./output_climwin_test/T_world.pdf')
# raster::plot(meanT_US[[1]])
# dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        World: SST                            ####
SST <- list.files('./data-raw/weather/SST', full.names = T)
daily_SST_single <- grep('day.mean', SST, value = TRUE)
all_SST <- raster::stack(lapply(1:length(daily_SST_single), FUN = function(x){
  raster::stack(daily_SST_single[x])}))


## read the prepared rotated SST stack
rot_SST <- raster::stack('./data-raw/weather/SST/SST_rotated.nc') ## does not preserve the names of the layers
names(rot_SST) <- names(all_SST)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        World: Precipit                            ####
precip_W <- list.files('./data-raw/weather/prec_world', full.names = T)
Precip_W_single <- grep('precip.', precip_W, value = T)
world_P <- raster::stack(lapply(1:length(Precip_W_single), FUN = function(x){
  raster::stack(Precip_W_single[x])}))

## read the prepared rotated world precip stack
rot_prec <- raster::stack('./data-raw/weather/prec_world/Prec_rotated.nc') ## does not preserve the names of the layers
names(rot_prec) <- names(world_P)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather Australia                           ####
# precip
precip_AU <- list.files('./data-raw/weather/Australia/Precipit_AU/', full.names = T)
totP_AU <- raster::stack(lapply(1:length(precip_AU), FUN = function(x){
  raster::stack(precip_AU[x])}))
totP_AU

## temp
temp_AU <- list.files('./data-raw/weather/Australia/Temp_AU/', full.names = T)

## minT - for names of the layers
min_T_AU <- grep('tmin', temp_AU, value = TRUE)
minT_AU <- raster::stack(lapply(1:length(min_T_AU), FUN = function(x){
  raster::stack(min_T_AU[x])}))

## read the prepared mean T
meanT_AU <- raster::stack('./data-raw/weather/Australia/Temp_AU/meanT_AU.nc') ## does not preserve the names of the layers
names(meanT_AU) <- names(minT_AU)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        BIOLOGICAL data                            ####

biol_dat <- read.csv('./data-raw/Data_1_02_2021.csv', header = T)
levels(biol_dat$Country)


## study by Moeller with non-comparable phenological measure, exclude
biol_dat <- droplevels(subset(biol_dat,  ID != 183))
## alaskan study with many NAs, exclude
biol_dat <- droplevels(subset(biol_dat, ! ID %in% c(228, 231, 232, 234)))


## still take care of dem_rate_Se == 0
unique(subset(biol_dat, Demog_rate_SE == 0)$Study_Authors) ## drop them??? or replace by NA?


biol_dat$Unit_trait[biol_dat$Unit_trait == 'KG'] <- 'kg'
biol_dat <- droplevels(biol_dat)

length(unique(biol_dat$ID))    ## 309

## some Demog_rate_Categ have white spaces... Correcting
biol_dat$Demog_rate_Categ <- trimws(biol_dat$Demog_rate_Categ)



## keep only European studies and US - for now (two separate datasets, because climwin are diff.)
biol_eu <- droplevels(subset(biol_dat, ! Country %in% c('Antarctica', 'Australia',
                                                        'Canada', 'Falkland Islands',
                                                        'Greenland', 'Mexico',
                                                        'New Zealand', 'South Africa',
                                                        'South Atlantic Ocean',
                                                        'South Georgia', 'Svalbard',
                                                        'Taiwan', 'USA', 'Venezuela')))

length(unique(biol_eu$ID))        ## 182
length(unique(biol_eu$Species))   ## 30
length(unique(biol_eu$Location))  ## 40
unique(biol_eu$Taxon)  # 3


# standardize units
levels(biol_eu$Unit_trait)

## use the funciton to converting the dates that are not Julian to Julian dates
biol_eu_stand <- convert_JulianDay(biol_data = biol_eu)

## US
biol_us <- droplevels(subset(biol_dat, Country %in% c('USA')))

length(unique(biol_us$ID))        ## 53
length(unique(biol_us$Species))   ## 10
length(unique(biol_us$Location))  ## 15
unique(biol_us$Taxon)             ## 4
table(biol_us$Taxon)   ## reptiles lead
table(biol_us$Trait_Categ)  ## mostly morphol

# standardize units
levels(biol_us$Unit_trait)  ## so convertion of dates to Julian not needed

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                               prepare data for analyses                                      ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. EUrope no sea birds
eu_noSea <- prep_subset(data = biol_eu_stand, Seabird = FALSE)
nrow(eu_noSea$Sel[[1]])   ## 121
length(unique(eu_noSea$subdata[[1]]$ID))   ## 140

max(eu_noSea$subdata[[1]]$Trait_mean[eu_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(eu_noSea$subdata[[1]]$Trait_mean[eu_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])


# 2. US no sea birds
## exclude sea birds
us_noSea <- prep_subset(data = biol_us, Seabird = FALSE)
nrow(us_noSea$Sel[[1]])  ## 47
length(unique(us_noSea$subdata[[1]]$ID))   ## 50

max(us_noSea$subdata[[1]]$Trait_mean[us_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(us_noSea$subdata[[1]]$Trait_mean[us_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])


## 3. Sea birds
## keep only sea birds
all_Sea <- prep_subset(data = biol_dat, Seabird = TRUE)
nrow(all_Sea$Sel[[1]])  ## 54
length(unique(all_Sea$subdata[[1]]$ID))    ## 68

max(all_Sea$subdata[[1]]$Trait_mean[all_Sea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(all_Sea$subdata[[1]]$Trait_mean[all_Sea$subdata[[1]]$Trait_Categ == 'Phenological'])

## 4. Other studies (no sea birds) for which we use world data - pay attention, this list has to be updated?
biol_w <- droplevels(subset(biol_dat, ID %in%
                              c(unique(subset(biol_dat, Country == 'Canada')$ID),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Wheelwright_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Cheng_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Nater_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Veiberg_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Jenouvrier_et_al']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Tarwater&Beissinger']),
                                unique(biol_dat$ID[biol_dat$Study_Authors == 'Veran&Beissinger']))))
length(unique(biol_w$ID))        ## 46
length(unique(biol_w$Species))   ## 22
length(unique(biol_w$Location))  ## 15
unique(biol_w$Taxon)             ## 3

# no need to standardize units
levels(biol_w$Unit_trait)

rest_w <- prep_subset(data = biol_w, Seabird = FALSE)
nrow(rest_w$Sel[[1]])                     ## 33
length(unique(rest_w$subdata[[1]]$ID))    ## 38

max(rest_w$subdata[[1]]$Trait_mean[rest_w$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(rest_w$subdata[[1]]$Trait_mean[rest_w$subdata[[1]]$Trait_Categ == 'Phenological'])


## 5. Australia
biol_au <- droplevels(subset(biol_dat, Country %in% c('Australia')))
biol_au_noSea <- prep_subset(biol_au, Seabird = FALSE)
nrow(biol_au_noSea$Sel[[1]])                      ## 11
length(unique(biol_au_noSea$subdata[[1]]$ID))     ## 13

max(biol_au_noSea$subdata[[1]]$Trait_mean[biol_au_noSea$subdata[[1]]$Trait_Categ == 'Phenological'], na.rm = T)
hist(biol_au_noSea$subdata[[1]]$Trait_mean[biol_au_noSea$subdata[[1]]$Trait_Categ == 'Phenological'])

