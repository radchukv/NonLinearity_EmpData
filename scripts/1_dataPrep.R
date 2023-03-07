## read in the temperature and biological study data

library(dplyr)
library(ggplot2)
library(magrittr)

# load functions
source('./R/convert_JulianDay.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        Weather world                        ####

temp_CPC <- list.files('./data-raw/weather/temp_world/', full.names = T)


## only read min_CPC, to assign the names of the layers for the mean stack
min_CPC <- grep('tmin', temp_CPC, value = TRUE)
minT <- raster::stack(lapply(1:length(min_CPC), FUN = function(x){
  raster::stack(min_CPC[x])}))


## read the prepared mean T stack
meanT_world <- raster::stack('./data-raw/weather/temp_world/meanT_world_CPC.nc', varname = 'variable') ## does not preserve the names of the layers
names(meanT_world) <- names(minT)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        World: SST                            ####
SST <- list.files('./data-raw/weather/SST', full.names = T)
daily_SST_single <- grep('day.mean', SST, value = TRUE)
all_SST <- raster::stack(lapply(1:length(daily_SST_single), FUN = function(x){
  raster::stack(daily_SST_single[x])}))


## read the prepared rotated SST stack
rot_SST <- raster::stack('./data-raw/weather/SST/SST_rotated.nc') ## check: this file is not in the folder
names(rot_SST) <- names(all_SST)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####                        BIOLOGICAL data                            ####

biol_dat <- read.csv('./data-raw/Data_1_02_2021.csv', header = T)
unique(biol_dat$Country)


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

## all characters as factors
biol_dat_fc <- biol_dat %>%
  mutate_if(is.character, as.factor)

## use the funciton to converting the dates that are not Julian to Julian dates
biol_dat_stand <- convert_JulianDay(biol_data = biol_dat_fc)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
####            Split into seabirds and others            ####

## splitting hte dataset into all who are not seabirds and the dataset with seabirds only
SeaBird <- droplevels(biol_dat_stand %>%
                                     dplyr::filter(., BirdType == 'Seabird'))
length(unique(SeaBird$ID)) # 68

others <- droplevels(biol_dat_stand %>%
                       dplyr::filter(., BirdType != 'Seabird'))

length(unique(others$ID))  ## 241


# quick check of the data in both datasets
max(SeaBird$Trait_mean[SeaBird$Trait_Categ == 'Phenological'], na.rm = T)
hist(SeaBird$Trait_mean[SeaBird$Trait_Categ == 'Phenological'])

unique(SeaBird$Unit_trait[SeaBird$Trait_Categ == 'Morphological'])
## so basically these would have to be looked at separatly for weight and for length....
max(SeaBird$Trait_mean[SeaBird$Trait_Categ == 'Morphological'], na.rm = T)
hist(SeaBird$Trait_mean[SeaBird$Trait_Categ == 'Morphological'])


max(others$Trait_mean[others$Trait_Categ == 'Phenological'], na.rm = T)
hist(others$Trait_mean[others$Trait_Categ == 'Phenological'])

unique(others$Unit_trait[others$Trait_Categ == 'Morphological']) ## here we have everything, kg, grams, cm, mm.....
# --> Standardize the trait prior to fitting the rleations between the yearly temperature and the trait

## so basically looking at histogrmas has to be done separatly for weight and for length + per unit....

