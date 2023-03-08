## calculating mean yearly temperatures for each study location

library(dplyr)
library(ggplot2)
library(magrittr)

## prepare the data
source('./scripts/1_dataPrep.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         mean yearly temperature                     ####

## checking if the function works
test <- climwin_proc(biol_data = others,
             clim_data = meanT_world,
             ID = 1, plot_check = TRUE,
             out_clim = './output/output_temp')

# All non-seabird studies
IDs_others <- unique(others$ID)
Oth_yearly <- lapply(IDs_others, FUN = function(x){
  temp_proc(biol_data = others,
               clim_data = meanT_world,
               ID = x, plot_check = TRUE,
               out_clim = './output/output_temp')})


# seabird studies
IDs_seab <- unique(SeaBird$ID)
SeaB_yearly <- lapply(IDs_seab, FUN = function(x){
  temp_proc(biol_data = SeaBird,
            clim_data = rot_SST,
            ID = x, plot_check = TRUE,
            out_clim = './output/output_temp_SeaB')})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`#`

## the rest can be trashed
## I actually do not have to do these next steps if I simply use world temperature
# with 0.5 by 0.5 resolution for all the studies in the dataset
## keep only European studies and US - for now (two separate datasets, because climwin are diff.)
biol_eu <- droplevels(subset(biol_dat_fc, ! Country %in% c('Antarctica', 'Australia',
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

