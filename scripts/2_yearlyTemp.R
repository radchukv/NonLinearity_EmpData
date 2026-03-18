## calculating mean yearly temperatures for each study location

library(dplyr)
library(ggplot2)
library(magrittr)

## prepare the data
source('./scripts/1_dataPrep.R')
source('./R/temp_proc.R')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####         mean yearly temperature                     ####

## checking if the function works
test <- temp_proc(biol_data = others,
             clim_data = meanT_world,
             ID = 1, plot_check = TRUE,
             out_clim = './output/output_temp',
             oneGrid = FALSE)

# All non-seabird studies
IDs_others <- unique(others$ID)
Oth_yearly <- lapply(IDs_others, FUN = function(x){
  temp_proc(biol_data = others,
            clim_data = meanT_world,
            ID = x, plot_check = TRUE,
            out_clim = './output/output_temp',
            oneGrid = FALSE)})


# seabird studies
IDs_seab <- unique(SeaBird$ID)
# check Seabird
test_sea <- temp_proc(biol_data = SeaBird,
                  clim_data = rot_SST,
                  ID = 196, plot_check = TRUE,  ## suspiciously for many studies fail to extract temp, check
                  out_clim = './output/output_temp',
                  oneGrid = FALSE) # check if the problem is with trying to extract the data from the cells that are too
# close to the mainland

# even with the neighbors, for 438 extracting the temperature fails, remove it
IDs_seab <- IDs_seab[-which(IDs_seab %in% c(438))]

# temp can't be extracted
# the pb was that some locarions are too close to the mainland and then for them no proper SST can be extracted, os
# using four neighbours - in that case, I will probably have to rerun the same also for the air temperature...
SeaB_yearly <- lapply(IDs_seab, FUN = function(x){
  temp_proc(biol_data = SeaBird,
            clim_data = rot_SST,
            ID = x, plot_check = TRUE,
            oneGrid = FALSE,
            out_clim = './output/output_temp_SeaB')})

