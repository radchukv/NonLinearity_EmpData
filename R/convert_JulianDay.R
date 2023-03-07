#' Convert phenological traits that are not specified
#' in Julian days to Julian days
#'
#' \code{convert_JulianDay} performs standardization of phenological
#' traits by converting those that are specified from other
#' starting day than the 1st of January to Julian days
#'
#' @param biol_data Data frame with trait data for a given population and species.
#'
#' @export
#'
#' @return A data frame with columns as those in the input data
#' (biol_data), but with any phenological traits that were specified
#' not in Julian Days converted to Julian Days.
#'
#' @examples
#' biol_dat <- read.csv('./data-raw/test_dateChange.csv')
#' out_biol <- convert_JulianDay(biol_data = biol_dat)
#'
convert_JulianDay <- function(biol_data){
  if (any(levels(biol_data$Trait_Categ) %in% c('Phenological'))){
    subs_phen <- biol_data[biol_data$Trait_Categ == 'Phenological', ]
    if (any(! levels(subs_phen$Unit_trait) %in%
            c('JulianDay', 'Time'))){
      biol_sub <- droplevels(subs_phen[! subs_phen$Unit_trait %in%
                                         c('JulianDay', 'Time'), ])
      Month <- unlist(lapply(1:nrow(biol_sub), FUN = function(x){
        substr(strsplit(strsplit(as.character(biol_sub$Unit_trait[x]),
                                 'DaySince')[[1]][2], '1')[[1]], 1, 3)
      }))
      biol_sub$DayFirst <- as.Date(paste(biol_sub$Year, Month, '1', sep = '/'),
                                   format = '%Y/%b/%d')
      biol_sub$Trait_mean <- biol_sub$Trait_mean +
        as.numeric(biol_sub$DayFirst - as.Date(paste(biol_sub$Year, 'Jan', '1', sep = '/'),
                                               format = '%Y/%b/%d'))
      biol_sub$Unit_trait <- rep('JulianDay', nrow(biol_sub))
      biol_sub$DayFirst <- NULL

      full_dat <- rbind(droplevels(biol_data[biol_data$Unit_trait %in%
                                               c('JulianDay', 'Time'), ]),
                        droplevels(biol_data[biol_data$Trait_Categ != 'Phenological', ]),
                        biol_sub)
      return(full_dat)

    } else {
      return(biol_data)
    }
  } else {
    return(biol_data)
  }
}
