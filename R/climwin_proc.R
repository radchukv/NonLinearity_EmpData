#' Conduct climate window analysis
#'
#' \code{climwin_proc} performs climate window analysis for
#' specified traits and climate data and returns the best
#' identified window together with results of randomisation
#' tests (if asked for)
#'
#' @param biol_data Data frame with trait data for a given population and species.
#' @param clim_data Either raster stack of climatic data for the region that encompases
#' the study location or a dataframe with the local weather from the closest weather station.
#' @param ID Numeric giving a unique ID of the current dataset for a given population and species.
#' @param randwin Logical (TRUE/FALSE). Should \code{\link[climwin]{randwin}}
#'  be run together with \code{\link[climwin]{slidingwin}}?
#' @param repeats If randwin is run, the number of times that data will be randomised.
#' @param plot_check Logical (TRUE/FALSE). Whether to display
#' one year (the first in the series) of the climatic data
#'  together with the study location.
#' @param cinterval A character specifying the resolution at which climate window
#' analysis will be conducted, defaults to 'week' (weekly resolution). For more details
#' see the same option in the function \code{\link[climwin]{slidingwin}}.
#' @param out_clim Character specifying the library on the path where
#' to save the results of climate window analysis
#' @param stat Character specifying which statistics to use for
#'  aggregating the climatic data (see also the same option in
#'  \code{\link[climwin]{slidingwin}}).
#' @param oneGrid Logical (TRUE/FALSE). Whether to extact climatic data from
#' a single grid cell where into which the study location falls or to use mean of
#' five cells: the focal one and four neighbours.
#' @param explanYear Logical (TRUE/FALSE). Whether to include year as an
#' explanatory variable in the baseline formula for window analysis.
#' @param endWindow A numeric specifying the furtherst number of time intervals
#' (set by cinterval) back from the reference day. See also the option range in
#' the function \code{\link[climwin]{slidingwin}}.
#' @param startWindow A numeric specifying the closest number of time intervals
#' (set by cinterval) back from the reference day. See also the option range in
#' the function \code{\link[climwin]{slidingwin}}.
#' @param RefMon A numeric specifying the month for the absolute window
#' to be used for all species. If NA is specified then a species-specific window is used for
#' each dataset using the available data on phenology and the period when morphology
#' was recorded.
#' @param weatherVar Character corresponding to the name of the climatic variable
#' in the data frame provided as clim_data.
#'
#' @export
#'
#' @return A tibble with six columns if randwin is set to FALSE
#' and seven columns if randwin is set TRUE. The columns are: "ID"
#' - study id, "Species" - study species, "Trait" - the trait considered,
#' "climwin_output" - a list returned by slidingwin(),
#' "randwin_output" - a list returned by randwin (if set to TRUE),
#' "clim_data" - a data frame with climate values,
#' "biol_data" - a data frame with traits, demographic rates
#' and population size.
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' meanT <- raster::stack('./data-raw/tg_ens_mean_0.1deg_reg_v18.0e.nc')
#' test_rand <- climwin_proc(biol_data = dat_birds,
#'                           clim_data = meanT, ID = 1,
#'                           randwin = FALSE, seednum = 1302,
#'                           repeats = 30, plot_check = FALSE,
#'                           out_clim = 'output_climwin',
#'                           cinterval = 'month',
#'                           stat = 'mean',
#'                           startWindow = 0, endWindow = 12,
#'                           oneGrid = FALSE, explanYear = TRUE,
#'                           RefMon = NA, weatherVar = NA)
#'
climwin_proc <- function(biol_data, clim_data,
                         ID, randwin = FALSE,
                         seednum = 1302, repeats = 200,
                         plot_check = FALSE,
                         cinterval = 'week',
                         out_clim = 'output_climwin_test',
                         stat = 'mean', oneGrid = TRUE,
                         explanYear = TRUE,
                         startWindow = 0, endWindow = 52,
                         RefMon = NA, weatherVar = NA){

  biol_data <- droplevels(biol_data[biol_data$ID == ID, ])
  # add Date to biol data (for slidingwin)
  biol_data$Date <- as.Date(paste('01', '06', biol_data$Year,
                                  sep = '/'), format = '%d/%m/%Y')

  ## imputing mean of the previous and next values for the mean of the traits, and
  ## median of all the values for the SE of the traits
  ## !!!! checking with no imputation for now   !!!!
  # biol_data <- impute_ma(data = biol_data, column = 'Trait_mean')
  # biol_data <- impute_median(data = biol_data, column = 'Trait_SE')
  # biol_data <- biol_data[! is.na(biol_data$Trait_mean), ] -  former approach, excluding NAs

  location_spatial <-  sp::SpatialPointsDataFrame(coords = biol_data[1, c('Longitude', 'Latitude')],
                                                  data = data.frame(ID = biol_data$ID[1]),
                                                  proj4string = sp::CRS(as.character('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')))

  ##  for a visual check
  if (plot_check){
    if(class(clim_data) != 'data.frame'){
      raster::plot(clim_data[[1]])
      points(location_spatial)
    } else {
      raster::plot(location_spatial)
    }
  }

  ####                          Extract temperature data                  ####
  # Determine date data for each layer of the raster (allows us to sort by each year).
  if(class(clim_data) == 'data.frame'){
    ## if weather comes from a single station
    Clim <- data.frame(Date = seq(as.Date(paste('01', '01', min(biol_data$Year) - 2, sep = '/'), format = '%d/%m/%Y'),
                                  as.Date(paste('01', '12', max(biol_data$Year) + 1, sep = '/'), format = '%d/%m/%Y'), 'day'))
    if(unique(biol_data$Study_Authors) %in% c('Tarwater&Beissinger')){
      clim_data$Date <- as.Date(clim_data$DATE, format = '%d/%m/%Y')
    }else{
    clim_data$Date <- as.Date(clim_data$DATE, format = '%Y-%m-%d')
    }
    subs_clim <- clim_data[, c('Date', weatherVar)]
    Clim <- merge(subs_clim, Clim, by = 'Date')
    Clim$Temp <- Clim[[weatherVar]]

  }else{ ## if weather is available as rasters or stacks
    temp_dates <- data.frame(Date = as.Date(substr(names(clim_data),
                                                   start = 1, stop = 11),
                                            format = 'X%Y.%m.%d'))
    temp_dates$Year <- lubridate::year(temp_dates$Date)

    ## extract data from Euro weather for all the necessary dates for the site
    ## (i.e. everything from the year before the first recorded nest until the year of the most recent brood).
    message(paste('Currently extracting temperature data for study',
                  biol_data$ID[1], 'for species',
                  biol_data$Species[1], 'in', biol_data$Location[1], 'for',
                  biol_data$Trait[1]))
    Clim <- data.frame(Date = seq(as.Date(paste('01', '01', min(biol_data$Year) - 2, sep = '/'), format = '%d/%m/%Y'),
                                  as.Date(paste('01', '12', max(biol_data$Year) + 1, sep = '/'), format = '%d/%m/%Y'), 'day'),  # why max should be this and not + 1 year?
                       Temp = NA)  ##longer end date for clim dtaa compared to biol data is needed in order for basewin()
    ## checks to not crush

    ptstart <- proc.time()

    ## using a single grid for Temp extraction
    if (oneGrid){
      Clim$Temp <- ifelse(is.na(Clim$Temp),
                          as.numeric(raster::extract(clim_data[[which(temp_dates$Date %in%
                                                                        Clim$Date)]], location_spatial)),
                          NA)
    }else{
      ## using multiple grids for Temp extraction

      ## marking the focal and the adjecent cells
      cellNum <- raster::cellFromXY(clim_data, location_spatial)
      FiveCell <- raster::adjacent(clim_data, cellNum, include = TRUE, pairs = FALSE)

      Clim$Temp <- ifelse(is.na(Clim$Temp),
                          colMeans(raster::extract(clim_data[[which(temp_dates$Date %in%
                                                                      Clim$Date)]], FiveCell), na.rm = T),
                          NA)
    }

    ptfinish <- proc.time() - ptstart
    cat('When extracting weather data from the raster the time elapsed is ',
        ptfinish[3], ';\n', 'time spend for ',
        names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
        'for ', names(ptfinish)[2], ' is ', ptfinish[2], '\n',
        sep = '')

  }


  # Check that data extracted properly!!
  if(all(is.na(Clim$Temp))){
    stop(paste0('Climate data failed to extract for ', biol_data$ID[1]))
  }

  # NOW THAT WE HAVE CLIMATE AND BIOLOGICAL DATA WE CAN RUN CLIMWIN!!
  set.seed(seednum)

  # biol_data$W <- 1 / biol_data$Trait_SE^2
  # get the refday for each species - the latest observed date (across years)
  # for phenological traits or the latest record made for morphological traits
  if (! is.na(RefMon)){
    refDay <-  as.Date(paste('01', RefMon, lubridate::year(Sys.Date()), sep = '/'),
                       format = '%d/%m/%Y')
    print(refDay)
  }else {
    if(unique(biol_data$Trait_Categ) == 'Phenological'){
      refDay <- max(biol_data$Trait_mean, na.rm = T)
      refDay <- as.Date(refDay, origin = as.Date(paste('01', '01', lubridate::year(Sys.Date()), sep = '/'),
                                                 format = '%d/%m/%Y'))
    }
    if(unique(biol_data$Trait_Categ) == 'Morphological'){
      if(length(grep('-', unique(biol_data$Record_date), value = T)) != 0){
        RecMonth <- strsplit( grep('-', unique(biol_data$Record_date), value = T), split = '-')[[1]][2]
        RecMonth <- substr(RecMonth, 1, 3)
      }else{
        RecMonth <- as.character(unique(biol_data$Record_date))
        RecMonth <- substr(RecMonth, 1, 3)
      }
      if (! is.na(RecMonth)){
      if(RecMonth == 'Feb'){
        refDay <- as.Date(paste('28', RecMonth, lubridate::year(Sys.Date()), sep = '/'),
                          format = '%d/%B/%Y')
      }else{
        refDay <- as.Date(paste('30', RecMonth, lubridate::year(Sys.Date()), sep = '/'),
                          format = '%d/%B/%Y')
      }
      } else {
        refDay <- NA
        }
      print(refDay)
    }
  }

  ## replacing the SEs for those studies where they are fully missing
  if (sum(is.na(biol_data$Trait_SE)) == nrow(biol_data)){
    biol_data$Trait_SE <- 1}

  ## have to exclude the rows with missing biological data before running the slidingwin
  biol_data_noNA <- biol_data %>%
    dplyr::filter(!is.na(Trait_mean) & !is.na(Trait_SE)) %>%
    # create weights here, otherwise slidingwin does not see them
    dplyr::mutate(W = 1/ Trait_SE^2)

  ## prepare the baseline formula
  if(explanYear){
    formBase <<- 'Trait_mean ~ Year'
  }else{
    formBase <<- 'Trait_mean ~ 1'
  }

  ptstart <- proc.time()
  climwin_output <- climwin::slidingwin(xvar = list(Temp = Clim$Temp),
                                        cdate = Clim$Date,
                                        bdate = biol_data_noNA$Date,
                                        baseline = stats::lm(stats::as.formula(formBase),
                                                             data = biol_data_noNA,
                                                             weights = W),
                                        range = c(endWindow, startWindow),
                                        cinterval = cinterval,
                                        stat = stat, func = 'lin',
                                        type = 'absolute',
                                        refday = c(lubridate::day(refDay),
                                                   lubridate::month(refDay)),
                                        cmissing = 'method1')
  ptfinish <- proc.time() - ptstart
  cat('When fitting slidingwin() the time elapsed is ',
      ptfinish[3], ';\n', 'time spend for ',
      names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
      'for ', names(ptfinish)[2], ' is ', ptfinish[2],
      '\n', sep = '')


  if(randwin){

    ptstart <- proc.time()
    randwin_output <- climwin::randwin(repeats = repeats,
                                       xvar = list(Temp = Clim$Temp),
                                       cdate = Clim$Date,
                                       bdate = biol_data_noNA$Date,
                                       baseline = lm(stats::as.formula(formBase),
                                                     data = biol_data_noNA,
                                                     weights = W),
                                       range = c(endWindow, startWindow),
                                       cinterval = cinterval,
                                       stat = stat, func = 'lin',
                                       type = 'absolute',
                                       refday = c(lubridate::day(refDay),
                                                  lubridate::month(refDay)),
                                       cmissing = 'method1')
    ptfinish <- proc.time() - ptstart
    cat('When running randwin() with ', repeats,
        'number of repeats, the time elapsed is ',
        ptfinish[3], ';\n', 'time spend for ',
        names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
        'for ', names(ptfinish)[2], ' is ', ptfinish[2],
        '\n', sep = '')

    # create a tibble to save all output together
    clim_out <- tibble::tibble(ID = biol_data$ID[1],
                               Species = biol_data$Species[1],
                               Trait = biol_data$Trait[1],
                               climwin_output = list(climwin_output[[1]]),
                               randwin_output = list(randwin_output[[1]]),
                               refDay = refDay,
                               clim_data = list(Clim),
                               biol_data = list(biol_data))
    saveRDS(object = clim_out,
            file = paste0('./', out_clim, '/', biol_data$ID[1], '_',
                          biol_data$Species[1], '_', biol_data$Location[1],
                          '_', biol_data$Trait[1], '_OneGrid_', oneGrid,
                          '_explYear_', explanYear, '_EndWindow_',
                          endWindow, '_RefMon_', RefMon, '_Rand',  '.RDS'))

  } else {
    # create a tibble to save all output together
    clim_out <- tibble::tibble(ID = biol_data$ID[1],
                               Species = biol_data$Species[1],
                               Trait = biol_data$Trait[1],
                               climwin_output = list(climwin_output[[1]]),
                               refDay = refDay,
                               clim_data = list(Clim),
                               biol_data = list(biol_data))
    saveRDS(object = clim_out,
            file = paste0('./', out_clim, '/', biol_data$ID[1], '_',
                          biol_data$Species[1], '_',
                          biol_data$Location[1], '_', biol_data$Trait[1],
                          '_OneGrid_', oneGrid, '_explYear_', explanYear,
                          '_EndWindow_', endWindow, '_RefMon_', RefMon,
                          '.RDS'))
  }
  return(clim_out)
}
