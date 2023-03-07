#' Prepare a subset of data focusing on a particular group of species
#'
#' \code{prep_subset} produces a subset of the supplied data
#' according to the requirements, i.e. whether the subset should be
#' for sea birds or not
#'
#' @param data Data frame to use for subsetting.
#' @param Seabird Logical (TRUE/FALSE) specifying whether to include
#' sea birds in the subset or not.
#'
#' @export
#'
#' @return Returns a tibble with three columns: "Sel" - a data frame where
#' each row represents a unique record per study ID for which climwin analyses
#' have to be conducted; "biol_NY" - a data frame with all the information as in
#' original study (longitudinal data per each study ID) and, additionally, NYears
#' column specifying the number of years per study; "subdata" -
#' a subset of biological data to be used in climwin analyses.
#'
#' @examples
#' biol_dat <- read.csv('./data-raw/Data_11_06_2019.csv')
#' biol_noSea <- prep_subset(data = biol_dat, Seabird = FALSE)

prep_subset <- function(data, Seabird = FALSE){

  ## choice on seabirds
  if(Seabird){
    sub_bird <- droplevels(data %>%
                             dplyr::filter(., BirdType == 'Seabird'))
  }else{
    sub_bird <- droplevels(data %>%
                             dplyr::filter(., BirdType != 'Seabird'))
  }

  ## have to check whether for each ID the trait is not being repeated
  biol_NY <- sub_bird %>%
    dplyr::group_by(., ID) %>%
    dplyr::mutate(., NYears = dplyr::n()) %>%
    dplyr::ungroup()

  nodupl <- biol_NY[! duplicated(biol_NY$ID), ]
  Sel <- droplevels(nodupl %>%
                      dplyr::distinct(., Study_Authors, Species, Location,
                                      Trait, NYears, .keep_all = T))

  return(tibble::tibble(Sel = list(Sel),
                        biol_NY = list(biol_NY),
                        subdata = list(sub_bird)))
}
