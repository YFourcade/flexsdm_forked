# data

#' A data set containing presences and absences of tree virtual species
#'
#' @format A tibble with 1150 rows and 3 variables:
#' \describe{
#'   \item{species}{virtual species names}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data('spp')
#' spp
#' }
"spp"

#' A data set containing environmental conditions of background points
#'
#' @format A tibble object with 5000 rows and 10 variables:
#' \describe{
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   \item{from column aet to landform}{columns with environmental variables}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data('backg')
#' spp
#' }
"backg"

#' A data set containing environmental condition of an Abies species
#'
#' @format A tibble object with 5000 rows and 10 variables:
#' \describe{
#'   \item{pr_ab}{presence and absences denoted by 1 and 0 respectively}
#'   \item{from column aet to landform}{columns with environmental variables}
#'   ...
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data('abies')
#' spp
#' }
"abies"