#' ZooMSS Functional Groups Data
#'
#' @title Default functional groups for the ZooMSS model
#' @description A dataset containing the biological parameters for different
#'   functional groups used in the ZooMSS size-structured marine ecosystem model.
#'   These represent various taxa from flagellates to large fish, each defined
#'   by their feeding behavior, size ranges, and physiological parameters.
#'
#' @format A data frame with 12 rows (functional groups) and 19 columns:
#' \describe{
#'   \item{Species}{Character. Name of the functional group/taxa}
#'   \item{Type}{Character. Broad category (Zooplankton or Fish)}
#'   \item{FeedType}{Character. Feeding strategy (Heterotroph, FilterFeeder, Omnivore, Carnivore)}
#'   \item{Prop}{Numeric. Initial proportion of total biomass}
#'   \item{W0}{Numeric. Log10 minimum body weight (g) for the group}
#'   \item{Wmax}{Numeric. Log10 maximum body weight (g) for the group}
#'   \item{Wmat}{Numeric. Log10 maturation body weight (g)}
#'   \item{SearchCoef}{Numeric. Search coefficient for predation interactions}
#'   \item{SearchExp}{Numeric. Search exponent for predation scaling}
#'   \item{PPMRscale}{Numeric. Predator-prey mass ratio scaling parameter}
#'   \item{PPMR}{Numeric. Predator-prey mass ratio (for fish groups)}
#'   \item{FeedWidth}{Numeric. Feeding kernel width parameter}
#'   \item{GrossGEscale}{Numeric. Gross growth efficiency scaling}
#'   \item{Carbon}{Numeric. Carbon content proportion}
#'   \item{Repro}{Numeric. Reproduction parameter}
#'   \item{Fmort}{Numeric. Fishing mortality rate}
#'   \item{Fmort_W0}{Numeric. Log10 minimum weight for fishing mortality}
#'   \item{Fmort_Wmax}{Numeric. Log10 maximum weight for fishing mortality}
#'   \item{PlotColour}{Character. Color code for plotting the functional group}
#' }
#'
#' @details The TestGroups dataset defines 12 functional groups spanning from
#'   small microzooplankton (flagellates, ciliates) through various mesozooplankton
#'   groups (copepods, euphausiids, chaetognaths) to gelatinous zooplankton (salps, jellyfish)
#'   and three fish size classes. Each group is characterized by:
#'   
#'   - **Size ranges**: W0 to Wmax define the body size spectrum
#'   - **Feeding behavior**: Different strategies for resource acquisition
#'   - **Interaction parameters**: Search rates and predator-prey relationships
#'   - **Physiological rates**: Growth efficiency and carbon content
#'   
#'   These parameters are based on marine ecological literature and represent
#'   typical values for temperate marine ecosystems.
#'
#' @source Marine ecological literature and ZooMSS model development
#' @family ZooMSS-data
#' @examples
#' data(TestGroups)
#' head(TestGroups)
#' 
#' # View size ranges across groups
#' plot(TestGroups$W0, TestGroups$Wmax, 
#'      col = TestGroups$PlotColour,
#'      xlab = "Log10 Min Weight", ylab = "Log10 Max Weight")
#' text(TestGroups$W0, TestGroups$Wmax, TestGroups$Species, pos = 3, cex = 0.7)
"TestGroups"
