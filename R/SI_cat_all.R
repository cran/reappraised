#' Example of 50 variables from different studies for categorical (cat_all_fn) analysis
#'
#' Sample from Sato/Iwamoto dataset of 31 studies with categorical data
#'
#' @format A data frame with 106 rows and 14 variables:
#' \describe{
#'   \item{study}{study ID}
#'   \item{var}{variable}
#'   \item{level}{level of variable}
#'   \item{recode}{value to recode level two if collapsing variable to two levels}
#'   \item{levels_no}{total number of levels for each variable}
#'   \item{group}{number of trial arms}
#'   \item{n1, n2, n3}{number of participants with characteristic in each group}
#'   \item{N1, N2, N3}{number of participants in each group}
#'   \item{p}{reported p-value}
#'   \item{stat}{reported statistical test used to calculate p-value}
#' }
"SI_cat_all"
