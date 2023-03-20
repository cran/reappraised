#' Example of 50 observations for p-value analysis
#'
#' Sample from Sato/Iwamoto dataset of 500 baseline variables in 41 trials\cr
#' for details see Bolland MJ, Avenell A, Gamble GD, Grey A. Systematic review and statistical analysis of the integrity of 33 randomized controlled trials. Neurology 2016;87:2391-2402.
#'
#'
#' @format A data frame with 50 rows and 15 variables:
#' \describe{
#'   \item{study}{study ID}
#'   \item{var}{variable name}
#'   \item{n1,n2,n3,n4}{size of group 1, 2, 3, 4}
#'   \item{m1,m2,m3,m4}{mean of variable for group 1, 2, 3, 4}
#'   \item{s1,s2,s3,s4}{sd of variable for group 1, 2, 3, 4}
#'   \item{p}{reported p-value}
#' }
"SI_pvals_cont"
