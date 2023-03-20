#' Compares differences between baseline means using Carlisle's montecarlo anova method
#'
#' Creates plots of distribution of p-values for differences in baseline means calculated using Carlisle's montecarlo anova method.\cr
#'
#' Method is from Carlisle JB, Loadsman JA. Evidence for non-random sampling in randomised, controlled trials by Yuhji Saitoh. Anaesthesia. 2017;72:17-27.\cr
#' R code is in appendix to paper. This function is adapted from that code.\cr
#' The function has two methods. The published code selects each variable from each study then generates
#' simulations for that variable using a row-wise approach with several loops. The adapted method is method = "orig".
#' The method = "alt" generates all the simulations at once and initially I thought was considerably faster, but in practice the time savings are small.\cr
#' The results from the two approaches will not be identical even if the same random number seed is used
#' because they use the generated random numbers in different orders but the p-values generated differ by about <0.1. Usually the
#' differences are close to 0.01 (although this depends on the number of simulations- more simulations = smaller differences).
#' The code that generates the p-value for each variable from the simulated means is essentially the same.
#'
#'
#' Returns a list containing 3 objects and (if verbose = TRUE) prints the plot anova_ecdf
#'
#' @param df dataframe generated from load_clean function
#' @param seed the seed to use for random number generation, default 0 = current date and time. Specify seed to make repeatable.
#' @param sims number of simulations, default -1 = function selects based on number of variables and sample size
#' @param btsp number of bootstrap repeats used to generate 95% confidence interval around AUC
#' @param method "orig" is adapted from original code; "alt" avoids using loops in the code (see details)
#' @param title optional title for plots
#' @param verbose TRUE or FALSE indicates whether progress bar and comments show and prints plot
#'
#' @return list containing 3 objects as described
#'\itemize{
#'   \item anova_ecdf = plot of cumulative distribution of calculated p-values compared to the expected uniform distribution
#'   \item anova_pvalues = plots of distribution of calculated p-values and AUC, as for pval_cont_fn()
#'   \item anova_all_results = list containing
#'      \itemize{
#'         \item anova_data = data frame of baseline data, with calculated p-values
#'         \item anova_pvals = plot of distribution of calculated p-values from anova_pvalues
#'         \item anova_auc = plot of AUC of calculated p-values from anova_pvalues
#'         }}
#'
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar
#' @importFrom dplyr select group_by left_join all_of summarise
#' @importFrom tidyr gather spread separate
#' @importFrom boot boot boot.ci
#' @importFrom ggpubr ggarrange
#' @importFrom data.table rbindlist as.data.table
#' @importFrom stats ks.test
#' @export anova_fn
#'
#'
#' @examples
#' # load example data
#' anova_data <- load_clean(import= "no", file.cont = "SI_pvals_cont",anova= "yes",
#' format.cont = "wide")$anova_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' anova_fn(seed=10, sims = 100, btsp = 100)$anova_ecdf
#'
#' # to import an excel spreadsheet (modify using local path,
#' # file and sheet name, range, and format):
#'
#' # get path for example files
#' path <- system.file("extdata", "reappraised_examples.xlsx", package = "reappraised",
#'                     mustWork = TRUE)
#' # delete file name from path
#' path <- sub("/[^/]+$", "", path)
#'
#' # load data
#' anova_data <- load_clean(import= "yes", anova = "yes", dir = path,
#'      file.name.cont = "reappraised_examples.xlsx", sheet.name.cont = "SI_pvals_cont",
#'      range.name.cont = "A:O", format.cont = "wide")$anova_data}
#'
#' @md


#function for anova pvalues ---------------------------------
anova_fn <- function (df= anova_data, method = "alt", seed= 0, sims = -1, btsp=500, title = "", verbose= TRUE) {

  if (! tolower(method) %in% c("orig", "alt")) {
    stop_mb("method must be 'orig', 'alt'")}

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  if (vb == "y") {pb = txtProgressBar(min=0, max=100, style = 2)}

  #Pb1
  if (vb == "y") {
    setTxtProgressBar(pb,10)
    cat("\nCalculating data for studies ...\n\n")
  }

  a <- df

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "anova") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
                   deparse(substitute(df)), "$data_type != 'anova'"))}

  #calculate dec places for mean for each variable
  #r can't handle 0s nor excel as numbers so take maximum dp (assuming that m1=x and m2 = x.1, m1 is actually x.0)
  nm.m <- colnames (a) [substring(colnames (a), 1, 1) == "m" ]
  nm.n <- colnames (a) [substring(colnames (a), 1, 1) == "n" ]
  nm.s <- colnames (a) [substring(colnames (a), 1, 1) == "s" & substring(colnames (a),2,2) %in% as.character(0:9)]

  #check if d exists
  d1 <- 0
  if ("d" %in% colnames (a)) {d1 <- 1}

  dp <- function (x) {
    ifelse(abs(x- round(x)) > .Machine$double.eps^0.5,
           nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(x))))), 0)}

  dps <- purrr::map_df (a [, nm.m], dp)
  a$dec <- apply(dps, 1, max, na.rm=TRUE)

  if (d1 ==1) {
    a$dec <- ifelse(is.na(a$d), a$dec, a$d)
    a$d <- NULL
  }

  rm(d1, dps, dp)

  #calculate stats for each trial variable
  a$study.mean <- rowSums(a[,nm.m]*a[,nm.n], na.rm = TRUE)/rowSums(a[,nm.n], na.rm=TRUE)
  a$study.sd <- (rowSums(a[,nm.s]^2*a[,nm.n], na.rm = TRUE)/rowSums(a[,nm.n], na.rm=TRUE))^.5
  a$study.nmean <- rowMeans(a[,nm.n],na.rm=TRUE)
  a$study.sem = a$study.sd/a$study.nmean^.5
  a$study.sumsqdiff <- rowSums((a[,nm.m]- a$study.mean)^2, na.rm = T)

  #create template for sims
  #suggested sims
  var <- c(20000, 10000, 5000, 1000, 500, 250, 100, 10)
  names (var) <- c(20, 50, 100, 500, 2000, 5000, 10000, 50000)
  sm <- names (var) [which.min(abs(var - nrow(a)))] %>% as.numeric()
  if (sims ==-1) {
    sims <- sm
    if (max(a [, nm.n], na.rm = TRUE) >500) {sm <- sims <- sims * .75}
  }

  #Pb2
  if(vb == "y") {
    setTxtProgressBar(pb,20)
    cat(paste0("\nSimulations: ",sims,", for ",nrow(a)," variables\n\n"))
    cat("Doing simulations...\nthis takes some time, usually at least 30-60s\n")
    cat("if cohorts have large n (eg n>500), this adds more time\n")
    cat(paste0("to speed up specify fewer simulations in the function call eg anova_fn (..., sims =",
               round(sims/2, 0)," ...)\n\n"))
  }

  if (tolower(method) == "alt") {

    #need popmean for each variable not each group, fastest to make single variable and merge 20% faster than
    #replicating a, generating the random#, then going to wide, which is 50% faster than going long and
    #manipulating from that point.

    #convert to long form them duplicate n=sim times
    sim <- a %>% dplyr::select(- dplyr::all_of(c(nm.m, nm.s))) %>%
      tidyr::gather(k,v, dplyr::all_of(c(nm.n)), na.rm = TRUE) %>%
      tidyr::separate(k, into = c("stat","grp"), sep=1) %>% tidyr::spread(stat, v) %>%
      dplyr::select(study, var, grp, n, study.mean, study.sem, study.sd, dec)

    sim_exp <- data.table::rbindlist(replicate(sims, sim, simplify = FALSE))
    sim_exp$sim <- rep(1:sims, each = nrow(sim))

    #create popmean
    sim1 <- a %>% dplyr::select(study, var, study.mean, study.sem)
    sim_exp1 <- data.table::rbindlist(replicate(sims, sim1, simplify = FALSE))
    sim_exp1$sim <- rep(1:sims, each = nrow(a))

    #generate random population mean for each study for each variable for each sim  (ie same mean for each group)
    if(seed ==0) {set.seed(Sys.time())
    } else {set.seed(seed)}

    sim_exp1$popmean <- rnorm(nrow(sim_exp1), sim_exp1$study.mean, sim_exp1$study.sem)
    sim_exp <- sim_exp1[, c("study", "var", "sim", "popmean")] [sim_exp, on=.(study, var, sim)]

    rm(sim_exp1)

    #take the mean of each sim
    m <- purrr::pmap(list(sim_exp$n, sim_exp$popmean, sim_exp$study.sd), rnorm)
    sim_exp$sim.mean <- vapply(m, function (x) sum (x) /length (x), FUN.VALUE = 1.0)
    sim_exp$round.sm <- round(sim_exp$sim.mean, sim_exp$dec)

    rm(m) #remove very large file

    sim1 <- sim_exp[, .(mean.sims = sum(round.sm * n)/sum (n)), by=list(study, var, sim)]
    sim_exp <- sim1[sim_exp, on=.(study, var, sim)]
    sim_exp$diff.simmean.meansims <- (sim_exp$round.sm - sim_exp$mean.sims)^2

    sim1 <-  sim_exp[, .(sum.sim.diffs = sum(diff.simmean.meansims)), by=list(study, var, sim)]
    sim2 <- as.data.table(a [, c("study", "var", "study.sumsqdiff")]) [sim1, on = .(study, var)]

    b <- sim2[, .(pvalue_anova = 1 -((sum(sum.sim.diffs < study.sumsqdiff) / sims
                                      + sum(sum.sim.diffs <= study.sumsqdiff) / sims)/2)), by=list(study, var)]
    setDF(b)
    }

  if (tolower(method) == "orig") {
    b <- a [, c("study", "var")]
    b$pvalue_anova <- NA
    if(seed ==0) {set.seed(Sys.time())
    } else {set.seed(seed)}
    for (r in 1:nrow(a)) {
      n <- a [r, nm.n] [!is.na( a [r, nm.n])]
      grp <- length(n)
      sim.mean <- sum.sim.diffs <- 0
      for (z in 1:sims){
        popmean <- rnorm(1, a$study.mean [r], a$study.sem [r])
        for(i in 1:grp){
          sim.mean [i] <- round(sum(rnorm(n [i], popmean, a$study.sd [r])/ n[i]), a$dec [r])
        }
        mean.sims <- sum(sim.mean * n) / sum (n)
        diff.simmean.meansims <- (sim.mean - mean.sims)^2
        sum.sim.diffs [z] = sum(diff.simmean.meansims)
      }
      b$pvalue_anova [r] = 1 -((sum(sum.sim.diffs < a$study.sumsqdiff [r]) / sims +
                                  sum(sum.sim.diffs <= a$study.sumsqdiff [r]) / sims)/2)

      #checked with faster method. With 1000 sims from sato/iwamoto Very similar p values; max difference 0.07
    }
  }

  #pb3
  if (vb == "y") {
    setTxtProgressBar(pb,75)
    cat("\nMaking plots ...")
    cat("\nFor AUC plot, Bootstrapping 95% CIs ...")
    cat("\nSpeed up by reducing value of btsp ....")
    cat("\nBut if btsp too small, warning message 'In norm. inter ...'\n")
  }

  #Plot ECDF with KS test
  #AUC functions
  #calculate AUC
  auc_calc <- function (x= auc, bt = "n") {
    auc <- x
    colnames(auc) <- "p" #p is x value
    auc$rank <- rank(auc$p, ties.method = "max")
    auc$y <- auc$rank/NROW(auc) # y is exp y value
    auc <- auc[order(auc$y),]
    auc <- rbind(c(0,0,0),auc)
    r1 <- 1:(NROW(auc)-1)
    r2 <- r1+1
    auc$h <- auc$w <- 0
    auc$w [2:(NROW(auc))] <- auc$p[r2] - auc$p[r1]
    auc$h [2:(NROW(auc))] <- (auc$y[r2] + auc$y[r1])/2
    auc$a <- auc$w *auc$h
    if (bt == "n") {
      auc$auc <- NA
      auc$auc [1] <- sum(auc$a)
      return (auc)
    } else {return (sum(auc$a))}
  }

  auc <- data.frame(b$pvalue_anova)
  auc <- auc_calc()
  ks <- suppressWarnings(ks.test(b$pvalue_anova, "punif")$p.value)
  #format for graph
  if (ks == 0) {pval <- paste("'p<", 2.2,"*'*","10^","-16")
  } else {
    if (ks <0.001) {p1 <- format(ks, scientific = TRUE, digit = 2)
    } else {p1 <- ks}
    p.e <- gregexpr("e",p1) [[1]]
    pval <- ifelse(p.e != -1,
                   paste("'p=",substring(as.character(p1),1,p.e-1),"*'*",
                         "10^",substring(as.character(p1),p.e+1,nchar(p1)),sep=""),
                   paste("'p='~'",signif(p1,2),"'", sep=""))
  }

  labels <- list(c("Carlisle Anova method"),
                 paste(length(unique(a$study)), "trials"),
                 paste(length(b$pvalue_anova), "variables"),
                 pval,"")


  ecdf <- auc_grph(auc, labels, round = "n")

  #do graph by decile and AUC from pval_cont function

  #now divide by deciles and plot graph
  b$p <- cut(b$pvalue_anova,c(seq(0,1,0.1)))
  p <- as.data.frame(table (b$p))

  #graph for calculated p-values
  labels <- list(c("Carlisle Anova method"),
                 paste(length(unique(a$study)), "trials"),
                 paste(length(b$pvalue_anova), "variables"))

  p.g <- pval_graph (x = p, round = "n", l= labels, t=title)

  pval <- p.g [[2]]
  p.g <- p.g [[1]]

  #AUC functions
  #AUC calculated above

  #helper function for boot
  auc_ci <- function(df, ind) {
    d <- data.frame(df [ind])
    sm <- auc_calc (x=d, bt = "y")
    return (sm)
  }

  R <- btsp # number of repeats
  cis <- boot::boot(b$pvalue_anova, auc_ci , R)
  ci <- boot::boot.ci (cis, type="perc")
  lci <- ci$percent [4]; uci <- ci$percent [5]

  #prepare values for AUC graph

  labels <- append(labels, list(pval, paste("AUC", round(auc$auc [1],2), "95% CI", round (lci,2), "-",round(uci,2))))

  auc.g <- auc_grph(auc, labels, round = "n")

  anova_p <- ggpubr::ggarrange(p.g, auc.g, ncol = 1, nrow= 2, align = "hv")

  results <- list (anova_ecdf = ecdf,
                   anova_pvalues = anova_p,
                   anova_all_results = list(anova_data = dplyr::left_join(a, b %>% dplyr::select(-p), by = c("study", "var")) %>%
                                              as.data.frame(), anova_pvals = p.g, anova_auc = auc.g))

  if (vb == "y") {
    print(ecdf)
    close(pb)}

  return (results)

}





