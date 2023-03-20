#' Compares observed and expected distribution of p-values for continuous variables
#'
#' Creates plots of calculated p-value distribution and AUC (area under curve)\cr
#'
#' Reference data is from (Carlisle 2017, Bolland 2021)\cr
#' Carlisle JB . Data fabrication and other reasons for non-random sampling in 5087 randomised, controlled trials in anaesthetic and general medical journals. Anaesthesia 2017;72:944–52 .2017\cr
#' Bolland MJ, Gamble GD, Grey A, Avenell A. Empirically generated reference proportions for baseline p values from rounded summary statistics. Anaesthesia 2020;75:1685-1687.\cr
#' See also Bolland MJ, Gamble GD, Avenell A, Grey A, Lumley T. Baseline P value distributions in randomized trials were uniform for continuous but not categorical variables. J Clin Epidemiol 2019;112:67-76.\cr
#' and Bolland MJ, Gamble GD, Avenell A, Grey A. Rounding, but not randomization method, non-normality, or correlation, affected baseline P-value distributions in randomized trials. J Clin Epidemiol 2019;110:50-62.
#'
#'
#' Returns a list containing 4 objects and (if verbose = TRUE) prints the plot pval_cont_calculated_pvalues
#'
#' @param df data frame generated from load_clean function
#' @param btsp number of bootstrap repeats used to generate 95% confidence interval around AUC
#' @param title optional title for plots
#' @param verbose TRUE or FALSE indicates whether progress bar and comments show and prints plot
#'
#' @return list containing 4 objects as described
#'\itemize{
#'   \item pval_cont_calculated_pvalues = plots of calculated p-value distribution and AUC
#'   \item pval_cont_reported_pvalues = plots of reported p-value distribution and AUC (if p-values were reported)
#'   \item pval_cont_ft_diff_calc_rep_p = flextable of distribution of differences in calculated and reported results
#'   \item all_results = list containing
#'      \itemize{
#'         \item pval_cont_baseline_pvalues_data = data frame of all results used in calculations
#'         \item pval_cont_diff_calc_rep_p = data frame of differences between calculated and reported p-values
#'         \item pval_cont_reported_pvalues= plot of reported p-value distribution
#'         \item pval_cont_auc_reported_pvalues = AUC of reported p-values
#'         \item pval_cont_calculated_pvalues = plot of calculated p-value distribution
#'         \item pval_cont_auc_calculated_pvalues= AUC of calculated p-values
#'         }}
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar
#' @importFrom dplyr select group_by left_join mutate count filter bind_cols rename full_join
#' @importFrom tidyr gather spread separate
#' @importFrom broom glance
#' @importFrom boot boot boot.ci
#' @importFrom ggpubr ggarrange
#' @importFrom flextable font fontsize valign set_table_properties
#' @export pval_cont_fn
#'
#'
#' @examples
#' # load example data
#' pval_cont_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", pval_cont= "yes",
#' format.cont = "wide")$pval_cont_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' pval_cont_fn(btsp=100)$pval_cont_calculated_pvalues
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
#' pval_cont_data <- load_clean(import= "yes", pval_cont = "yes", dir = path,
#'      file.name.cont = "reappraised_examples.xlsx", sheet.name.cont = "SI_pvals_cont",
#'      range.name.cont = "A1:O51", format.cont = "wide")$pval_cont_data}
#'
#' @md

#function for baseline p-values and AUC ---------------------------------
pval_cont_fn <- function (df= pval_cont_data, btsp = 500, title = "", verbose = TRUE) {

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  if (vb == "y") {pb = txtProgressBar(min=0, max=100, style = 2)}

  #convert to long format
  a <- df

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "pval_cont") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
                       deparse(substitute(df)), "$data_type != 'pval_cont'"))}

  #Pb1
  if (vb == "y") {
  setTxtProgressBar(pb,10)
  cat("Preparing data ...")
  }

  #add separator
  names <- colnames(a) [substr(colnames(a),1,1) %in% c("n","m") |
                          (substr(colnames(a),1,1) == "s" & colnames (a) != "study")]
  colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1),"_", gsub("\\D", "", names))

  b <- a %>%
    dplyr::select(-c("p","p_num","data_type")) %>%
    tidyr::gather(v, value, -c("study","var")) %>%
    tidyr::separate(v,c("var1", "group")) %>%
    tidyr::spread(var1, value)

  #delete missing data
  b <- b [! (is.na(b$m) | is.na(b$s) | is.na(b$n)),]

  #now generate dataset for analysis
  #uses anova for summary data
  xis <-  b$m + sqrt((b$s**2)/b$n)
  xns <-  b$n*b$m - (b$n-1)*xis;
  x1 <- seq(1,length(xis)*2,by=2)
  x2 <- seq(2,length(xis)*2,by=2)
  y <- f <- study <- var <- grp <- vector(mode="numeric", length = length(xis)*2)
  y[x1] <- xis
  y[x2] <- xns
  f[x1] <- b$n-1
  f[x2] <- 1
  study [x1] <- b$study
  study [x2] <- b$study
  var [x1] <- b$var
  var [x2] <- b$var
  grp [x1] <- b$group
  grp [x2] <- b$group
  z <- data.frame(y,f,study,var, grp)
  z.expand <- as.data.frame(lapply(z, rep, z$f))

  #Pb2
  if (vb == "y") {
  setTxtProgressBar(pb,25)
  cat("Calculating p-values ... this usually takes the longest time")
  }

  #lm for each variable
  fitted_models = z.expand %>% dplyr::group_by(study,var) %>%
    dplyr::group_modify(~ broom::glance(lm(y ~ as.factor(grp), data = .)))

  #now divide by deciles and plot graph
  fitted_models$p <- cut(fitted_models$p.value,c(seq(0,1,0.1)))
  p <- as.data.frame(table (fitted_models$p))

  #graph for calculated p-values
  labels <- list(c("Calculated from summary data"),
                 paste(length(unique(b$study)), "trials"),
                 paste(length(fitted_models$p.value), "variables"))

  p.g <- pval_graph (x = p, round = "y", l= labels, t=title)

  pval <- p.g [[2]]
  p.g <- p.g [[1]]

  # graph for original p
  if(sum(!is.na(a$p))==0) {p.orig.g <- "There are no reported p-values"
  } else {

    a$p_num <- ifelse(substring(a$p,1,1) == "<", 0.0001, a$p_num)
    a$p1 <- cut(a$p_num,breaks =seq(0,1,0.1))

    p.orig <- as.data.frame(table (a$p1))

    labels.orig <- list(c("Reported p-values"),
                        paste(length(unique(a$study[!is.na(a$p_num)])), "trials"),
                        paste(sum(!is.na(a$p_num)), "variables"))
    p.orig.g <- pval_graph (x=p.orig, round= "n", l = labels.orig, t= title)

    pval.orig <- p.orig.g [[2]]
    p.orig.g <- p.orig.g [[1]]
  }

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

  #helper function for boot
  auc_ci <- function(df, ind) {
    d <- data.frame(df [ind])
    sm <- auc_calc (x=d, bt = "y")
    return (sm)
  }

  auc <- data.frame(fitted_models$p.value)
  auc <- auc_calc()

  #pb3
  if (vb == "y") {
  setTxtProgressBar(pb,75)
  cat("\nBootstrapping 95% CIs ...")
  cat("\nSpeed up by reducing value of btsp ....")
  cat("\nBut if btsp too small, warning message 'In norm. inter ...'\n")
  }

  R <- btsp # number of repeats
  cis <- boot::boot(fitted_models$p.value, auc_ci , R)
  ci <- boot::boot.ci (cis, type="perc")
  lci <- ci$percent [4]; uci <- ci$percent [5]

  #prepare values for AUC graph

  labels <- append(labels, list(pval, paste("AUC", round(auc$auc [1],2), "95% CI", round (lci,2), "-",round(uci,2))))

  auc.g <- auc_grph(auc, labels, round = "y")

  #AUC original
  if(sum(!is.na(a$p))==0) {auc.orig.g <- "There are no reported p-values"
  } else {

    auc.orig <- data.frame(a$p_num [!is.na(a$p_num)])
    auc.orig <- auc_calc(auc.orig)

    cis <- boot::boot(a$p_num [!is.na(a$p_num)], auc_ci , R)
    ci <- boot::boot.ci (cis, type="perc")
    lci <- ci$percent [4]; uci <- ci$percent [5]

    labels.orig <- append(labels.orig, list(pval.orig,
                                            paste("AUC", round(auc.orig$auc [1],2), "95% CI", round (lci,2), "-", round(uci,2))))

    auc.orig.g <- auc_grph(auc.orig, labels.orig, round = "n")
  }

  #compare reported vs unreported p-values
  a <- dplyr::left_join (a, fitted_models[,c("study","var","p.value")], by = c("study", "var")) %>% as.data.frame()
  a$p.value <- round (a$p.value,3)

  if (sum (!is.na(a$p)) >0) {
    a$df <- round(abs(a$p.value- a$p_num),3)
    a$df.groups <- cut(a$df,seq(0,1,0.1))
    diff <- as.data.frame(table(a$df.groups))
    diff$prop <- diff$Freq/ sum(diff$Freq)
    diff1 <- diff
    diff[,1] <- paste(seq(0,.9,0.1)," - ", seq(0.1,1,0.1), "")
    diff[,2] <- paste(diff$Freq," (",round(diff$prop*100,0),")","")
    diff$prop <- NULL
    colnames (diff) <- c("Difference between \nreported and\n calculated p-values", "n (%)")
    colnames (diff1) [1] <- "difference_rep_calc_p-values"

    #make a flextable
    f <- flextable::flextable(diff)
    f <- flextable::font(f, fontname = "Arial")
    f <- flextable::fontsize(f, size = 8)
    f <- flextable::fontsize(f, part = "header", size = 8)
    f <- flextable::valign(f, i = 1, j = 2, valign = "bottom", part = "head")
    f <- flextable::set_table_properties(f, width = .3, layout = "autofit")

    #change names
    colnames (a) [colnames(a) == "df"] <- "difference_rep_calc_pvalues"
    a$df.groups <- NULL

    #reported pval graph
    rep_p <- ggpubr::ggarrange(p.orig.g, auc.orig.g, ncol = 1, nrow= 2, align = "hv")

  } else {
    diff1 <- rep_p <- f <- "There are no reported p-values"
  }

  #change names
  colnames (a) [colnames(a) == "p"] <- "reported_p-value"
  colnames (a) [colnames(a) == "p_num"] <- "reported_p-value_numeric"
  colnames (a) [colnames(a) == "p.value"] <- "p-value_calc_from_summary_stats"
  colnames (a) [colnames(a) == "var"] <- "variable"
  a$p1 <- NULL;

  calc_p <- ggpubr::ggarrange(p.g, auc.g, ncol = 1, nrow= 2, align = "hv")

  results <- list (pval_cont_calculated_pvalues = calc_p,
                   pval_cont_reported_pvalues = rep_p,
                   pval_cont_ft_diff_calc_rep_p = f,
                   all_results = list(pval_cont_baseline_pvalues_data = a, pval_cont_diff_calc_rep_p = diff1,
                                      pval_cont_reported_pvalues= p.orig.g, pval_cont_auc_reported_pvalues = auc.orig.g,
                                      pval_cont_calculated_pvalues = p.g, pval_cont_auc_calculated_pvalues= auc.g))

  if (vb == "y") {
    print(calc_p)
    close(pb)}

  return (results)

}

