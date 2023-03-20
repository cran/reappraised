#' Compares observed and expected distribution of p-values for categorical variables
#'
#' Creates plots of calculated p-value distribution and AUC (area under curve)\cr
#'
#' See also Bolland MJ, Gamble GD, Avenell A, Grey A, Lumley T. Baseline P value distributions in randomized trials were uniform for continuous but not categorical variables. J Clin Epidemiol 2019;112:67-76.\cr
#'
#' Returns a list containing 3 objects and (if verbose = TRUE) prints the plot pval_cat_calculated_pvalues
#'
#' @param df data frame generated from load_clean function
#' @param seed the seed to use for random number generation, default 0 = current date and time. Specify seed to make repeatable.
#' @param sims number of simulations, default -1 = function selects based on number of variables.
#' @param btsp number of bootstrap repeats used to generate 95% confidence interval around AUC
#' @param title optional title for plots
#' @param stat statistical test to be used 'chisq' or 'fisher'
#' @param stat.overide if 'yes' then test specified in stat will be used rather than values for stat in data frame
#' @param method 'sm', 'mix', or 'ind'. 'ind' does test on individual data, 'sm' summarises data and then does test on summary data, 'mix' does 'ind' for fisher and 'sm' for chisq. Duration varies with size of studies, test, and number of simulations. Experiment before running large simulations.
#' @param verbose TRUE or FALSE indicates whether progress bar and comments show and prints plot
#'
#'
#' @return list containing 3 objects as described
#'\itemize{
#'   \item pval_cat_calculated_pvalues = plots of calculated p-value distribution and AUC
#'   \item pval_cat_reported_pvalues = plots of reported p-value distribution and AUC (if p-values were reported)
#'   \item all_results = list containing
#'      \itemize{
#'         \item pval_cat_baseline_pvalues_data = data frame of all results used in calculations
#'         \item pval_cat_reported_pvalues= plot of reported p-value distribution
#'         \item pval_cat_auc_reported_pvalues = AUC of reported p-values
#'         \item pval_cat_calculated_pvalues = plot of calculated p-value distribution
#'         \item pval_cat_auc_calculated_pvalues= AUC of calculated p-values
#'         }}
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar
#' @importFrom dplyr select group_by left_join mutate summarise left_join filter slice rename starts_with
#' @importFrom tidyr gather spread separate unite
#' @importFrom boot boot boot.ci
#' @importFrom ggpubr ggarrange
#' @importFrom data.table rbindlist setDT setDF .N
#' @export pval_cat_fn
#'
#'
#' @examples
#' # load example data
#' pval_cat_data <- load_clean(import= "no", file.cat = "SI_cat_all", pval_cat= "yes",
#' format.cont = "wide")$pval_cat_data
#'
#' \donttest{
#' # run function (takes a few seconds)
#' pval_cat_fn(seed=10, sims = 50, btsp = 100)$pval_cat_calculated_pvalues
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
#' pval_cat_data <- load_clean(import= "yes", pval_cat = "yes", dir = path,
#'      file.name.cat = "reappraised_examples.xlsx", sheet.name.cat = "SI_cat_all",
#'      range.name.cat = "A:n", format.cat = "wide")$pval_cat_data}
#'
#' @md

#function for baseline p-values and AUC ---------------------------------
pval_cat_fn <- function (df = pval_cat_data, seed= 0, sims = -1, btsp = 500, title = "", stat= "chisq", stat.overide = "no",
                         method = "mix", verbose = TRUE) {

  #import dataset
  a <- df

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "pval_cat") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
                   deparse(substitute(df)), "$data_type != 'pval_cat'"))}

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  if (vb == "y") {pb = txtProgressBar(min=0, max=100, style = 2)}

  stat <- tolower(stat)

  if (!stat %in% c("chisq", "fisher")) {
    stop_mb ("Only statistical test options are chisq and fisher")
  }

  if (!method %in% c("mix", "sm", "ind")) {
    stop_mb ("Only method options are 'ind', 'sm', 'mix'")
  }

  if (tolower(substr(stat.overide,1,1)) == "y") {a$stat <- stat
  } else if ("stat" %in% colnames(a)) {a$stat <- ifelse(is.na(a$stat),stat, tolower(a$stat))
  } else {a$stat <- stat}

  if (length(setdiff(a$stat, c("chisq", "fisher")))>0) {
    if (vb == "y") {cat("\nChanging statistical tests other than chisq and fisher to chisq\n\n")}
    a$stat <- ifelse(!a$stat %in% c("chisq", "fisher"), "chisq", a$stat)
  }

  if (vb == "y") {
    #Pb1
    setTxtProgressBar(pb,5)
    cat("\nGenerating reference data ...\n\n")}

  #get reference p-values and proportions
  a$study_no <- as.numeric(factor(a$study))
  a$var_no <- as.numeric(factor(a$var))

  b <- a %>% tidyr::gather(k, v, starts_with("n"), na.rm= TRUE) %>%
    dplyr::mutate(k = paste0(gsub("[0-9]+","", k), "_", gsub("\\D", "", k))) %>%
    tidyr::separate(k, into = c("k1","grp")) %>% tidyr::spread (k1, v) %>% dplyr::mutate(grp = as.numeric(grp)) %>%
    dplyr::select(study_no, var_no, group, grp, level, levels_no, n, N, stat)

  b.exp <- b[rep(seq.int(1,nrow(b)), b$n), ]

  #rbindlist much faster than bind_rows
  if (sims ==-1) {sims <- round(10000/(nrow(a %>% filter(level ==1))),-1)} #choose sims so about 10,000 tests

  if (vb == "y") {
    if(sims * nrow(b.exp) > 1000000) {cat('There are more than 1 million simulated observations, and',
                                          sims*nrow(a %>% filter(level ==1)), 'Statistical tests to be done\n')
      cat("\nThis may take some time....\n\n")}
  }

  sim_exp <- data.table::rbindlist(replicate(sims, b.exp, simplify = FALSE)) #make simulations
  sim_exp$sim <- rep(1:sims, each = nrow(b.exp))

  if(seed == 0) {set.seed(Sys.time())
  } else {set.seed(seed)}
  sim_exp$rand <- runif(nrow(sim_exp))
  sim_exp$g1 <- ceiling(sim_exp$rand/ (1/sim_exp$group))

  #calculating from individual data
  if (method == "ind") {
    if (vb == "y" ) {cat("\nCalculating tests from individual data....\n\n")}

    #split into fisher or chisq
    if (vb == "y" & nrow(a %>% filter(level ==1 & stat == "fisher"))>0) {
      #Pb2
      setTxtProgressBar(pb,15)
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "fisher")), 'fisher`s exact tests ...\n\n')}

    f <- ch <- 0
    if ("fisher" %in% a$stat) {f <- 1
    sim_exp.f <- sim_exp [sim_exp$stat =="fisher",]
    r.f <- sim_exp.f %>%
      dplyr::group_by(sim, study_no, var_no) %>%
      dplyr::summarise(p=fisher.test(g1,level, workspace = 2e8)$p.value, .groups = "drop") %>% as.data.frame()}

    if (vb == "y" & nrow(a %>% filter(level ==1 & stat == "chisq"))>0) {
      #Pb3
      setTxtProgressBar(pb,30)
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "chisq")), 'Chi-square tests ...\n\n')}

    if ("chisq" %in% a$stat) {ch <- 1
    sim_exp.ch <- sim_exp [sim_exp$stat =="chisq",]
    r.ch <- sim_exp.ch %>%
      dplyr::group_by(sim, study_no, var_no) %>%
      dplyr::summarise(p=suppressWarnings(chisq.test(g1,level, correct= FALSE)$p.value), .groups = "drop") %>% as.data.frame()}

    if (f == 1 & ch ==1) {
      r <- rbind(cbind(r.ch, stat= "Chisq"), cbind(r.f, stat = "Fisher"))
      rm(r.ch, r.f)
    } else if (f ==1) {
      r <- r.f; rm (r.f); r$stat <- "Fisher"
    } else if (ch ==1) {
      r <- r.ch; rm (r.ch); r$stat <- "Chisq"
    }

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat)
  }

  #calculating from summary statistics
  if (method == "sm") {
    if (vb == "y" ) {cat("\nCalculating tests from summary statistics....\n\n")}

    #faster than using dplyr
    #sim_exp %>% count(sim, study_no, var_no, g1, level) etc

    data.table::setDT(sim_exp)
    z <- sim_exp[, .N, by=list(sim, study_no, var_no, levels_no, g1, level, group, stat)]
    data.table::setDF(z)
    z <- z %>% tidyr::unite(g, g1, level) %>% tidyr::spread(g, N) %>% replace(is.na(.), 0)

    #split into fisher or chisq
    if ("fisher" %in% a$stat) {f <- 1

    if (vb == "y" & nrow(a %>% dplyr::filter(level ==1 & stat == "fisher"))>0) {
      #Pb2
      setTxtProgressBar(pb,15)
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "fisher")), 'fisher`s exact tests ...\n\n')}

    z.f <- z [z$stat =="fisher",]
    nm <- colnames(z.f) [!colnames(z.f) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]
    for (i in 1:nrow(z.f)) {
      x <- as.numeric(z.f[i, nm [sub("_.*", "", nm) <= z.f$group [i] & sub(".*_", "", nm) <= z.f$levels_no [i]]])
      x <- matrix(x, byrow = TRUE, ncol = z.f$levels_no [i])
      z.f$p [i] <- fisher.test(x, workspace = 2e8)$p.value}
    }

    if ("chisq" %in% a$stat) {ch <- 1

    if (vb == "y" & nrow(a %>% dplyr::filter(level ==1 & stat == "chisq"))>0) {
      #Pb3
      setTxtProgressBar(pb,30)
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "chisq")), 'Chi-square tests ...\n\n')}

    z.ch <- z [z$stat =="chisq",]
    nm <- colnames(z.ch) [!colnames(z.ch) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]
    for (i in 1:nrow(z.ch)) {
      x <- as.numeric(z.ch[i, nm [sub("_.*", "", nm) <= z.ch$group [i] & sub(".*_", "", nm) <= z.ch$levels_no [i]]])
      x <- matrix(x, byrow = TRUE, ncol = z.ch$levels_no [i])
      z.ch$p [i] <- suppressWarnings(chisq.test(x, correct= FALSE)$p.value)}
    }

    if (f == 1 & ch ==1) {
      r <- rbind(z.ch, z.f)
    } else if (f ==1) {r <- z.f
    } else if (ch ==1) {r <- z.ch}

    r$stat <- paste0(toupper(substr(r$stat, 1, 1)), substr(r$stat, 2, nchar(r$stat)))

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat)
  }

  #calculating from mixed data
  if (method == "mix") {
    if (vb == "y" ) {cat("\nCalculating fisher from individual data and chisq from summary statistics....\n\n")}

    #split into fisher or chisq
    if (vb == "y" & nrow(a %>% dplyr::filter(level ==1 & stat == "fisher"))>0) {
      #Pb2
      setTxtProgressBar(pb,15)
      cat('\nDoing',sims*nrow(a %>% dplyr::filter(level ==1 & stat == "fisher")), 'fisher`s exact tests ...\n\n')}

    f <- ch <- 0
    if ("fisher" %in% a$stat) {f <- 1
    sim_exp.f <- sim_exp [sim_exp$stat =="fisher",]
    r.f <- sim_exp.f %>%
      dplyr::group_by(sim, study_no, var_no) %>%
      dplyr::summarise(p=fisher.test(g1,level, workspace = 2e8)$p.value, .groups = "drop") %>% as.data.frame()
    r.f$stat <- "Fisher"}

    if ("chisq" %in% a$stat) {ch <- 1
    z1 <- sim_exp [sim_exp$stat =="chisq",]
    data.table::setDT(z1)
    z <- z1[, .N, by=list(sim, study_no, var_no, levels_no, g1, level, group, stat)]
    data.table::setDF(z)
    z <- z %>% tidyr::unite(g, g1, level) %>% tidyr::spread(g, N) %>% replace(is.na(.), 0)

    if (vb == "y" & nrow(a %>% filter(level ==1 & stat == "chisq"))>0) {
      #Pb3
      setTxtProgressBar(pb,30)
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "chisq")), 'Chi-square tests ...\n\n')}

    z.ch <- z [z$stat =="chisq",]
    nm <- colnames(z.ch) [!colnames(z.ch) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]
    for (i in 1:nrow(z.ch)) {
      x <- as.numeric(z.ch[i, nm [sub("_.*", "", nm) <= z.ch$group [i] & sub(".*_", "", nm) <= z.ch$levels_no [i]]])
      x <- matrix(x, byrow = TRUE, ncol = z.ch$levels_no [i])
      z.ch$p [i] <- suppressWarnings(chisq.test(x, correct= FALSE)$p.value)}
    }

    if (f == 1 & ch ==1) {
      r <- rbind(z.ch %>% dplyr::select(sim, study_no, var_no, p, stat), r.f)
    } else if (f ==1) {r <- r.f
    } else if (ch ==1) {r <- z.ch}

    r$stat <- paste0(toupper(substr(r$stat, 1, 1)), substr(r$stat, 2, nchar(r$stat)))

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat)
  }


  if (vb == "y") {
    #Pb4
    setTxtProgressBar(pb,70)
    cat('\nDoing calculations and making graphs\n')}

  #now divide by deciles to get ref_p
  r$pval <- cut(r$p,c(seq(0,1,0.1)))

  p_ref <- as.data.frame(table (r$pval))
  p_ref$prop <- p_ref$Freq/sum(p_ref$Freq)

  #now get calculated p-values
  #calculate p-values
  #convert to wide format for level
  b <- a %>% dplyr::select(study, var, level, dplyr::starts_with("n", ignore.case = TRUE), stat) %>%
    tidyr::gather (k, v, dplyr::starts_with("n", ignore.case = TRUE)) %>%
    dplyr::mutate(k = paste0(gsub("[0-9]+","", k), "_", gsub("\\D", "", k))) %>%
    tidyr::separate(k, into = c("k1", "k2")) %>%
    tidyr::unite(k, k1, level, sep = "_") %>% tidyr::unite(k, k, k2, sep = "") %>% tidyr::spread (k, v)

  b <- dplyr::left_join(b, a %>% dplyr::select(study, var, levels_no ,group) %>% dplyr::group_by(study, var) %>%
                          dplyr::slice (1), by = c("study", "var"))
  b <- b %>% dplyr::select(-dplyr::starts_with("N", ignore.case = FALSE))

  #do as loop
  nm <- which(substr(colnames(b),1,1) == "n")
  b$p <- NA

  #calculate p by different methods
  for (i in 1:nrow(b)) {
    x <- na.omit(as.numeric(b [i, nm]))
    x <- matrix(x, byrow = TRUE, nrow = b$levels_no [i])
    if (b$stat [i] == "chisq") {b$p [i] <- suppressWarnings(chisq.test(x, correct = FALSE)$p.value)}
    if (b$stat [i] == "fisher") {b$p [i] <-  fisher.test(x, workspace = 2e8)$p.value}
  }

  p_calc <- cut(b$p, c(seq(0,1,0.1)))

  p_calc <- as.data.frame(table (p_calc))
  p_calc$prop <- p_calc$Freq/sum(p_calc$Freq)


  #graph for calculated p-values
  labels <- list(c("Calculated from summary data"),
                 paste(length(unique(a$study)), "trials"),
                 paste(nrow(a %>% filter(level ==1)), "variables"))

  p.g <- pval_graph (x = p_calc, round = "n", l= labels, t= title, ref = "y", rp = p_ref)

  pval <- p.g [[2]]
  p.g <- p.g [[1]]

  # graph for original p
  if(sum(!is.na(a$p))==0) {p.orig.g <- "There are no reported p-values"
  } else {

    #now get reported p-values
    a$p_num <- suppressWarnings(as.numeric(a$p))
    a$p_num <- ifelse(substring(a$p,1,1) == "<", 0.0001, a$p_num) #assume that < = <.05 or lower value so all decile 1

    p.orig <- as.data.frame(table(cut(a$p_num [a$level == 1],c(seq(0,1,0.1)))))

    labels.orig <- list(c("Reported p-values"),
                        paste(length(unique(a$study[!is.na(a$p_num)])), "trials"),
                        paste(nrow(a %>% filter (!is.na(p_num) & level ==1)), "variables"))

    p.orig.g <- pval_graph (x=p.orig, round= "n", l = labels.orig, t=title, ref = "y", rp = p_ref)

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

  #calculate ref data
  auc <- data.frame(r$p)
  auc_ref <- auc_calc()

  #calculated data
  auc <- data.frame(b$p)
  auc <- auc_calc()

  if (vb == "y") {
    #Pb4
    setTxtProgressBar(pb,75)
    cat("\nBoostrapping 95% CIs ...")
    cat("\nSpeed up by reducing value of btsp ....")
    cat("\nBut if btsp too small, warning message 'In norm. inter ...\n\n'")}

  R <- btsp # number of repeats
  cis <- boot::boot(b$p, auc_ci , R)
  ci <- boot::boot.ci (cis, type="perc")
  lci <- ci$percent [4]; uci <- ci$percent [5]

  #prepare values for AUC graph

  labels <- append(labels, list(pval, paste("AUC", round(auc$auc [1],2), "95% CI", round (lci,2), "-",round(uci,2))))

  auc.g <- auc_grph(auc, labels, round = "n", ref = "y", ra = auc_ref)

  #AUC original
  if(sum(!is.na(a$p))==0) {auc.orig.g <- "There are no reported p-values"
  } else {

    auc.orig <- data.frame(a$p_num [!is.na(a$p_num) & a$level == 1])
    auc.orig <- auc_calc(auc.orig)

    cis <- boot::boot(a$p_num [!is.na(a$p_num) & a$level == 1], auc_ci , R)
    ci <- boot::boot.ci (cis, type="perc")
    lci <- ci$percent [4]; uci <- ci$percent [5]

    labels.orig <- append(labels.orig, list(pval.orig,
                                            paste("AUC", round(auc.orig$auc [1],2), "95% CI", round (lci,2), "-", round(uci,2))))

    auc.orig.g <- auc_grph(auc.orig, labels.orig, round = "n", ref = "y", ra = auc_ref)
  }

  #reported pval graph
  if(sum(!is.na(a$p))==0) {rep_p <- "There are no reported p-values"
  } else {rep_p <- ggpubr::ggarrange(p.orig.g, auc.orig.g, ncol = 1, nrow= 2, align = "hv")}

  calc_p <- ggpubr::ggarrange(p.g, auc.g, ncol = 1, nrow= 2, align = "hv")

  a$p_num <- NULL

  a <- a %>% dplyr::select(-data_type, -study_no, -var_no) %>%
    dplyr::left_join(b %>% dplyr::select(study, var, p) %>% dplyr::rename (`p-value_calc_from_summary_stats` = p), by = c("study", "var")) %>%
    as.data.frame()

  results <- list (pval_cat_calculated_pvalues = calc_p,
                   pval_cat_reported_pvalues = rep_p,
                   all_results = list(pval_cat_baseline_pvalues_data = a,
                                      pval_cat_reported_pvalues= p.orig.g, pval_cat_auc_reported_pvalues = auc.orig.g,
                                      pval_cat_calculated_pvalues = p.g, pval_cat_auc_calculated_pvalues= auc.g))



  if (vb == "y") {
    print(calc_p)
    close(pb)}

  return(results)
}
