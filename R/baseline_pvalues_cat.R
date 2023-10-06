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
#' @param stat statistical test to be used 'chisq', 'fisher', 'midp' or 'midp.epitools' (from epitools package), 'midp.sas' (as calculated in SAS), or combinations -if chisq is not appropriate because expected cells<5, use second test: 'chi_fish', 'chi_midp' or 'chi_midp.epi','chi_midp.sas'
#' @param stat.override if 'yes' then test specified in stat will be used rather than values for stat in data frame
#' @param fisher.sim "yes" or "no" indicator whether to allow fisher test to simulate p-values for >2*2 tables
#' @param fish.n.sims number of simulations to use in Fisher test, default 10,000
#' @param method 'sm', 'mix', or 'ind'. 'ind' does test on individual data, 'sm' summarises data and then does test on summary data, 'mix' does 'ind' for fisher and 'sm' for others. Duration varies with size of studies, test, and number of simulations. Experiment before running large simulations.
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
#' @importFrom dplyr select group_by left_join mutate summarise left_join filter slice rename starts_with arrange
#' @importFrom tidyr gather spread separate unite
#' @importFrom boot boot boot.ci
#' @importFrom ggpubr ggarrange
#' @importFrom data.table rbindlist setDT setDF .N
#' @importFrom purrr possibly
#' @importFrom epitools riskratio.wald
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
pval_cat_fn <- function (df = pval_cat_data, seed= 0, sims = -1, btsp = 500, title = "", stat= "chi_midp", stat.override = "no",
                         fisher.sim = "y", fish.n.sims = 10000, method = "mix", verbose = TRUE) {

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

  if (!method %in% c("mix", "sm", "ind")) {
    stop_mb ("Only method options are 'ind', 'sm', 'mix'")
  }

  if (tolower(substr(stat.override,1,1)) == "y") {a$stat <- stat
  } else if ("stat" %in% colnames(a)) {a$stat <- ifelse(is.na(a$stat),stat, tolower(a$stat))
  } else {a$stat <- stat}

  if (!stat %in% c("chisq", "fisher", "chi_fish", "chi_midp", "chi_midp.epi","chi_midp.sas","midp",
                   "midp.epitools", "midp.sas")) {
    stop_mb ("Only statistical test options are chisq, fisher, midp or midp.epitools, midp.sas, \nand combinations (chi_fish, chi_midp,, chi_midp.epi, chi_midp.sas)")
  }

  if (length(setdiff(a$stat, c("chisq", "fisher", "chi_fish", "chi_midp", "chi_midp.epi","chi_midp.sas","midp",
                               "midp.epitools", "midp.sas")))>0) {
    if (vb == "y") {cat("\nChanging statistical tests other than permitted tests to default of chi_midp\n\n")}
    a$stat <- ifelse(!a$stat %in% c("chisq", "fisher", "chi_fish", "chi_midp", "chi_midp.epi","chi_midp.sas","midp",
                                    "midp.epitools", "midp.sas"), "chi_midp", a$stat)
  }

  if (sum(a$levels_no*a$group > 4 & a$stat %in% c("chi_midp", "chi_midp.epi","chi_midp.sas",
          "midp", "midp.epitools", "midp.sas"))>0) {
    if (vb == "y") {cat("\nMid-p tests only available for 2*2 tables, if levels or groups >2, changing to chi_fish\n\n")}
    a$stat <- ifelse(a$levels_no*a$group > 4 & a$stat %in% c("chi_midp", "chi_midp.epi","chi_midp.sas",
                    "midp", "midp.epitools", "midp.sas"),"chi_fish", a$stat)
  }

  if (method == "ind" & sum (a$stat %in% c("chi_midp", "chi_midp.epi","chi_midp.sas",
                                           "midp", "midp.epitools", "midp.sas"))>0) {
    if (vb == "y") {cat("\nMid-p tests only available for summary data, changing to chi_fish\n\n")}
    a$stat <- ifelse(a$stat %in% c("chi_midp", "chi_midp.epi","chi_midp.sas",
                                   "midp", "midp.epitools", "midp.sas"),"chi_fish", a$stat)
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
    cat('\nNumber of sims: ', sims, '\n\n')
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
    f <- ch <- 0

    if (vb == "y" & nrow(a %>% dplyr::filter(level ==1 & stat != "fisher"))>0) {
      #Pb2
      setTxtProgressBar(pb,30)
      cat('\nDoing',sims*nrow(a %>% dplyr::filter(level ==1 & stat != "fisher")), 'Chi-square tests ...\n\n')}

    if ("chisq" %in% sim_exp$stat | "chi_fish" %in% sim_exp$stat) {ch <- 1
    sim_exp.ch <- sim_exp [sim_exp$stat =="chisq" | sim_exp$stat =="chi_fish", ]
    ch_f <- function (x,y) {
      chi <- suppressWarnings(chisq.test(x, y, correct= FALSE))
      chi1 <- ifelse(length(chi$expected) == 4 & sum (chi$expected <5) >0 |
                       length(chi$expected) > 4 & (sum(chi$expected <1) >0 |
                                                     sum(chi$expected <5) / length(chi$expected) > 0.2),
                     "exp < 5", "")
      res <- list(chi$p.value, chi1) }

    r.ch <- sim_exp.ch %>%
      dplyr::group_by(sim, study_no, var_no, stat) %>%
      dplyr::summarise(p1= suppressWarnings(list(purrr::possibly(ch_f, otherwise = list(NaN, ""))(g1,level))), .groups = "drop") %>%
      as.data.frame()
    r.ch$warn <- r.ch$p <- NA
    for (i in 1:nrow(r.ch)) {
      r.ch$p [i] <- r.ch$p1 [[i]] [[1]]
      r.ch$warn [i] <- ifelse(r.ch$p1 [[i]] [[2]] == "exp < 5" ,"exp < 5", NA)}
    r.ch$p1 <- NULL

    #select warning if needed
    warn <- r.ch [!is.na(r.ch$warn) & r.ch$stat == "chi_fish", ]
    if (nrow(warn) >0) {
      sim_exp <- dplyr::left_join(sim_exp, warn %>% dplyr::select(-p, -stat), by = c("sim", "study_no", "var_no")) %>%
        dplyr::mutate(stat = ifelse(!is.na(warn), "fisher", stat))
      r.ch <-  r.ch [!(!is.na(r.ch$warn) & r.ch$stat == "chi_fish"), ]}
    }

    if (vb == "y") {
      nfish <- nrow(sim_exp %>% dplyr::filter(level ==1 & grp==1 & stat == "fisher") %>% dplyr::group_by(sim, study_no, var_no) %>%
                      dplyr::filter(row_number()==1))
      if(nfish > 0) {
      #Pb3
      setTxtProgressBar(pb,15)
      cat('\nDoing', nfish, 'fisher`s exact tests ...\n\n')}
    }

    if ("fisher" %in% sim_exp$stat) {f <- 1
    sim_exp.f <- sim_exp [sim_exp$stat =="fisher",]

    if (tolower(substr(fisher.sim,1,1)) == "y") {
      fn_f <- function (x,y) {fisher.test(x, y, workspace = 2e8, simulate.p.value = TRUE, B= fish.n.sims)$p.value}
      suppressWarnings(r.f <- sim_exp.f %>%
                         dplyr::group_by(sim, study_no, var_no, stat, warn) %>%
                         dplyr::summarise(p= purrr::possibly(fn_f, otherwise = 1)(g1, level), .groups = "drop") %>%
                         as.data.frame())
    } else {
      fn_f <- function (x,y) {fisher.test(x, y, workspace = 2e8)$p.value}
      suppressWarnings(r.f <- sim_exp.f %>%
                         dplyr::group_by(sim, study_no, var_no, stat, warn) %>%
                         dplyr::summarise(p= purrr::possibly(fn_f, otherwise = 1)(g1, level), .groups = "drop") %>% as.data.frame())
    }
    }

    if (f == 1 & ch ==1) {
      r <- rbind(r.ch, r.f)
     } else if (f ==1) {
      r <- r.f
    } else if (ch ==1) {
      r <- r.ch;
    }

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat, warn)

    r$stat <- paste0(toupper(substr(r$stat, 1, 1)), substr(r$stat, 2, nchar(r$stat)))
  }

  #calculating from summary statistics
  #note in this case it is group_levels
  if (method == "sm") {
    if (vb == "y" ) {cat("\nCalculating tests from summary statistics....\n\n")}

    #faster than using dplyr
    #sim_exp %>% count(sim, study_no, var_no, g1, level) etc

    data.table::setDT(sim_exp)
    z <- sim_exp[, .N, by=list(sim, study_no, var_no, levels_no, g1, level, group, stat)]
    data.table::setDF(z)
    z <- z %>% arrange(level, g1) %>% tidyr::unite(g, g1, level) %>%
      dplyr::mutate(g = factor(g, levels = unique(g))) %>% tidyr::spread(g, N) %>% as.data.frame()

    g_ord <- colnames(b %>% dplyr::arrange (level, grp) %>% tidyr::unite(g, grp, level) %>%
          dplyr::mutate(g = factor(g, levels = unique(g))) %>% tidyr::spread(g, N) %>%
          dplyr::select(-(study_no:stat), -n))

    z [, setdiff(g_ord, colnames(z))] <- as.integer(NA)
    z <- z [, c("sim","study_no","var_no","levels_no", "group", "stat", g_ord)]

    #%>% replace(is.na(.), 0) doesn't work here, replace 0s only missing from group data

    nm <- colnames(z) [!colnames(z) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]

    zy <- rowSums (!is.na(z [, nm])) != z$group * z$levels_no
    for (i in 1:nrow(z)) {
      if (zy [i]) {
        zz <- nm [which(as.numeric(sub("_.*", "", nm)) <= z$group [i] & as.numeric(sub(".*_", "", nm)) <= z$levels_no [i])]
        z [i, zz] <- replace (z [i, zz], is.na(z [i, zz]), 0)}
    }
    rm (zy, zz)

    if (vb == "y") {
      #Pb2
      setTxtProgressBar(pb,15)
      if (nrow(a %>% dplyr::filter(level ==1 & stat == "fisher"))>0) {
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat == "fisher")), 'fisher`s exact tests ...\n\n')}

      if (nrow(a %>% dplyr::filter(level ==1 & stat %in% c("chisq","chi_fish","chi_midp",
                                                           "chi_midp.epi","chi_midp.sas")))>0) {
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat %in% c("chisq","chi_fish", "chi_midp",
                              "chi_midp.epi","chi_midp.sas" ))), 'Chi-square tests ...\n\n')}

      if (nrow(a %>% dplyr::filter(level ==1 & stat %in% c("midp","midp.epitools", "midp.sas")))>0) {
        cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat %in% c("midp","midp.epitools", "midp.sas" ))),
                                'Mid-p tests ...\n\n')}
    }

    nm <- colnames(z) [!colnames(z) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]
    z$warn <- NA; zz <- 0
    for (i in 1:nrow(z)) {
      x <- as.numeric(z [i, nm [as.numeric(sub("_.*", "", nm)) <= z$group [i] & as.numeric(sub(".*_", "", nm)) <= z$levels_no [i]]])
      x <- matrix(x, byrow = TRUE, nrow = z$levels_no [i])
      #do chisquare first
      if(z$stat [i] %in% c("chisq","chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
        chi <- suppressWarnings(chisq.test(x, correct = FALSE))
        z$p [i] <- chi$p.value
        #if warning, record and change stat if needed
        if(length(chi$expected) == 4 & sum (chi$expected <5) >0 |
           length(chi$expected) > 4 & (sum(chi$expected <1) >0 | sum(chi$expected <5) / length(chi$expected) > 0.2)) {
          z$warn [i] <- "exp < 5"
          if(z$stat [i] %in% c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
            z$stat [i] <- c("fisher", "midp", "midp", "midp.sas") [which (z$stat[i] == c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas"))]
          }
        }
      }

        if(z$stat [i] == "fisher") {
          if (sum(rowSums(x)==0) >0) {z$p [i] <- 1
          } else if (tolower(substr(fisher.sim,1,1)) == "y") {
            z$p [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8, simulate.p.value = TRUE, B= fish.n.sims)$p.value)
          } else {z$p [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8)$p.value)}
        }

        if(z$stat [i] %in% c("midp", "midp.epitools")) {
          #only for 2*2 tables
          z$p [i] <- suppressWarnings(epitools::riskratio.wald(x)$p.value [2] [1])}

        if(z$stat [i] == "midp.sas") {
          #sas midpvalues only for 2*2 tables also
          d1 = x [[1]]; d2 = x [[2]]; d3= x[[3]]; d4 = x[[4]];
          d5 = d1 +d2 +d3 + d4
          #multiply loop
          num = c(seq(max(d1+d2,1),1,-1), seq(max(d3+d4,1),1,-1), seq(max(d1+d3,1),1,-1), seq(max(d2+d4,1),1,-1))
          denom = c(seq(max(d5,1),1,-1), seq(max(d1,1),1,-1), seq(max(d2,1),1,-1), seq(max(d3,1),1,-1) ,seq(max(d4,1),1,-1))
          len = length(num) -length (denom)
          if (len > 0) {denom <- c(denom, rep(1, len)) }
          if (len < 0) {num <- c(num, rep(1, abs(len))) }
          p <- 100000
          num <- sort(num, decreasing = TRUE)
          denom <- sort(denom, decreasing = TRUE)
          for (j in 1:length(denom)) {
            p <- p * num [j] / denom [j]
          }
          p <- p / 100000
          z$p [i] <- round(suppressWarnings(fisher.test(x, workspace = 2e8)$p.value) - 0.5* p,9)
        }
      if (vb== "y" & zz==0 & i/nrow(z) >.33) {setTxtProgressBar(pb,30); zz <- 1}
      if (vb== "y" & zz==1 & i/nrow(z) >.66) {setTxtProgressBar(pb,50); zz <- 2}
    }

    r <- z

    r$stat <- paste0(toupper(substr(r$stat, 1, 1)), substr(r$stat, 2, nchar(r$stat)))

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat, warn)
  }

  if (method == "mix") {
    if (vb == "y" ) {cat("\nCalculating fisher from individual data and others from summary statistics....\n\n")}

    #split into fisher or others
    f <- ch <- 0
    if (sum(sim_exp$stat %in% c("chisq","chi_fish","chi_midp","chi_midp.epi","chi_midp.sas",
                            "midp","midp.epitools", "midp.sas")) >0) {ch <- 1

    z1 <- sim_exp [sim_exp$stat !="fisher",]
    data.table::setDT(z1)
    z <- z1[, .N, by=list(sim, study_no, var_no, levels_no, g1, level, group, stat)]
    data.table::setDF(z)
    z <- z %>% arrange(level, g1) %>% tidyr::unite(g, g1, level) %>%
      dplyr::mutate(g = factor(g, levels = unique(g))) %>% tidyr::spread(g, N) %>% as.data.frame()

    g_ord <- colnames(b %>% dplyr::arrange (level, grp) %>% tidyr::unite(g, grp, level) %>%
                        dplyr::mutate(g = factor(g, levels = unique(g))) %>% tidyr::spread(g, N) %>%
                        dplyr::select(-(study_no:stat), -n))

    z [, setdiff(g_ord, colnames(z))] <- as.integer(NA)
    z <- z [, c("sim","study_no","var_no","levels_no", "group", "stat", g_ord)]

    #%>% replace(is.na(.), 0) doesn't work here, replace 0s only missing from group data

    nm <- colnames(z) [!colnames(z) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]

    zy <- rowSums (!is.na(z [, nm])) != z$group * z$levels_no
    for (i in 1:nrow(z)) {
      if (zy [i]) {
        zz <- nm [which(as.numeric(sub("_.*", "", nm)) <= z$group [i] & as.numeric(sub(".*_", "", nm)) <= z$levels_no [i])]
        z [i, zz] <- replace (z [i, zz], is.na(z [i, zz]), 0)}
    }
    rm (zy, zz)

    if (vb == "y") {
      #Pb2
      setTxtProgressBar(pb,15)
    if (nrow(a %>% dplyr::filter(level ==1 & stat %in% c("chisq","chi_fish","chi_midp",
                                                         "chi_midp.epi","chi_midp.sas")))>0) {
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat %in% c("chisq","chi_fish", "chi_midp",
                                                                   "chi_midp.epi","chi_midp.sas" ))), 'Chi-square tests ...\n\n')}

    if (nrow(a %>% dplyr::filter(level ==1 & stat %in% c("midp","midp.epitools", "midp.sas")))>0) {
      cat('\nDoing',sims*nrow(a %>% filter(level ==1 & stat %in% c("midp","midp.epitools", "midp.sas" ))),
          'Mid-p tests ...\n\n')}
    }
    nm <- colnames(z) [!colnames(z) %in% c("sim","study_no","var_no", "stat","levels_no", "group")]
    z$warn <- NA; zz <- 0

    for (i in 1:nrow(z)) {
      x <- as.numeric(z [i, nm [as.numeric(sub("_.*", "", nm)) <= z$group [i] & as.numeric(sub(".*_", "", nm)) <= z$levels_no [i]]])
      x <- matrix(x, byrow = TRUE, nrow = z$levels_no [i])
      #do chisquare first
      if(z$stat [i] %in% c("chisq","chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
        chi <- suppressWarnings(chisq.test(x, correct = FALSE))
        z$p [i] <- chi$p.value
        #if warning, record and change stat if needed
        if(length(chi$expected) == 4 & sum (chi$expected <5) >0 |
           length(chi$expected) > 4 & (sum(chi$expected <1) >0 | sum(chi$expected <5) / length(chi$expected) > 0.2)) {
          z$warn [i] <- "exp < 5"
          if(z$stat [i] %in% c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
            z$stat [i] <- c("fisher", "midp", "midp", "midp.sas") [which (z$stat[i] == c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas"))]
          }
        }
      }

      if(z$stat [i] %in% c("midp", "midp.epitools")) {
        #only for 2*2 tables
        z$p [i] <- suppressWarnings(epitools::riskratio.wald(x)$p.value [2] [1])}

      if(z$stat [i] == "midp.sas") {
        #sas midpvalues only for 2*2 tables also
        d1 = x [[1]]; d2 = x [[2]]; d3= x[[3]]; d4 = x[[4]];
        d5 = d1 +d2 +d3 + d4
        #multiply loop
        num = c(seq(max(d1+d2,1),1,-1), seq(max(d3+d4,1),1,-1), seq(max(d1+d3,1),1,-1), seq(max(d2+d4,1),1,-1))
        denom = c(seq(max(d5,1),1,-1), seq(max(d1,1),1,-1), seq(max(d2,1),1,-1), seq(max(d3,1),1,-1) ,seq(max(d4,1),1,-1))
        len = length(num) -length (denom)
        if (len > 0) {denom <- c(denom, rep(1, len)) }
        if (len < 0) {num <- c(num, rep(1, abs(len))) }
        p <- 100000
        num <- sort(num, decreasing = TRUE)
        denom <- sort(denom, decreasing = TRUE)
        for (j in 1:length(denom)) {
          p <- p * num [j] / denom [j]
        }
        p <- p / 100000
        z$p [i] <- round(suppressWarnings(fisher.test(x, workspace = 2e8)$p.value) - 0.5* p,9)
      }
      if (vb== "y" & zz==0 & i/nrow(z) >.33) {setTxtProgressBar(pb,30); zz <- 1}
      if (vb== "y" & zz==1 & i/nrow(z) >.66) {setTxtProgressBar(pb,45); zz <- 2}
    }

    #select warning if needed
    warn <- z [z$stat == "fisher", ] %>% dplyr::select(sim, study_no, var_no, warn)
    if (nrow(warn) >0) {
    sim_exp <- dplyr::left_join(sim_exp, warn, by = c("sim", "study_no", "var_no")) %>%
      dplyr::mutate(stat = ifelse(!is.na(warn), "fisher", stat))
    z <-  z [z$stat != "fisher", ]}
  }

  if (vb == "y") {
    nfish <- nrow(sim_exp %>% dplyr::filter(level ==1 & grp==1 & stat == "fisher") %>% dplyr::group_by(sim, study_no, var_no) %>%
                    dplyr::filter(row_number()==1))
    if(nfish > 0) {
      #Pb3
      setTxtProgressBar(pb,60)
      cat('\nDoing', nfish, 'fisher`s exact tests ...\n\n')}
  }

  if ("fisher" %in% sim_exp$stat) {f <- 1
    sim_exp.f <- sim_exp [sim_exp$stat =="fisher",]

    if (tolower(substr(fisher.sim,1,1)) == "y") {
      fn_f <- function (x,y) {fisher.test(x, y, workspace = 2e8, simulate.p.value = TRUE, B= fish.n.sims)$p.value}
      suppressWarnings(r.f <- sim_exp.f %>%
                         dplyr::group_by(sim, study_no, var_no) %>%
                         dplyr::summarise(p= purrr::possibly(fn_f, otherwise = 1)(g1, level), .groups = "drop") %>%
                         as.data.frame())
    } else {
      fn_f <- function (x,y) {fisher.test(x, y, workspace = 2e8)$p.value}
      suppressWarnings(r.f <- sim_exp.f %>%
                         dplyr::group_by(sim, study_no, var_no) %>%
                         dplyr::summarise(p= purrr::possibly(fn_f, otherwise = 1)(g1, level), .groups = "drop") %>% as.data.frame())
    }
    r.f$warn <- NA}

 if (f == 1 & ch ==1) {
      r <- rbind(z %>% dplyr::select(sim, study_no, var_no, p, stat, warn), cbind(r.f, stat = "Fisher"))
      } else if (f ==1) {
      r <- r.f; r$stat <- "Fisher"
    } else if (ch ==1) {
      r <- z %>% dplyr::select(sim, study_no, var_no, p, stat, warn)}

    r$stat <- paste0(toupper(substr(r$stat, 1, 1)), substr(r$stat, 2, nchar(r$stat)))

    #merge back
    r <- dplyr::left_join(r, a %>% dplyr::filter (level ==1) %>% dplyr::select(study, var, study_no, var_no), by= c("study_no","var_no")) %>%
      dplyr::select(sim, study, var, p, stat, warn)
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
  #unlike sims above this is level * group as with cat_all, however order of variables is same, but suffixes arent

  b <- a %>% dplyr::select(study, var, level, dplyr::starts_with("n", ignore.case = TRUE), stat) %>%
    tidyr::gather (k, v, dplyr::starts_with("n", ignore.case = TRUE)) %>%
    dplyr::mutate(k = paste0(gsub("[0-9]+","", k), "_", gsub("\\D", "", k))) %>%
    tidyr::separate(k, into = c("k1", "k2")) %>%
    dplyr::arrange(level, k2) %>%  tidyr::unite(k, k1, level, sep = "_") %>%
    tidyr::unite(k, k, k2, sep = "_") %>% dplyr::mutate(k = factor(k, levels = unique(k))) %>%
    tidyr::spread (k, v)

  b <- dplyr::left_join(b, a %>% dplyr::select(study, var, levels_no ,group) %>% dplyr::group_by(study, var) %>%
                          dplyr::slice (1), by = c("study", "var"))
  b <- b %>% dplyr::select(-dplyr::starts_with("N", ignore.case = FALSE))

  #do as loop
  nm <- which(substr(colnames(b),1,1) == "n")
  b$p <- NA

  #calculate p by different methods
  b$Chisquare.warn <- NA
  for (i in 1:nrow(b)) {
    x <- na.omit(as.numeric(b [i, nm]))
    x <- matrix(x, byrow = TRUE, nrow = b$levels_no [i])
    if (b$stat [i] %in% c("chisq","chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
      chi <- suppressWarnings(chisq.test(x, correct = FALSE))
      b$p [i] <- chi$p.value
      if(length(chi$expected) == 4 & sum (chi$expected <5) >0 |
         length(chi$expected) > 4 & (sum(chi$expected <1) >0 | sum(chi$expected <5) / length(chi$expected) > 0.2)) {
        b$Chisquare.warn [i] <- "exp < 5"
        if(b$stat [i] %in% c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas")) {
          b$stat [i] <- c("fisher", "midp", "midp", "midp.sas") [which (b$stat[i] == c("chi_fish","chi_midp", "chi_midp.epi","chi_midp.sas"))]
        }
      }
    }

    if (b$stat [i] == "fisher") {
      if (sum(rowSums(x)==0) >0) {b$p [i] <- 1
      } else if (tolower(substr(fisher.sim,1,1)) == "y") {
        b$p [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8, simulate.p.value = TRUE, B= fish.n.sims)$p.value)
      } else {b$p [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8)$p.value)}
    }

    if(b$stat [i] %in% c("midp", "midp.epitools")) {
      #only for 2*2 tables
      b$p [i] <- suppressWarnings(epitools::riskratio.wald(x)$p.value [2] [1])}

    if(b$stat [i] == "midp.sas") {
      #sas midpvalues only for 2*2 tables also
      d1 = x [[1]]; d2 = x [[2]]; d3= x[[3]]; d4 = x[[4]];
      d5 = d1 +d2 +d3 + d4
      #multiply loop
      num = c(seq(max(d1+d2,1),1,-1), seq(max(d3+d4,1),1,-1), seq(max(d1+d3,1),1,-1), seq(max(d2+d4,1),1,-1))
      denom = c(seq(max(d5,1),1,-1), seq(max(d1,1),1,-1), seq(max(d2,1),1,-1), seq(max(d3,1),1,-1) ,seq(max(d4,1),1,-1))
      len = length(num) -length (denom)
      if (len > 0) {denom <- c(denom, rep(1, len)) }
      if (len < 0) {num <- c(num, rep(1, abs(len))) }
      p <- 100000
      num <- sort(num, decreasing = TRUE)
      denom <- sort(denom, decreasing = TRUE)
      for (j in 1:length(denom)) {
        p <- p * num [j] / denom [j]
      }
      p <- p / 100000
      b$p [i] <- round(suppressWarnings(fisher.test(x, workspace = 2e8)$p.value) - 0.5* p,9)
    }
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
    cat("\nBut if btsp too small, warning message 'In norm. inter ...'\n\n")}

  R <- btsp # number of repeats
  cis <- boot::boot(b$p [!is.na(b$p)], auc_ci , R)
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
    dplyr::left_join(b %>% dplyr::select(study, var, p, Chisquare.warn) %>% dplyr::rename (`p-value_calc_from_summary_stats` = p), by = c("study", "var")) %>%
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
