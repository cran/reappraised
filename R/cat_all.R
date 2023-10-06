#' Compares observed and expected distribution of all categorical (binomial) variables
#'
#' Creates plots of observed to expected numbers and ratios for the binomial variables and/or compares reported and calculated p-values for the variables\cr
#' Reference: Bolland MJ, Gamble GD, Avenell A, Cooper DJ, Grey A. Distributions of baseline categorical variables were different from the expected distributions in randomized trials with integrity concerns. J Clin Epidemiol. 2023;154:117-124
#'
#'
#' Returns a list containing objects described below and (if verbose = TRUE) prints the flextable cat_all_diff_calc_rep_ft and/or graph cat_all_graph depending on options chosen
#'
#' @param df data frame generated from load_clean function
#' @param comp.pvals "yes" or "no" indicator whether reported and calculated p-values should be compared
#' @param fisher.sim "yes" or "no" indicator whether to allow fisher test to simulate p-values for >2*2 tables
#' @param fish.n.sims number of simulations to use in Fisher test, default 10,000
#' @param binom "yes" or "no" indicator whether observed to expected distributions of binomial variables should be calculated
#' @param two_levels "yes" or "no" indicator whether variables with more than 2 levels should be collapsed to 2 levels
#' @param del.disparate if yes, data in which the absolute difference between group sizes is >20% are deleted
#' @param excl.level "yes" or "no" indicator whether one level of a variable should be deleted. Deleted level is chosen randomly using seed parameter.
#' @param seed seed for random number generator, default 0 = current date and time. Specify seed to make repeatable.
#' @param title title name for plots (optional)
#' @param verbose TRUE or FALSE indicates whether progress bar and comments show and flextable or plot or both are printed
#'
#' @return list containing objects as described\cr
#'
#' if p-value comparison used:
#'\itemize{
#'   \item cat_all_pvals = data frame of data for comparison of reported and calculated p-values
#'   \item cat_all_diff_calc_rep_ft = flextable of comparison of reported and calculated p-values
#'   \item cat_all_diff_calc_rep_data = data frame used to make flextable
#'   \item cat_all_diff_thresh_ft = flextable of comparison of reported and calculated p-values when only threshold given
#'   \item cat_all_diff_thresh_data = data frame used to make flextable for p-value thresholds
#'   }
#'
#' if comparing categorical variables used
#' \itemize{
#'   \item cat_all_graph = plot of observed to expected numbers and differences between groups, top panels are the absolute numbers, bottom panels are the differences between trial arms in two arm studies
#'   \item cat_all_graph_pc = plot of observed to expected numbers expressed as percentages and differences between groups, top panels are the percentages, bottom panels are the differences between trial arms in two arm studies
#'   \item cat_all_data_abs = data frame of data for absolute numbers
#'   \item cat_all_data_df = data frame of data for difference between groups in two arm studies
#'   \item cat_all_dataset_abs = data frame of dataset used for all trials
#'   \item cat_all_dataset_df = data frame of dataset used for two arm trials
#'   \item cat_all_all_graphs list containing
#'      \itemize{
#'         \item abs = plot for absolute numbers only
#'         \item df = plot for difference between groups in two arm studies only
#'         \item pc = plot for percentages only
#'         \item all_pc = composite plot of percentages and absolute numbers
#'         \item individual_graphs list of 6 individual plots making up composite figures
#'     }
#'   }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate left_join group_by slice arrange count filter rename ungroup across summarise bind_cols
#' @importFrom tidyr spread gather unite separate
#' @importFrom ggpubr ggarrange
#' @importFrom data.table setorder as.data.table
#' @importFrom purrr map map_int
#' @importFrom vcd assocstats
#' @importFrom vcdExtra CMHtest
#' @importFrom epitools riskratio.wald
#' @importFrom flextable flextable font fontsize valign align add_header_row merge_h set_table_properties
#' @importFrom stats runif
#' @export cat_all_fn
#'
#'
#' @examples
#' # load example data
#' cat_all_data <- load_clean(import= "no", file.cat = "SI_cat_all", cat_all= "yes",
#' format.cat = "wide")$cat_all_data
#'
#' \donttest{
#' # run function comparing p-values only (takes only a few seconds)
#' cat_all_fn (comp.pvals = "yes")$cat_all_diff_calc_rep_ft
#'
#' # run function comparing distribution of binomial variables only
#'
#' # to speed example up limit to 12 2-arm trials with 20 variables
#' # (takes close to 5 secs)
#'
#' cat_all_data <- cat_all_data [1:41, c(1:8,10:11,13:15)]
#'
#' cat_all_fn (binom = "yes", two_levels = "yes", del.disparate = "yes",
#' excl.level = "yes", seed = 10)$cat_all_graph
#'
#'
#' # to import an excel spreadsheet (modify using local path,
#' # file and sheet name, range, and format):
#'
#' # get path for example files
#' path <- system.file("extdata", "reappraised_examples.xlsx", package = "reappraised",
#'                    mustWork = TRUE)
#' # delete file name from path
#' path <- sub("/[^/]+$", "", path)
#'
#' # load data
#' cat_all_data <- load_clean(import= "yes", cat_all = "yes", dir = path,
#'    file.name.cat = "reappraised_examples.xlsx", sheet.name.cat = "SI_cat_all",
#'    range.name.cat = "A:N", format.cat = "wide")$cat_all_data}
#'
#' @md


cat_all_fn <- function (df = cat_all_data, comp.pvals = "no", fisher.sim = "y", fish.n.sims =10000, binom = "no", two_levels= "no", del.disparate = "yes",
                        excl.level= "yes", seed = 0, title = "", verbose = TRUE) {

  #import dataset
  a <- df
  results <- list()
  results2 <- list()

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "cat_all") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
               deparse(substitute(df)), "$data_type != 'cat_all'"))}

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  if(vb == "y") {
    pb = txtProgressBar(min=0, max=100, style = 2)
    if (tolower(substr(comp.pvals,1,1)) == "y" & tolower(substr(binom,1,1)) == "y") {inc <- 33
    } else if(tolower(substr(comp.pvals,1,1)) == "y" | tolower(substr(binom,1,1)) == "n") {inc <- 50
    } else if(tolower(substr(comp.pvals,1,1)) == "n" | tolower(substr(binom,1,1)) == "y") {inc <- 50}
  }

  #compare- pvalues
  if (tolower(substr(comp.pvals,1,1)) == "y") {

    a <- a %>% dplyr::select (-data_type)

    if(!"p" %in% colnames(a)) {
      b1 <- "There were no reported p-values"
    } else {

      #calculate p-values
      #convert to wide format for level
      #matrix is levels = row & groups = column = n_levels_group
      #make sure variables are in right order by arrange and factor
      b <- a %>% dplyr::select(study, var, level, starts_with("n", ignore.case = TRUE)) %>%
        tidyr::gather (k, v, starts_with("n", ignore.case = TRUE)) %>%
        dplyr::mutate(k = paste0(gsub("[0-9]+","", k), "_", gsub("\\D", "", k))) %>%
        tidyr::separate(k, into = c("k1", "k2")) %>%
        dplyr::arrange(level, k2) %>%  tidyr::unite(k, k1, level, sep = "_") %>%
        tidyr::unite(k, k, k2, sep = "_") %>%
        dplyr::mutate(k = factor(k, levels = unique(k))) %>% tidyr::spread (k, v)

      b <- dplyr::left_join(b, a %>% dplyr::select(study, var, levels_no ,group) %>%
                              dplyr::group_by(study, var) %>% dplyr::slice (1), by = c("study", "var"))
      b <- b %>% dplyr::select(-starts_with("N", ignore.case = FALSE))

      #Pb1
      if(vb == "y") {setTxtProgressBar(pb, getTxtProgressBar (pb)+ inc)
        cat("\nCalculating p-values ...\n\n")}

      #do as loop
      nm <- which(substr(colnames(b),1,1) == "n")
      nm_p <- c("pc", "pcc", "pf", "plr", "pchm", "p.m", "p.f", "p.c", "pmdpsas")
      b [, nm_p] <- NA

      b$Chisquare.warn <- NA

      #calculate p by different methods
      #matrix is levels = row & groups = column = n_levels_group
      for (i in 1:nrow(b)) {
        x <- na.omit(as.numeric(b [i, nm]))
        x <- matrix(x, byrow = TRUE, nrow = b$levels_no [i])
        chi <- suppressWarnings(chisq.test(x, correct = FALSE))
        b$pc [i] <- chi$p.value
        if(length(chi$expected) == 4 & sum (chi$expected <5) >0 |
              length(chi$expected) > 4 & (sum(chi$expected <1) >0 | sum(chi$expected <5) / length(chi$expected) > 0.2)) {
          b$Chisquare.warn [i] <- "exp <5"}
        b$pcc [i] <- suppressWarnings(chisq.test(x, correct = TRUE)$p.value)
        if (tolower(substr(fisher.sim,1,1)) == "y") {
          b$pf [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8, simulate.p.value = TRUE, B= fish.n.sims)$p.value)
        } else {b$pf [i] <-  suppressWarnings(fisher.test(x, workspace = 2e8)$p.value)}
        b$plr [i] = vcd::assocstats(x)$chisq_tests [1,3]
        b$pchm [i] = vcdExtra::CMHtest(x)$table [1,3]
        if (b$group [i] == 2 & b$levels_no [i] == 2) { #epitools only for 2*2
          pm <- list(suppressWarnings(epitools::riskratio.wald(x)$p.value [2,]))
          b$p.m [i] = pm [[1]] [1]
          b$p.f [i] = pm [[1]] [2]
          b$p.c [i] = pm [[1]] [3]
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
          b$pmdpsas [i] <- round(b$p.f [i] - 0.5* p,9)
        }
      }

      #now merge original p/test back in
      if (!"stat" %in% colnames(a)) {a$stat <- NA}

      a$stat <- tolower(a$stat)

      a$stat1 <- a$stat #keep originals

      if (sum(! a$stat %in% c("chisq", "chisqc", "fisher", "midp", "lr", "mh", NA))>0) {
        if(vb == "y") {cat("\nValues for 'stat' are 'chisq', 'chisqc', 'fisher', 'midp', 'lr', 'mh'\nOther values have been removed\n\n")
        }
        a$stat [!a$stat %in% c("chisq", "chisqc", "fisher", "midp", "lr", "mh")] <- NA
      }

      if (sum(a$levels_no*a$group >4 & a$stat %in% "midp") >0) {
        if(vb == "y") {cat("\nMid-p test not available for variables with group >2 or levels >2\n'midp' for these cases have been removed\n\n")
        }
        a$stat [a$levels_no*a$group >4 & a$stat %in% "midp"] <- NA
      }

      if (!"p" %in% colnames(a)) {a$p <- NA}
      b <- dplyr::left_join(b,
                            a %>% dplyr::select(study, var, p, stat, stat1) %>% dplyr::arrange(study, var, p, stat) %>%
                              dplyr::group_by(study, var) %>% dplyr::slice(1),
                            by = c("study", "var"))

      #get pvalues into usuable form
      b$p_n <- suppressWarnings(as.numeric(b$p))
      b <- b %>% dplyr::mutate(p_g = suppressWarnings(ifelse(grepl(">", p), as.numeric(sub('.*>', '', p)), NA)),
                               p_l = suppressWarnings(ifelse(grepl("<", p), as.numeric(sub('.*<', '', p)), NA))) %>%
        dplyr::mutate(p_g = suppressWarnings(ifelse(grepl("ns", tolower(p)), 0.05, p_g)))
      b$p <- tolower(b$p)
      b$p <- ifelse(is.na(b$p_n),b$p, as.character(b$p_n))

      #add digits
      b$p_d <- ifelse(!is.na(b$p_n),
                      ifelse(abs(b$p_n - round(b$p_n)) > .Machine$double.eps^0.5,
                             nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(b$p_n))))), 0), NA)

      #round all p-values to # digits
      b [, nm_p] <- purrr::map(b [,nm_p], ~ifelse(!is.na(b$p_n), round(.x, b$p_d), round(.x, max(b$p_d,3, na.rm= TRUE))), b$p_n, b$p_d)

      #compare p results
      nm_p <- c("pc", "pf", "pcc", "plr", "pchm", "p.m", "p.f", "p.c", "pmdpsas")
      names (nm_p) <- c("chisq", "fisher", "chisqc", "lr", "mh", "midp","na", "na2", "na3")

      for (i in 1:NROW(b)) {b$p_match [i] <- ifelse(!is.na(nm_p [b$stat [i]]), b [i, nm_p [b$stat [i]]], NA)}
      b$p_match2 <- ifelse(!is.na(b$stat) & b$stat == "midp", b$pmdpsas, NA)

      #thresholds
      b$m_t <- NA
      b$m_t_comment <- NA

      #quite complex logic so do in loop- works through options of test vs no test, and thresholds
      for (i in 1:NROW(b)) {
        if (is.na(b$p_g [i]) & is.na(b$p_l [i])) {b$m_t [i] <- NA; next} #if no thresholds no match

        if (is.na(b$stat [i])) { #no test result
          if (b$group [i] * b$levels_no [i] > 4) {
            if (!is.na(b$p_g [i])) {
              #find if p_calc > p_rep
              #no match then match, greater than
              if (sum (b [i, nm_p] > b$p_g [i], na.rm = TRUE) ==0) {
                b$m_t [i] <- "No match"; b$m_t_comment [i] <- "No reported test, >2*2"
              } else {
                b$m_t [i] <- "Match"; b$m_t_comment [i] <- paste0("No reported test, >2*2, test=",
                                                names(nm_p [which.max (b [i, nm_p] > b$p_g [i])]))
              }
          }

          if (!is.na(b$p_l [i])) {
            #lower than
            if (sum (b [i, nm_p] < b$p_l [i], na.rm = TRUE) ==0) {
              b$m_t [i] <- "No match"; b$m_t_comment [i] <- "No reported test, >2*2"
            } else {
              b$m_t [i] <- "Match"; b$m_t_comment [i] <- paste0("No reported test, >2*2, test=",
                                                                names(nm_p [which.max (b [i, nm_p] < b$p_l [i])]))
            }
          }
        } #end of >2*2

          if (b$group [i] * b$levels_no [i] <= 4) {
            if (!is.na(b$p_g [i])) {
              if (sum (b [i, nm_p] > b$p_g [i], na.rm = TRUE) ==0) {
                b$m_t [i] <- "No match"; b$m_t_comment [i] <- "No reported test, 2*2"
              } else {
                b$m_t [i] <- "Match"; b$m_t_comment [i] <- paste0("No reported test, 2*2, test=",
                                                                  names(nm_p [which.max (b [i, nm_p] > b$p_g [i])]))
              }
            }
            if (!is.na(b$p_l [i])) {
              if (sum (b [i, nm_p] < b$p_l [i], na.rm = TRUE) ==0) {
                b$m_t [i] <- "No match"; b$m_t_comment [i] <- "No reported test, 2*2"
              } else {
                b$m_t [i] <- "Match"; b$m_t_comment [i] <- paste0("No reported test, 2*2, test=",
                                                                  names(nm_p [which.max (b [i, nm_p] < b$p_l [i])]))
              }
            }
          }#end of 2*2
        } #end of no test result

        if (!is.na(b$stat [i])) {  #test result
          if (b$stat [i] == "midp") {
            if (!is.na(b$p_g [i])) {
              b$m_t [i] <- ifelse(b$p_match [i] > b$p_g [i] & b$p_match2 [i] > b$p_g [i], "Match", "No match")} #> and midp
            if (!is.na(b$p_l [i])) {
              b$m_t [i] <- ifelse(b$p_match [i] < b$p_l [i] & b$p_match2 [i] < b$p_l [i], "Match", "No match")} #< and >2*2
            b$m_t_comment [i] <- "Midp SAS and R"}

          if (b$stat [i] != "midp") {
            if (!is.na(b$p_g [i])) {
              b$m_t [i] <- ifelse(b$p_match [i] > b$p_g [i], "Match", "No match")} #> and other
            if (!is.na(b$p_l [i])) {
              b$m_t [i] <- ifelse(b$p_match [i] < b$p_l [i], "Match", "No match")} #< and 2*2
            b$m_t_comment [i] <- b$stat [i]}
        } #end of test result
      }
      b$m_t_comment <- ifelse(!is.na(b$m_t_comment),
                              paste0(toupper(substr(b$m_t_comment, 1, 1)), substr(b$m_t_comment, 2, nchar(b$m_t_comment))), NA)


      b$p_df_stat <- b$p_df <- NA
      nm_ph <- c("pc","pf","p.m","pmdpsas","pcc","plr","pchm","p.f","p.c")
      #difference between reported values
      for (i in 1:NROW(b)) {
        if (is.na(b$p_n [i])) {next}
        if (is.na(b$stat [i])) {#no stat result
          #are there any matches df =0 or lowest dfs
          b$p_df [i] <- ifelse(!is.na(b$p_n [i]), min(abs(b$p_n [i] - b [i, nm_p]), na.rm = TRUE), NA)
          b$p_df_stat [i]  <-  names(which(nm_p == nm_ph [which.max(b$p_n [i] == b [i, nm_ph] + b$p_df [i])]))
        }
        if (!is.na(b$stat[i])) {
          b$p_df_stat[i] <- b$stat[i]
          if (!b$stat [i] == "midp") {b$p_df[i] <- abs(b$p_n[i] - b[i, nm_p[b$stat[i]]])
          } else {
             if (abs(b$p_n [i] - b$p.m [i]) < abs(b$p_n [i] - b$pmdpsas [i])) {
               b$p_df_stat [i] <- "midp.epitools"
               b$p_df [i] <-abs(b$p_n [i] - b$p.m [i])
             } else {
               b$p_df_stat [i] <-  "midp.SAS"
               b$p_df [i] <- abs(b$p_n [i] - b$pmdpsas [i])
             }
          }
        }
      }
      b$p_df_stat <- ifelse(!is.na(b$p_df_stat),
                            paste0(toupper(substr(b$p_df_stat, 1, 1)), substr(b$p_df_stat, 2, nchar(b$p_df_stat))), NA)

      #change midp titles
      b$p_df_stat <- ifelse(b$p_df_stat %in% "Na3", "Midp.sas",
                            ifelse(b$p_df_stat %in% "Midp", "Midp.epitools", b$p_df_stat))

      #output
      #compare reported vs unreported p-values
      b$df.groupsp <- cut(b$p_n,c(-1e-99,1e-99,seq(0.1,1.0,0.1)))
      b$df.groupsdf <- cut(b$p_df,c(-1e-99,1e-99,seq(0.1,1,0.1)))
      nm <- data.frame(diff_in_p_value= c(0, paste(seq(0,.9,0.1),"-", seq(0.1,1,0.1), "")), name = levels(b$df.groupsdf),
                       stringsAsFactors = FALSE)

      diff <- dplyr::left_join(nm %>% dplyr::select(diff_in_p_value),
                               dplyr::left_join (b %>% dplyr::filter(!is.na(p_n)) %>% dplyr::count(df.groupsdf) %>%
                                                   dplyr::rename(name =df.groupsdf), nm, by = "name") %>%
                                 dplyr::select (-name),
                               by = "diff_in_p_value")
      d1 <- ncol(diff)
      if (sum(!is.na(b$stat))>0) {
        diff <- dplyr::left_join(diff,
                                 dplyr::left_join (b %>% dplyr::filter(!is.na(p_n)) %>%
                                                     dplyr::mutate(stat1 = ifelse(is.na(stat1), "No test reported", stat1)) %>%
                                                     dplyr::count(stat1, df.groupsdf) %>% dplyr::rename(name =df.groupsdf),
                                                   nm, by = "name") %>% dplyr::select (-name) %>%
                                   dplyr::mutate (stat1 = paste0(toupper(substr(stat1, 1, 1)), substr(stat1, 2, nchar(stat1)))) %>%
                                   tidyr::spread (stat1, n),
                                 by = "diff_in_p_value")}
      d2 <- ncol(diff) - d1
      diff <- dplyr::left_join(diff,
                               dplyr::left_join(
                                 dplyr::left_join (b %>% dplyr::filter(!is.na(p_n)) %>% dplyr::count(df.groupsp, df.groupsdf) %>%
                                                     dplyr::rename(name =df.groupsdf), nm, by = "name") %>%
                                   dplyr::select (-name) %>% dplyr::rename(name = df.groupsp),
                                 nm %>% dplyr::rename (bpv = diff_in_p_value), by = "name") %>% dplyr::select (-name) %>%
                                 tidyr::spread (bpv, n),
                               by = "diff_in_p_value")
      d3 <- ncol(diff) - d1 -d2
      diff[is.na(diff)] <- ""

      diff <- diff %>% dplyr::mutate(n = ifelse(n=="", "",
                                                paste0 (as.numeric(n), " (",
                                                        round(as.numeric(n)/sum(as.numeric(n), na.rm = TRUE)*100,0),")")))
      colnames (diff) [1:2] <- c("Difference between \nreported and\ncalculated p-values", "n (%)")

      f <- flextable::flextable(diff)
      f <- flextable::font(f, fontname = "Arial")
      f <- flextable::fontsize(f, size = 8)
      f <- flextable::valign(f, i = 1, j = 2:ncol(diff), valign = "bottom", part = "head")
      f <- flextable::align(f, j = 2:ncol(diff), align = "center", part = "all")
      h <- c(rep("", d1), rep ("Test",d2), rep ("Baseline p-value", d3))
      f <- flextable::add_header_row(f, values = h, top = TRUE)
      f <- flextable::merge_h(f, part = "header")
      f <- flextable::align(f, i = 1, align = "center", part = "header")
      f <- flextable::fontsize(f, part = "header", size = 8)
      f <- flextable::set_table_properties(f, width = 1, layout = "autofit")

      #compare thresholds
      if (sum (!is.na(b$m_t)) == 0) {diff1 <- f1 <- "There were no threshold p-values reported"
      } else {
        b$tp <- factor(ifelse(rowSums(!is.na(cbind(b$p_g,b$p_l)), na.rm = TRUE) >0, b$p, NA))
        diff1 <- dplyr::left_join(b %>% dplyr::filter(!is.na(tp)) %>% dplyr::count(tp),
                                  b %>% dplyr::filter(!is.na(tp)) %>% dplyr::count(tp, m_t) %>% tidyr::spread (m_t, n),
                                  by = "tp")
        diff1[is.na(diff1)] <- ""
        if (! "Match" %in% colnames(diff1)) {diff1$Match <- 0}
        if (! "No match" %in% colnames(diff1)) {diff1$`No match` <- 0}
        diff1 <- diff1 %>%
          dplyr::mutate(Match = ifelse(Match=="", "",
                                       paste0 (as.numeric(Match), " (", round(as.numeric(Match)/n*100,0),")")),
                        `No match` = ifelse(`No match`=="", "", paste0 (as.numeric(`No match`), " (",
                                                                        round(as.numeric(`No match`)/n*100,0),")")),
                        n = ifelse(n=="", "",
                                   paste0 (as.numeric(n), " (",
                                           round(as.numeric(n)/sum(as.numeric(n), na.rm = TRUE)*100,0),")")))
        colnames (diff1) [1:2] <- c("Threshold\np-value", "n (%)")

        f1 <- flextable::flextable(diff1)
        f1 <- flextable::font(f1, fontname = "Arial")
        f1 <- flextable::fontsize(f1, size = 8)
        f1 <- flextable::valign(f1, i = 1, j = 2:ncol(diff1), valign = "bottom", part = "head")
        f1 <- flextable::align(f1, j = 2:ncol(diff1), align = "center", part = "all")
        f1 <- flextable::add_header_row(f1, values = c("","", "Thresholds vs calculated p value", "Thresholds vs calculated p value"), top = TRUE)
        f1 <- flextable::merge_h(f1, part = "header")
        f1 <- flextable::align(f1, i = 1, align = "center", part = "header")
        f1 <- flextable::fontsize(f1, part = "header", size = 8)
        f1 <- flextable::set_table_properties(f1, width = 0.4, layout = "autofit")
      }

      #return data
      b$warn <- ifelse(b$p_df_stat %in% "Chisq", b$Chisquare.warn, NA)
      b1 <- b %>% dplyr::select(study,var, p, stat1, p_df, p_df_stat, warn, m_t, m_t_comment, pc, Chisquare.warn, pf, p.m, pcc, plr, pchm, pmdpsas, group, levels_no) %>%
        dplyr::ungroup () %>% as.data.frame ()
      colnames (b1) <-  c("study", "variable", "reported_p_value", "test_reported", "difference_reported_calculated_p_values",
                          "test_calculated", "warn.message", "reported_p_matches_threshold", "threshold_comment", "chisquare", "chisquare.warning","fisher",
                          "midp_r_epitools","chisquare_continuty", "likelihood_ratio", "mantel_haenszel", "midp_SAS_calculation",
                          "groups", "levels_of_variable")

      results <- list(cat_all_pvals = b1, cat_all_diff_calc_rep_ft =f, cat_all_diff_calc_rep_data = diff,
                      cat_all_diff_thresh_ft = f1, cat_all_diff_thresh_data = diff1)
      if (vb == "y") {print (f)}
    }

    if(all(b1 == "There were no reported p-values")) {results <- list(cat_all_pvals = b1, cat_all_diff_calc_rep_ft =b1,
                                                                      cat_all_diff_calc_rep_data = b1, cat_all_diff_thresh_ft = b1, cat_all_diff_thresh_data = b1)}
  }


  #binomial function # note this is essentially cat_fn but uses n and N not v and n
  if (tolower(substr(binom,1,1)) == "y") {

    #collapse to 2 levels and delete one level
    if (tolower(substr(two_levels,1,1)) == "n") {

      #choose a random level to remove
      if(seed ==0) {set.seed(Sys.time())
      } else {set.seed(seed)}

      df <- df %>% dplyr::mutate(r1 = runif(nrow(df))) %>% dplyr::group_by (study, var) %>% dplyr::arrange (study, var, r1) %>%
        dplyr::mutate (l1 = dplyr::row_number()) %>% dplyr::select(-r1) %>% dplyr::ungroup () %>% as.data.frame()

      if (tolower(substr(excl.level,1,1)) == "y") {
        a <- df %>% dplyr::filter (!l1 == 1) %>% dplyr::select (study,var, recode, group, starts_with("n")) #exclude 1 level
      } else {a <- df %>% dplyr::select (study,var, recode, group, starts_with("n"))} # don't exclude

    } else {

      a <- df %>% dplyr::select(-level, -levels_no)
      if ("p" %in% colnames(a)) {a <- a %>% dplyr::select(-p)}

      nm1 <- colnames(a) [substr(colnames(a),1,1) == "n"]
      nm2 <- colnames(a) [substr(colnames(a),1,1) == "N"]

      a <- dplyr::left_join (a %>% dplyr::group_by(study, var, recode) %>%
                               dplyr::summarise (dplyr::across(c((all_of(nm1))), ~ sum(.x, na.rm = FALSE)), .groups = "keep"),
                             a %>% dplyr::select(study, var, group, all_of(nm2)) %>% dplyr::group_by (study, var) %>% dplyr::slice(1),
                             by = c("study", "var")) %>%
        dplyr::select(study:recode, group, all_of(nm1), all_of(nm2)) %>% dplyr::ungroup() %>% as.data.frame()

      if (tolower(substr(excl.level,1,1)) == "y") {
        if(seed ==0) {set.seed(Sys.time())
        } else {set.seed(seed)}

        a <- a %>% dplyr::mutate(r1 = runif(nrow(a))) %>% dplyr::group_by (study, var) %>% dplyr::arrange (study, var, r1) %>%
          dplyr::mutate (recode1 = dplyr::row_number()) %>% dplyr::select(-r1) %>% dplyr::ungroup () %>% as.data.frame() %>%
          filter(recode1 ==1) %>% select(-recode1)
      }
    }

    #Pb1/2
    if(vb == "y") {setTxtProgressBar(pb, getTxtProgressBar (pb)+ inc)
      cat("\nCalculating data for raw numbers and making graphs ...\n\n")}

    g <- max(a$group)#largest number of groups
    nm.n <- colnames(a) [substr(colnames(a),1,1) == "n"]
    nm.N <- colnames(a) [substr(colnames(a),1,1) == "N"]

    #calculate N, n and p # p = prop of people with trait
    a$p <- rowSums(a[, nm.n], na.rm = TRUE) / rowSums(a [, nm.N], na.rm = TRUE)

    #get largest group and add 1 because there can be 0 people with trait
    n <- max(a[, nm.N], na.rm = TRUE)+1

    #starting row for each group
    rw_g <- list()
    for (i in 1:g) {rw_g [[i]] <-  (i-1)*n +1}

    #get expected p
    exp_p <- expect_p(a, n, g, nm.N, rw_g)

    #sum the expected probs
    sum_exp <- 0
    for (i in 1:g) {sum_exp <- sum_exp + exp_p [rw_g [[i]]: (rw_g [[i]]+n-1),]}

    #sum the observed
    obs <- a %>% dplyr::select(all_of(nm.n)) %>% tidyr::gather (k, v) %>% dplyr::group_by (v) %>% dplyr::count () %>%
      dplyr::filter (!is.na(v)) %>% dplyr::rename (obs = n) %>% dplyr::ungroup()

    b <- dplyr::left_join(data.frame(num= 0: (n-1), exp= rowSums(sum_exp)), obs, by = c("num" = "v")) %>% replace(is.na(.), 0)

    #if large number of groups then compress
    n_gp <- sum(!is.na(b$obs))
    for (sz in seq(5,ceiling(n_gp/5)*5,5)) {
      b_temp <- compress(b, "obs", "exp", size=sz)
      if (sz == 5) {
        b <- dplyr::bind_cols(b, b_temp)
        if (sum(!is.na(b_temp [ , "obsm__5"])) < 22) {break}
      }
      if (sum(!is.na(b_temp [ , paste0("obsm__",sz)])) < 22) {
        b <- dplyr::bind_cols(b, b_temp)
        break}
    }

    graph <- list()
    graph [1:2] <- cat_graph (gph= b,  xtitle = "Number with characteristic per group", ytitle = "Number of trial groups*variables",
                              size = sz, sfx = "", text = list(length(unique(a$study)), paste(sum(b$obs), "groups*variables")),
                              fn = "cat_all", ti = "Y", top = "Y", t= title)

    rr_g <- list()
    rr_g [1:2] <- rr_graph(rr = b, xtitle =  "Number with characteristic per group", ytitle = "Observed/Expected ratio",
                           size =sz, sfx= "", text = list(length(unique(a$study)), sum(b$obs)),
                           fn = "cat_all", ti = "Y", top = "yes", t=title)

    #percent
    colnames(exp_p) <- 1:ncol(exp_p)
    exp_pc <- dplyr::bind_cols(num = rep(0:(n-1),g), g= rep(1:g, each= n), exp_p)
    exn <- as.data.frame(t(a [, nm.N]))
    exn_exp <- dplyr::bind_cols(num = rep(0:(n-1),g), g= rep(1:g, each= n), as.data.frame(lapply(exn, rep, rep(n,g))))
    colnames(exn_exp) <- c("num","g", 1:ncol(exn))
    ex1 <- dplyr::left_join(exp_pc %>% tidyr::gather(k, v, -num, -g) %>% dplyr::mutate(k = as.numeric(k)),
                            exn_exp %>% tidyr::gather(k, v1, -num, -g) %>% dplyr::mutate(k = as.numeric(k)), by = c("num", "g", "k"))
    ex1 <- ex1 %>% dplyr::filter(!is.na(v1)) %>% dplyr:: mutate(num1 = round(num/v1*100,0)) %>% dplyr::filter(num1 <101) %>%
      dplyr::group_by(k,g, num1) %>% dplyr::summarise (v = sum(v), .groups= "drop") %>% tidyr::spread(k, v)
    ex2 <- dplyr::left_join(data.frame(num1 = rep(0:100, g), g= rep(1:g, each=101)), ex1, by = c("num1","g")) %>%
      replace(is.na(.), 0)

    sum_exp_pc <- 0
    for (i in 1:g) {sum_exp_pc <- sum_exp_pc + ex2 [ex2$g == i, 3:ncol(ex2)]}

    obs_pc <- dplyr::bind_cols(a %>% dplyr::select(all_of(nm.n)) %>% tidyr::gather (k, v), a %>%
                                 dplyr::select(all_of(nm.N)) %>% tidyr::gather (k, v1) %>% dplyr::select(v1)) %>%
      dplyr::mutate(v = round(v/v1*100,0)) %>% dplyr::group_by (v) %>% dplyr::count () %>% dplyr::filter(!is.na(v)) %>%
      dplyr::rename(obs = n) %>% dplyr::ungroup ()

    b_pc= dplyr::left_join(data.frame(num= 0:100, exp = rowSums(sum_exp_pc)), obs_pc, by = c("num"= "v")) %>%
      replace(is.na(.), 0)

    n_gp <- sum(!is.na(b_pc$obs))
    for (sz in seq(5,ceiling(n_gp/5)*5,5)) {
      b_temp <- compress(b_pc, "obs", "exp", size=sz)
      if (sz == 5) {
        b_pc <- dplyr::bind_cols(b_pc, b_temp)
        if (sum(!is.na(b_temp [ , "obsm__5"])) < 22) {break}
      }
      if (sum(!is.na(b_temp [ , paste0("obsm__",sz)])) < 22) {
        b_pc <- dplyr::bind_cols(b_pc, b_temp)
        break}
    }

    graph [5:6] <- cat_graph (gph= b_pc,  xtitle = "Percentage with characteristic per group", ytitle = "Number of trial groups*variables",
                              size = sz, sfx = "", text = list(length(unique(a$study)), paste(sum(b_pc$obs), "groups*variables")),
                              fn = "cat_all", ti = "Y", top = "Y", t= title)

    rr_g [5:6] <- rr_graph(rr = b_pc, xtitle =  "Percentage with characteristic per group", ytitle = "Observed/Expected ratio",
                           size =sz, sfx= "", text = list(length(unique(a$study)), sum(b_pc$obs)),
                           fn = "cat_all", ti = "Y", top = "yes", t=title)

    #Pb2/3
    if(vb == "y") {setTxtProgressBar(pb, getTxtProgressBar (pb)+ inc)
      cat("\nCalculating data for differences betweeen groups and making graphs ...\n")
      cat("\nThis usually takes the most time ...\n\n")}

    #repeat for 2 arm studies only;
    a.2 <- a %>% dplyr::filter (group == 2) %>% dplyr::select(study, var, all_of(nm.n) [1:2], all_of(nm.N) [1:2])
    exp_p2 <- exp_p [, a$group == 2]

    #delete disparate rows
    if (tolower(substr(del.disparate,1,1)) == "y") {
      dd_f <- with(a.2, !(N1/N2 < 0.8 | N1/N2 >1.2 | N2/N1 < 0.8 | N2/N1 > 1.2))
      a.2 <- a.2 %>% dplyr::filter(dd_f)
      exp_p2 <- exp_p2 [, dd_f, drop= FALSE]
      rm (dd_f)}

    #expected differences
    #calculate probability of group numbers differing by 0, 1, 2 etc# this usually takes the longest amount of time.
    exp_p3 <- expect_p_diff (a.2, exp_p2, n, nm.N)

    #count observed diffs
    a.2$obs = abs(a.2 [, nm.n [[1]]] - a.2[ , nm.n [[2]]])
    obs <- a.2 %>% dplyr::count(obs) %>% dplyr::rename (num= obs, obs= n)
    b_df <- dplyr::left_join(data.frame(num= 0: (n-1), exp= rowSums(exp_p3)), obs, by = "num") %>% replace(is.na(.), 0)

    n_gp_df <- sum(!is.na(b_df$obs))
    for (sz_df in seq(5,ceiling(n_gp_df/5)*5,5)) {
      b_df_temp <- compress(b_df, "obs", "exp", size=sz_df)
      if (sz_df == 5) {
        b_df <- dplyr::bind_cols (b_df, b_df_temp)
        if (sum(!is.na(b_df_temp [ , "obsm__5"])) < 22) {break}
      }
      if (sum(!is.na(b_df_temp [ , paste0("obsm__",sz_df)])) < 22) {
        b_df <- dplyr::bind_cols (b_df, b_df_temp)
        break}
    }

    graph [3:4] <-  cat_graph (gph= b_df,  xtitle = "Differences between trial groups",
                               ytitle = "Number of variables", size = sz_df, sfx = "",
                               text = list(length(unique(a.2$study)), paste(sum(b_df$obs), "variables")),
                               fn = "cat_all", ti = "n", top = "n", t= title)

    rr_g [3:4] <- rr_graph(rr = b_df ,xtitle = "Differences between trial groups",
                           ytitle = "Observed/Expected ratio", size =sz_df, sfx= "",
                           text = list(length(unique(a.2$study)), sum(b_df$obs)),
                           fn = "cat_all", ti = "n", top = "n", t= title)

    graphs <- ggpubr::ggarrange(graph [[2]], graph [[4]], ncol = 1, nrow= 3, align = "hv")
    rr_graphs <- ggpubr::ggarrange(rr_g [[2]], rr_g [[4]], ncol = 1, nrow= 2, align = "hv")

    both_graphs <- list ()
    both_graphs [[1]] <- ggarrange(graph [[2]], rr_g [[2]], graph [[4]], rr_g [[4]],
                                   ncol = 2, nrow= 3, align = "hv")
    both_graphs [[2]] <- ggarrange(graph [[1]], rr_g [[2]], ncol = 2, nrow= 1, align = "hv")
    both_graphs [[3]] <- ggarrange(graph [[3]], rr_g [[4]], ncol = 2, nrow= 1, align = "hv")
    both_graphs [[4]] <- ggarrange(graph [[5]], rr_g [[6]], ncol = 2, nrow= 1, align = "hv")
    both_graphs [[5]] <- ggarrange(graph [[1]], rr_g [[2]], graph [[5]], rr_g [[6]], ncol = 2, nrow= 2, align = "hv")
    both_graphs [[6]] <- ggarrange(graph [[5]], rr_g [[6]], graph [[4]], rr_g [[4]], ncol = 2, nrow= 2, align = "hv")


    results2 <- list(cat_all_graph= both_graphs [[1]], cat_all_graph_pc= both_graphs [[6]],
                     cat_all_data_abs = b, cat_all_data_df = b_df, cat_all_dataset_abs = a, cat_all_dataset_df = a.2,
                     cat_all_all_graphs =  list(abs = both_graphs [[2]], df = both_graphs [[3]], pc = both_graphs [[4]],
                                                all_pc = both_graphs [[5]], ind_graphs = c(graph[c(1,3,5)], rr_g [c(1,3,5)])))

   if (vb == "y") {print(both_graphs [[1]])}
  }


results3 <- c(results, results2)

if (vb == "y") {close(pb)}

return (results3)

}
