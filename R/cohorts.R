#' Compares proportions of matching summary statistics in different cohorts
#'
#' Creates flextable of probability of matching mean, SD, and mean and SD for each variable in different cohorts in the
#' specified number of simulations\cr
#'
#' Reference data is from Bolland 2021\cr
#' Bolland MJ, Gamble GD, Avenell A, Grey A. Identical summary statistics were uncommon in randomized trials and cohort studies. J Clin Epidemiol 2021;136:180-188.

#'
#' Returns a list containing 6 objects and (if verbose = TRUE) prints the flextable cohort_ft
#'
#' @param df data frame generated from load_clean function
#' @param seed the seed to use for random number generation, default 0 = current date and time. Specify seed to make repeatable.
#' @param sims number of simulations, default -1 = function selects based on number of variables and sample size.
#' @param n_vars restrict analyses to variables in at least (>=) this number of cohorts, default = 10 (ie variable has mean in 10 or more cohorts).
#' @param popn if dataset contains studies in different sub-populations, code this in cohort_data$population and studies are subsetted if match in this variable. 'All' overrides this and uses all data regardless of information in this variable.
#' @param title title name for plots (optional)
#' @param verbose TRUE or FALSE indicates whether progress bar and comments show and flextable is printed
#'
#' @return list containing 6 objects as described
#'\itemize{
#'   \item cohort_ft = flextable of results
#'   \item cohort_graph = plot of observed to expected numbers of matches per cohort for mean; SD; and mean and SD
#'   \item all_graphs = list containing
#'      \itemize{
#'         \item all_graphs = all plots on single plot
#'         \item both_graphs = list of 3 plots row by row used to form all_graphs
#'         \item individual_graphs= list of 6 individual plots used to form all_graphs
#'         }
#'   \item cohort_cohort_data = data frame used to generate results data
#'   \item cohort_prob_data = data frame used to make flextable
#'   \item cohort_oe_data= data frame used to make observed to expected plots
#'   }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select group_by ungroup filter mutate_all arrange summarise full_join left_join count mutate_at rename bind_rows bind_cols
#' @importFrom tidyr gather spread separate
#' @importFrom officer fp_border
#' @importFrom purrr pmap_dfr pmap map2_dbl modify pmap_dfc map2_dfc map
#' @importFrom flextable set_header_labels add_header_row merge_at align border_remove hline_bottom hline footnote font fontsize set_table_properties as_paragraph
#' @importFrom ggpubr ggarrange
#' @importFrom data.table setDF setDT setorderv as.data.table rbindlist
#' @export cohort_fn
#'
#'
#' @examples
#' # load example data
#' cohort_data <- load_clean(import= "no", file.cont = "SI_cohort", cohort= "yes",
#' format.cont = "long")$cohort_data
#'
#' \donttest{
#' # run function (takes close to 5 seconds)
#' cohort_fn(seed=10, sims = 100)$cohort_ft
#'
#' # to import an excel spreadsheet (modify using local path,
#' # file and sheet name, range, and format):
#'
#' # get path for example files
#' path <- system.file("extdata", "reappraised_examples.xlsx", package = "reappraised",
#'                      mustWork = TRUE)
#' # delete file name from path
#' path <- sub("/[^/]+$", "", path)
#'
#' # load data
#' cohort_data <- load_clean(import= "yes", cohort = "yes", dir = path,
#'      file.name.cont = "reappraised_examples.xlsx", sheet.name.cont = "SI_cohort",
#'      range.name.cont = "A1:F101", format.cont = "long")$cohort_data}
#'
#' @md

cohort_fn <- function (df= cohort_data, seed= 0 , sims = -1, n_vars = 10, popn = "", title= "", verbose = TRUE) {

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  if(vb == "y") {pb = txtProgressBar(min=0, max=100, style = 2)}

  #Pb1
  if(vb == "y") {setTxtProgressBar(pb,10)
  cat("\nCounting observed matches ...\n\n")}

  #import dataset
  a <- df

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "cohort") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
           deparse(substitute(df)), "$data_type != 'cohort'"))}

  #filter by population if present
  if (tolower(popn) == "all") {
    popn <- ""
    a$population <- NULL}

  if (popn != "") {
    if (!"population" %in% colnames (a)) {
      stop_mb ("variable 'population' is missing from dataframe")}
    if (! popn %in% a$population) {
      stop_mb (paste0(popn, "is missing from the variable 'population'"))}
    a <- a %>% dplyr::filter (population == popn)
    cat(paste0("\nSelecting population as ",popn,"\n"))
  }

  if("population" %in% colnames (a)){
    if (length(unique(a$population)) >1) {
      stop_mb (paste0("variable 'population' has more than one level but has not been selected ","\n",
                      "delete population to run without this variable"))}
  }

  a$population <- ifelse("population" %in% colnames (a),a$population, "")

  #count decimal places
  nm <- colnames (a) [substring(colnames (a), 1, 1) == "m" |
                        (substring(colnames (a), 1, 1) == "s" & substring(colnames (a),2,2) %in% as.character(0:9))]

  dp <- function (x) {
    ifelse(abs(x- round(x)) > .Machine$double.eps^0.5,
           nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(x))))), 0)}

  a <- cbind(a, purrr::map(a %>% dplyr::select(all_of(nm)) %>% setNames (paste0("dp",nm)), dp))

  #r and excel don't include final 0s (as numbers) so take maximum dp (assuming that if m1=x and m2 = x.1, m1 is actually x.0)
  nm.m <- paste0("dp", nm [substr(nm, 1,1) == "m"])
  nm.s <- paste0("dp", nm [substr(nm, 1,1) == "s"])

  a[ , nm.m] <- purrr::pmap_df(a %>% dplyr::select(all_of(nm.m)), ~{x = c(...)
  ifelse(is.na(x), NA, max(x, na.rm = TRUE)) })
  a[ , nm.s] <- purrr::pmap_df(a %>% dplyr::select(all_of(nm.s)), ~{x = c(...)
  ifelse(is.na(x), NA, max(x, na.rm = TRUE)) })

  #make into long format
  b <- a %>% dplyr::select(-data_type) %>%  tidyr::gather(k, value, -c("study","var","population")) %>%
    dplyr::mutate(k = paste0(gsub("[0-9]+","", k), "_", gsub("\\D", "", k))) %>%
    tidyr::separate(k, into = c("v", "group"), sep="_") %>%
    tidyr::spread(v, value) %>% dplyr::filter (!is.na(n))

  #nvar = number of vars, min_dp minimum dps per variable, nmax = study size,
  b <- suppressWarnings( #add this because if NA for SD get warning for infinity)
    b %>% dplyr::group_by(var) %>% dplyr::mutate(min_dpm = min (dpm, na.rm = TRUE), max_dpm = max(dpm, na.rm = TRUE),
                                                 min_dps = min(dps, na.rm = TRUE), max_dps = max(dps, na.rm= TRUE), n_var = dplyr::n()) %>%
      dplyr::group_by(study, group) %>% dplyr::mutate (n_max = max(n)) %>% dplyr::ungroup () %>%
      dplyr::mutate(ms = paste(b$m,b$s,sep="_")) %>%
      dplyr::filter (n_var >= n_vars) %>% # remove if number of variables is less than nvars (default 10)
      dplyr::mutate_all(., .funs = function(x){ifelse(is.infinite(x),NA,x)})) #remove any infinites

  if (NROW(b) == 0) {
    stop_mb (paste0("After excluding variables wih less than ",n_vars," observations, there were no observations"))}

  #observed values
  #identify matches
  b <- b %>% dplyr::group_by(var, m) %>% dplyr::mutate(mn_m = sum(!is.na(m))) %>% dplyr::group_by(var, s) %>%
    dplyr::mutate(mn_s = sum(!is.na(s))) %>% dplyr::group_by(var, ms) %>%
    dplyr::mutate(mn_ms = sum(!is.na(m)&!is.na(s))) %>% dplyr::ungroup ()

  #count number of matches for each variable
  match_count <- function (df1,g, v, v_m) {
    x <- {{df1}} %>% dplyr::group_by({{g}}) %>% dplyr::filter (!is.na({{v}})|!grepl("NA",{{v}})) %>%
      dplyr::select({{g}}, {{v_m}}) %>% dplyr::arrange({{g}}, desc({{v_m}})) %>%
      dplyr::filter ({{v_m}} > 1) %>% dplyr::summarise ("{{v_m}}" := dplyr::n())
    return(x)}

  b_all <- data.frame("var"=unique(b$var)) %>%
    dplyr::full_join(match_count(b, var, m, mn_m), by = "var") %>%
    dplyr::full_join(match_count(b, var, s, mn_s), by = "var") %>%
    dplyr::full_join(match_count(b, var, ms, mn_ms), by = "var") %>% replace(is.na(.), 0)

  #make table of matches by cohort
  obs <- data.frame(freq = 0:(max(b_all %>% dplyr::select(-var)))) %>%
    dplyr::left_join(b_all %>% tidyr::gather(k, freq, -var) %>% dplyr::count(k,freq) %>%
                       tidyr::spread(k, n), by= "freq") %>%
    dplyr::select(freq:mn_m, mn_s, mn_ms) %>% replace(is.na(.), 0) %>% setNames(c("freq", "ob_m", "ob_s","ob_ms"))

  #expected values
  #get modes
  mode <- function(x,n) {
    n <- n [!is.na(x)]
    x <- x[!is.na(x)]
    ux <- unique(x)
    if(!anyDuplicated(x)){
      m <- sum(x*n)/sum(n)
    } else {
      tbl <-   tabulate(match(x, ux))
      m <- mean(ux[tbl==max(tbl)])
    }
    return (m)
  }

  #for m_s variable
  modems <- function(x) {
    x <- x[!grepl("NA",x)]
    ux <- unique(x)
    if(!anyDuplicated(x)){
      m <- NA
    } else {
      tbl <-   tabulate(match(x, ux))
      m <- ifelse(sum(tbl==max(tbl))>1, NA, ux[tbl == max(tbl)])
    }
    return (m)
  }

  #calculate mode, or if no mode, weighted mean
  c <- b %>% dplyr::group_by(var) %>%
    dplyr::summarise(modem = mode(m, n), modes= mode(s,n), modems = modems(ms), dm = max(max_dpm), ds = max(max_dps),
                     m_wt = sum(m*n, na.rm = TRUE)/sum(n, na.rm = TRUE),
                     s_wt = sum(s*n, na.rm = TRUE)/sum(n, na.rm = TRUE))  %>%
    tidyr::separate(modems,into = c("modemn", "modesd"),sep = "_", remove = FALSE) %>%
    dplyr::mutate_at(c('modemn','modesd'),as.numeric) %>%
    dplyr::mutate(m_wt = round(m_wt, dm), s_wt = round(s_wt, ds)) %>%
    dplyr::mutate (modemn = ifelse(is.na(modemn), round(m_wt,dm), modemn)) %>% #if no mode take weighted mean
    dplyr::mutate (modesd = ifelse(is.na(modesd), round(s_wt,ds), modesd))

  #suggested sims
  n_v <- NROW(b)
  n_max <- max(b$n)
  var <- c(10000, 4000, 2000, 1000, 500, 200, 100)
  names (var) <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  sm <- names (var) [which.min(abs(var - n_v))] %>% as.numeric()
  if (sims ==-1) {
    sims <- sm
    if (n_max >500) {sm <- sims <- sims * .75}
  }

  #Pb2
  if(vb == "y") {
    setTxtProgressBar(pb,20)
    cat(paste0("\nSimulations: ",sims,", for ",n_v," variables\n\n"))
    cat("Doing simulations...\nthis takes some time, usually at least 30-60s\n")
    cat("if cohorts have large n (eg n>500), this adds more time\n")
    cat(paste0("to speed up specify fewer simulations in the function call eg cohort_fn (..., sims =",
               round(sims/2, 0)," ...)\n\n"))
  }

  sim <- dplyr::left_join(b %>% dplyr::arrange(study, var, group) %>%
                            dplyr::select(study, var, group, n_max, dpm, max_dpm, dps, max_dps),
                          c %>% dplyr::select(var, modemn, modesd), by = "var") %>%
    dplyr::rename(n = n_max, m= modemn, s= modesd)

  #do the simulations
  if(seed ==0) {set.seed(Sys.time())
  } else {set.seed(seed)}
  sim$row <- 1:n_v
  #sim_exp <- dplyr::bind_rows(replicate(sims, sim, simplify = FALSE))
  sim_exp <- data.table::rbindlist(replicate(sims, sim, simplify = FALSE)) #slightly faster
  sim_exp$sim <- rep(1:sims, each = n_v)

  m <- purrr::pmap(list(sim_exp$n, sim_exp$m, sim_exp$s), rnorm)
  sim_exp$mn <- vapply(m, function (x) sum (x) /length (x), FUN.VALUE = 1.0)
  #vapply faster than sapply, sum/length faster than mean
  sim_exp$sd <- vapply(m, function(z) sqrt((sum(z^2)/length (z) - (sum(z)/length (z))^2) *length (z)/(length (z)-1)), FUN.VALUE = 1.0)
  # as with mean, this is substantially faster than using the SD function

  rm(m) #remove very large file

  #identify matches

  #Pb3
  if(vb == "y") {setTxtProgressBar(pb,50)
    cat("\nCounting expected matches ...\nThis takes the most time, especially if large number of simulations.\n\n")}

  d <- sim_exp %>% dplyr::mutate(mn = round(mn, dpm), sd = round(sd, dps), msd = paste0(mn, sd))
  data.table::setDT(d)
  d <- d[,  mn_m:=.N, by=.(var, sim, mn)]
  d <- d[,  mn_s:=.N, by=.(var, sim, sd)]
  d <- d[,  mn_ms:=.N, by=.(var, sim, msd)]
  data.table::setDF(d)

  #count number of matches for each variable- based on earlier function but uses datatable for speed
  match_count1 <- function (v, v_m) {
    x <- d %>% dplyr::select(var, sim, {{v}}, {{v_m}}) %>% dplyr::mutate("{{v_m}}" := ifelse({{v_m}} < 2, 0, {{v_m}}))
    ord <- c("var", "sim", deparse(substitute(v_m)), deparse(substitute(v)))
    x <- data.table::as.data.table(x)
    data.table::setorderv(x, ord, c(1,1,-1,1))
    x <- unique(x)
    x <- data.table::setDF(x) %>% dplyr::group_by(var, sim) %>%
      dplyr::summarise ("{{v_m}}" := sum({{v_m}}), .groups = "keep") %>% dplyr::ungroup ()
    return(x)}

  d_all <- data.frame("var"= rep(unique(b$var), each= sims), "sim"= rep(1:sims, times = length(unique(b$var)))) %>%
    dplyr::left_join(match_count1(mn, mn_m), by = c("var","sim")) %>%
    dplyr::left_join(match_count1(sd, mn_s), by = c("var","sim")) %>%
    dplyr::left_join(match_count1(msd, mn_ms), by = c("var","sim")) %>% replace(is.na(.), 0) %>% dplyr::arrange(var, sim)

  #make table of matches by cohort
  mx <- d_all %>% dplyr::select(mn_m:mn_ms) %>% max(., na.rm = TRUE)
  exp <- data.frame(freq = 0:mx) %>%
    dplyr::left_join(d_all %>% tidyr::gather(k, freq, -var) %>% dplyr::count(k,freq) %>%
                       dplyr::mutate(n = n/sims) %>% tidyr::spread(k, n), by= "freq") %>%
    dplyr::select(freq:mn_m, mn_s, mn_ms) %>% replace(is.na(.), 0) %>% setNames(c("freq", "exp_m", "exp_s","exp_ms"))

  rm(sim_exp,d)

  #Pb4
  if(vb == "y") {setTxtProgressBar(pb,90)
    cat("\nAlmost done... making table and graphs.\n")}

  #probability
  #expected data
  prob <- dplyr::bind_rows(
    data.frame(matrix(ncol=mx+3, nrow = 0)) %>% setNames(., c("var", "stat", 0:mx)) %>% dplyr::mutate_all(as.numeric) %>%
      dplyr::mutate(across(1:2, as.character)),
    d_all %>% dplyr::select (var, mn_m) %>% dplyr::group_by(var, mn_m) %>%
      dplyr::summarise (n = dplyr::n(), .groups = "keep") %>%
      tidyr::spread(mn_m, n) %>% dplyr::mutate (stat = "mean") %>% dplyr::select(var, stat, everything()),
    d_all %>% dplyr::select (var, mn_s) %>% dplyr::group_by(var, mn_s) %>%
      dplyr::summarise (n = dplyr::n(), .groups = "keep") %>%
      tidyr::spread(mn_s, n) %>% dplyr::mutate (stat = "SD") %>% dplyr::select(var, stat, everything()),
    d_all %>% dplyr::select (var, mn_ms) %>% dplyr::group_by(var, mn_ms) %>%
      dplyr::summarise (n = dplyr::n(), .groups = "keep") %>%
      tidyr::spread(mn_ms, n) %>% dplyr::mutate (stat = "mean_SD") %>% dplyr::select(var, stat, everything())) %>%
    replace(is.na(.), 0) %>% dplyr::select(- '1') %>% dplyr::rename ('1' = '0')

  colnames(prob) [1] <- "variable"

  # add observed data
  prob1 <- b_all %>%setNames(., c("variable", "mean", "SD", "mean_SD")) %>% tidyr::gather(stat, val, -variable) %>%
    dplyr::left_join (prob, by = c("variable", "stat")) %>% replace(is.na(.), 0)

  #find median row for each one
  df1 <- t(prob1[, 4:NCOL(prob1)]) %>% as.data.frame
  prob1$med <- apply(df1, 2, function (x) {which.min(abs(cumsum(x) -sims/2))})
  prob1$val_1 <- ifelse(prob1$val ==0, 1, prob1$val) #0 matches = 1 occurence

  #calculate number of matches for observed value or more extreme
  #observed > median - 0 to observed, obs < median, obs to max, obs = med smallest (ie more extreme) of val1/val2
  prob1 <- prob1 %>% dplyr::mutate(
    med1 = purrr::map2_dbl(df1, med, ~sum(.x [1:.y], na.rm= TRUE)), #sum 0 to med
    med2 = purrr::map2_dbl(df1, med, ~sum(.x [.y:length(.x)], na.rm= TRUE)), #sum med to max
    val1 = purrr::map2_dbl(df1, val, ~sum(.x [.y:length(.x)], na.rm= TRUE)/sims), #sum 0 to observed
    val2 = purrr::map2_dbl(df1, val, ~sum(.x [1:.y], na.rm= TRUE)/sims), #sum observed to max
    p = ifelse(val_1 > med | (val_1 == med & med1 > med2), val1,
               ifelse(val_1 < med | (val_1 == med & med1 < med2), val2, NA)))

  #make into suitable table
  #pvalues
  prob <- prob1 %>% dplyr::select(variable:stat, p, -val) %>%  tidyr::spread(stat,p) %>%
    dplyr::select(variable, mean, SD, mean_SD) %>% dplyr::arrange(desc(mean_SD)) %>%
    dplyr::mutate(m1 = ifelse(mean==0, 1/sims, mean), s1 = ifelse(SD==0, 1/sims, SD),
                  ms1 = ifelse(mean_SD==0, 1/sims, mean_SD))
  #add number of cohorts
  prob <- prob %>% dplyr::left_join (b %>% dplyr::select(var, m,s, ms) %>% dplyr::group_by(var) %>%
                                       dplyr::summarise(across(everything (), ~sum(!is.na(.)))),
                                     by = c("variable" = "var")) %>%
    dplyr::rename(mean_n = m, sd_n = s, msd_n = ms)
  # matches
  prob <- prob %>% dplyr::left_join(b_all, by = c("variable" = "var")) %>%
    dplyr::rename (mean_m = mn_m, sd_m = mn_s, msd_m = mn_ms)
  #max matches
  prob <- prob %>% dplyr::left_join(b %>% dplyr::select(var, mn_m, mn_s, mn_ms) %>% dplyr::group_by(var) %>%
                                      dplyr::summarise(across(everything(), ~max(., na.rm = TRUE))),
                                    by = c("variable" = "var")) %>%
    dplyr::rename(mean_mx = mn_m, sd_mx = mn_s, msd_mx = mn_ms)
  #make any 0s <- 1
  prob[, c("mean_mx", "sd_mx", "msd_mx")]  <- purrr::modify(prob[, c("mean_mx", "sd_mx", "msd_mx")], function (x) ifelse(x ==1, 0, x))
  #add mode
  prob <- prob %>% dplyr::left_join(c %>% dplyr::mutate (mean_mode = paste0(modemn, " (", modesd,")")) %>%
                                      dplyr::select(var, mean_mode),  by = c("variable" = "var"))

  #improve format
  rnd <- function(x,y,z) {setNames(data.frame(paste0(prob [, y], " (", round(prob [, y] / prob [ , z]*100,0),")")), x)}

  prob <- prob %>% dplyr::bind_cols(
    purrr::pmap_dfc(list(c("mn_m","s_m","ms_m","mean_mx_m","sd_mx_m","msd_mx_m"),
                         c("mean_m","sd_m","msd_m","mean_mx","sd_mx","msd_mx"),
                         c("mean_n","sd_n","msd_n","mean_n","sd_n","msd_n")), rnd))

  rnd1 <- function(x,y) {setNames(data.frame(
    ifelse(prob [, y]<0.001, "<0.001",
           ifelse(prob [, y]<0.01, as.character(round(prob [, y],3)),
                  ifelse(prob [, y]>0.99, ">0.99", as.character(round(prob [, y],2)))))),x)}

  prob <- prob %>%  dplyr::bind_cols(purrr::map2_dfc(c("m2", "s2", "ms2"), c("m1", "s1", "ms1"), rnd1))

  prob_f <- prob %>%
    dplyr::select(variable, mean_n, mn_m, mean_mx_m, m2, sd_n, s_m, sd_mx_m, s2, mean_mode, ms_m, msd_mx_m, ms2) %>%
    dplyr::arrange (desc(mean_n))
  colnames(prob_f) <- c("Variable","Mean \n(N)", "Mean Match \nN (%)", "Largest Number Mean Matches\nN (%)", "P Mean Match",
                        "SD \n(N)","SD Match\nN (%)", "Largest Number SD Matches\nN (%)","P SD Match",
                        "Mean (SD) \nMode or Average", "Mean and SD Match \nN (%)", "Largest Number Mean and SD Matches\nN (%)",
                        "P Mean and SD Match")

  #make flextable
  f <- flextable::flextable(prob_f)
  f <- flextable::add_header_row(f,
                                 values = c("", "Mean", "", "", "", "SD", "", "", "", "Both Mean and SD", "", "", ""), top = TRUE )
  f <- flextable::merge_at(f, i=1, j=2:5, part = "header")
  f <- flextable::merge_at(f, i=1, j=6:9, part = "header")
  f <- flextable::merge_at(f, i=1, j=10:13, part = "header")
  f <- flextable::align(f, i=1, j = 2:13, align = "center", part = "header") #note spelling of centre
  f <- flextable::align(f, j = 2:13, align = "center", part = "all")
  f <- flextable::border_remove(f)
  f <- flextable::hline_bottom(f, border= officer::fp_border(color="black", width = 2), part = "header")
  f <- flextable::hline_bottom(f, border= officer::fp_border(color="black", width = 2), part = "body")
  f <- flextable::hline(f, i = 1, j=2:13, border= officer::fp_border(color="black", width = 1), part = "header")
  #make footnotes
  f <- flextable::footnote(f, i = 2, j = c(2:4,6:8, 12),
                           value = flextable::as_paragraph(
                             paste0("The columns represent the number of reported means and SD; ",
                                    "the number of means, SD, and both mean and SD that match amongst the different cohorts; ",
                                    "and the largest number of means, SD, or both matching mean and SD in a single cohort ",
                                    "that had a matching value in a different cohort")),
                           ref_symbols = c("a"),
                           part = "header")
  f <- flextable::footnote(f, i = 2, j = 10:11,
                           value = flextable::as_paragraph(
                             paste0("where there was no matching mean (SD), the average mean and SD in the cohorts weighted by participant ",
                                    "are reported, otherwise the most common mean (SD), ie the mode, is reported together with the number of ",
                                    "occurrences of the mode")),
                           ref_symbols = c("b"),
                           part = "header")
  f <- flextable::footnote(f, i = 2, j = c(5,9,13),
                           value = flextable::as_paragraph(
                             paste0("P refers to the probability of the reported number of matching means, matching SDs ",
                                    "or both matching mean and matching SD for each variable in ", sims, " simulations")),
                           ref_symbols = c("c"),
                           part = "header")
  f <- flextable::font(f, fontname = "Arial", part = "all")
  f <- flextable::font(f, fontname = "Times", part = "footer")
  f <- flextable::fontsize(f, size = 8, part= "all")
  f <- flextable::set_table_properties(f, width = 1, layout = "autofit")

  #obs to expected
  oe <- dplyr::full_join(obs, exp, by = "freq") %>% replace(is.na(.), 0) %>%
    dplyr::mutate (p_XS_m = NA, p_XS_s = NA, p_XS_ms = NA)

  cs <- function(x,y) {suppressWarnings(chisq.test(x, p= y/sum(y)) [[3]])}
  oe$p_XS_m [1] <-cs(oe$ob_m, oe$exp_m)
  oe$p_XS_s [1] <-cs(oe$ob_s, oe$exp_s)
  oe$p_XS_ms [1] <-cs(oe$ob_ms, oe$exp_ms)

  #compress
  n_gp <- sum(oe$ob_m)
  for (sz in seq(5,ceiling(n_gp/5)*5,5)) {
    oe_temp <- compress(oe, "ob_m", "exp_m", sfx= "m", size= sz) %>%
      dplyr::bind_cols(compress(oe, "ob_s", "exp_s", sfx= "s", size= sz)) %>%
      dplyr::bind_cols(compress(oe, "ob_ms", "exp_ms", sfx= "ms", size= sz))
    if (sz == 5) {
      oe <- dplyr::bind_cols(oe, oe_temp)
      if (sum(!is.na(oe [ , "obsm_m_5"])) < 22) {break}
    }
    if (sum(!is.na(oe_temp [ , paste0("obsm_m_",sz)])) < 22) {
      oe <- dplyr::bind_cols(oe, oe_temp)
      break}
  }

  #make a graph
  graph <- list()
  graph [1:2] <- cat_graph(gph = oe, xtitle = "Number of mean matches per cohort", ytitle = "Number of cohorts",
                           size =5, sfx= "m", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                           fn = "cohort", ti = "Y", top = "yes", t= title)
  graph [3:4] <- cat_graph(gph = oe, xtitle = "Number of SD matches per cohort", ytitle = "Number of cohorts",
                           size =5, sfx= "s", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                           fn = "cohort", ti = "Y", top = "no", t= title)
  graph [5:6] <- cat_graph(gph = oe, xtitle = "Number of mean and SD matches per cohort", ytitle = "Number of cohorts",
                           size =5, sfx= "ms", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                           fn = "cohort", ti = "Y", top = "no", t= title)

  graphs <- ggpubr::ggarrange(graph [[2]], graph [[4]], graph [[6]], ncol = 1, nrow= 3, align = "hv")

  #obs to exp ratio
  rr_g <- list()
  rr_g [1:2] <- rr_graph(rr = oe, xtitle = "Number of mean matches per cohort", ytitle = "Number of cohorts",
                         size =5, sfx= "m", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                         fn = "cohort", ti = "Y", top = "yes", t= title)
  rr_g [3:4] <- rr_graph(rr = oe, xtitle = "Number of SD matches per cohort", ytitle = "Number of cohorts",
                         size =5, sfx= "s", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                         fn = "cohort", ti = "Y", top = "no", t=title)
  rr_g [5:6] <- rr_graph(rr = oe, xtitle = "Number of mean and SD matches per cohort", ytitle = "Number of cohorts",
                         size =5, sfx= "ms", text = list(length(unique(b$study)),length(unique(b$var)), sims),
                         fn = "cohort", ti = "Y", top = "no", t=title)

  rr_graphs <- ggpubr::ggarrange(rr_g [[2]], rr_g [[4]], rr_g [[6]], ncol = 1, nrow= 3, align = "hv")

  both_graphs <- list ()
  both_graphs [[1]] <- ggpubr::ggarrange(graph [[2]], rr_g [[2]], graph [[4]], rr_g [[4]],
                                         graph [[6]], rr_g [[6]], ncol = 2, nrow= 3, align = "hv")
  both_graphs [[2]] <- ggpubr::ggarrange(graph [[1]], rr_g [[2]], ncol = 2, nrow= 1, align = "hv")
  both_graphs [[3]] <- ggpubr::ggarrange(graph [[3]], rr_g [[4]], ncol = 2, nrow= 1, align = "hv")
  both_graphs [[4]] <- ggpubr::ggarrange(graph [[5]], rr_g [[6]], ncol = 2, nrow= 1, align = "hv")

  results <- list(cohort_ft = f, cohort_graph = graphs,
                  cohort_all_graphs =  list(all_graph = both_graphs [[1]], both_graphs = both_graphs [c(2,3,4)],
                                            individual_graphs = c(graph[c(1,3,5)], rr_g [c(1,3,5)])),
                  cohort_cohort_data= b, cohort_prob_data =prob_f, cohort_oe_data = oe)

  if(vb == "y") {
    print(f)
    close(pb)}

  return(results)
}

