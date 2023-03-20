#' Compares observed and expected distribution of a categorical (binomial) variable
#'
#' Creates plots of observed to expected numbers and ratios for the specified binomial variable \cr
#'
#' An example is for trial withdrawls in Bolland 2021\cr
#' Bolland MJ, Gamble GD, Avenell A, Cooper DJ, Grey A. Participant withdrawals were unusually distributed in randomized trials with integrity concerns: a statistical investigation. J Clin Epidemiol 2021;131:22-29.
#'
#' Returns a list containing 4 objects and (if verbose = TRUE) prints the plot cat_graph
#'
#' @param df data frame generated from load_clean function
#' @param x_title name of the variable for use on the x-axis
#' @param prefix letter for variable columns in data frame
#' @param del.disparate if yes, data in which the absolute difference between group sizes is >20% are deleted
#' @param title title name for plots (optional)
#' @param verbose TRUE or FALSE indicates whether to print plot
#'
#' @return list containing 4 objects as described
#'
#'\itemize{
#'   \item cat_graph = plot of observed to expected numbers and differences between groups, top panels are the absolute numbers, bottom panels are the differences between trial arms in two arm studies
#'   \item cat_data_abs = data frame of data for absolute numbers
#'   \item cat_data_df = data frame of data for difference between groups in two arm studies
#'   \item cat_all_graphs = list containing
#'      \itemize{
#'         \item abs = plot for absolute numbers only
#'         \item df = plot for difference between groups in two arm studies only
#'         \item individual_graphs list of 4 individual plots making up composite figures
#'     }
#'   }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr select count group_by filter rename left_join bind_cols
#' @importFrom tidyr gather
#' @importFrom ggpubr ggarrange
#' @importFrom data.table setorder as.data.table
#' @export cat_fn
#'
#'
#' @examples
#' # load example data
#' cat_data <- load_clean(import= "no", file.cat = "SI_cat", cat= "yes",
#' format.cat = "wide", cat.names = c("n", "w"))$cat_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' cat_fn(x_title= "withdrawals", prefix="w", del.disparate = "yes")$cat_graph
#'
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
#' cat_data <- load_clean(import= "yes", cat = "yes", dir = path,
#'     file.name.cat = "reappraised_examples.xlsx", sheet.name.cat = "SI_cat",
#'     range.name.cat = "A:G", cat.names = c("n", "w"), format.cat = "wide")$cat_data}
#'
#' @md

cat_fn <- function (df = cat_data, x_title = "", prefix="", del.disparate = "yes", title="", verbose = TRUE) {

  #import dataset
  a <- df

  if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type [1]) | a$data_type [1] != "cat") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
                   deparse(substitute(df)), "$data_type != 'cat'"))
  } else if(nrow(a) ==1) {stop_mb("There is only one row of study data- the code will not work")}

  g <- max(a$group) #largest number of groups
  nm.n <- colnames(a) [substr(colnames(a),1,1) == "n"]
  nm.v <- colnames(a) [substr(colnames(a),1,1) == prefix]

  #calculate p
  a$p <- rowSums(a[nm.v],na.rm =TRUE)/rowSums(a[nm.n],na.rm = TRUE)

  #get largest group and add 1 because there can be 0 withdrawals
  n <- max(a[nm.n], na.rm = TRUE)+1

  #starting row for each group
  rw_g <- list()
  for (i in 1:g) {rw_g [[i]] <-  (i-1)*n +1}

  #get expected p
  exp_p <- expect_p(a, n, g, nm.n, rw_g)

  #sum the expected probs
  sum_exp <- 0
  for (i in 1:g) {sum_exp <- sum_exp + exp_p [rw_g [[i]]: (rw_g [[i]]+n-1),]}

  #sum the observed
  obs <- a %>% dplyr::select(all_of(nm.v)) %>% tidyr::gather (k, v) %>% dplyr::group_by (v) %>%
    dplyr::count () %>% dplyr::filter (!is.na(v)) %>% dplyr::rename (obs = n) %>% dplyr::ungroup()

  b <- dplyr::left_join(data.frame(num= 0: (n-1), exp= rowSums(sum_exp)), obs, by = c("num" = "v")) %>%
    replace(is.na(.), 0)

  #if large number of groups then reduce
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
  graph [1:2] <- cat_graph (gph= b,  xtitle = paste0("Number of ", x_title," per group"), ytitle = "Number of trial groups",
                            size = sz, sfx = "", text = list(length(a$study), sum(b$obs)),
                            fn = "cat", ti = "Y", top = "Y", t= title)

  rr_g <- list()
  rr_g [1:2] <- rr_graph(rr = b, xtitle = paste0("Number of ", x_title," per group"), ytitle = "Observed/Expected ratio",
                         size =sz, sfx= "", text = list(length(a$study), sum(b$obs)), fn = "cat", ti = "Y", top = "yes", t=title)

  #repeat for 2 arm studies only;
  a.2 <- a %>% dplyr::filter (group == 2) %>% dplyr::select(study, all_of(nm.n) [1:2], all_of(nm.v) [1:2])
  exp_p2 <- exp_p [, a$group == 2]

  #delete disparate rows
  if (tolower(substr(del.disparate,1,1)) == "y") {
    dd_f <- with(a.2, !(n1/n2 < 0.8 | n1/n2 >1.2 | n2/n1 < 0.8 | n2/n1 > 1.2))
    a.2 <- a.2 %>% dplyr::filter(dd_f)
    exp_p2 <- exp_p2 [, dd_f, drop= FALSE]
    rm (dd_f)}

  #expected differences
  #calculate probability of group numbers differing by 0, 1, 2 etc# this usually takes the longest amount of time.
  exp_p3 <- expect_p_diff (a.2, exp_p2, n, nm.n)

  #count observed diffs
  a.2$obs = abs(a.2 [, nm.v [[1]]] - a.2[ , nm.v [[2]]])
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

  graph [3:4] <-  cat_graph (gph= b_df,  xtitle = paste0("Differences in ", x_title," between trial groups"),
                             ytitle = "Number of trials", size = sz_df, sfx = "", text = list(sum(b_df$obs)),
                             fn = "cat", ti = "n", top = "n", t= title)

  rr_g [3:4] <- rr_graph(rr = b_df ,xtitle = paste0("Differences in ", x_title," between trial groups"),
                         ytitle = "Observed/Expected ratio", size =sz_df, sfx= "", text = list(sum(b_df$obs)),
                         fn = "cat", ti = "n", top = "n", t= title)

  #graphs <- ggpubr::ggarrange(graph [[2]], graph [[4]], ncol = 1, nrow= 3, align = "hv")
  #rr_graphs <- ggpubr::ggarrange(rr_g [[2]], rr_g [[4]], ncol = 1, nrow= 2, align = "hv")

  both_graphs <- list ()
  both_graphs [[1]] <- ggpubr::ggarrange(graph [[2]], rr_g [[2]], graph [[4]], rr_g [[4]],
                                         ncol = 2, nrow= 2, align = "hv")
  both_graphs [[2]] <- ggpubr::ggarrange(graph [[1]], rr_g [[2]], ncol = 2, nrow= 1, align = "hv")
  both_graphs [[3]] <- ggpubr::ggarrange(graph [[3]], rr_g [[4]], ncol = 2, nrow= 1, align = "hv")

  results <- list(cat_graph= both_graphs [[1]], cat_data_abs = b, cat_data_df = b_df,
                  cat_all_graphs =  list(abs = both_graphs [[2]], df = both_graphs [[3]],
                                         individual_graphs = c(graph[c(1,3)], rr_g [c(1,3)])))

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {print(both_graphs [[1]])}

  return (results)
}

