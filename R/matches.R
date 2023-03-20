#' Compares proportions of matching summary statistics within two-arm randomised trials
#'
#' Creates flextable of matching summary statistics by significant figures with Reference data\cr
#'
#' Reference data is from Bolland 2021\cr
#' Bolland MJ, Gamble GD, Avenell A, Grey A. Identical summary statistics were uncommon in randomized trials and cohort studies. J Clin Epidemiol 2021;136:180-188.
#'
#' Returns a list containing 6 objects and (if verbose = TRUE ) prints the flextable match_ft_all
#'
#' @param df data frame generated from load_clean function
#' @param verbose TRUE or FALSE indicates whether to print flextable
#'
#' @return list containing 6 objects as described
#'\itemize{
#'   \item match_ft_all = flextable of matches with reference data
#'   \item match_ft = flextable of matches (no reference data)
#'   \item ref_match_ft = flextable of reference data
#'   \item match_match_data = data frame of results used in calculations
#'   \item match_table = data frame of matches used to make flextable
#'   \item ref_table = data frame of reference data used to make flextable
#'   }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select count filter left_join
#' @importFrom tidyr gather spread
#' @importFrom officer fp_border
#' @importFrom purrr map modify2 pmap
#' @importFrom flextable merge_at font fontsize align bold border_remove hline_bottom hline set_table_properties
#' @export match_fn
#'
#'
#' @examples
#' # load example data
#' match_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", match= "yes",
#' format.cont = "wide")$match_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' match_fn()$match_ft_all
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
#' match_data <- load_clean(import= "yes", match = "yes", dir = path,
#'      file.name.cont = "reappraised_examples.xlsx", sheet.name.cont = "SI_pvals_cont",
#'      range.name.cont = "A:O", format.cont = "wide")$match_data}
#'
#' @md



# function for baseline matching vars within an RCT ---------------------------------
match_fn <- function(df = match_data, verbose = TRUE) {

  # import dataset
  a <- df

  if (!"data_type" %in% colnames(a)) {a$data_type <- 0}
  if (is.na(a$data_type[1]) | a$data_type[1] != "match") {
    stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
               deparse(substitute(df)), "$data_type != 'match'"))}

  # get sig figs
  sigfigs <- function(j) {
    x <- format(a [, j], scientific = FALSE, trim = TRUE)
    x <- sub("\\.", "", x)
    x <- gsub("(^0+|0+$)", "", x)
    return(nchar(x))
  }

  sf <- c("m1", "m2", "s1", "s2")
  sf1 <- paste0(sf, "_sf")

  a[, sf1] <- purrr::map(sf, sigfigs)

  a <- a %>%
    dplyr::mutate(sfm = pmax(m1_sf, m2_sf, na.rm = TRUE), sfs = pmax(s1_sf, s2_sf, na.rm = TRUE)) %>%
    dplyr::select(-all_of(sf1))

  # now do analysis
  a <- within(a, {
    mnm <- ifelse(is.na(m1) | is.na(m2), NA, ifelse(m1 == m2, 1, 0)) # mean match
    sdm <- ifelse(is.na(s1) | is.na(s2), NA, ifelse(s1 == s2, 1, 0)) # sd match
    msm <- ifelse(is.na(mnm) | is.na(sdm), NA, ifelse(mnm == 1 & sdm == 1, 1, 0)) # both match
    sfsf <- ifelse(sfs > 5, 5, sfs)
    sfmf <- ifelse(sfm > 5, 5, sfm)
    sfmn <- pmin(sfsf, sfmf, na.rm = TRUE)
    sfmx <- pmax(sfsf, sfmf, na.rm = TRUE)
  })

  max_sf <- max(a$sfmx, na.rm = TRUE)

  # create a table
  ressf <- c("", "Number of variables, n (%)", "Mean", "SD", "Proportion of matching summary statistics in both treatment groups (%)",
             "Means match", "SDs match", "Both Mean and SD matches")

  ressf <- data.frame(ressf, stringsAsFactors = FALSE) # rownames
  ressf[1, 1] <- paste(length(unique(a$study)), "trials") # trials
  ressf[, 2:(max_sf + 2)] <- "" # columns for SFs
  ressf[3, 2] <- sum(!is.na(a$mnm)) # number of mean, SD for whole cohort
  ressf[4, 2] <- sum(!is.na(a$sdm))
  ressf[6, 2] <- ifelse(is.na(table(a$mnm)[2]), 0, round(table(a$mnm)[2] / as.numeric(ressf[3, 2]) * 100, 1)) # proportion of mean matches
  ressf[7, 2] <- ifelse(is.na(table(a$sdm)[2]), 0, round(table(a$sdm)[2] / as.numeric(ressf[4, 2]) * 100, 1)) # SD
  ressf[8, 2] <- ifelse(is.na(table(a$msm)[2]), 0, round(table(a$msm)[2] / min(as.numeric(ressf[3, 2]), as.numeric(ressf[4, 2])) * 100, 1)) # both

  # count variables
  ct <- function(v, v1 = "", f1 = "") {
    if (v1 != "") {
      t <- a %>% dplyr::count(.data[[v]], .data[[v1]]) %>% tidyr::spread(.data[[v1]], n)
    } else if (f1 != "") {
      t <- a %>% dplyr::filter(.data[[f1]] == 1) %>% dplyr::count(.data[[v]])
    } else {
      t <- a %>% dplyr::count(.data[[v]])}

    if (v1 == "") {colnames(t) <- c("Var1", paste0(v, v1, f1))
    } else {colnames(t)[1] <- "Var1"}

    t <- dplyr::left_join(tb, t, by = "Var1") %>% dplyr::select(-Var1)
    return(t)
  }

  tb <- data.frame("Var1" = as.numeric(1:max_sf))

  # generates counts
  tb <- dplyr::bind_cols(tb, purrr::pmap_dfc(list(c("sfmf", "sfsf", "sfmf", "sfsf", "sfmn", "sfmn"),
                                                  c(rep("", 5), "mnm"), c(rep("", 2), "mnm", "sdm", "msm", "")), ct))

  # format data
  if (!"1" %in% colnames(tb)) {tb$"1" <- 0}
  tb[is.na(tb)] <- 0
  tb[, 4:5] <- purrr::modify2(tb[, 4:5], tb[, 2:3], function(x, y) sprintf("%.1f", round(x / y * 100, 1)))
  tb[, 2:3] <- purrr::modify2(tb[, 2:3], as.data.frame(t(ressf[3:4, 2])[rep(1, nrow(tb)), ]),
                              function(x, y) paste0(x, " (", sprintf("%.1f", round(x / as.numeric(y) * 100, 1)), ")"))
  tb <- cbind(tb[, 1:3], "blank" = rep("", nrow(tb)), tb[, 4:8])
  tb$sfmnmsm <- sprintf("%.1f", round(tb$sfmnmsm / rowSums(tb[, c("0", "1")]) * 100, 1))

  # import into results
  ressf[3:8, 3:(2 + max_sf)] <- t(tb[, 2:7])
  ressf[ressf == "NaN"] <- ""

  # format dataset
  a <- a %>% dplyr::select(study:s2, sfmf:sfsf, mnm, sdm, msm, data_type, -c(sfmx:sfmn), -c(sfm:sfs))

  colnames(a)[colnames(a) == "mnm"] <- "means_match"
  colnames(a)[colnames(a) == "sdm"] <- "sd_match"
  colnames(a)[colnames(a) == "msm"] <- "both_mean_sd_match"
  colnames(a)[colnames(a) == "sfsf"] <- "sig_dig_SD"
  colnames(a)[colnames(a) == "sfmf"] <- "sig_dig_mean"

  # format results
  # if less than 5 sig figs add blank columns
  colnames(ressf) <- c("Variable", "All variables", paste(1:max_sf, " Sig. dig.", sep = ""))

  zz <- NCOL(ressf) - 2
  if (zz == 1) {ressf <- cbind(ressf, "2 Sig. dig." = "", "3 Sig. dig." = "", "4 Sig. dig." = "", "5 Sig. dig." = "", stringsAsFactors = FALSE)}
  if (zz == 2) {ressf <- cbind(ressf, "3 Sig. dig." = "", "4 Sig. dig." = "", "5 Sig. dig." = "", stringsAsFactors = FALSE)}
  if (zz == 3) {ressf <- cbind(ressf, "4 Sig. dig." = "", "5 Sig. dig." = "", stringsAsFactors = FALSE)}
  if (zz == 4) {ressf <- cbind(ressf, "5 Sig. dig." = "", stringsAsFactors = FALSE)}

  # make flextables
  f.sf <- flextable::flextable(ressf)
  f.sf <- flextable::merge_at(f.sf, i = 5, j = 1:3)
  f.sf <- flextable::font(f.sf, fontname = "Arial")
  f.sf <- flextable::fontsize(f.sf, size = 8)
  f.sf <- flextable::fontsize(f.sf, part = "header", size = 10)
  f.sf <- flextable::align(f.sf, j = 2:5, align = "center", part = "all")
  f.sf <- flextable::bold(f.sf, j = 1, part = "body")
  f.sf <- flextable::bold(f.sf, part = "header")
  f.sf <- flextable::border_remove(f.sf)
  f.sf <- flextable::hline_bottom(f.sf, border = officer::fp_border(color = "black", width = 2), part = "header")
  f.sf <- flextable::hline(f.sf, i = 2, j = 1, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.sf <- flextable::hline(f.sf, i = 5, j = 1:3, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.sf <- flextable::set_table_properties(f.sf, width = 1, layout = "autofit")

  # reference data
  ref <- data.frame(
    "Variable" = c("3585 trials", "Number of variables, n (%)", "Mean", "SD", "Proportion of matching summary statistics in both treatment groups (%)", "Means match", "SDs match", "Both Mean and SD matches"),
    "All variables" = c("", "", "21948", "21145", "", "13.4", "14.8", "5.1"),
    "1 Sig. dig." = c("", "", "607 (2.8)", "4320 (20.4)", "", "66.4", "40.5", "16.9"),
    "2 Sig. dig." = c("", "", "8149 (37.1)", "10682 (50.5)", "", "20.6", "11.6", "2.7"),
    "3 Sig. dig." = c("", "", "10978 (50.0)", "5618 (26.6)", "", "7.5", "2.7", "0.2"),
    "4 Sig. dig." = c("", "", "2082 (9.5)", "552 (2.6)", "", "1.6", "0.4", "0.0"),
    "5 Sig. dig." = c("", "", "132 (0.6)", "36 (0.2)", "", "0.8", "0.0", "0.0"))

  colnames(ref) <- c("Variable", "All variables", "1 Sig. dig.", "2 Sig. dig.", "3 Sig. dig.", "4 Sig. dig.", "5 Sig. dig.")

  f.ref <- flextable::flextable(ref)
  f.ref <- flextable::merge_at(f.ref, i = 5, j = 1:3)
  f.ref <- flextable::font(f.ref, fontname = "Arial")
  f.ref <- flextable::fontsize(f.ref, size = 8)
  f.ref <- flextable::fontsize(f.ref, part = "header", size = 10)
  f.ref <- flextable::align(f.ref, j = 2:5, align = "center", part = "all")
  f.ref <- flextable::bold(f.ref, j = 1, part = "body")
  f.ref <- flextable::bold(f.ref, part = "header")
  f.ref <- flextable::border_remove(f.ref)
  f.ref <- flextable::hline_bottom(f.ref, border = officer::fp_border(color = "black", width = 2), part = "header")
  f.ref <- flextable::hline(f.ref, i = 2, j = 1, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.ref <- flextable::hline(f.ref, i = 5, j = 1:3, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.ref <- flextable::set_table_properties(f.ref, width = 1, layout = "autofit")

  # merge flextables

  f.all <- flextable::flextable(rbind(ressf, c(rep("", 7)), c("Reference data", rep("", 6)), ref))
  f.all <- flextable::merge_at(f.all, i = 5, j = 1:3)
  f.all <- flextable::merge_at(f.all, i = 10, j = 1:6)
  f.all <- flextable::merge_at(f.all, i = 15, j = 1:3)
  f.all <- flextable::font(f.all, fontname = "Arial")
  f.all <- flextable::fontsize(f.all, size = 8)
  f.all <- flextable::fontsize(f.all, part = "header", size = 10)
  f.all <- flextable::fontsize(f.all, i = 10, j = 1, part = "body", size = 10)
  f.all <- flextable::align(f.all, j = 2:5, align = "center", part = "all")
  f.all <- flextable::bold(f.all, j = 1, part = "body")
  f.all <- flextable::bold(f.all, part = "header")
  f.all <- flextable::border_remove(f.all)
  f.all <- flextable::hline_bottom(f.all, border = officer::fp_border(color = "black", width = 2), part = "header")
  f.all <- flextable::hline(f.all, i = 2, j = 1, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.all <- flextable::hline(f.all, i = 5, j = 1:3, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.all <- flextable::hline(f.all, i = 12, j = 1, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.all <- flextable::hline(f.all, i = 15, j = 1:3, border = officer::fp_border(color = "black", width = 1), part = "body")
  f.all <- flextable::hline(f.all, i = 10, part = "body", border = officer::fp_border(color = "black", style = "solid", width = 2))
  f.all <- flextable::set_table_properties(f.all, width = 1, layout = "autofit")

  results <- list(match_ft_all = f.all, match_ft = f.sf, ref_match_ft = f.ref, match_match_data = a, match_table = ressf, ref_table = ref)

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {print(f.all)}

  return (results)

}
