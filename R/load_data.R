#' Load data then clean and format it
#'
#' Function loads and cleans data for the nine functions
#'
#' Function can load continuous or categorical data.
#' Continuous data can be used for comparison of baseline p-values (pval_cont_fn),
#' matching summary stats within a trial (match_fn), matching summary stats in different cohorts (cohort_fn),
#' or comparing means of baseline p-values (anova_fn).
#' Categorical data can be used for comparisons of observed with expected distributions for single variable (cat_fn),
#' for group numbers in trials using simple randomisation (sr_fn), for all variables (cat_all_fn), and for comparison
#' of baseline p-values (pval_cat_fn).\cr
#'
#' There is one function in development that allows assessment of proportion of final digits in summary statistics (final_digit_fn).
#' This function works using summary statistics but could be adapted to use on raw continuous or categorical data.
#'
#' Only 1 continuous and/or 1 categorical data set allowed per load to avoid clashes
#'
#' Data can be imported from a file (import = "yes") or taken from an existing data frame, import = "no"
#'
#' If loading from an existing data use file.cont and file.cat
#'
#' If loading from common directory or file, can use dir and file.name rather than more specific dir.cont, dir.cat,
#' file.name.cont, or file.name.cat.
#'
#' **Comments about each indicator:**
#' *pval_cont*\cr
#' loads continuous data for pval_cont_fn, outputs as list of 1 containing named data frame pval_cont_data.\cr
#'
#' format should be study, variable or var, n, m, s, p. Can be in any order. n = sample size, m = mean, s = standard deviation,
#' p = baseline p value (can omit if not reported)\cr
#'
#' can be in wide or long format\cr
#' wide: study, var, n1, n2, n3 ..., m1, m2, m3 ... s1, s2, s3..., p\cr
#' long: study, var, group, m , s,  n , p\cr
#'
#' group or g or grp required for long format\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *match*\cr
#' loads continuous data for match_fn, outputs as list of 1 containing named data frame match_data\cr
#'
#' remainder is same as for pval_cont above.\cr
#' only difference between pval_cont and match is that match allows for missing mean or SD whereas pval_cont does not\cr
#'
#' format should be study, variable or var, n, m, s. Can be in any order. n = sample size, m = mean, s = standard deviation\cr
#'
#' can be in wide or long format\cr
#' wide: study, var, n1, n2, m1, m2, s1, s2, p\cr
#' long: study, var, group, m , s,  n \cr
#'
#' group or g or grp required for long format\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *cohort*\cr
#' loads continuous data for cohort_fn, outputs as list of 1 containing named data frame cohort_data\cr
#'
#' same as pval_cont but allows a lookup variable for variable names\cr
#'
#' format should be study, variable or var, n, m, s, p. Can be in any order. n = sample size, m = mean, s = standard deviation\cr
#'
#' can be in wide or long format\cr
#' wide: study, var, n1, n2, n3 ..., m1, m2, m3 ... s1, s2, s3... \cr
#' long: study, var, group, m , s,  n \cr
#'
#' group or g or grp required for long format\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#' lookup table is var_name_final, var_name_orig and allows you to specify a list of all variables names (var_name_orig)
#' from all studies and a lookup table of standardised names (var_name_final) allowing different names in different studies to
#' be standardised\cr
#'
#' has optional variable 'population' which can be used to subset the data if trials in different populations are reported\cr
#'
#'
#' *anova*\cr
#' loads continuous data for anova_fn, outputs as list of 1 containing named data frame anova_data\cr
#'
#' same as for pval_cont above but allows for optional value for decimal place\cr
#'
#' format should be study, variable or var, n, m, s, p. Can be in any order. n = sample size, m = mean, s = standard deviation,
#' d= decimal place of mean (if omitted, this is calculated automatically in anova_fn)\cr
#'
#' can be in wide or long format\cr
#' wide: study, var, n1, n2, n3 ..., m1, m2, m3 ... s1, s2, s3..., d\cr
#' long: study, var, group, m , s,  n , d\cr
#'
#' group or g or grp required for long format\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *cat*\cr
#' loads categorical data for cat_fn, outputs as list of 1 containing named data frame cat_data\cr
#'
#' format should be study, n, v. Can be in any order, n= group size, v= number with characteristic\cr
#'
#' can be in wide or long format\cr
#' wide: study, n1, n2, n3 ..., v1, v2, v3...\cr
#' long: study, group, n, v\cr
#'
#' group or g or grp required for long format\cr
#' use cat.names to name variable eg c("n", "v") , c("n", "g") ...\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *sr*\cr
#' loads categorial data for sr_fn, outputs as list of 1 containing named data frame sr_data\cr
#'
#' as for cat but only requires study and n
#'
#' format should be study, n. n= group size\cr
#'
#' can be in wide or long format\cr
#' wide: study, n1, n2, n3 ...\cr
#' long: study, group, n\cr
#'
#' group or g or grp required for long format\cr
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *cat_all*\cr
#' loads categorical data for cat_all_fn, outputs as list of 1 containing named data frame cat_all_data\cr
#'
#' format should be study, var or variable, n, N, level, stat, recode, p. Can be in any order, n = number with characteristic, N = group size,
#' p = baseline p value (can omit if not reported), can use "ns" for not significant or "<" or ">" to indicate threshold (eg "<0.05") \cr
#'
#' optional level - number for level of variable (eg y/n =1,2; high/med/low =1,2,3)\cr
#' optional recode- for variables with >2 levels to tell how to recode into 2 groups\cr
#' optional stat: statistical test used for p-value : chisq - Chisquare, chisqc- Chisquare with correction,
#' fisher- Fisher's exact, midp - midp -calculated using two different methods, lr- likelihood ratio,
#' mh - Mantel-Haenszel test\cr
#'
#' can be in wide or long format\cr
#' wide study, var, n1, n2, n3, ... N1, N2, N3... p, stat, level, recode\cr
#' long study, var, group, n, N, p, stat, level, recode\cr
#'
#' group or g or grp required for long format\cr
#'
#' if variable has 2 levels, only 1 required, other will be calculated.\cr
#'
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *pval_cat*\cr
#' loads categorical data for pval_cat_fn, outputs as list of 1 containing named data frame pval_cat_data\cr
#'
#' as for cat_all but recode variable is not generated
#'
#' format should be study, var or variable, n, N, p. Can be in any order, n = number with characteristic, N = group size,
#' p = baseline p value (can omit if not reported), can use "ns" for not significant or "<" or ">" to indicate threshold (eg "<0.05") \cr
#'
#' optional level - number for level of variable (eg y/n =1,2; high/med/low =1,2,3)\cr
#' optional stat: statistical test used for p-value : chisq - Chisquare, fisher- Fisher's exact\cr
#'
#' can be in wide or long format\cr
#' wide study, var, n1, n2, n3, ... N1, N2, N3... p, stat, level\cr
#' long study, var, group, n, N, p, stat, level\cr
#'
#' group or g or grp required for long format\cr
#'
#' if variable has 2 levels, only 1 required, other will be calculated.\cr
#'
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#'
#'
#' *generic*\cr
#' loads data for use generic use, outputs as list of 1 containing named data frame generic_data\cr
#'
#' use cont suffixes for file details: dir.cont (or dir), file.name.cont (or file.name), sheet.name,cont, range.name.cont)\cr
#'
#' format should be study, var or variable, variable names\cr
#'
#' optional gen.vars.keep = vector of variables to keep\cr
#' optional gen.vars.del = vector of variables to delete
#'
#' can be in wide or long format\cr
#' wide study, var, a1, a2..., b1, b2 ... \cr
#' long study, var, group, a, b, ....\cr
#'
#' group or g or grp required for long format\cr
#'
#' separators (eg n1 n_1 n.1) are stripped and replaced\cr
#' no data checking or other transformations take place\cr
#'
#'
#' @param import 'yes' indicates import excel file. 'no' indicates takes dataset already loaded into R as data frame
#' @param file.cont If import = 'no', name of data frame containing continuous data
#' @param file.cat If import = 'no', name of data frame containing categorical data
#' @param dir If import = 'yes', path to location of excel file for continuous and categorical data
#' @param file.name If import = 'yes', file name of excel file containing continuous and categorical data
#' @param pval_cont 'yes'/'no' indicating if data will be used for pval_cont_fn. Only data for 1 continuous data function can be loaded with each run of this function.
#' @param match 'yes'/'no' indicating if data will be used for match_fn. Only data for 1 continuous data function can be loaded with each run of this function.
#' @param cohort 'yes'/'no' indicating if data will be used for cohort_fn. Only data for 1 continuous data function can be loaded with each run of this function.
#' @param anova 'yes'/'no' indicating if data will be used for anova_fn. Only data for 1 continuous data function can be loaded with each run of this function.
#' @param dir.cont If import = 'yes', path to location of excel file for continuous data
#' @param file.name.cont If import = 'yes', file name of excel file containing continuous data
#' @param sheet.name.cont Sheet name containing continuous data
#' @param range.name.cont Range of cells containing continuous data. Can be in format 'a1:b20' or 'a:b'
#' @param format.cont 'wide'/'long' indicating continuous data is in wide or long format
#' @param cat 'yes'/'no' indicating if data will be used for cat_fn. Only data for 1 categorical data function can be loaded with each run of this function.
#' @param sr 'yes'/'no' indicating if data will be used for sr_fn. Only data for 1 categorical data function can be loaded with each run of this function.
#' @param cat_all 'yes'/'no' indicating if data will be used for cat_all_fn. Only data for 1 categorical data function can be loaded with each run of this function.
#' @param pval_cat 'yes'/'no' indicating if data will be used for cat_all_fn. Only data for 1 categorical data function can be loaded with each run of this function.
#' @param cat.names names of variables to be used in cat_fn and sr_fn
#' @param dir.cat If import = 'yes', path to location of excel file for categorical data
#' @param file.name.cat If import = 'yes', file name of excel file containing categorical data
#' @param sheet.name.cat Sheet name containing categorical data
#' @param range.name.cat Range of cells containing categorical. Can be in format 'a1:b20' or 'a:b'
#' @param format.cat 'wide'/'long' indicating categorical data is in wide or long format
#' @param generic 'yes'/'no' indicating if data to be loaded for generic use
#' @param gen.vars.keep Vector of variables in data to keep
#' @param gen.vars.del Vector of variables in data to delete
#' @param verbose TRUE/FALSE TRUE indicates comments will be printed during loading
#'
#' @return list containing a named data frame containing data in suitable format for appropriate function as described in Details
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom stats chisq.test fisher.test na.omit pbinom pnorm qnorm setNames
#' @importFrom utils getTxtProgressBar setTxtProgressBar txtProgressBar
#' @importFrom dplyr across all_of desc everything starts_with select mutate arrange count filter left_join mutate_all group_by summarise ungroup slice bind_rows bind_cols row_number rename
#' @importFrom tidyr spread gather unite pivot_wider separate
#' @importFrom readxl read_excel cell_cols
#' @importFrom purrr map map_int
#' @importFrom data.table as.data.table setorder
#'
#' @export load_clean
#' @examples
#' # examples of loading data for each function are given in the individual functions.
#' # Here is one- for pval_cont_fn():
#'
#' pval_cont_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", pval_cont= "yes",
#' format.cont = "wide")$pval_cont_data
#'
#' \donttest{
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

load_clean <- function (
    #can only load one continuous and one categorical file at a time
    #if loading a dataset already in R- import= "no" and file.cont/file.cat = name of data frame
    #if loading a excel file- import = "yes"
    import = "yes",
    file.cont = "",
    file.cat = "",
    #if common directory/file name can use dir/file.name, if different use specific names
    dir = "",
    file.name = "",
    #continuous data function to be used
    pval_cont = "no",
    match = "no",
    cohort = "no",
    anova = "no",
    #details of excel spreadsheet holding continuous data
    dir.cont ="",
    file.name.cont = "",
    sheet.name.cont = "Sheet1",
    range.name.cont = "",
    #format long or wide
    format.cont = "wide",
    #categorical data function to be used
    cat= "no",
    sr = "no",
    cat_all = "no",
    pval_cat = "no",
    #names of variable for categorical function- must have n eg c(n, v)
    cat.names = c("n"),
    #details of excel spreadsheet holding continuous data
    dir.cat ="",
    file.name.cat = "",
    sheet.name.cat = "Sheet1",
    range.name.cat = "",
    #format long or wide
    format.cat = "wide",
    #generic option
    generic = "",
    gen.vars.keep = "",
    gen.vars.del ="",
    verbose = TRUE
    ) {

# functions first ---------------------------------------------------------


  load_file <- function (dir1, range1, file, sheet1) {
    #check type of path
    dir1 <- ifelse("/" %in% dir1, dir1, gsub("\\\\","/",dir1))

    #check format of range
    r <- strsplit(range1, ":")
    r1 <- substring(r [[1]], nchar(r [[1]]), nchar(r [[1]])) #gets last char of each cell reference
    if (grepl("\\d", r1 [1]) != grepl("\\d", r1 [2])) {
      stop_mb("Incorrect range format")
    }
    if (grepl("\\d", r1 [2])) {
      a <- as.data.frame(readxl::read_excel (paste0(dir1,"/",file), sheet = sheet1 , range = range1),
                         stringsasfactors = FALSE)
    } else {
      a <- as.data.frame(readxl::read_excel (paste0(dir1,"/",file), sheet = sheet1, range = readxl::cell_cols(range1)),
                         stringsasfactors = FALSE)
    }
    return (a)
  }

  #pval
  pval_import <- function (x= "") {

    if (cohort == "y" | anova == "y") {
      a <- x
    } else {
      if (tolower(substr(import,1,1)) == "n") {a <- get(file.cont)
      } else {a <- load_file(dir= dir.cont, range1 = range.name.cont, file = file.name.cont, sheet1 = sheet.name.cont)}
    }

    #get names;
    colnames(a) <- lapply(colnames (a), tolower)
    colnames (a) [colnames(a) == "variable"] <- "var" #don't need an if statement here for this to work
    colnames (a) [colnames(a) %in% c("g", "group", "grp")] <- "group"

    #long
    if (tolower(format.cont) == "long") {
      a [, c("m", "n", "s")] <-   a [, c("m", "n", "s")] %>%  dplyr::mutate_all(as.numeric)

      #check if a p column;
      if (!"p" %in% colnames(a)) {a$p <- NA}

      a <- a %>% dplyr::select(study, var, group, m, n, s,p)

      #make wide format
      a_wide <- a %>% tidyr::pivot_wider(names_from = group, values_from = c(m, n, s, p), names_sep = "")

      #shift p values left
      names.p <- colnames(a_wide) [substring (colnames(a_wide), 1, 1) == "p"]
      p.check <- as.data.frame(t(apply(a_wide[, names.p], 1, function(x) {
        return(c(x[!is.na(x)], x[is.na(x)]))
      })))
      a_wide[, names.p] <- NULL
      if (sum(p.check, na.rm = TRUE) == 0) {a_wide$p <- NA
      } else {a_wide$p <- p.check[, 1]}

      a <- a_wide
      rm (a_wide, p.check, names.p)
    }

    #remove separators
    names <- colnames(a) [substr(colnames(a),1,1) %in% c("n","m","s") & colnames(a) != "study"]
    colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1), gsub("\\D", "", names))

    names <- colnames (a)
    names.m <- names [substring(names, 1, 1) == "m"]
    names.n <- names [substring(names, 1, 1) == "n"]
    names.s <- names [substring(names, 1, 1) == "s" & substring(names,2,2) %in% as.character(0:9)]


    #check groups
    g.m1 <- length(names.m)
    g.n1 <- length(names.n)
    g.s1 <- length(names.s)
    #stop if not matching#
    if (g.n1 != g.m1 | g.m1 != g.s1) {
      stop_mb ("Cont data error- numbers of groups don't match")}
    g <- g.n1 #number of groups

    #check equal and more than one group
    a <- within (a, {
      nm <- rowSums(!is.na(a[, names.m]))
      ns <- rowSums(!is.na(a[, names.s]))
      nn <- rowSums(!is.na(a[, names.n]))
      check_eq <- ifelse(nm == ns & ns == nn, 0, 1) # check equal groups
      check_n <- ifelse(nm == 1, 1, 0) #check more than one group
    })

    #for cohort ignore both of these errors, for match only check more than one group
    if (pval_cont == "y" ) {
      if (sum(a$check_eq > 0)) {
        stop_mb ("Cont data error: some groups do not match")}
      if (sum(a$check_n > 0)) {
        stop_mb ("Cont data error: some groups only 1 value")}
    }

    if (match == "y" ) {
      if (sum(a$check_n > 0)) {
        stop_mb ("Cont data error: some groups only 1 value")}
    }

    a <- a %>% dplyr::select (- c(check_n : nm))

    #rearrange
    names <- c("study", "var", paste0("n",1:g), paste0 ("m",1:g), paste0 ("s", 1:g))

    #add p column if there is one
    if ("p" %in% colnames(a)) {
      if (is.numeric(a$p)) {names <- c(names,"p")
      } else {
        suppressWarnings(a$p_num <- as.numeric(a$p))
        names <- c(names, "p", "p_num")
      }
    } else {a$p <- NA ; names <- c(names, "p")}

    #add p _num column
    if (!"p_num" %in% colnames(a)) {
      a$p_num <- a$p; names <- c(names, "p_num")}

    #select names
    a <- a %>% dplyr::select(all_of(names))

    #drop rows if all n, m, s are missing
    if (match != "yes") {
    a$del <- ifelse(rowSums(is.na(a[, c(names.n, names.m, names.s)])) != ncol(a[,c(names.n, names.m, names.s)]) |
                      rowSums(is.na(a[, names.m])) != ncol(a[, names.m]) |
                      rowSums(is.na(a[, names.n])) != ncol(a[, names.n]) |
                      rowSums(is.na(a[, names.s])) != ncol(a[, names.s]), 0,1)
    a <- a [a$del == 0, ]
    a$del <- NULL}

    #add rearranged columns
    if(g >1) {
    n.check <- as.data.frame(t(apply(a[,names.n], 1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )), stringsAsFactors = FALSE)
    m.check <- as.data.frame(t(apply(a[,names.m], 1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )), stringsAsFactors = FALSE)
    s.check <- as.data.frame(t(apply(a[,names.s], 1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} )), stringsAsFactors = FALSE)

    a[,names.m] <- m.check
    a[,names.n] <- n.check
    a[,names.s] <- s.check
    }

    a [, c(names.n, names.m, names.s)] <- dplyr::mutate_all(a [, c(names.n, names.m, names.s)], as.numeric)

    a$data_type <- NA
    a$data_type [1] <- c("pval_cont", "match", "cohort", "anova") [c(pval_cont,match,cohort,anova) == "y"]

    return(a)
  }

  #match
  match_import <- function () {
    a <- pval_import()
    a$p_num <- a$p <- NULL

    #drop if groups !=2
    a$grp <- rowSums(!is.na(a [, colnames(a) [substr(colnames(a),1,1) == "n"]]))
    a <- dplyr::left_join(a, a %>% dplyr::group_by(study) %>% dplyr:: select(study, grp) %>%
                            dplyr::summarise (gm = max(grp, na.rm = TRUE)),
                   by = "study")

    if (max(a$gm >2, na.rm=TRUE)) {
      if (vb == "y") {cat("Matching data- deleting any studies that have more than 2 groups\n")}
      a <- a %>% dplyr::filter (gm == 2)}

    if (1 %in% a$grp) {
      if (vb == "y") {cat("Matching data- deleting any variables that have only one group\n")}
      a <- a %>% dplyr::filter (!grp == 1)}

    #mark as NA if s1 or s2 == 0 or missing
    a$del <- ifelse(is.na(a$s1)|is.na(a$s2), 0, ifelse(a$s1 <= 0 | a$s2 <=0, 1, 0))

    if(1 %in% a$del) {
      if (vb == "y") {cat("Matching data- mark any variables that have a SD of 0 or less as NA")}
      a$s1 <- ifelse(is.na(a$s1) | a$s1 <=0, NA, a$s1)
      a$s2 <- ifelse(is.na(a$s2) | a$s2 <=0, NA, a$s2)}

    names <- colnames(a)
    names.m <- names [substring(names, 1, 1) == "m"]
    names.n <- names [substring(names, 1, 1) == "n"]
    names.s <- names [substring(names, 1, 1) == "s" & substring(names,2,2) %in% as.character(0:9)]
    names_del <- c(names.m [3:length(names.m)], names.n [3:length(names.n)], names.s [3:length(names.s)])

    if (max(a$gm == 2)) {a <- a %>% dplyr::select(-gm, -grp, -del)
    } else if (max(a$gm >2)) {a <- a %>% dplyr::select(-gm, -grp, -del, - all_of(names_del))}

    return (a)
  }

  #cohort
  cohort_import <- function () {
    if (tolower(substr(import,1,1)) == "n") {a <- get(file.cont)
    } else {a <- load_file(dir= dir.cont, range1 = range.name.cont, file = file.name.cont, sheet1 = sheet.name.cont)}

    #get names;
    colnames(a) <- lapply(colnames (a), tolower)

    colnames (a) [colnames(a) == "variable"] <- "var"
    colnames (a) [colnames(a) %in% c("g", "group", "grp")] <- "group"

    #apply lookup table if there
    if ("var_name_original" %in% colnames(a)) {
      if (!"var_name_final" %in% colnames(a)) {
        stop_mb("Cohort data: there is an original variable name (var_name_original) but no final variable name (var_name_final)")
      }
      var1 <- a$var_name_final[!is.na(a$var_name_final)]
      var1 <- gsub(" ","_",var1)
      var2 <- tolower(a$var_name_original[!is.na(a$var_name_original)])
      a <- dplyr::left_join(a %>% dplyr::mutate (var = tolower(var)), data.frame(var1 = var1,var = var2), by = "var")
      a$var <- a$var1
      a <- a %>% dplyr::select(-c(var_name_final, var_name_original, var1))
    }

    #keep population aside for later
    pop_i <-  0
    if ("population" %in% colnames(a)) {
      pop_i <- 1
      pop <- a %>% dplyr::select(study, var, population)
      if (tolower(format.cont) == "long") {
        pop <- pop %>% dplyr::group_by (study, var) %>% dplyr::slice (1) %>% dplyr::ungroup ()
      }
      pop$population [is.na(pop$population)] <- "NA"
    }
    #tidy up using pval_import
    a <- pval_import(x = a)
    a$p <- a$p_num <-  NULL #unwanted

    #add back in population variable
    if (pop_i == 1) {a <- dplyr::left_join (a, pop , by = c("study", "var"))
    #reorder final two variables (data.type and pop)
    nc <- ncol(a)
    a <- a [, c(1:(nc-2), nc, nc-1)]
    }
    return (a)
  }

  #anova
  anova_import <- function () {
    #put aside "d" for later
    if (tolower(substr(import,1,1)) == "n") {a <- get(file.cont)
    } else {a <- load_file(dir= dir.cont, range1 = range.name.cont, file = file.name.cont, sheet1 = sheet.name.cont)}

    #get names;
    colnames(a) <- lapply(colnames (a), tolower)
    colnames(a) [colnames(a) == "variable"] <- "var" #don't need an if statement here for this to work
    d_i = 0
    if("d" %in% colnames(a)) {
      d_i <- 1
      d <- a %>% dplyr::select(study, var, d)
      if (tolower(format.cont) == "long") {
        d <- d %>% dplyr::group_by (study, var, d) %>% dplyr::mutate(r = dplyr::row_number()) %>%
          dplyr::filter(r==1) %>% dplyr::ungroup ()
        d$r <- NULL
        if(nrow(d) > nrow(unique(cbind(a$study,a$var)))) {
            stop_mb("Anova data: some decimal places differ by group- they should be the same")
        }
      }
    }
      a <- pval_import(x = a)
      a$p <- a$p_num <-  NULL #unwanted

      #add back in d variable
      if (d_i == 1) {a <- dplyr::left_join (a, d , by = c("study", "var"))
      #reorder final two variables (data.type and pop)
      nc <- ncol(a)
      a <- a [, c(1:(nc-2), nc, nc-1)]
      }
      return (a)
    }


  #cat
  cat_import <- function () {

    if (substr(tolower(import),1,1) == "n") {a <- get(file.cat)
    } else {a <- load_file(dir1= dir.cat, range1 = range.name.cat, file = file.name.cat, sheet1 = sheet.name.cat)}

    colnames(a) <- lapply(colnames (a), tolower)
    colnames (a) [colnames(a) %in% c("g", "group", "grp")] <- "group"

    #long
    if (tolower(format.cat) == "long") {

      #check n column
      if (!"n" %in% colnames(a)) {stop_mb("For cat='yes', format='long' data, no 'n' column")}

      #get names
      a [, cat.names] <- a [, cat.names] %>% as.data.frame() %>%  dplyr::mutate_all(as.numeric)
      a <- a %>% dplyr::select(study, group, all_of(cat.names))

      #check rows correct
      if (length(cat.names) == 1 & nrow(a %>% dplyr::group_by(study, group) %>% dplyr::slice (1)) != nrow(a)) {
        stop_mb("For cat='yes', format='long' data, more than 1 n per study/group")}

      #make data wider, if only n - n_, else n_, v_
      a <- a %>% tidyr::pivot_wider(names_from = group, values_from = all_of(cat.names),
                             names_prefix = ifelse(length(cat.names) == 1, paste0(cat.names,"_"),""))
    }

    #get names;
    #remove separators
    names <- colnames(a) [substr(colnames(a),1,1) %in% cat.names]
    colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1), gsub("\\D", "", names))

    a <- a %>% dplyr::select(study, dplyr::all_of(names))

    #check equal number of groups
    names <- colnames(a %>% dplyr::select(-study))
    #avoid loop with map
    names.n <- purrr::map_int(cat.names, function (x) sum (x == substr(names,1,1)))

    #stop if not matching#
    if (length(unique(names.n)) != 1) {
      stop_mb ("Cat data error- some data for groups are missing")}

    g <- names.n [[1]] #number of groups

    #rearrange
    #get into standard format ie n1 ... nx; w1 ... wx; make missing data go to end;
    n <- purrr::map(cat.names, function (x) paste0(x,1:g))

    #drop rows if all N or all v are missing, or do not match
    #need to use recursive loop rather than map/walk- if use walk plus function need to << to env. and slower
    for (i in n) {
      chk <- rowSums(is.na(a[,i])) != ncol(a[,i])
      a <- a %>% dplyr::filter (chk)
    }

    #add rearranged columns
    n.check <- NULL
    for (i in 1:length(cat.names)) {
      n.1 <- a[, n [[i]]]
      n.check [[i]] <- t(apply(n.1,1, function(x) { return(c(x[!is.na(x)],x[is.na(x)]) )} ))  #this shifts non-missing cells to left
    }
    suppressMessages(n.check1 <- dplyr::bind_cols(n.check))
    colnames(n.check1) <- unlist (n)
    a <- data.frame ("study" = a$study, n.check1)

    #check n and cat columns match
    n.c1 <- rowSums(!is.na(a[, n [[1]]]))

    if( length(cat.names) >1) {
      for (i in 2:length(cat.names)) {
        cat.c1 <- rowSums(!is.na(a[, n [[i]]]))
        if (sum(cat.c1 != n.c1) >0 ) {
          stop_mb ("Cat data error- number of cat variables and groups for at least one study don't match")}
      }
    }

    #check number of grps match;
    a$group <- n.c1
    if (g != max(a$group)) {
      stop_mb ("Cat data error- problem with number of groups")}

    a$data_type <- NA
    a$data_type [1] <- "cat"

    return (a)
  }

  #sr
  sr_import <- function () {
    a <- cat_import()
    a$data_type <- NA
    a$data_type [1] <- "sr"
    return (a)
  }

  # cat_all
  cat_all_import <- function () {
    if (substr(tolower(import),1,1) == "n") {a <- get(file.cat)
    } else {a <- load_file(dir1= dir.cat, range1 = range.name.cat, file = file.name.cat, sheet1 = sheet.name.cat)}

    colnames (a) [tolower(colnames(a)) == "variable"] <- "var"
    colnames (a) [tolower(colnames(a)) %in% c("g", "group", "grp")] <- "group"
    colnames (a) [substr(colnames (a),1,1) == "N"] <- paste0("t",
                                                             substr(colnames (a) [substr(colnames (a),1,1) == "N"],2,nchar(colnames(a))))
    colnames(a) <- lapply(colnames (a), tolower)

    #keep stat till the end
    stat = 0

    if ("stat" %in% colnames(a)) {
      nm <- a %>% dplyr::select(study, var, stat) %>% dplyr::group_by (study, var) %>%
        dplyr::arrange (study, var) %>% dplyr::mutate (l1 = dplyr::row_number()) %>%
        tidyr::gather(k, v, stat) %>% tidyr::spread(l1,v) %>% dplyr::ungroup () %>% dplyr::select(-study, -var, -k)
      nm$n <- apply(nm, 1, function (x) length(unique(x [!is.na(x)]))>1)
      if (sum(nm$n)> 0 ) {
        stop_mb ("Cat_all data error- some stat tests don't match for different levels")}
    }

    if ("stat" %in% colnames(a)) {stat_df <- a %>% dplyr::select(study, var, stat) %>% dplyr::group_by (study,var) %>%
      dplyr::arrange (study, var) %>% slice (1)
    } else {stat <- 1}
    a$stat <- NULL

    #long
    if (tolower(format.cat) == "long") {

      #get names
      a [, c("n","t")] <-   a [, c("n","t")] %>%  dplyr::mutate_all(as.numeric)

      #make data wider, if only n - n_, else n_, v_
      a <- a %>% tidyr::pivot_wider(names_from = group, values_from = c("n","t"))
    }

    a[, substr(colnames(a),1,1) %in% c("n","t")] <- a[, substr(colnames(a),1,1) %in% c("n","t")] %>%
      dplyr::mutate_all(as.numeric)

    #remove and add temporary separator
    names <- colnames(a) [substr(colnames(a),1,1) %in% c("n","t")]
    colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1),"_", gsub("\\D", "", names))

    #add p if none
    if (!"p" %in% colnames (a)) {a$p <- NA}

    #get names;
    nm1 <- colnames(a) [substr(colnames(a),1,1) == "n"]
    nm2 <- colnames(a) [substr(colnames(a),1,1) == "t"]

    #check groups
    if (!"group" %in% colnames(a)) {a$group <- rowSums(!is.na( a [,nm1]))}

    #stop if not matching#
    if (length(nm1) != length(nm2) | sum(a$group != rowSums(!is.na( a [,nm2]))) >0 ) {
      stop_mb ("Cat_all data error- some data for groups are missing or group numbers don't match")}

    #number of levels
    if ("level" %in% colnames(a)) {

      #check missing levels
      if (sum(is.na(a$level)) > 0) {
        #if one level is missing, delete all, doesn't matter for 2 level variables.
        a <- dplyr::left_join (a,
                        a %>% dplyr::filter(is.na(level)) %>% dplyr::select(study, var) %>% dplyr::group_by(study, var) %>%
                          dplyr::slice (1) %>% dplyr::mutate (l1 = 99), by = c("study", "var")) %>% dplyr::ungroup ()
        a <- a %>% dplyr::mutate(level = ifelse(!is.na(l1), NA, level)) %>% dplyr::select (-l1)
        # now give them a new level
        a <- a %>% dplyr::mutate(x= rowSums(dplyr::across(all_of(nm1)), na.rm = TRUE)) %>%
          dplyr::arrange(study, var, desc(x)) %>%
          dplyr::group_by (study,var) %>% dplyr::mutate(l1 = dplyr::row_number()) %>%
          dplyr::mutate (level = ifelse(is.na(level), l1, level)) %>%
          dplyr::select(-x, -l1) %>% dplyr::ungroup ()
      }

    } else {
      a <- a %>% dplyr::mutate(x= rowSums(dplyr::across(all_of(nm1)), na.rm = TRUE)) %>%
        dplyr::arrange(study, var, desc(x)) %>%
        dplyr::group_by (study,var) %>% dplyr::mutate(level = dplyr::row_number()) %>%
        dplyr::select (-x) %>% dplyr::ungroup ()
    }

    a <- dplyr::left_join(a, a %>% dplyr::count(study, var) %>%  dplyr::rename(n_all = n), by = c("study", "var"))

    #create 2nd row for studies with only 1 line eg Sato
    if (sum(a$n_all ==1) >0) {
      nm11 <- paste0(substr(nm1,1,2),9,substr(nm1,3, nchar(nm1)))
      nm21 <- paste0(substr(nm2,1,2),9,substr(nm2,3, nchar(nm2)))
      chk <- a %>% dplyr::filter (n_all == 1)  %>% dplyr::select(study, var, starts_with("n"), starts_with("t"), level, - n_all)
      chk [, nm11] <- chk [, nm2] - chk [, nm1]
      chk [, nm21] <- chk [, nm2]

      chk <- chk %>%  tidyr::gather (k, v, -study, -var, -level) %>% dplyr::filter (!is.na(v)) %>%
        tidyr::separate (k, into = c("key1", "key2")) %>%
        dplyr::mutate(level = ifelse(as.numeric(key2)>90, 2, level)) %>% #if new var then level =2
        dplyr::mutate(key2 = ifelse(level ==2, substr(key2,2,nchar(key2)), key2)) %>%  #make group number 2
        tidyr::unite (k, key1, key2, sep ="_") %>%
        tidyr::spread (k,v) %>% dplyr::mutate (n_all =2, recode = level)

      chk <- dplyr::left_join (chk, a %>% dplyr::filter (n_all == 1) %>% dplyr::select(study, var, group, p), by = c("study","var"))
      #add relevant cols

      a <- dplyr::bind_rows(a %>% dplyr::filter (n_all > 1), chk)

    } else {chk <- a}

    #check data matches
    #calculate study/variable size
    chk_long <- a %>% dplyr::select(study, var, level, starts_with("n"), - n_all) %>%
      tidyr::gather(k,v,-study,-var, -level) %>% dplyr::group_by (study, var, k) %>%
      dplyr::summarise (x = sum(v, na.rm = TRUE), .groups = "keep") %>%
      dplyr::filter (x >0) %>% dplyr::filter (substr(k,1,1) == "n") %>%
      dplyr::mutate (k = paste0("x",k)) %>% tidyr::spread (k,x) %>% dplyr::ungroup ()

    chk1 <- dplyr::left_join(a, chk_long, by = c("study", "var")) %>%
      dplyr::select (study, var, level, starts_with("t", ignore.case = FALSE), starts_with("x", ignore.case = FALSE)) %>%
      tidyr::gather(k,v,-study,-var,-level) %>%
      tidyr::separate(k, into = c("k1", "k2"), sep = "_") %>% tidyr::spread (k1,v)

    if (!isTRUE(all.equal(chk1$t, chk1$xn))) {
      stop_mb ("Cat_all error- some data for levels are missing or don't match internally calculated values,
               or some values are not correct ie n1+n2 != N1 or N2")}

    #recodes
    if ("recode" %in% colnames(a)) {
      if (max(a$recode, na.rm = TRUE) >2 ) {
        stop_mb ("Cat_all error- recode levels are greater than 2- please recode to 1 or 2")}
    }

    #do recodes
    a <- a %>% dplyr::mutate (x = rowSums(dplyr::across(all_of(nm1)), na.rm = TRUE)) %>%
      dplyr::arrange (study, var, desc(x)) %>%
      dplyr::group_by (study, var) %>% dplyr::mutate (x1 = dplyr::row_number()) %>% dplyr::ungroup()

    a$t1 <- a$t2 <- a$ r <- 0
    for (i in 1:NROW(a)) {
      if (a$x1 [i] == 1) {a$t1 [i] <- a$x [i]; a$t2 [i] <- 0; a$r [i] <- 1}
      if (a$x1 [i] == 2) {a$t2 [i] <- a$x [i]; a$t1 [i] <- a$t1 [i-1]; a$r [i] <- 2}
      if (a$x1 [i] > 2) {
        if (a$t1 [i-1] > a$t2 [i-1]) {a$r [i] <- 2; a$t2 [i] <- a$t2 [i-1] + a$x [i]; a$t1 [i] <- a$t1 [i-1]
        } else {a$r [i] <- 1; a$t1 [i] <- a$t1 [i-1] + a$x [i]; a$t2 [i] <- a$t2 [i-1]}
      }
    }
    a <- a %>% dplyr::select(-x, -x1, -t1, -t2)

    if (! "recode" %in% colnames(a)) {a <- a %>% dplyr::rename(recode = r) }

    #check for missing recodes
    if (sum(is.na(a$recode)) > 0) {
      #if one level is missing, delete all, doesn't matter for 2 level variables.
      a <- dplyr::left_join (a,
                      a %>% dplyr::filter(is.na(recode)) %>% dplyr::select(study, var) %>% dplyr::group_by(study, var) %>%
                        dplyr::slice (1) %>% dplyr::mutate (r1 = 99), by = c("study", "var"))
      # now give them a new level
      a <- a %>% dplyr::mutate(recode = ifelse(!is.na(r1), r, recode)) %>% dplyr::select (-r1, -r)
    }
    a$r <- NULL

    #now arrange into standard format

    #remove temporary p-value
    if (sum(!is.na(a$p)) == 0 ) {a$p <- NULL}

    nm <- c("study","var","level","recode","n_all","group", nm1, nm2)
    if ("p" %in% colnames(a)) {nm <- c(nm, "p")}
    if ("stat" %in% colnames(a)) {nm <- c(nm, "stat")}
    a <- a [,nm]
    colnames(a) [5] <- "levels_no"
    names <- colnames(a) [substr(colnames(a),1,1) == "n"]
    colnames (a) [colnames (a) %in% names] <- paste0("n", gsub("\\D", "", names))
    names <- colnames(a) [substr(colnames(a),1,1) == "t"]
    colnames (a) [colnames (a) %in% names] <- paste0("N", gsub("\\D", "", names))


    #check p-vals/tests are correct
    if ("p" %in% colnames(a)) {
      nm <- a %>% dplyr::select(study:level, p) %>%  tidyr::gather(k, v, p) %>% tidyr::spread(level,v) %>%
        dplyr::ungroup () %>% dplyr::select(-study, -var, -k)
      nm$n <- apply(nm, 1, function (x) length(unique(x [!is.na(x)]))>1)

      if (sum(nm$n)> 0 ) {
        stop_mb ("Cat_all data error- some p-values don't match for different levels")}
    }

    #add p-value to each level
    if ("p" %in% colnames(a)) {
      nm <- a %>% dplyr::select(study, var, p) %>%  dplyr::arrange (study, var, p) %>% dplyr::group_by (study, var) %>%
        dplyr::slice(1) %>% dplyr::ungroup()
      a$p <- NULL
      a <- left_join(a, nm, by = c("study", "var"))}

    #add in stat
    if (stat == 0) {
      a <- left_join(a, stat_df, by = c("study", "var"))
    }

    a$data_type <- NA
    a$data_type [1] <- "cat_all"

    return (a)
  }

  #pval_cat
  pval_cat_import <- function () {
    a <- cat_all_import()
    a$data_type <- NA
    a$data_type [1] <- "pval_cat"
    a <- a %>% dplyr::select(-recode)
    return (a)
  }

  #gen
  gen_import <- function () {
    if (tolower(substr(import,1,1)) == "n") {a <- get(file.cont)
    } else {a <- load_file(dir= dir.cont, range1 = range.name.cont, file = file.name.cont, sheet1 = sheet.name.cont)}

    colnames(a) <- lapply(colnames (a), tolower)
    gen.vars.keep <- tolower(gen.vars.keep)
    gen.vars.del <- tolower(gen.vars.del)
    colnames (a) [colnames(a) == "variable"] <- "var"
    colnames (a) [colnames(a) %in% c("g", "group", "grp")] <- "group"

    #delete and keep names;
    if (sum(gen.vars.keep == "") == 0 & length(nchar(setdiff(gen.vars.keep, colnames(a)))) > 0) {
      stop_mb ("Error in gen.var.keep: some names not in file")}
    if (sum(gen.vars.del == "") == 0 & length(nchar(setdiff(gen.vars.del, colnames(a)))) > 0) {
      stop_mb ("Error in gen.var.del: some names not in file")}

    if (sum(gen.vars.keep == "") == 0)  {a <- a %>% dplyr::select(all_of(gen.vars.keep))}
    if (sum(gen.vars.del == "") == 0) {a <- a %>% dplyr::select (-all_of(gen.vars.del))}

    if (! "study" %in% colnames(a) | ! "var" %in% colnames(a)) {
      stop_mb ("file must have variables named study and variable or var")}

    #long
    if (tolower(format.cont) == "long") {

      #make wide format
      a$group <- as.numeric(a$group)
      a <- a %>% tidyr::gather(k, v, -study, -var, -group) %>% tidyr::unite( k1, k, group) %>% tidyr::spread (k1, v)
    }

    #remove separators
    names <- colnames(a) [! (colnames(a) == "study" | colnames(a) =="var")]
    colnames (a) [colnames (a) %in% names] <- paste0(substr(names,1,1), gsub("\\D", "", names))

    names <- setdiff(colnames (a), c("study", "var"))
    prefix <- gsub("\\d", "", names)
    pref <- unique(prefix)

    nm <- list()
    for (i in 1:length(pref)) {
      nm [[i]] <- paste0(pref[i],1:(sum(prefix == pref [i])))
    }
    nm <- c("study", "var", unlist(nm))

    #select names
    a <- a %>% dplyr::select(all_of(nm))

    return (a)
  }


# Code here ---------------------------------------------------------------

   #simplify indicators
   pval_cont <- substr(tolower(pval_cont),1,1)
   match <- substr(tolower(match),1,1)
   cohort <- substr(tolower(cohort),1,1)
   anova <- substr(tolower(anova),1,1)
   cat <- substr(tolower(cat),1,1)
   sr <-  substr(tolower(sr),1,1)
   cat_all <-  substr(tolower(cat_all),1,1)
   pval_cat <- substr(tolower(pval_cat),1,1)
   generic <-  substr(tolower(generic),1,1)

   if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
   } else {vb <- "n"}

   #check for only 1 cat and cont dataset
  if(sum(c(pval_cont, match, cohort, generic, anova) == "y") >1) {
    stop_mb("There can be only one cont data set loaded each run")}

   if(sum(c(cat, sr, cat_all, pval_cat) == "y") >1) {
    stop_mb("There can be only one cat data set loaded each run")}

   #check if only 1 dataset so don't need separate loadings
  if (dir !="") {dir.cont <- dir.cat <- dir}
  if (file.name != "") {file.name.cont <- file.name.cat <- file.name}

  #pval
  if(pval_cont == "y") {cont.data <- pval_import()}

  #match
  if(match == "y") {cont.data <- match_import()}

  #cohort
  if(cohort == "y") {cont.data <- cohort_import()}

  #cohort
  if(anova == "y") {cont.data <- anova_import()}

  #cat
  if (cat == "y" & length (cat.names) == 1) {
    stop_mb("Cat function error- there is only one variable: n and one other variable needed")}
  if(cat == "y") {cat.data <- cat_import()}

  #sr
  if (sr == "y" & !("n" %in% tolower(cat.names))) {
    stop_mb("SR function error- cat.names is not 'n'") }
  if(sr == "y") {cat.data <- sr_import()}

  #cat_all
  if (cat_all == "y") {cat.data <- cat_all_import()}

  #gen
  if (generic == "y") {generic.data <- gen_import()}

  #pval_cat
  if (pval_cat == "y") {cat.data <- pval_cat_import()}


  #return data
   results <- list()
   if (pval_cont == "y") {results = append(results, list(pval_cont_data = as.data.frame(cont.data)))}
   if (match == "y"){results = append(results, list(match_data = as.data.frame(cont.data)))}
   if (cohort == "y"){results = append(results, list(cohort_data = as.data.frame(cont.data)))}
   if (anova == "y"){results = append(results, list(anova_data = as.data.frame(cont.data)))}
   if (cat == "y") {results = append(results, list(cat_data = as.data.frame(cat.data)))}
   if (cat_all == "y") {results = append(results, list (cat_all_data = as.data.frame(cat.data)))}
   if (sr == "y") {results = append(results, list (sr_data = as.data.frame(cat.data)))}
   if (pval_cat == "y") {results = append(results, list(pval_cat_data = as.data.frame(cat.data)))}
   if (generic == "y") {results = append(results, list (generic_data = as.data.frame(generic.data)))}

   return (results)
  }

#global variables
utils::globalVariables(c(".", ".N", "Match", "No match", "SD", "Var1", "bpv", "data_type", "del", "df.groupsdf",
                         "df.groupsp", "diff_in_p_value", "dm", "dpm", "dps", "ds", "freq", "gm", "group", "grp",
                         "k", "k1", "k2", "key1", "key2", "l1", "leg", "level", "levels_no", "m", "m1", "m1_sf",
                         "m2", "m2_sf", "m_t", "m_t_comment", "m_wt", "max_dpm", "max_dps", "mean_SD", "mean_mode",
                         "mean_mx_m", "mean_n", "med", "med1", "med2", "mn", "mn_m", "mn_ms", "mn_s", "modemn",
                         "modesd", "ms", "ms2", "ms_m", "msd", "msd_mx_m", "n", "n1", "n_all", "n_var",
                         "name", "num", "num1", "oe", "p", "p.m", "p_df", "p_df_stat", "p_g", "p_n", "pc",
                         "pcc", "pchm", "pf", "plr", "pmdpsas", "population", "prob", "prop", "pval_cont_data",
                         "r", "r1", "recode", "recode1", "ref", "s", "s1", "s1_sf", "s2", "s2_sf", "s3", "s4", "s5", "s_m",
                         "s_wt", "sd", "sd_mx_m", "sd_n", "sfm", "sfs", "stat", "study", "t1", "t2", "title",
                         "tp", "v", "v1", "val", "val1", "val2", "val_1", "value", "var", "var1", "var_name_final",
                         "var_name_original", "variable", "x", "x1", "y", "cat_all_data", "cat_data", "cohort_data",
                         "rnorm", "match_data", "sr_data", "final", "generic_data", "pval_cat_data",
                         "N", "g", "g1", "p_num", "sim", "study_no", "var_no", "anova_data", "study.mean",
                         "study.sem", "study.sd", "dec", "round.sm", "study.sumsqdiff", "stat1",
                         "warn", "Chisquare.warn"))
