#' Compares proportions of final digits from summary statistics
#'
#' Creates graph of proportion of final digits for summary statistics of specified variables\cr
#'
#' This approach is still in development and needs validation and discussion about its place in integrity assessment.\cr
#'
#' Requires data frame containing columns for study, variable (named var), summary statistic(s) (named with single letter eg m or s),
#'  and optional columns for decimal places for each statistic (named dp_* eg dp_m, dp_s). Data can be imported using the generic option of load_clean function\cr
#'
#'
#' Returns a list containing 5 objects and prints the plot digit_graph
#'
#' @param df data frame generated from load_clean function
#' @param vars vector of the summary statistics to be used
#' @param dec.pl "yes" or "no" indicator whether columns for decimal places are included (yes) or should be calculated (no)
#' @param dec.pl.vars vector of the names of the columns for decimal places for each statistics
#' @param title title name for plots (optional)
#' @param verbose TRUE or FALSE indicates whether print plot
#'
#'
#' @return list containing 5 objects as described
#'
#'\itemize{
#'   \item digit_graph = plot of proportions of final digits
#'   \item digit_ft = flextable of results
#'   \item digit_table = data frame of results
#'   \item digit_dataset = data frame of data set used to generate results data
#'   \item digit_data = results of analyses used to generate results data
#'   }
#'
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select arrange full_join rename
#' @importFrom purrr map_df
#' @importFrom flextable align border_remove hline_bottom font fontsize set_table_properties bold
#' @export final_digit_fn
#'
#'
#' @examples
#' # load example data
#' generic_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", generic= "yes",
#' gen.vars.del = c("p"), format.cont = "wide")$generic_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' final_digit_fn(vars = c("m","s"), dec.pl = "n")$digit_graph
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
#' generic_data <- load_clean(import= "yes", generic = "yes", dir = path,
#'     file.name.cont = "reappraised_examples.xlsx", sheet.name.cont = "SI_pvals_cont",
#'     range.name.cont = "A1:O51", gen.vars.del = c("p"),
#'     format.cont = "wide")$generic_data}
#' @md


#function for distribution of final digits ---------------------------------
final_digit_fn <- function (df= generic_data, vars = "", dec.pl = "no", dec.pl.vars ="", title="", verbose = TRUE) {

  a <- df

  #install_libs ()

  if(substr(tolower(dec.pl),1,1) == "y") {
    if (!all.equal(dec.pl.vars, paste0("dp_", vars)) | sum (dec.pl.vars %in% colnames(a)) != length(dec.pl.vars)) {
      stop_mb("The supplied dec.pl.vars does not match variable vars. Format should be 'dp_m', 'dp_s' etc in same order")}
    }

  names <- colnames(a) [substr(colnames(a),1,1) %in% vars & !colnames(a) %in% c("study", "var")]

  #select variables
  if(substr(tolower(dec.pl),1,1) == "y") {a <- a %>% dplyr::select(study, var, all_of(names), all_of(dec.pl.vars))
  } else {a <- a %>% dplyr::select(study, var, all_of(names))}

  if(substr(tolower(dec.pl),1,1) == "n") {dec.pl.vars =""}

  names <- setdiff(colnames (a), c("study", "var", dec.pl.vars))
  prefix <- gsub("\\d", "", names)

  #if no decimal places, calculate them
  if (substr(tolower(dec.pl),1,1) != "y") {

    #add digits
    dec <- function (x) {
      ifelse(!is.na(x),
             ifelse(abs(x - round(x)) > .Machine$double.eps^0.5,
                    nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(abs(x))))), 0)
             , NA)}

    dp <- purrr::map_df(a[, names], dec) %>% as.data.frame()
    a [, paste0("dp_",vars)] <- NA

    for (i in 1:length(vars)) {
      a [ , paste0("dp_",vars [i])] <- apply(dp [, substr(colnames(dp),1,1) == vars [i]], 1,
                                             function(x) ifelse(!all(is.na(x)), max(x, na.rm=TRUE), NA))}
  }

  #convert to long and check final digits match decimal places
  a [, vars] <- a [, paste0("dp_",vars)]

  rm(dp, i, dec)

  b <- a %>% dplyr::select (- all_of(paste0("dp_",vars))) %>%  gather(k, val, -study, -var, -all_of(vars), na.rm = TRUE) %>%
    dplyr::mutate(k1 = gsub("\\d", "", k))
  b$match = match(b$k1, names(b))
  b$dp <- as.numeric(b [cbind(seq_along(b$match), b$match)])
  b$digits = ifelse(b$dp == 0, b$val, sub(".*\\.","", b$val))
  b$digits0 <- ifelse(b$dp ==0, b$digits, ifelse(nchar(b$digits) != b$dp, paste0(b$digits,"0"), b$digits))
  b$final <- as.numeric(substr(b$digits0, nchar(b$digits0), nchar(b$digits0)))
  b$diff <- ifelse(b$dp ==0, "", ifelse(nchar(b$digits)!= b$dp, "Zero", ""))

  p <- as.data.frame(table (b$final))
  p$n <- as.numeric(as.character(p$Var1))
  p <- dplyr::full_join(p, data.frame (n = 0:9), by = "n") %>% dplyr::arrange(n) %>% replace(is.na(.), 0) %>% as.data.frame()
  p$Var1 <- factor(p$n)

  #graph for calculated p-values
  labels <- list(paste(length(unique(b$study)), "trials"),
                 paste(length(b$val), "values"))

  p.g <- digit_g (x = p, l= labels, t= title)

  pval <- p.g [[2]]
  p.g <- p.g [[1]]

  #create a flextable
  p$Percent = round(p$Freq/sum(p$Freq)*100,1)
  colnames(p) [2:3] <- c("Number", "Final digit")

  p1 <- p [, c(3,2,4)]

  #make flextables
  f <- flextable::flextable(p1)
  f <- flextable::font(f, fontname = "Arial")
  f <- flextable::fontsize(f, size = 8)
  f <- flextable::fontsize(f, part = "header", size = 10)
  f <- flextable::align(f, align = "center", part = "all")
  f <- flextable::bold(f, j=1, part= "body")
  f <- flextable::bold(f, part= "header")
  f <- flextable::border_remove(f)
  f <- flextable::hline_bottom(f, border= fp_border(color="black", width = 2), part = "header")
  f <- flextable::set_table_properties(f, width = 1, layout = "autofit")

  b <-  b %>% dplyr::rename(vars_grp = k, final_digit = final) %>% dplyr::select("study", "var", "vars_grp", "val", "dp", "final_digit")

  results <- list( digit_graph= p.g, digit_ft = f, digit_table = p1, digit_dataset = a, digit_data = b)

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {print (p.g)}

  return (results)

}





