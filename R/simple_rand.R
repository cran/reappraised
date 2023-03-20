#' Compares observed and expected distribution of difference in numbers of participants between groups in two-arm randomised trials
#'
#' Creates plot of observed to expected numbers and ratios for differences in numbers of participants between trial groups\cr
#'
#' An example is for Sato and Iwamoto trials in Bolland 2016\cr
#' Bolland MJ, Avenell A, Gamble GD, Grey A. Systematic review and statistical analysis of the integrity of 33 randomized controlled trials. Neurology 2016;87:2391-2402.
#'
#' Returns a list containing 4 objects and (if verbose = TRUE) prints the plot sr_graph
#'
#' @param df data frame generated from load_clean function
#' @param br block randomisation: 'yes' or 'no'. If 'no' runs function as if all trials used simple randomisation. If 'yes' performs simple calculations as if block randomised
#' @param block an additional option for studies using block randomisation. block is a data frame containing columns named study (for study id); fbsz (for the final block size); n_fb (for number of participants in the final block); df (the difference between groups)
#' @param title title name for plots (optional)
#' @param verbose TRUE or FALSE indicates whether progress bar and comments in block randomisation function show and whether to print plot
#'
#' @return list containing 4 objects as described
#'
#'\itemize{
#'   \item sr_graph = plot of observed to expected numbers for differences between numbers of participants in trial groups
#'   \item sr_graph = plot of observed to expected numbers and ratios for differences between numbers of participants in trial groups
#'   \item sr_individual_graphs = list containing 2 plots making up composite figure
#'   \item sr_data = data frame containing data for plots
#'   }
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select count filter left_join bind_cols rename full_join
#' @importFrom ggpubr ggarrange
#' @importFrom purrr pmap map2
#' @export sr_fn
#
#' @examples
#' # load example data
#' sr_data <- load_clean(import= "no", file.cat = "SI_cat", sr= "yes",
#' format.cat = "wide")$sr_data
#'
#' \donttest{
#' # run function (takes only a few seconds)
#' sr_fn()$sr_graph
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
#' sr_data <- load_clean(import= "yes", sr = "yes", dir = path,
#'    file.name.cat = "reappraised_examples.xlsx", sheet.name.cat = "SI_cat",
#'    range.name.cat = "A:D", format.cat = "wide")$sr_data
#'
#' # function has an additional option for block randomisation.
#' #  If studies are block randomised and the final block size is known,
#' #  the number of participants in the final block can be determined.
#' #  The distribution of differences between groups for the final block
#' #  can be compared to the expected distribution
#' #
#' #  Few studies provide all these details so it seems unlikely this function
#' #  would get used often
#'
#' # Example takes only a few seconds to run
#'
#' sr_fn(br = "yes", block = data.frame(study = c(1,2,3,4,5,6,7,8,9,10),
#'          fb_sz= c(2,4,6,8,10,12,8,8,6,14), n_fb = c(1,1,4,5,7,8,4,6,2,10),
#'          df=c(1,1,0,1,3,4,2,2,0,0)))$sr_graph}
#'
#' @md

sr_fn <- function (df=sr_data, br = "no", block = data.frame(study="", fbsz="", n_fb= "", df="") , title= "", verbose = TRUE) {

  if (isTRUE(verbose) | tolower(substr(verbose,1,1)) == "t") {vb <- "y"
  } else {vb <- "n"}

  #block_rand function
  if (tolower(substr(br,1,1)) == "y") {
    if(max(block$fb_sz) > 20) {
      stop_mb("Block sizes >20 take >30s using this simple code, so upper threshold for block size is 20")}

    if(min(block$fb_sz) <1 | min (block$n_fb <1)) {
      stop_mb("Block sizes or number in final block <1")}

    results <- data.frame(df = 0:max(block$fb_sz))



    if (vb == "y") {pb = txtProgressBar(min=0, max=100, style = 2)}

    #create dataset of all combinations, delete unequal, calc probabilities
    block_exp_p <- function (study, fb_sz, n_fb) {
      perms <- expand.grid(as.data.frame(matrix(c("a","b"), 2, fb_sz)),stringsAsFactors = FALSE)
      #delete unequal
      perms$n1 <- apply(perms, 1, function(x) length(which(x=="a")))
      perms$n2 <- apply(perms %>% dplyr::select(-n1), 1, function(x) length(which(x=="b")))
      perms <- perms [perms$n1 == perms$n2, 1:n_fb , drop = FALSE]
      perms$n1 <- apply(perms, 1, function(x) length(which(x=="a")))
      perms$n2 <- apply(perms %>% dplyr::select(-n1), 1, function(x) length(which(x=="b")))
      perms$df <- abs(perms$n1-perms$n2)

      p <- perms %>% dplyr::count(df) %>% dplyr::mutate(exp = n/ sum(n)) %>% dplyr::select(-n)
      colnames (p) [colnames(p) == "exp"] <- paste0("exp_",study)
      p <- dplyr::left_join(results, p, by = "df") %>% replace(is.na(.), 0) %>% dplyr::select(-df)
      if (vb == "y") {setTxtProgressBar(pb, getTxtProgressBar (pb) + 100/nrow(block))}
      return(p)
    }

    results <- dplyr::bind_cols(results, purrr::pmap_dfc(block %>% dplyr::select (-df), block_exp_p))
    if (vb == "y") {close(pb)}

    results$exp <- rowSums (results %>% dplyr::select(-df))
    colnames(results) [1] <- "num"

    obs <- block %>% dplyr::count (df) %>% dplyr::rename (num = df, obs = n)

    b <- dplyr::left_join(results %>% dplyr::select(num, exp), obs, by = "num") %>% replace(is.na(.), 0)

    a <- block %>% dplyr::select(study)

    title <-paste0("Block rand.: ", title)
    }

    #import dataset
  if (tolower(substr(br,1,1)) != "y") {
    a <- df

    if (! "data_type" %in% colnames(a)) {a$data_type <- 0}
    if (is.na(a$data_type [1]) | a$data_type [1] != "sr") {
      stop_mb(paste0(deparse(substitute(df)), " is either incorrect data for this function or ",
                     deparse(substitute(df)), "$data_type != 'sr'"))
    } else if(nrow(a) ==1) {stop_mb("There is only one row of study data- the code will not work")}


    if (max(a$group) > 2) {
      if (vb == "y") {cat("Deleting studies with >2 groups")}
      a <- a %>% dplyr::filter(a$group < 3)
      nm.n <- colnames(a) [substr(colnames(a),1,1) == "n"]
      a [, nm.n [as.numeric(gsub("\\D","", nm.n)) >2]] <- NULL
    }

    #split into odd and even
    nm.n <- colnames(a) [substr(colnames(a),1,1) == "n"]

    a$type <- ifelse(rowSums(a[,nm.n])%%2 == 1, "Odd", "Even")
    a$s <- ifelse(a$type == "Odd", (rowSums(a [, nm.n])-1)/2, (rowSums(a [,nm.n]))/2)
    a$t <- ifelse(a$type == "Odd", rowSums(a [,nm.n]), rowSums(a [,nm.n]))

    odd <- a[a$type == "Odd",]
    even <- a[a$type == "Even",]

    #odd numbers
    if (NROW(odd) >0) {
      op <- purrr::map2(odd$s, odd$t, ~dbinom(.x - (seq(0,.x,by=1)), .y, 0.5)*2) %>% setNames(odd$study)
      op <- sapply(op, '[', 1:(max(odd$s)+1)) %>% as.data.frame ()
      op$n <- seq(1,max(odd$s)*2+1, by=2)
      op <- op %>% dplyr::select(n, everything())
    }

    #even numbers
    if (NROW(even)>0) {
      ep <- purrr::map2(even$s, even$t, ~c(dbinom(.x, .y , 0.5), dbinom(.x - (seq(1, .x, by=1)), .y, 0.5)*2)) %>%
        setNames(even$study)
      ep <- sapply(ep, '[', 1:(max(even$s)+1)) %>% as.data.frame ()
      ep$n <- seq(0,max(even$s)*2, by=2)
      ep <- ep %>% dplyr::select(n, everything())
    }

    #join odds and evens
    if (NROW(even) == 0) {exp <- op
    } else if (NROW(odd)== 0) {exp <- ep
    } else {exp <- dplyr::full_join(ep, op, by="n")}

    #sum expected
    exp <- dplyr::left_join(data.frame("n" = 0:max(exp$n)), exp, by= "n")
    expsum <- rowSums(exp[2:NCOL(exp)], na.rm = TRUE)

    #sum the observed differences
    a$obs = abs(a [, nm.n [[1]]] - a[ , nm.n [[2]]])
    obs <- a %>% dplyr::count(obs) %>% dplyr::rename (num= obs, obs= n)
    b <- dplyr::left_join(data.frame(num= exp$n, exp= expsum), obs, by = "num") %>% replace(is.na(.), 0)
  }

  #this is where block rejoins

#compress down
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
graph [1:2] <-  cat_graph (gph= b,  xtitle = "Differences between trial groups", ytitle = "Number of trials",
                           size = sz, sfx = "", text = list(length(unique(a$study))),
                           fn = "cat", ti = "y", top = "y", t= title)
rr_g <- list()
rr_g [1:2] <- rr_graph(rr = b ,xtitle = "Differences between trial groups", ytitle = "Observed/Expected ratio",
                       size =sz, sfx= "", text = list(length(unique(a$study))),
                       fn = "cat", ti = "y", top = "n", t= title)

#results
both_graphs <- ggpubr::ggarrange(graph [[1]], rr_g [[1]], ncol = 1, nrow= 2, align = "hv")

results <- list(sr_graph = graph [[1]], sr_both_graph = both_graphs,
                sr_individual_graphs= list(graph [[1]], rr_g [[1]]), sr_data = b)

if (vb == "y") {print(graph [[1]])}

return(results)
}






