## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, , echo = FALSE, results = "hide", message = FALSE-----------------
library(reappraised)

## -----------------------------------------------------------------------------
 pval_cont_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", pval_cont= "yes", format.cont = "wide")$pval_cont_data

## ---- eval = FALSE------------------------------------------------------------
#  # get path for example files
#  path <- system.file("extdata", "reappraised_examples.xlsx",
#               package = "reappraised", mustWork = TRUE)
#  # delete file name from path
#  path <- sub("/[^/]+$", "", path)
#  
#  # load data
#  pval_cont_data <- load_clean(import= "yes", pval_cont = "yes", dir = path,
#          file.name.cont = "reappraised_examples.xlsx",
#          sheet.name.cont = "SI_pvals_cont",
#          range.name.cont = "A1:O51", format.cont = "wide")$pval_cont_data

## ---- fig.width=4, fig.align='center', fig.asp = 0.7--------------------------
base_pval_cont <- pval_cont_fn(pval_cont_data, btsp= 500, verbose = FALSE)
base_pval_cont$all_results$pval_cont_calculated_pvalues

## ---- fig.width=4, fig.align='center', fig.asp = 0.7--------------------------
base_pval_cont$all_results$pval_cont_auc_calculated_pvalues

## ---- fig.width=4, fig.align='center', fig.asp = 0.7--------------------------
base_pval_cont$all_results$pval_cont_reported_pvalues

## ---- fig.width=4, fig.align='center', fig.asp = 0.7--------------------------
base_pval_cont$all_results$pval_cont_auc_reported_pvalues

## -----------------------------------------------------------------------------
match_data <- load_clean(import= "no", file.cont = "SI_pvals_cont", match= "yes", format.cont = "wide", verbose = FALSE)$match_data
match_fn(match_data, verbose= FALSE)$match_ft_all

## -----------------------------------------------------------------------------
cohort_data <- load_clean(import= "no", file.cont = "SI_cohort", cohort= "yes", format.cont = "long", verbose= FALSE)$cohort_data
cohort <- cohort_fn(cohort_data, seed = 10, sims = 100, verbose= FALSE)
cohort$cohort_ft

## ---- fig.width=9, fig.align='center', fig.height=10--------------------------
cohort$cohort_all_graphs$all_graph

## ---- fig.width=4, fig.align='center', fig.height=4, warning = FALSE----------

anova_data <- load_clean(import= "no", file.cont = "SI_pvals_cont",anova= "yes", format.cont = "wide")$anova_data
anova_fn(anova_data, seed = 10, sims = 100, btsp = 100, verbose = FALSE)$anova_ecdf

## ---- fig.width=7, fig.align='center', fig.height=8---------------------------
cat_data <- load_clean(import= "no", file.cat = "SI_cat", cat= "yes", format.cat = "wide", cat.names = c("n", "w"))$cat_data
cat_fn(cat_data, x_title="withdrawals", prefix="w", del.disparate="yes", verbose = FALSE)$cat_graph

## ---- fig.width=4, fig.align='center', fig.asp=.67, warning = FALSE-----------
sr_data <- load_clean(import= "no", file.cat = "SI_cat", sr= "yes", format.cat = "wide")$sr_data
sr_fn(sr_data, verbose = FALSE)$sr_graph

## ---- fig.width=4, fig.align='center', fig.asp=.67, warning = FALSE-----------
#create data set
block = data.frame(study = c(1,2,3,4,5,6,7,8,9,10), #study numbers
                   fb_sz= c(2,4,6,8,10,12,8,8,6,14), #final block size
                   n_fb = c(1,1,4,5,7,8,4,6,2,10), #number in final block
                   df=c(1,1,0,1,3,4,2,2,0,0)) #difference between groups


sr_fn(br = "yes", block = block, verbose = FALSE)$sr_graph

## ---- fig.width=7, fig.align='center', fig.height=7, warning = FALSE----------

cat_all_data <- load_clean(import= "No", file.cat = "SI_cat_all", cat_all= "yes", format.cat = "wide")$cat_all_data

cat_all_fn (cat_all_data, binom = "yes", two_levels = "y", del.disparate = "y", excl.level = "yes", seed = 10, verbose = FALSE)$cat_all_graph_pc

## -----------------------------------------------------------------------------
cat_all_fn (cat_all_data, comp.pvals = "yes", verbose= FALSE)$cat_all_diff_calc_rep_ft

## ---- fig.width=4, fig.align='center', fig.height=6, warning = FALSE----------

pval_cat_data <- load_clean(import= "no", file.cat = "SI_cat_all", pval_cat= "yes", format.cont = "wide", verbose = FALSE)$pval_cat_data
pval_cat_fn(pval_cat_data, sims = 10, seed=10, btsp = 100, verbose = FALSE)$pval_cat_calculated_pvalues

## ---- fig.width=4, fig.align='center', fig.asp=.7, warning = FALSE------------

generic_data <- load_clean(import= "No", file.cont = "SI_pvals_cont", generic= "yes", gen.vars.del = c("p"), format.cont = "wide")$generic_data
final_digit_fn(generic_data, vars = c("m","s"), dec.pl = "n", verbose = FALSE)$digit_graph 

