# Scalable Record Linkage with FlexRL

``` r

library(FlexRL)
```

`FlexRL` is a package proposing flexible probabilistic **R**ecord
**L**inkage with diagnostic tools for downstream inference on linked
data. The record linkage algorithm
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/stEM.md) uses a
**St**ochastic **E**xpectation **M**aximisation approach to combine 2
data sources and outputs the set of records referring to the same
entities.

More details on the [record linkage
methodology](https://doi.org/10.1093/jrsssc/qlaf016); on the [false
discovery proportion estimation](https://doi.org/10.1002/sim.70292); and
on the potential risks of downstream inference on linked data (Robach et
al., in preparation).

Several other algorithms for record linkage are available in R. The
`link_with_*()` wrappers of FlexRL run them through one common
interface, so that the same diagnostics and false discovery proportion
(FDP) estimators apply to all of them:
[`link_with_FlexRL()`](https://robachowyk.github.io/FlexRL/reference/link_with_FlexRL.md),
[`link_with_BRL()`](https://robachowyk.github.io/FlexRL/reference/link_with_BRL.md),
[`link_with_fastLink()`](https://robachowyk.github.io/FlexRL/reference/link_with_fastLink.md),
[`link_with_multilink()`](https://robachowyk.github.io/FlexRL/reference/link_with_multilink.md),
[`link_with_reclin2()`](https://robachowyk.github.io/FlexRL/reference/link_with_reclin2.md),
[`link_with_diyar()`](https://robachowyk.github.io/FlexRL/reference/link_with_diyar.md),
[`link_with_fedmatch()`](https://robachowyk.github.io/FlexRL/reference/link_with_fedmatch.md).
Each wrapper takes the two encoded sources and a named list of arguments
forwarded to the underlying package, and returns a list with `idxA`,
`idxB` (row indices of the linked records in dataA and dataB) and
`LinkScore` (the linkage scores, `NULL` when the package does not return
any).

This vignette runs these wrappers on the same pipeline for: - a regional
subset of the
[NLTCS](https://www.icpsr.umich.edu/web/NACDA/studies/9681/versions/V5)
data (1982 and 1994 waves) - the
[SHIW](https://www.bancaditalia.it/statistiche/tematiche/indagini-famiglie-imprese/bilanci-famiglie/distribuzione-microdati/index.html)
data (2012 and 2014 waves) data with and without regional blocking

## One pipeline for all the wrappers

``` r

peakRAMNew <- function(...){
  # Capture R expressions or function calls
  args <- c(as.list(match.call(expand.dots = FALSE)$`...`), NULL)
  Function <- sapply(args, function(e) paste0(deparse(e), collapse = ""))
  Function <- sapply(Function, function(e) gsub(" ", "", e))
  # Initialize containers for output
  numfunc <- length(args)
  totaltime <- vector("numeric", numfunc)
  RAMused <- vector("numeric", numfunc)
  RAMpeak <- vector("numeric", numfunc)
  i <- 1
  for(arg in args){
    # Reset garbage collector and save baseline
    start <- gc(verbose = FALSE, reset = TRUE)
    # Evaluate regular and anonymous functions
    evalTime <- system.time(result <- eval.parent(arg))
    if(inherits(result, "function")){
      evalTime <- system.time(output <- result())
      rm(result)
    }else{
      output <- result
      rm(result)
    }
    # Call garbage collector and save post-eval
    end <- gc(verbose = FALSE, reset = FALSE)
    rm(output)
    # Calculate total and peak RAM used
    totaltime[i] <- as.numeric(evalTime["elapsed"])
    get_column_mb = function(object, before_column) {
      ind = which(colnames(object) == before_column)
      stopifnot(colnames(object)[ind+1] == "(Mb)")
      if (length(ind) == 0) {
        stop("No value found for that column")
      }
      # need to ref index and not colname because duplicated colnames for (Mb)
      object[2, ind + 1]
    }
    RAMused[i] <- get_column_mb(end, "used")- get_column_mb(start, "used")
    RAMpeak[i] <- get_column_mb(end, "max used") - get_column_mb(start, "max used")
    i <- i + 1
  }
  data.frame("Function_Call" = Function,
             "Elapsed_Time_sec" = totaltime,
             "Total_RAM_Used_MiB" = RAMused,
             "Peak_RAM_Used_MiB" = RAMpeak,
             row.names = 1:numfunc)
}

wrapper_args <- function(method, prep_data, ...) {
  PIVs <- names(prep_data$PIVs_config)
  args <- c(...)
    switch(method,
         FlexRL    = c(list(data = prep_data),
                       args[intersect(names(args), names(formals(FlexRL::StEM)))]),
         BRL       = c(list(flds = PIVs, types = rep("bi", length(PIVs)), nIter = 1000),
                       args[intersect(names(args), c(names(formals(BRL::compareRecords)), names(formals(BRL::bipartiteGibbs))))]),
         fastLink  = c(list(varnames = PIVs, tol.em = 1e-06, threshold.match = 0.5, return.all = TRUE, n.cores = 1),
                       args[intersect(names(args), names(formals(fastLink::fastLink)))]),
         multilink = c(list(types = rep("bi", length(PIVs)), breaks = stats::setNames(as.list(rep(NA, length(PIVs))), PIVs),
                            duplicates = c(0, 0), n_iter = 1000, max_cc_size = 50),
                       args[intersect(names(args), c(names(formals(multilink::create_comparison_data)),
                                                     names(formals(multilink::reduce_comparison_data)),
                                                     names(formals(multilink::specify_prior)),
                                                     names(formals(multilink::gibbs_sampler)),
                                                     names(formals(multilink::find_bayes_estimate)),
                                                     names(formals(multilink::relabel_bayes_estimate))))]),
         reclin2   = c(list(on = PIVs, formula = stats::as.formula(paste("~", paste(PIVs, collapse = "+"))),
                            type = "mpost", add = TRUE, variable = "selected", score = "mpost", threshold = 0.5),
                       args[intersect(names(args), c(names(formals(reclin2::pair)),
                                                     names(formals(reclin2::compare_pairs)),
                                                     names(formals(reclin2::problink_em)),
                                                     names(formals(stats::predict)),
                                                     names(formals(reclin2::select_threshold))))]),
         diyar     = c(list(probabilistic = TRUE, return_weights = TRUE),
                       args[intersect(names(args), c(names(formals(diyar::prob_score_range)),
                                                     names(formals(diyar::links_wf_probabilistic))))]),
         fedmatch  = c(list(by = PIVs, match_type = "multivar", unique_key_1 = "unique_key_A",
                            unique_key_2 = "unique_key_B", multivar_settings = fedmatch::build_multivar_settings(
                              compare_type = rep("indicator", length(PIVs)),
                              wgts = rep(1 / length(PIVs), length(PIVs)))),
                       args[intersect(names(args), names(formals(fedmatch::merge_plus)))])
  )
}

run_wrapper <- function(method, prep_data, arguments) {
  PIVs <- names(prep_data$PIVs_config)
  link_with <- get(paste0("link_with_", method), envir = asNamespace("FlexRL"))
  if (method == "diyar") {
    prep_data$encodedA <- prep_data$encodedA[, c(PIVs, "local_id", "source")]
    prep_data$encodedB <- prep_data$encodedB[, c(PIVs, "local_id", "source")]
  }
  tryCatch(link_with(prep_data$encodedA, prep_data$encodedB, arguments),
           error = function(e) {
             if (grepl("cannot allocate|memory|bad_alloc|too large|long vectors", conditionMessage(e), ignore.case = TRUE)) {
               message(method, ": memory error (", conditionMessage(e), "), results reported as NA")
               return(NULL)
             }
             stop(e)
           })
}

evaluate_linkage <- function(fit, true_pairs, threshold = 0.5) {
  if (is.null(fit)) return(rep(NA_real_, 5))
  delta_result <- data.frame(idxA = fit$idxA, idxB = fit$idxB)
  if (!is.null(fit$LinkScore)) delta_result <- delta_result[fit$LinkScore > threshold, , drop = FALSE]
  linked_pairs   <- do.call(paste, c(delta_result, list(sep = "_")))
  true_positive  <- length(intersect(linked_pairs, true_pairs))
  false_positive <- length(setdiff(linked_pairs, true_pairs))
  false_negative <- length(setdiff(true_pairs, linked_pairs))
  c(TP = true_positive, FP = false_positive, FN = false_negative,
    sensitivity = true_positive / (true_positive + false_negative),
    FDP = false_positive / (true_positive + false_positive))
}

results_rows <- c("TP               ", "FP               ", "FN               ", 
                  "sensitivity      ", "FDP              ", "est. FDP score   ", 
                  "est. FDP* score  ", "est. FDP* synth  ", "max SMD          ", 
                  "min support IoU  ", "max MMD          ", "time (secs)      ", 
                  "memory (Mb)      ")

run_all_wrappers <- function(prep_data, PIVs_type, methods, FDP_estimation = TRUE,
                             maxIter4CV = 1, n_repeats = 2, ...) {
  args <- list(...)
  op <- options()
  on.exit(options(op))
  PIVs <- names(prep_data$PIVs_config)
  true_pairs <- do.call(paste, c(prep_data$true_pairs, list(sep = "_")))
  df_results <- data.frame(matrix(NA, nrow = length(results_rows), ncol = 0))
  rownames(df_results) <- results_rows
  
  t0 <- Sys.time()
  MeM <- peakRAMNew(naive_fit <- naive_linkage(PIVs = PIVs, encodedA = prep_data$encodedA, encodedB = prep_data$encodedB))
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  df_results[1:5, "Naive"] <- evaluate_linkage(list(idxA = naive_fit[[1]], idxB = naive_fit[[2]]), true_pairs)
  df_results[12,"Naive"] <- elapsed
  df_results[13, "Naive"] <- MeM$Peak_RAM_Used_MiB

  for (method in methods) {
    arguments <- wrapper_args(method, prep_data, args)
    
    t0 <- Sys.time()
    MeM <- peakRAMNew(fit <- run_wrapper(method, prep_data, arguments))
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    df_results[1:5, method] <- evaluate_linkage(fit, true_pairs)
    if (is.null(fit) || length(fit$idxA) == 0) next
    
    encodedA <- prep_data$encodedA
    encodedB <- prep_data$encodedB
    if (method == "diyar") {
      encodedA <- encodedA[, c(PIVs, "local_id", "source")]
      encodedB <- encodedB[, c(PIVs, "local_id", "source")]
    }
    
    diag <- do.call(RL_diagnostics, c(list(fit, encodedA, encodedB, PIVs, PIVs_type, PIVs,
                                           true_pairs = prep_data$true_pairs, FDP_estimation = FDP_estimation,
                                           RL_method = method, data = prep_data, maxIter4CV = maxIter4CV, n_repeats = n_repeats),
                                      arguments))
    if (FDP_estimation) {
      df_results[6, method] <- diag$FDP_measures$curve$RL_FDP_score[1]
      df_results[7, method] <- diag$FDP_measures$curve$augmRL_FDP_score[1]
      df_results[8, method] <- diag$FDP_measures$curve$augmRL_FDP_synth[1]
    }
    df_results[9,  method] <- max(diag$discrepancy_measures[1, grep("smd", names(diag$discrepancy_measures))])
    df_results[10, method] <- min(diag$discrepancy_measures[1, grep("iou", names(diag$discrepancy_measures))])
    df_results[11, method] <- max(diag$discrepancy_measures[1, grep("mmd", names(diag$discrepancy_measures))])
    
    df_results[12, method] <- elapsed
    df_results[13, method] <- MeM$Peak_RAM_Used_MiB
  }
  df_results
}

show_results <- function(df_results) {
  cat("Record Linkage\n");      print(df_results[1:5, ],  2)
  cat("\nFDP estimation\n");    print(df_results[6:8, ],  2)
  cat("\nLinked divergence\n"); print(df_results[9:11, ], 2)
  cat("\nComputational time and memory\n"); print(df_results[12:13, ], 2)
}

methods <- c("BRL", "multilink", "fastLink", "reclin2", "diyar", "fedmatch", "FlexRL")
```

We will fill a table of results throughout the vignette, gathering for
each record linkage method: linkage performance, data quality
diagnostics: FDP estimation, measures of divergence of linked data from
sources.

## Regional subset of the NLTCS data

We use 4 **P**artially **I**dentifying **V**ariable**s** to link the
records: sex, birth month and birth year are `stable`, the state code
may change over time but we cannot model these changes so we set this
variable as `flexible`.

``` r

df1982 <- read.csv("ex-reg-NLTCS-82.csv", row.names = 1)
df1994 <- read.csv("ex-reg-NLTCS-94.csv", row.names = 1)

PIVs_config <- list(
  SEX       = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  DOB_MONTH = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  DOB_YEAR  = list(dynamics = "stable",   bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  STATE     = list(dynamics = "flexible", bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA)))

PIVs_type <- list(SEX = FALSE, DOB_MONTH = TRUE, DOB_YEAR = TRUE, STATE = FALSE)
```

[`prepare_data()`](https://robachowyk.github.io/FlexRL/reference/prepare_data.md)
removes the records for which the linking variables are outside of the
common support, checks that PIVs are not strongly associated with each
other (which may degrade record linkage performance), maps PIV values to
natural numbers and encodes missing values as 0. The StEM algorithm
denotes the smallest source as `A` and the largest as `B`. The true
pairs are known from the identifier ID (they are indexed on the rows of
`encodedA` and `encodedB`, after the support restriction).

``` r

prep_data <- prepare_data(df1982, df1994, "1982", "1994", PIVs_config, same_mistakes = TRUE, uniq_id = "ID")
#> '1982' is the larger source, saved as B; '1994' saved as A.

c(n_A = nrow(prep_data$encodedA), n_B = nrow(prep_data$encodedB), n_links = nrow(prep_data$true_pairs))
#>     n_A     n_B n_links 
#>     491     901     332

prep_data$n_values
#>       SEX DOB_MONTH  DOB_YEAR     STATE 
#>         2        12        28         9
```

``` r

df_results_NLTCS <- run_all_wrappers(prep_data, PIVs_type, methods)
```

``` r

show_results(df_results_NLTCS)
#> Record Linkage
#>                    Naive     BRL multilink fastLink reclin2   diyar fedmatch
#> TP                320.00 210.000    221.00   265.00  320.00 126.000   266.00
#> FP                163.00  18.000     26.00    78.00  163.00   5.000   223.00
#> FN                 12.00 122.000    111.00    67.00   12.00 206.000    66.00
#> sensitivity         0.96   0.633      0.67     0.80    0.96   0.380     0.80
#> FDP                 0.34   0.079      0.11     0.23    0.34   0.038     0.46
#>                   FlexRL
#> TP                229.00
#> FP                 27.00
#> FN                103.00
#> sensitivity         0.69
#> FDP                 0.11
#> 
#> FDP estimation
#>                   Naive  BRL multilink fastLink reclin2 diyar fedmatch FlexRL
#> est. FDP score       NA 0.15       NaN     0.30    0.49   NaN    0.075   0.16
#> est. FDP* score      NA 0.17       NaN     0.31      NA   NaN    0.073   0.18
#> est. FDP* synth      NA 0.12      0.17     0.38      NA     0      NaN   0.21
#> 
#> Linked divergence
#>                   Naive   BRL multilink fastLink reclin2 diyar fedmatch FlexRL
#> max SMD              NA 0.146     0.133    0.141   0.177  0.43   0.1125  0.170
#> min support IoU      NA 0.875     0.875    0.875   0.875  0.88   0.8750  0.875
#> max MMD              NA 0.014     0.012    0.016   0.025  0.03   0.0093  0.015
#> 
#> Computational time and memory
#>                   Naive   BRL multilink fastLink reclin2 diyar fedmatch FlexRL
#> time (secs)         0.7   5.6       168        4     0.8   3.6      3.3    189
#> memory (Mb)         2.3 414.6       825      834   239.7 834.6    157.8    834
```

## Full SHIW data

We use 5 **P**artially **I**dentifying **V**ariable**s** to link the
records: birth year, sex, birth region are `stable`, marital status and
occupation are `flexible`. The occupation is built by combining the
working status and the branch of activity together as one PIV. We create
such variable because the high association level between working status
and branch of activity may damage the record linkage performance.

``` r

df2012 <- comp[comp$ANNO == 2012, var_of_interest]
df2014 <- comp[comp$ANNO == 2014, var_of_interest]
df2012 <- df2012[!df2012$ID %in% df2012$ID[duplicated(df2012$ID)], ]
df2014 <- df2014[!df2014$ID %in% df2014$ID[duplicated(df2014$ID)], ]
rownames(df2012) <- 1:nrow(df2012)
rownames(df2014) <- 1:nrow(df2014)

PIVs_config <- list(
  ANASC          = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  SESSO          = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  STACIV         = list(dynamics = "flexible", bound_mistakes = c(NA, NA), fix_mistakes = c(NA, NA)),
  NASCREG        = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  SETTP11QUALP10 = list(dynamics = "flexible", bound_mistakes = c(NA, NA), fix_mistakes = c(NA, NA)))

PIVs_type <- list(ANASCI = TRUE, SESSO = FALSE, STACIV = FALSE, NASCREG = FALSE, SETTP11QUALP10 = FALSE)

PIVs_config_full <- c(PIVs_config,
  list(IREG = list(dynamics = "flexible", bound_mistakes = c(NA, NA), fix_mistakes = c(NA, NA))))
PIVs_type_full <- c(PIVs_type, IREG = FALSE)

prep_data_full <- prepare_data(df2012, df2014, "2012", "2014", PIVs_config_full, same_mistakes = TRUE, uniq_id = "ID")
#> '2012' is the larger source, saved as B; '2014' saved as A.
true_pairs_full <- do.call(paste, c(prep_data_full$true_pairs, list(sep = "_")))
```

``` r

df_results_full <- data.frame(matrix(NA, nrow = 6, ncol = 0))
rownames(df_results_full) <- c("TP", "FP", "FN", "sensitivity", "FDP", "minutes")
for (method in methods) {
  t0 <- Sys.time()       
  args <- list(StEM_iter = 10, StEM_burnin = 5, gibbs_iter = 10, gibbs_burnin = 5, 
               n_post_sample = 10, n_iter = 100, max_cc_size = 50, n.cores = 4)
  fit <- run_wrapper(method, prep_data_full, wrapper_args(method, prep_data_full, args))
  df_results_full[1:5, method] <- evaluate_linkage(fit, true_pairs_full)
  df_results_full[6, method] <- if (is.null(fit)) NA else round(as.numeric(Sys.time() - t0, units = "mins"))
}

# on Apple M1 Max ARM chip, 8-core CPU, 32 GB of memory

# we reduce the number of iterations of FlexRL to gain in computational time
```

``` r

cat("Record Linkage performance on large files")
#> Record Linkage performance on large files
knitr::kable(df_results_full, digits = 2, align=rep('c', ncol(df_results_full)))
```

|             | BRL | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL  |
|:------------|:---:|:---------:|:--------:|:-------:|:-----:|:--------:|:-------:|
| TP          | NA  |    NA     | 5082.00  |   NA    |  NA   | 5057.00  | 4598.00 |
| FP          | NA  |    NA     | 11154.00 |   NA    |  NA   | 14154.00 | 8080.00 |
| FN          | NA  |    NA     | 5486.00  |   NA    |  NA   | 5511.00  | 5970.00 |
| sensitivity | NA  |    NA     |   0.48   |   NA    |  NA   |   0.48   |  0.44   |
| FDP         | NA  |    NA     |   0.69   |   NA    |  NA   |   0.74   |  0.64   |
| minutes     | NA  |    NA     |   9.00   |   NA    |  NA   |   6.00   |  52.00  |

Note that when methods do not ensure the one-to-one assignment
constraint the number of False Positives may be high due to many-to-many
falsely linked pairs of records.

## Loop over the regions of the SHIW data with blocking

Blocking on the region code `IREG` splits the task into 20 smaller ones.
Within a region, `IREG` is constant and the linkage uses the 4 other
PIVs. Region 2 is the regional subset used above.

``` r

regions <- sort(unique(c(df2012$IREG, df2014$IREG)))
results_IREG <- list()
for (r in regions) {
  block2012 <- df2012[df2012$IREG == r, ]
  block2014 <- df2014[df2014$IREG == r, ]
  prep_data_r <- prepare_data(block2012, block2014, "2012", "2014", PIVs_config, same_mistakes = TRUE, uniq_id = "ID")
  true_pairs_r <- do.call(paste, c(prep_data_r$true_pairs, list(sep = "_")))
  for (method in methods) {
    t0 <- Sys.time()
    fit <- run_wrapper(method, prep_data_r, wrapper_args(method, prep_data_r))
    results_IREG[[length(results_IREG) + 1]] <- data.frame(
      IREG = r, n_A = nrow(prep_data_r$encodedA),
      n_B = nrow(prep_data_r$encodedB), n_links = length(true_pairs_r),
      method = method, sprintf("%.2f",t(evaluate_linkage(fit, true_pairs_r))),
      seconds = round(as.numeric(Sys.time() - t0, units = "secs")))
  }
}

results_IREG <- do.call(rbind, results_IREG)
```

``` r

cat("Record Linkage on regional subsets")
#> Record Linkage on regional subsets
knitr::kable(results_IREG, digits = 2, align=rep('c', ncol(results_IREG)))
```

| IREG | n_A | n_B | n_links | method | sprintf….2f…t.evaluate_linkage.fit..true_pairs_r… | seconds |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 1 | 1360 | 1495 | 645 | BRL | 63.00 | 17 |
| 1 | 1360 | 1495 | 645 | BRL | 24.00 | 17 |
| 1 | 1360 | 1495 | 645 | BRL | 582.00 | 17 |
| 1 | 1360 | 1495 | 645 | BRL | 0.10 | 17 |
| 1 | 1360 | 1495 | 645 | BRL | 0.28 | 17 |
| 1 | 1360 | 1495 | 645 | multilink | 153.00 | 1108 |
| 1 | 1360 | 1495 | 645 | multilink | 399.00 | 1108 |
| 1 | 1360 | 1495 | 645 | multilink | 492.00 | 1108 |
| 1 | 1360 | 1495 | 645 | multilink | 0.24 | 1108 |
| 1 | 1360 | 1495 | 645 | multilink | 0.72 | 1108 |
| 1 | 1360 | 1495 | 645 | fastLink | 321.00 | 2 |
| 1 | 1360 | 1495 | 645 | fastLink | 812.00 | 2 |
| 1 | 1360 | 1495 | 645 | fastLink | 324.00 | 2 |
| 1 | 1360 | 1495 | 645 | fastLink | 0.50 | 2 |
| 1 | 1360 | 1495 | 645 | fastLink | 0.72 | 2 |
| 1 | 1360 | 1495 | 645 | reclin2 | 632.00 | 0 |
| 1 | 1360 | 1495 | 645 | reclin2 | 25520.00 | 0 |
| 1 | 1360 | 1495 | 645 | reclin2 | 13.00 | 0 |
| 1 | 1360 | 1495 | 645 | reclin2 | 0.98 | 0 |
| 1 | 1360 | 1495 | 645 | reclin2 | 0.98 | 0 |
| 1 | 1360 | 1495 | 645 | diyar | 47.00 | 10 |
| 1 | 1360 | 1495 | 645 | diyar | 5.00 | 10 |
| 1 | 1360 | 1495 | 645 | diyar | 598.00 | 10 |
| 1 | 1360 | 1495 | 645 | diyar | 0.07 | 10 |
| 1 | 1360 | 1495 | 645 | diyar | 0.10 | 10 |
| 1 | 1360 | 1495 | 645 | fedmatch | 402.00 | 10 |
| 1 | 1360 | 1495 | 645 | fedmatch | 958.00 | 10 |
| 1 | 1360 | 1495 | 645 | fedmatch | 243.00 | 10 |
| 1 | 1360 | 1495 | 645 | fedmatch | 0.62 | 10 |
| 1 | 1360 | 1495 | 645 | fedmatch | 0.70 | 10 |
| 1 | 1360 | 1495 | 645 | FlexRL | 197.00 | 2305 |
| 1 | 1360 | 1495 | 645 | FlexRL | 119.00 | 2305 |
| 1 | 1360 | 1495 | 645 | FlexRL | 448.00 | 2305 |
| 1 | 1360 | 1495 | 645 | FlexRL | 0.31 | 2305 |
| 1 | 1360 | 1495 | 645 | FlexRL | 0.38 | 2305 |
| 2 | 44 | 53 | 24 | BRL | 21.00 | 0 |
| 2 | 44 | 53 | 24 | BRL | 11.00 | 0 |
| 2 | 44 | 53 | 24 | BRL | 3.00 | 0 |
| 2 | 44 | 53 | 24 | BRL | 0.88 | 0 |
| 2 | 44 | 53 | 24 | BRL | 0.34 | 0 |
| 2 | 44 | 53 | 24 | multilink | 23.00 | 1 |
| 2 | 44 | 53 | 24 | multilink | 10.00 | 1 |
| 2 | 44 | 53 | 24 | multilink | 1.00 | 1 |
| 2 | 44 | 53 | 24 | multilink | 0.96 | 1 |
| 2 | 44 | 53 | 24 | multilink | 0.30 | 1 |
| 2 | 44 | 53 | 24 | fastLink | 21.00 | 1 |
| 2 | 44 | 53 | 24 | fastLink | 14.00 | 1 |
| 2 | 44 | 53 | 24 | fastLink | 3.00 | 1 |
| 2 | 44 | 53 | 24 | fastLink | 0.88 | 1 |
| 2 | 44 | 53 | 24 | fastLink | 0.40 | 1 |
| 2 | 44 | 53 | 24 | reclin2 | 23.00 | 0 |
| 2 | 44 | 53 | 24 | reclin2 | 50.00 | 0 |
| 2 | 44 | 53 | 24 | reclin2 | 1.00 | 0 |
| 2 | 44 | 53 | 24 | reclin2 | 0.96 | 0 |
| 2 | 44 | 53 | 24 | reclin2 | 0.68 | 0 |
| 2 | 44 | 53 | 24 | diyar | 10.00 | 0 |
| 2 | 44 | 53 | 24 | diyar | 0.00 | 0 |
| 2 | 44 | 53 | 24 | diyar | 14.00 | 0 |
| 2 | 44 | 53 | 24 | diyar | 0.42 | 0 |
| 2 | 44 | 53 | 24 | diyar | 0.00 | 0 |
| 2 | 44 | 53 | 24 | fedmatch | 22.00 | 0 |
| 2 | 44 | 53 | 24 | fedmatch | 22.00 | 0 |
| 2 | 44 | 53 | 24 | fedmatch | 2.00 | 0 |
| 2 | 44 | 53 | 24 | fedmatch | 0.92 | 0 |
| 2 | 44 | 53 | 24 | fedmatch | 0.50 | 0 |
| 2 | 44 | 53 | 24 | FlexRL | 22.00 | 22 |
| 2 | 44 | 53 | 24 | FlexRL | 3.00 | 22 |
| 2 | 44 | 53 | 24 | FlexRL | 2.00 | 22 |
| 2 | 44 | 53 | 24 | FlexRL | 0.92 | 22 |
| 2 | 44 | 53 | 24 | FlexRL | 0.12 | 22 |
| 3 | 2103 | 2140 | 918 | BRL | 42.00 | 35 |
| 3 | 2103 | 2140 | 918 | BRL | 13.00 | 35 |
| 3 | 2103 | 2140 | 918 | BRL | 876.00 | 35 |
| 3 | 2103 | 2140 | 918 | BRL | 0.05 | 35 |
| 3 | 2103 | 2140 | 918 | BRL | 0.24 | 35 |
| 3 | 2103 | 2140 | 918 | multilink | 121.00 | 3069 |
| 3 | 2103 | 2140 | 918 | multilink | 585.00 | 3069 |
| 3 | 2103 | 2140 | 918 | multilink | 797.00 | 3069 |
| 3 | 2103 | 2140 | 918 | multilink | 0.13 | 3069 |
| 3 | 2103 | 2140 | 918 | multilink | 0.83 | 3069 |
| 3 | 2103 | 2140 | 918 | fastLink | 450.00 | 4 |
| 3 | 2103 | 2140 | 918 | fastLink | 1383.00 | 4 |
| 3 | 2103 | 2140 | 918 | fastLink | 468.00 | 4 |
| 3 | 2103 | 2140 | 918 | fastLink | 0.49 | 4 |
| 3 | 2103 | 2140 | 918 | fastLink | 0.75 | 4 |
| 3 | 2103 | 2140 | 918 | reclin2 | 911.00 | 1 |
| 3 | 2103 | 2140 | 918 | reclin2 | 54782.00 | 1 |
| 3 | 2103 | 2140 | 918 | reclin2 | 7.00 | 1 |
| 3 | 2103 | 2140 | 918 | reclin2 | 0.99 | 1 |
| 3 | 2103 | 2140 | 918 | reclin2 | 0.98 | 1 |
| 3 | 2103 | 2140 | 918 | diyar | 54.00 | 18 |
| 3 | 2103 | 2140 | 918 | diyar | 1.00 | 18 |
| 3 | 2103 | 2140 | 918 | diyar | 864.00 | 18 |
| 3 | 2103 | 2140 | 918 | diyar | 0.06 | 18 |
| 3 | 2103 | 2140 | 918 | diyar | 0.02 | 18 |
| 3 | 2103 | 2140 | 918 | fedmatch | 545.00 | 19 |
| 3 | 2103 | 2140 | 918 | fedmatch | 1558.00 | 19 |
| 3 | 2103 | 2140 | 918 | fedmatch | 373.00 | 19 |
| 3 | 2103 | 2140 | 918 | fedmatch | 0.59 | 19 |
| 3 | 2103 | 2140 | 918 | fedmatch | 0.74 | 19 |
| 3 | 2103 | 2140 | 918 | FlexRL | 280.00 | 3910 |
| 3 | 2103 | 2140 | 918 | FlexRL | 130.00 | 3910 |
| 3 | 2103 | 2140 | 918 | FlexRL | 638.00 | 3910 |
| 3 | 2103 | 2140 | 918 | FlexRL | 0.31 | 3910 |
| 3 | 2103 | 2140 | 918 | FlexRL | 0.32 | 3910 |
| 4 | 437 | 541 | 229 | BRL | 84.00 | 2 |
| 4 | 437 | 541 | 229 | BRL | 29.00 | 2 |
| 4 | 437 | 541 | 229 | BRL | 145.00 | 2 |
| 4 | 437 | 541 | 229 | BRL | 0.37 | 2 |
| 4 | 437 | 541 | 229 | BRL | 0.26 | 2 |
| 4 | 437 | 541 | 229 | multilink | 96.00 | 74 |
| 4 | 437 | 541 | 229 | multilink | 83.00 | 74 |
| 4 | 437 | 541 | 229 | multilink | 133.00 | 74 |
| 4 | 437 | 541 | 229 | multilink | 0.42 | 74 |
| 4 | 437 | 541 | 229 | multilink | 0.46 | 74 |
| 4 | 437 | 541 | 229 | fastLink | 141.00 | 3 |
| 4 | 437 | 541 | 229 | fastLink | 233.00 | 3 |
| 4 | 437 | 541 | 229 | fastLink | 88.00 | 3 |
| 4 | 437 | 541 | 229 | fastLink | 0.62 | 3 |
| 4 | 437 | 541 | 229 | fastLink | 0.62 | 3 |
| 4 | 437 | 541 | 229 | reclin2 | 228.00 | 0 |
| 4 | 437 | 541 | 229 | reclin2 | 16275.00 | 0 |
| 4 | 437 | 541 | 229 | reclin2 | 1.00 | 0 |
| 4 | 437 | 541 | 229 | reclin2 | 1.00 | 0 |
| 4 | 437 | 541 | 229 | reclin2 | 0.99 | 0 |
| 4 | 437 | 541 | 229 | diyar | 10.00 | 2 |
| 4 | 437 | 541 | 229 | diyar | 0.00 | 2 |
| 4 | 437 | 541 | 229 | diyar | 219.00 | 2 |
| 4 | 437 | 541 | 229 | diyar | 0.04 | 2 |
| 4 | 437 | 541 | 229 | diyar | 0.00 | 2 |
| 4 | 437 | 541 | 229 | fedmatch | 157.00 | 2 |
| 4 | 437 | 541 | 229 | fedmatch | 280.00 | 2 |
| 4 | 437 | 541 | 229 | fedmatch | 72.00 | 2 |
| 4 | 437 | 541 | 229 | fedmatch | 0.69 | 2 |
| 4 | 437 | 541 | 229 | fedmatch | 0.64 | 2 |
| 4 | 437 | 541 | 229 | FlexRL | 81.00 | 599 |
| 4 | 437 | 541 | 229 | FlexRL | 29.00 | 599 |
| 4 | 437 | 541 | 229 | FlexRL | 148.00 | 599 |
| 4 | 437 | 541 | 229 | FlexRL | 0.35 | 599 |
| 4 | 437 | 541 | 229 | FlexRL | 0.26 | 599 |
| 5 | 1181 | 1319 | 750 | BRL | 146.00 | 13 |
| 5 | 1181 | 1319 | 750 | BRL | 29.00 | 13 |
| 5 | 1181 | 1319 | 750 | BRL | 604.00 | 13 |
| 5 | 1181 | 1319 | 750 | BRL | 0.19 | 13 |
| 5 | 1181 | 1319 | 750 | BRL | 0.17 | 13 |
| 5 | 1181 | 1319 | 750 | multilink | 266.00 | 606 |
| 5 | 1181 | 1319 | 750 | multilink | 410.00 | 606 |
| 5 | 1181 | 1319 | 750 | multilink | 484.00 | 606 |
| 5 | 1181 | 1319 | 750 | multilink | 0.35 | 606 |
| 5 | 1181 | 1319 | 750 | multilink | 0.61 | 606 |
| 5 | 1181 | 1319 | 750 | fastLink | 365.00 | 22 |
| 5 | 1181 | 1319 | 750 | fastLink | 710.00 | 22 |
| 5 | 1181 | 1319 | 750 | fastLink | 385.00 | 22 |
| 5 | 1181 | 1319 | 750 | fastLink | 0.49 | 22 |
| 5 | 1181 | 1319 | 750 | fastLink | 0.66 | 22 |
| 5 | 1181 | 1319 | 750 | reclin2 | 598.00 | 1 |
| 5 | 1181 | 1319 | 750 | reclin2 | 341135.00 | 1 |
| 5 | 1181 | 1319 | 750 | reclin2 | 152.00 | 1 |
| 5 | 1181 | 1319 | 750 | reclin2 | 0.80 | 1 |
| 5 | 1181 | 1319 | 750 | reclin2 | 1.00 | 1 |
| 5 | 1181 | 1319 | 750 | diyar | 4.00 | 9 |
| 5 | 1181 | 1319 | 750 | diyar | 0.00 | 9 |
| 5 | 1181 | 1319 | 750 | diyar | 746.00 | 9 |
| 5 | 1181 | 1319 | 750 | diyar | 0.01 | 9 |
| 5 | 1181 | 1319 | 750 | diyar | 0.00 | 9 |
| 5 | 1181 | 1319 | 750 | fedmatch | 349.00 | 8 |
| 5 | 1181 | 1319 | 750 | fedmatch | 832.00 | 8 |
| 5 | 1181 | 1319 | 750 | fedmatch | 401.00 | 8 |
| 5 | 1181 | 1319 | 750 | fedmatch | 0.47 | 8 |
| 5 | 1181 | 1319 | 750 | fedmatch | 0.70 | 8 |
| 5 | 1181 | 1319 | 750 | FlexRL | 199.00 | 2309 |
| 5 | 1181 | 1319 | 750 | FlexRL | 86.00 | 2309 |
| 5 | 1181 | 1319 | 750 | FlexRL | 551.00 | 2309 |
| 5 | 1181 | 1319 | 750 | FlexRL | 0.27 | 2309 |
| 5 | 1181 | 1319 | 750 | FlexRL | 0.30 | 2309 |
| 6 | 459 | 513 | 238 | BRL | 96.00 | 2 |
| 6 | 459 | 513 | 238 | BRL | 36.00 | 2 |
| 6 | 459 | 513 | 238 | BRL | 142.00 | 2 |
| 6 | 459 | 513 | 238 | BRL | 0.40 | 2 |
| 6 | 459 | 513 | 238 | BRL | 0.27 | 2 |
| 6 | 459 | 513 | 238 | multilink | 122.00 | 70 |
| 6 | 459 | 513 | 238 | multilink | 84.00 | 70 |
| 6 | 459 | 513 | 238 | multilink | 116.00 | 70 |
| 6 | 459 | 513 | 238 | multilink | 0.51 | 70 |
| 6 | 459 | 513 | 238 | multilink | 0.41 | 70 |
| 6 | 459 | 513 | 238 | fastLink | 153.00 | 2 |
| 6 | 459 | 513 | 238 | fastLink | 202.00 | 2 |
| 6 | 459 | 513 | 238 | fastLink | 85.00 | 2 |
| 6 | 459 | 513 | 238 | fastLink | 0.64 | 2 |
| 6 | 459 | 513 | 238 | fastLink | 0.57 | 2 |
| 6 | 459 | 513 | 238 | reclin2 | 236.00 | 0 |
| 6 | 459 | 513 | 238 | reclin2 | 2752.00 | 0 |
| 6 | 459 | 513 | 238 | reclin2 | 2.00 | 0 |
| 6 | 459 | 513 | 238 | reclin2 | 0.99 | 0 |
| 6 | 459 | 513 | 238 | reclin2 | 0.92 | 0 |
| 6 | 459 | 513 | 238 | diyar | 6.00 | 2 |
| 6 | 459 | 513 | 238 | diyar | 0.00 | 2 |
| 6 | 459 | 513 | 238 | diyar | 232.00 | 2 |
| 6 | 459 | 513 | 238 | diyar | 0.03 | 2 |
| 6 | 459 | 513 | 238 | diyar | 0.00 | 2 |
| 6 | 459 | 513 | 238 | fedmatch | 169.00 | 3 |
| 6 | 459 | 513 | 238 | fedmatch | 289.00 | 3 |
| 6 | 459 | 513 | 238 | fedmatch | 69.00 | 3 |
| 6 | 459 | 513 | 238 | fedmatch | 0.71 | 3 |
| 6 | 459 | 513 | 238 | fedmatch | 0.63 | 3 |
| 6 | 459 | 513 | 238 | FlexRL | 93.00 | 602 |
| 6 | 459 | 513 | 238 | FlexRL | 52.00 | 602 |
| 6 | 459 | 513 | 238 | FlexRL | 145.00 | 602 |
| 6 | 459 | 513 | 238 | FlexRL | 0.39 | 602 |
| 6 | 459 | 513 | 238 | FlexRL | 0.36 | 602 |
| 7 | 670 | 710 | 321 | BRL | 79.00 | 4 |
| 7 | 670 | 710 | 321 | BRL | 31.00 | 4 |
| 7 | 670 | 710 | 321 | BRL | 242.00 | 4 |
| 7 | 670 | 710 | 321 | BRL | 0.25 | 4 |
| 7 | 670 | 710 | 321 | BRL | 0.28 | 4 |
| 7 | 670 | 710 | 321 | multilink | 132.00 | 175 |
| 7 | 670 | 710 | 321 | multilink | 151.00 | 175 |
| 7 | 670 | 710 | 321 | multilink | 189.00 | 175 |
| 7 | 670 | 710 | 321 | multilink | 0.41 | 175 |
| 7 | 670 | 710 | 321 | multilink | 0.53 | 175 |
| 7 | 670 | 710 | 321 | fastLink | 210.00 | 2 |
| 7 | 670 | 710 | 321 | fastLink | 337.00 | 2 |
| 7 | 670 | 710 | 321 | fastLink | 111.00 | 2 |
| 7 | 670 | 710 | 321 | fastLink | 0.65 | 2 |
| 7 | 670 | 710 | 321 | fastLink | 0.62 | 2 |
| 7 | 670 | 710 | 321 | reclin2 | 316.00 | 0 |
| 7 | 670 | 710 | 321 | reclin2 | 6320.00 | 0 |
| 7 | 670 | 710 | 321 | reclin2 | 5.00 | 0 |
| 7 | 670 | 710 | 321 | reclin2 | 0.98 | 0 |
| 7 | 670 | 710 | 321 | reclin2 | 0.95 | 0 |
| 7 | 670 | 710 | 321 | diyar | 34.00 | 4 |
| 7 | 670 | 710 | 321 | diyar | 2.00 | 4 |
| 7 | 670 | 710 | 321 | diyar | 287.00 | 4 |
| 7 | 670 | 710 | 321 | diyar | 0.11 | 4 |
| 7 | 670 | 710 | 321 | diyar | 0.06 | 4 |
| 7 | 670 | 710 | 321 | fedmatch | 227.00 | 3 |
| 7 | 670 | 710 | 321 | fedmatch | 443.00 | 3 |
| 7 | 670 | 710 | 321 | fedmatch | 94.00 | 3 |
| 7 | 670 | 710 | 321 | fedmatch | 0.71 | 3 |
| 7 | 670 | 710 | 321 | fedmatch | 0.66 | 3 |
| 7 | 670 | 710 | 321 | FlexRL | 120.00 | 895 |
| 7 | 670 | 710 | 321 | FlexRL | 59.00 | 895 |
| 7 | 670 | 710 | 321 | FlexRL | 201.00 | 895 |
| 7 | 670 | 710 | 321 | FlexRL | 0.37 | 895 |
| 7 | 670 | 710 | 321 | FlexRL | 0.33 | 895 |
| 8 | 1479 | 1524 | 898 | BRL | 144.00 | 17 |
| 8 | 1479 | 1524 | 898 | BRL | 25.00 | 17 |
| 8 | 1479 | 1524 | 898 | BRL | 754.00 | 17 |
| 8 | 1479 | 1524 | 898 | BRL | 0.16 | 17 |
| 8 | 1479 | 1524 | 898 | BRL | 0.15 | 17 |
| 8 | 1479 | 1524 | 898 | multilink | 290.00 | 1127 |
| 8 | 1479 | 1524 | 898 | multilink | 421.00 | 1127 |
| 8 | 1479 | 1524 | 898 | multilink | 608.00 | 1127 |
| 8 | 1479 | 1524 | 898 | multilink | 0.32 | 1127 |
| 8 | 1479 | 1524 | 898 | multilink | 0.59 | 1127 |
| 8 | 1479 | 1524 | 898 | fastLink | 452.00 | 34 |
| 8 | 1479 | 1524 | 898 | fastLink | 865.00 | 34 |
| 8 | 1479 | 1524 | 898 | fastLink | 446.00 | 34 |
| 8 | 1479 | 1524 | 898 | fastLink | 0.50 | 34 |
| 8 | 1479 | 1524 | 898 | fastLink | 0.66 | 34 |
| 8 | 1479 | 1524 | 898 | reclin2 | 721.00 | 1 |
| 8 | 1479 | 1524 | 898 | reclin2 | 339232.00 | 1 |
| 8 | 1479 | 1524 | 898 | reclin2 | 177.00 | 1 |
| 8 | 1479 | 1524 | 898 | reclin2 | 0.80 | 1 |
| 8 | 1479 | 1524 | 898 | reclin2 | 1.00 | 1 |
| 8 | 1479 | 1524 | 898 | diyar | 20.00 | 11 |
| 8 | 1479 | 1524 | 898 | diyar | 0.00 | 11 |
| 8 | 1479 | 1524 | 898 | diyar | 878.00 | 11 |
| 8 | 1479 | 1524 | 898 | diyar | 0.02 | 11 |
| 8 | 1479 | 1524 | 898 | diyar | 0.00 | 11 |
| 8 | 1479 | 1524 | 898 | fedmatch | 510.00 | 9 |
| 8 | 1479 | 1524 | 898 | fedmatch | 969.00 | 9 |
| 8 | 1479 | 1524 | 898 | fedmatch | 388.00 | 9 |
| 8 | 1479 | 1524 | 898 | fedmatch | 0.57 | 9 |
| 8 | 1479 | 1524 | 898 | fedmatch | 0.66 | 9 |
| 8 | 1479 | 1524 | 898 | FlexRL | 293.00 | 2778 |
| 8 | 1479 | 1524 | 898 | FlexRL | 104.00 | 2778 |
| 8 | 1479 | 1524 | 898 | FlexRL | 605.00 | 2778 |
| 8 | 1479 | 1524 | 898 | FlexRL | 0.33 | 2778 |
| 8 | 1479 | 1524 | 898 | FlexRL | 0.26 | 2778 |
| 9 | 1363 | 1500 | 711 | BRL | 43.00 | 18 |
| 9 | 1363 | 1500 | 711 | BRL | 13.00 | 18 |
| 9 | 1363 | 1500 | 711 | BRL | 668.00 | 18 |
| 9 | 1363 | 1500 | 711 | BRL | 0.06 | 18 |
| 9 | 1363 | 1500 | 711 | BRL | 0.23 | 18 |
| 9 | 1363 | 1500 | 711 | multilink | 107.00 | 1065 |
| 9 | 1363 | 1500 | 711 | multilink | 371.00 | 1065 |
| 9 | 1363 | 1500 | 711 | multilink | 604.00 | 1065 |
| 9 | 1363 | 1500 | 711 | multilink | 0.15 | 1065 |
| 9 | 1363 | 1500 | 711 | multilink | 0.78 | 1065 |
| 9 | 1363 | 1500 | 711 | fastLink | 343.00 | 34 |
| 9 | 1363 | 1500 | 711 | fastLink | 857.00 | 34 |
| 9 | 1363 | 1500 | 711 | fastLink | 368.00 | 34 |
| 9 | 1363 | 1500 | 711 | fastLink | 0.48 | 34 |
| 9 | 1363 | 1500 | 711 | fastLink | 0.71 | 34 |
| 9 | 1363 | 1500 | 711 | reclin2 | 8.00 | 1 |
| 9 | 1363 | 1500 | 711 | reclin2 | 188790.00 | 1 |
| 9 | 1363 | 1500 | 711 | reclin2 | 703.00 | 1 |
| 9 | 1363 | 1500 | 711 | reclin2 | 0.01 | 1 |
| 9 | 1363 | 1500 | 711 | reclin2 | 1.00 | 1 |
| 9 | 1363 | 1500 | 711 | diyar | 27.00 | 12 |
| 9 | 1363 | 1500 | 711 | diyar | 0.00 | 12 |
| 9 | 1363 | 1500 | 711 | diyar | 684.00 | 12 |
| 9 | 1363 | 1500 | 711 | diyar | 0.04 | 12 |
| 9 | 1363 | 1500 | 711 | diyar | 0.00 | 12 |
| 9 | 1363 | 1500 | 711 | fedmatch | 338.00 | 11 |
| 9 | 1363 | 1500 | 711 | fedmatch | 1025.00 | 11 |
| 9 | 1363 | 1500 | 711 | fedmatch | 373.00 | 11 |
| 9 | 1363 | 1500 | 711 | fedmatch | 0.48 | 11 |
| 9 | 1363 | 1500 | 711 | fedmatch | 0.75 | 11 |
| 9 | 1363 | 1500 | 711 | FlexRL | 194.00 | 2578 |
| 9 | 1363 | 1500 | 711 | FlexRL | 84.00 | 2578 |
| 9 | 1363 | 1500 | 711 | FlexRL | 517.00 | 2578 |
| 9 | 1363 | 1500 | 711 | FlexRL | 0.27 | 2578 |
| 9 | 1363 | 1500 | 711 | FlexRL | 0.30 | 2578 |
| 10 | 658 | 687 | 384 | BRL | 72.00 | 4 |
| 10 | 658 | 687 | 384 | BRL | 11.00 | 4 |
| 10 | 658 | 687 | 384 | BRL | 312.00 | 4 |
| 10 | 658 | 687 | 384 | BRL | 0.19 | 4 |
| 10 | 658 | 687 | 384 | BRL | 0.13 | 4 |
| 10 | 658 | 687 | 384 | multilink | 119.00 | 159 |
| 10 | 658 | 687 | 384 | multilink | 139.00 | 159 |
| 10 | 658 | 687 | 384 | multilink | 265.00 | 159 |
| 10 | 658 | 687 | 384 | multilink | 0.31 | 159 |
| 10 | 658 | 687 | 384 | multilink | 0.54 | 159 |
| 10 | 658 | 687 | 384 | fastLink | 206.00 | 6 |
| 10 | 658 | 687 | 384 | fastLink | 359.00 | 6 |
| 10 | 658 | 687 | 384 | fastLink | 178.00 | 6 |
| 10 | 658 | 687 | 384 | fastLink | 0.54 | 6 |
| 10 | 658 | 687 | 384 | fastLink | 0.64 | 6 |
| 10 | 658 | 687 | 384 | reclin2 | 321.00 | 0 |
| 10 | 658 | 687 | 384 | reclin2 | 102152.00 | 0 |
| 10 | 658 | 687 | 384 | reclin2 | 63.00 | 0 |
| 10 | 658 | 687 | 384 | reclin2 | 0.84 | 0 |
| 10 | 658 | 687 | 384 | reclin2 | 1.00 | 0 |
| 10 | 658 | 687 | 384 | diyar | 2.00 | 4 |
| 10 | 658 | 687 | 384 | diyar | 0.00 | 4 |
| 10 | 658 | 687 | 384 | diyar | 382.00 | 4 |
| 10 | 658 | 687 | 384 | diyar | 0.01 | 4 |
| 10 | 658 | 687 | 384 | diyar | 0.00 | 4 |
| 10 | 658 | 687 | 384 | fedmatch | 238.00 | 4 |
| 10 | 658 | 687 | 384 | fedmatch | 420.00 | 4 |
| 10 | 658 | 687 | 384 | fedmatch | 146.00 | 4 |
| 10 | 658 | 687 | 384 | fedmatch | 0.62 | 4 |
| 10 | 658 | 687 | 384 | fedmatch | 0.64 | 4 |
| 10 | 658 | 687 | 384 | FlexRL | 113.00 | 1004 |
| 10 | 658 | 687 | 384 | FlexRL | 47.00 | 1004 |
| 10 | 658 | 687 | 384 | FlexRL | 271.00 | 1004 |
| 10 | 658 | 687 | 384 | FlexRL | 0.29 | 1004 |
| 10 | 658 | 687 | 384 | FlexRL | 0.29 | 1004 |
| 11 | 762 | 872 | 486 | BRL | 113.00 | 5 |
| 11 | 762 | 872 | 486 | BRL | 33.00 | 5 |
| 11 | 762 | 872 | 486 | BRL | 373.00 | 5 |
| 11 | 762 | 872 | 486 | BRL | 0.23 | 5 |
| 11 | 762 | 872 | 486 | BRL | 0.23 | 5 |
| 11 | 762 | 872 | 486 | multilink | 182.00 | 216 |
| 11 | 762 | 872 | 486 | multilink | 213.00 | 216 |
| 11 | 762 | 872 | 486 | multilink | 304.00 | 216 |
| 11 | 762 | 872 | 486 | multilink | 0.37 | 216 |
| 11 | 762 | 872 | 486 | multilink | 0.54 | 216 |
| 11 | 762 | 872 | 486 | fastLink | 249.00 | 8 |
| 11 | 762 | 872 | 486 | fastLink | 422.00 | 8 |
| 11 | 762 | 872 | 486 | fastLink | 237.00 | 8 |
| 11 | 762 | 872 | 486 | fastLink | 0.51 | 8 |
| 11 | 762 | 872 | 486 | fastLink | 0.63 | 8 |
| 11 | 762 | 872 | 486 | reclin2 | 394.00 | 0 |
| 11 | 762 | 872 | 486 | reclin2 | 135933.00 | 0 |
| 11 | 762 | 872 | 486 | reclin2 | 92.00 | 0 |
| 11 | 762 | 872 | 486 | reclin2 | 0.81 | 0 |
| 11 | 762 | 872 | 486 | reclin2 | 1.00 | 0 |
| 11 | 762 | 872 | 486 | diyar | 4.00 | 5 |
| 11 | 762 | 872 | 486 | diyar | 0.00 | 5 |
| 11 | 762 | 872 | 486 | diyar | 482.00 | 5 |
| 11 | 762 | 872 | 486 | diyar | 0.01 | 5 |
| 11 | 762 | 872 | 486 | diyar | 0.00 | 5 |
| 11 | 762 | 872 | 486 | fedmatch | 251.00 | 5 |
| 11 | 762 | 872 | 486 | fedmatch | 511.00 | 5 |
| 11 | 762 | 872 | 486 | fedmatch | 235.00 | 5 |
| 11 | 762 | 872 | 486 | fedmatch | 0.52 | 5 |
| 11 | 762 | 872 | 486 | fedmatch | 0.67 | 5 |
| 11 | 762 | 872 | 486 | FlexRL | 149.00 | 1330 |
| 11 | 762 | 872 | 486 | FlexRL | 62.00 | 1330 |
| 11 | 762 | 872 | 486 | FlexRL | 337.00 | 1330 |
| 11 | 762 | 872 | 486 | FlexRL | 0.31 | 1330 |
| 11 | 762 | 872 | 486 | FlexRL | 0.29 | 1330 |
| 12 | 977 | 1142 | 455 | BRL | 39.00 | 10 |
| 12 | 977 | 1142 | 455 | BRL | 15.00 | 10 |
| 12 | 977 | 1142 | 455 | BRL | 416.00 | 10 |
| 12 | 977 | 1142 | 455 | BRL | 0.09 | 10 |
| 12 | 977 | 1142 | 455 | BRL | 0.28 | 10 |
| 12 | 977 | 1142 | 455 | multilink | 125.00 | 510 |
| 12 | 977 | 1142 | 455 | multilink | 235.00 | 510 |
| 12 | 977 | 1142 | 455 | multilink | 330.00 | 510 |
| 12 | 977 | 1142 | 455 | multilink | 0.27 | 510 |
| 12 | 977 | 1142 | 455 | multilink | 0.65 | 510 |
| 12 | 977 | 1142 | 455 | fastLink | 240.00 | 16 |
| 12 | 977 | 1142 | 455 | fastLink | 567.00 | 16 |
| 12 | 977 | 1142 | 455 | fastLink | 215.00 | 16 |
| 12 | 977 | 1142 | 455 | fastLink | 0.53 | 16 |
| 12 | 977 | 1142 | 455 | fastLink | 0.70 | 16 |
| 12 | 977 | 1142 | 455 | reclin2 | 375.00 | 1 |
| 12 | 977 | 1142 | 455 | reclin2 | 3911.00 | 1 |
| 12 | 977 | 1142 | 455 | reclin2 | 80.00 | 1 |
| 12 | 977 | 1142 | 455 | reclin2 | 0.82 | 1 |
| 12 | 977 | 1142 | 455 | reclin2 | 0.91 | 1 |
| 12 | 977 | 1142 | 455 | diyar | 6.00 | 8 |
| 12 | 977 | 1142 | 455 | diyar | 0.00 | 8 |
| 12 | 977 | 1142 | 455 | diyar | 449.00 | 8 |
| 12 | 977 | 1142 | 455 | diyar | 0.01 | 8 |
| 12 | 977 | 1142 | 455 | diyar | 0.00 | 8 |
| 12 | 977 | 1142 | 455 | fedmatch | 254.00 | 9 |
| 12 | 977 | 1142 | 455 | fedmatch | 722.00 | 9 |
| 12 | 977 | 1142 | 455 | fedmatch | 201.00 | 9 |
| 12 | 977 | 1142 | 455 | fedmatch | 0.56 | 9 |
| 12 | 977 | 1142 | 455 | fedmatch | 0.74 | 9 |
| 12 | 977 | 1142 | 455 | FlexRL | 148.00 | 1580 |
| 12 | 977 | 1142 | 455 | FlexRL | 61.00 | 1580 |
| 12 | 977 | 1142 | 455 | FlexRL | 307.00 | 1580 |
| 12 | 977 | 1142 | 455 | FlexRL | 0.33 | 1580 |
| 12 | 977 | 1142 | 455 | FlexRL | 0.29 | 1580 |
| 13 | 454 | 493 | 311 | BRL | 137.00 | 2 |
| 13 | 454 | 493 | 311 | BRL | 24.00 | 2 |
| 13 | 454 | 493 | 311 | BRL | 174.00 | 2 |
| 13 | 454 | 493 | 311 | BRL | 0.44 | 2 |
| 13 | 454 | 493 | 311 | BRL | 0.15 | 2 |
| 13 | 454 | 493 | 311 | multilink | 164.00 | 56 |
| 13 | 454 | 493 | 311 | multilink | 95.00 | 56 |
| 13 | 454 | 493 | 311 | multilink | 147.00 | 56 |
| 13 | 454 | 493 | 311 | multilink | 0.53 | 56 |
| 13 | 454 | 493 | 311 | multilink | 0.37 | 56 |
| 13 | 454 | 493 | 311 | fastLink | 197.00 | 3 |
| 13 | 454 | 493 | 311 | fastLink | 203.00 | 3 |
| 13 | 454 | 493 | 311 | fastLink | 114.00 | 3 |
| 13 | 454 | 493 | 311 | fastLink | 0.63 | 3 |
| 13 | 454 | 493 | 311 | fastLink | 0.51 | 3 |
| 13 | 454 | 493 | 311 | reclin2 | 310.00 | 0 |
| 13 | 454 | 493 | 311 | reclin2 | 17730.00 | 0 |
| 13 | 454 | 493 | 311 | reclin2 | 1.00 | 0 |
| 13 | 454 | 493 | 311 | reclin2 | 1.00 | 0 |
| 13 | 454 | 493 | 311 | reclin2 | 0.98 | 0 |
| 13 | 454 | 493 | 311 | diyar | 39.00 | 2 |
| 13 | 454 | 493 | 311 | diyar | 0.00 | 2 |
| 13 | 454 | 493 | 311 | diyar | 272.00 | 2 |
| 13 | 454 | 493 | 311 | diyar | 0.13 | 2 |
| 13 | 454 | 493 | 311 | diyar | 0.00 | 2 |
| 13 | 454 | 493 | 311 | fedmatch | 198.00 | 3 |
| 13 | 454 | 493 | 311 | fedmatch | 255.00 | 3 |
| 13 | 454 | 493 | 311 | fedmatch | 113.00 | 3 |
| 13 | 454 | 493 | 311 | fedmatch | 0.64 | 3 |
| 13 | 454 | 493 | 311 | fedmatch | 0.56 | 3 |
| 13 | 454 | 493 | 311 | FlexRL | 146.00 | 689 |
| 13 | 454 | 493 | 311 | FlexRL | 36.00 | 689 |
| 13 | 454 | 493 | 311 | FlexRL | 165.00 | 689 |
| 13 | 454 | 493 | 311 | FlexRL | 0.47 | 689 |
| 13 | 454 | 493 | 311 | FlexRL | 0.20 | 689 |
| 14 | 280 | 299 | 183 | BRL | 100.00 | 1 |
| 14 | 280 | 299 | 183 | BRL | 26.00 | 1 |
| 14 | 280 | 299 | 183 | BRL | 83.00 | 1 |
| 14 | 280 | 299 | 183 | BRL | 0.55 | 1 |
| 14 | 280 | 299 | 183 | BRL | 0.21 | 1 |
| 14 | 280 | 299 | 183 | multilink | 113.00 | 18 |
| 14 | 280 | 299 | 183 | multilink | 44.00 | 18 |
| 14 | 280 | 299 | 183 | multilink | 70.00 | 18 |
| 14 | 280 | 299 | 183 | multilink | 0.62 | 18 |
| 14 | 280 | 299 | 183 | multilink | 0.28 | 18 |
| 14 | 280 | 299 | 183 | fastLink | 131.00 | 1 |
| 14 | 280 | 299 | 183 | fastLink | 106.00 | 1 |
| 14 | 280 | 299 | 183 | fastLink | 52.00 | 1 |
| 14 | 280 | 299 | 183 | fastLink | 0.72 | 1 |
| 14 | 280 | 299 | 183 | fastLink | 0.45 | 1 |
| 14 | 280 | 299 | 183 | reclin2 | 182.00 | 0 |
| 14 | 280 | 299 | 183 | reclin2 | 7282.00 | 0 |
| 14 | 280 | 299 | 183 | reclin2 | 1.00 | 0 |
| 14 | 280 | 299 | 183 | reclin2 | 0.99 | 0 |
| 14 | 280 | 299 | 183 | reclin2 | 0.98 | 0 |
| 14 | 280 | 299 | 183 | diyar | 2.00 | 1 |
| 14 | 280 | 299 | 183 | diyar | 0.00 | 1 |
| 14 | 280 | 299 | 183 | diyar | 181.00 | 1 |
| 14 | 280 | 299 | 183 | diyar | 0.01 | 1 |
| 14 | 280 | 299 | 183 | diyar | 0.00 | 1 |
| 14 | 280 | 299 | 183 | fedmatch | 137.00 | 2 |
| 14 | 280 | 299 | 183 | fedmatch | 143.00 | 2 |
| 14 | 280 | 299 | 183 | fedmatch | 46.00 | 2 |
| 14 | 280 | 299 | 183 | fedmatch | 0.75 | 2 |
| 14 | 280 | 299 | 183 | fedmatch | 0.51 | 2 |
| 14 | 280 | 299 | 183 | FlexRL | 105.00 | 344 |
| 14 | 280 | 299 | 183 | FlexRL | 23.00 | 344 |
| 14 | 280 | 299 | 183 | FlexRL | 78.00 | 344 |
| 14 | 280 | 299 | 183 | FlexRL | 0.57 | 344 |
| 14 | 280 | 299 | 183 | FlexRL | 0.18 | 344 |
| 15 | 1939 | 2160 | 1069 | BRL | 23.00 | 35 |
| 15 | 1939 | 2160 | 1069 | BRL | 21.00 | 35 |
| 15 | 1939 | 2160 | 1069 | BRL | 1046.00 | 35 |
| 15 | 1939 | 2160 | 1069 | BRL | 0.02 | 35 |
| 15 | 1939 | 2160 | 1069 | BRL | 0.48 | 35 |
| 15 | 1939 | 2160 | 1069 | multilink | 69.00 | 2239 |
| 15 | 1939 | 2160 | 1069 | multilink | 648.00 | 2239 |
| 15 | 1939 | 2160 | 1069 | multilink | 1000.00 | 2239 |
| 15 | 1939 | 2160 | 1069 | multilink | 0.06 | 2239 |
| 15 | 1939 | 2160 | 1069 | multilink | 0.90 | 2239 |
| 15 | 1939 | 2160 | 1069 | fastLink | 276.00 | 65 |
| 15 | 1939 | 2160 | 1069 | fastLink | 1526.00 | 65 |
| 15 | 1939 | 2160 | 1069 | fastLink | 793.00 | 65 |
| 15 | 1939 | 2160 | 1069 | fastLink | 0.26 | 65 |
| 15 | 1939 | 2160 | 1069 | fastLink | 0.85 | 65 |
| 15 | 1939 | 2160 | 1069 | reclin2 | 1069.00 | 2 |
| 15 | 1939 | 2160 | 1069 | reclin2 | 1572964.00 | 2 |
| 15 | 1939 | 2160 | 1069 | reclin2 | 0.00 | 2 |
| 15 | 1939 | 2160 | 1069 | reclin2 | 1.00 | 2 |
| 15 | 1939 | 2160 | 1069 | reclin2 | 1.00 | 2 |
| 15 | 1939 | 2160 | 1069 | diyar | 1.00 | 19 |
| 15 | 1939 | 2160 | 1069 | diyar | 0.00 | 19 |
| 15 | 1939 | 2160 | 1069 | diyar | 1068.00 | 19 |
| 15 | 1939 | 2160 | 1069 | diyar | 0.00 | 19 |
| 15 | 1939 | 2160 | 1069 | diyar | 0.00 | 19 |
| 15 | 1939 | 2160 | 1069 | fedmatch | 280.00 | 20 |
| 15 | 1939 | 2160 | 1069 | fedmatch | 1659.00 | 20 |
| 15 | 1939 | 2160 | 1069 | fedmatch | 789.00 | 20 |
| 15 | 1939 | 2160 | 1069 | fedmatch | 0.26 | 20 |
| 15 | 1939 | 2160 | 1069 | fedmatch | 0.86 | 20 |
| 15 | 1939 | 2160 | 1069 | FlexRL | 113.00 | 4093 |
| 15 | 1939 | 2160 | 1069 | FlexRL | 124.00 | 4093 |
| 15 | 1939 | 2160 | 1069 | FlexRL | 956.00 | 4093 |
| 15 | 1939 | 2160 | 1069 | FlexRL | 0.11 | 4093 |
| 15 | 1939 | 2160 | 1069 | FlexRL | 0.52 | 4093 |
| 16 | 1130 | 1354 | 683 | BRL | 46.00 | 13 |
| 16 | 1130 | 1354 | 683 | BRL | 12.00 | 13 |
| 16 | 1130 | 1354 | 683 | BRL | 637.00 | 13 |
| 16 | 1130 | 1354 | 683 | BRL | 0.07 | 13 |
| 16 | 1130 | 1354 | 683 | BRL | 0.21 | 13 |
| 16 | 1130 | 1354 | 683 | multilink | 112.00 | 527 |
| 16 | 1130 | 1354 | 683 | multilink | 470.00 | 527 |
| 16 | 1130 | 1354 | 683 | multilink | 571.00 | 527 |
| 16 | 1130 | 1354 | 683 | multilink | 0.16 | 527 |
| 16 | 1130 | 1354 | 683 | multilink | 0.81 | 527 |
| 16 | 1130 | 1354 | 683 | fastLink | 266.00 | 23 |
| 16 | 1130 | 1354 | 683 | fastLink | 774.00 | 23 |
| 16 | 1130 | 1354 | 683 | fastLink | 417.00 | 23 |
| 16 | 1130 | 1354 | 683 | fastLink | 0.39 | 23 |
| 16 | 1130 | 1354 | 683 | fastLink | 0.74 | 23 |
| 16 | 1130 | 1354 | 683 | reclin2 | 681.00 | 1 |
| 16 | 1130 | 1354 | 683 | reclin2 | 404295.00 | 1 |
| 16 | 1130 | 1354 | 683 | reclin2 | 2.00 | 1 |
| 16 | 1130 | 1354 | 683 | reclin2 | 1.00 | 1 |
| 16 | 1130 | 1354 | 683 | reclin2 | 1.00 | 1 |
| 16 | 1130 | 1354 | 683 | diyar | 9.00 | 8 |
| 16 | 1130 | 1354 | 683 | diyar | 1.00 | 8 |
| 16 | 1130 | 1354 | 683 | diyar | 674.00 | 8 |
| 16 | 1130 | 1354 | 683 | diyar | 0.01 | 8 |
| 16 | 1130 | 1354 | 683 | diyar | 0.10 | 8 |
| 16 | 1130 | 1354 | 683 | fedmatch | 290.00 | 8 |
| 16 | 1130 | 1354 | 683 | fedmatch | 840.00 | 8 |
| 16 | 1130 | 1354 | 683 | fedmatch | 393.00 | 8 |
| 16 | 1130 | 1354 | 683 | fedmatch | 0.42 | 8 |
| 16 | 1130 | 1354 | 683 | fedmatch | 0.74 | 8 |
| 16 | 1130 | 1354 | 683 | FlexRL | 149.00 | 2293 |
| 16 | 1130 | 1354 | 683 | FlexRL | 68.00 | 2293 |
| 16 | 1130 | 1354 | 683 | FlexRL | 534.00 | 2293 |
| 16 | 1130 | 1354 | 683 | FlexRL | 0.22 | 2293 |
| 16 | 1130 | 1354 | 683 | FlexRL | 0.31 | 2293 |
| 17 | 264 | 276 | 193 | BRL | 123.00 | 1 |
| 17 | 264 | 276 | 193 | BRL | 32.00 | 1 |
| 17 | 264 | 276 | 193 | BRL | 70.00 | 1 |
| 17 | 264 | 276 | 193 | BRL | 0.64 | 1 |
| 17 | 264 | 276 | 193 | BRL | 0.21 | 1 |
| 17 | 264 | 276 | 193 | multilink | 122.00 | 13 |
| 17 | 264 | 276 | 193 | multilink | 78.00 | 13 |
| 17 | 264 | 276 | 193 | multilink | 71.00 | 13 |
| 17 | 264 | 276 | 193 | multilink | 0.63 | 13 |
| 17 | 264 | 276 | 193 | multilink | 0.39 | 13 |
| 17 | 264 | 276 | 193 | fastLink | 130.00 | 1 |
| 17 | 264 | 276 | 193 | fastLink | 100.00 | 1 |
| 17 | 264 | 276 | 193 | fastLink | 63.00 | 1 |
| 17 | 264 | 276 | 193 | fastLink | 0.67 | 1 |
| 17 | 264 | 276 | 193 | fastLink | 0.43 | 1 |
| 17 | 264 | 276 | 193 | reclin2 | 193.00 | 0 |
| 17 | 264 | 276 | 193 | reclin2 | 7266.00 | 0 |
| 17 | 264 | 276 | 193 | reclin2 | 0.00 | 0 |
| 17 | 264 | 276 | 193 | reclin2 | 1.00 | 0 |
| 17 | 264 | 276 | 193 | reclin2 | 0.97 | 0 |
| 17 | 264 | 276 | 193 | diyar | 16.00 | 1 |
| 17 | 264 | 276 | 193 | diyar | 1.00 | 1 |
| 17 | 264 | 276 | 193 | diyar | 177.00 | 1 |
| 17 | 264 | 276 | 193 | diyar | 0.08 | 1 |
| 17 | 264 | 276 | 193 | diyar | 0.06 | 1 |
| 17 | 264 | 276 | 193 | fedmatch | 127.00 | 2 |
| 17 | 264 | 276 | 193 | fedmatch | 137.00 | 2 |
| 17 | 264 | 276 | 193 | fedmatch | 66.00 | 2 |
| 17 | 264 | 276 | 193 | fedmatch | 0.66 | 2 |
| 17 | 264 | 276 | 193 | fedmatch | 0.52 | 2 |
| 17 | 264 | 276 | 193 | FlexRL | 122.00 | 376 |
| 17 | 264 | 276 | 193 | FlexRL | 44.00 | 376 |
| 17 | 264 | 276 | 193 | FlexRL | 71.00 | 376 |
| 17 | 264 | 276 | 193 | FlexRL | 0.63 | 376 |
| 17 | 264 | 276 | 193 | FlexRL | 0.27 | 376 |
| 18 | 501 | 554 | 362 | BRL | 90.00 | 2 |
| 18 | 501 | 554 | 362 | BRL | 29.00 | 2 |
| 18 | 501 | 554 | 362 | BRL | 272.00 | 2 |
| 18 | 501 | 554 | 362 | BRL | 0.25 | 2 |
| 18 | 501 | 554 | 362 | BRL | 0.24 | 2 |
| 18 | 501 | 554 | 362 | multilink | 151.00 | 60 |
| 18 | 501 | 554 | 362 | multilink | 131.00 | 60 |
| 18 | 501 | 554 | 362 | multilink | 211.00 | 60 |
| 18 | 501 | 554 | 362 | multilink | 0.42 | 60 |
| 18 | 501 | 554 | 362 | multilink | 0.46 | 60 |
| 18 | 501 | 554 | 362 | fastLink | 205.00 | 4 |
| 18 | 501 | 554 | 362 | fastLink | 246.00 | 4 |
| 18 | 501 | 554 | 362 | fastLink | 157.00 | 4 |
| 18 | 501 | 554 | 362 | fastLink | 0.57 | 4 |
| 18 | 501 | 554 | 362 | fastLink | 0.55 | 4 |
| 18 | 501 | 554 | 362 | reclin2 | 311.00 | 0 |
| 18 | 501 | 554 | 362 | reclin2 | 118067.00 | 0 |
| 18 | 501 | 554 | 362 | reclin2 | 51.00 | 0 |
| 18 | 501 | 554 | 362 | reclin2 | 0.86 | 0 |
| 18 | 501 | 554 | 362 | reclin2 | 1.00 | 0 |
| 18 | 501 | 554 | 362 | diyar | 7.00 | 3 |
| 18 | 501 | 554 | 362 | diyar | 0.00 | 3 |
| 18 | 501 | 554 | 362 | diyar | 355.00 | 3 |
| 18 | 501 | 554 | 362 | diyar | 0.02 | 3 |
| 18 | 501 | 554 | 362 | diyar | 0.00 | 3 |
| 18 | 501 | 554 | 362 | fedmatch | 194.00 | 3 |
| 18 | 501 | 554 | 362 | fedmatch | 307.00 | 3 |
| 18 | 501 | 554 | 362 | fedmatch | 168.00 | 3 |
| 18 | 501 | 554 | 362 | fedmatch | 0.54 | 3 |
| 18 | 501 | 554 | 362 | fedmatch | 0.61 | 3 |
| 18 | 501 | 554 | 362 | FlexRL | 117.00 | 912 |
| 18 | 501 | 554 | 362 | FlexRL | 61.00 | 912 |
| 18 | 501 | 554 | 362 | FlexRL | 245.00 | 912 |
| 18 | 501 | 554 | 362 | FlexRL | 0.32 | 912 |
| 18 | 501 | 554 | 362 | FlexRL | 0.34 | 912 |
| 19 | 1524 | 1625 | 1010 | BRL | 56.00 | 18 |
| 19 | 1524 | 1625 | 1010 | BRL | 8.00 | 18 |
| 19 | 1524 | 1625 | 1010 | BRL | 954.00 | 18 |
| 19 | 1524 | 1625 | 1010 | BRL | 0.06 | 18 |
| 19 | 1524 | 1625 | 1010 | BRL | 0.12 | 18 |
| 19 | 1524 | 1625 | 1010 | multilink | 138.00 | 1006 |
| 19 | 1524 | 1625 | 1010 | multilink | 574.00 | 1006 |
| 19 | 1524 | 1625 | 1010 | multilink | 872.00 | 1006 |
| 19 | 1524 | 1625 | 1010 | multilink | 0.14 | 1006 |
| 19 | 1524 | 1625 | 1010 | multilink | 0.81 | 1006 |
| 19 | 1524 | 1625 | 1010 | fastLink | 354.00 | 36 |
| 19 | 1524 | 1625 | 1010 | fastLink | 1057.00 | 36 |
| 19 | 1524 | 1625 | 1010 | fastLink | 656.00 | 36 |
| 19 | 1524 | 1625 | 1010 | fastLink | 0.35 | 36 |
| 19 | 1524 | 1625 | 1010 | fastLink | 0.75 | 36 |
| 19 | 1524 | 1625 | 1010 | reclin2 | 1006.00 | 1 |
| 19 | 1524 | 1625 | 1010 | reclin2 | 677113.00 | 1 |
| 19 | 1524 | 1625 | 1010 | reclin2 | 4.00 | 1 |
| 19 | 1524 | 1625 | 1010 | reclin2 | 1.00 | 1 |
| 19 | 1524 | 1625 | 1010 | reclin2 | 1.00 | 1 |
| 19 | 1524 | 1625 | 1010 | diyar | 5.00 | 14 |
| 19 | 1524 | 1625 | 1010 | diyar | 0.00 | 14 |
| 19 | 1524 | 1625 | 1010 | diyar | 1005.00 | 14 |
| 19 | 1524 | 1625 | 1010 | diyar | 0.00 | 14 |
| 19 | 1524 | 1625 | 1010 | diyar | 0.00 | 14 |
| 19 | 1524 | 1625 | 1010 | fedmatch | 351.00 | 11 |
| 19 | 1524 | 1625 | 1010 | fedmatch | 1173.00 | 11 |
| 19 | 1524 | 1625 | 1010 | fedmatch | 659.00 | 11 |
| 19 | 1524 | 1625 | 1010 | fedmatch | 0.35 | 11 |
| 19 | 1524 | 1625 | 1010 | fedmatch | 0.77 | 11 |
| 19 | 1524 | 1625 | 1010 | FlexRL | 162.00 | 3182 |
| 19 | 1524 | 1625 | 1010 | FlexRL | 68.00 | 3182 |
| 19 | 1524 | 1625 | 1010 | FlexRL | 848.00 | 3182 |
| 19 | 1524 | 1625 | 1010 | FlexRL | 0.16 | 3182 |
| 19 | 1524 | 1625 | 1010 | FlexRL | 0.30 | 3182 |
| 20 | 826 | 827 | 478 | BRL | 70.00 | 5 |
| 20 | 826 | 827 | 478 | BRL | 28.00 | 5 |
| 20 | 826 | 827 | 478 | BRL | 408.00 | 5 |
| 20 | 826 | 827 | 478 | BRL | 0.15 | 5 |
| 20 | 826 | 827 | 478 | BRL | 0.29 | 5 |
| 20 | 826 | 827 | 478 | multilink | 142.00 | 241 |
| 20 | 826 | 827 | 478 | multilink | 273.00 | 241 |
| 20 | 826 | 827 | 478 | multilink | 336.00 | 241 |
| 20 | 826 | 827 | 478 | multilink | 0.30 | 241 |
| 20 | 826 | 827 | 478 | multilink | 0.66 | 241 |
| 20 | 826 | 827 | 478 | fastLink | 232.00 | 9 |
| 20 | 826 | 827 | 478 | fastLink | 491.00 | 9 |
| 20 | 826 | 827 | 478 | fastLink | 246.00 | 9 |
| 20 | 826 | 827 | 478 | fastLink | 0.49 | 9 |
| 20 | 826 | 827 | 478 | fastLink | 0.68 | 9 |
| 20 | 826 | 827 | 478 | reclin2 | 478.00 | 0 |
| 20 | 826 | 827 | 478 | reclin2 | 71205.00 | 0 |
| 20 | 826 | 827 | 478 | reclin2 | 0.00 | 0 |
| 20 | 826 | 827 | 478 | reclin2 | 1.00 | 0 |
| 20 | 826 | 827 | 478 | reclin2 | 0.99 | 0 |
| 20 | 826 | 827 | 478 | diyar | 13.00 | 6 |
| 20 | 826 | 827 | 478 | diyar | 2.00 | 6 |
| 20 | 826 | 827 | 478 | diyar | 465.00 | 6 |
| 20 | 826 | 827 | 478 | diyar | 0.03 | 6 |
| 20 | 826 | 827 | 478 | diyar | 0.13 | 6 |
| 20 | 826 | 827 | 478 | fedmatch | 245.00 | 5 |
| 20 | 826 | 827 | 478 | fedmatch | 580.00 | 5 |
| 20 | 826 | 827 | 478 | fedmatch | 233.00 | 5 |
| 20 | 826 | 827 | 478 | fedmatch | 0.51 | 5 |
| 20 | 826 | 827 | 478 | fedmatch | 0.70 | 5 |
| 20 | 826 | 827 | 478 | FlexRL | 134.00 | 1455 |
| 20 | 826 | 827 | 478 | FlexRL | 61.00 | 1455 |
| 20 | 826 | 827 | 478 | FlexRL | 344.00 | 1455 |
| 20 | 826 | 827 | 478 | FlexRL | 0.28 | 1455 |
| 20 | 826 | 827 | 478 | FlexRL | 0.31 | 1455 |

Averaged over the regions and total computing time:

``` r

summary_IREG <- do.call(rbind, lapply(split(results_IREG, results_IREG$method), function(d) {
  data.frame(TP = mean(d$TP, na.rm=TRUE), FP = mean(d$FP, na.rm=TRUE), 
             FN = mean(d$FN, na.rm=TRUE), sensitivity = mean(d$sensitivity, na.rm=TRUE),
             FDP = mean(d$FDP, na.rm=TRUE), method = d$method[1],
             minutes = sum(d$seconds) / 60)
}))
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$TP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$FN, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
#> Warning in mean.default(d$sensitivity, na.rm = TRUE): argument is not numeric
#> or logical: returning NA
#> Warning in mean.default(d$FDP, na.rm = TRUE): argument is not numeric or
#> logical: returning NA
rownames(summary_IREG) <- summary_IREG$method
```

``` r

newdf <- as.data.frame(apply(t(summary_IREG)[c("TP","FP","FN","sensitivity","FDP","minutes"),methods], 2, function(x)  sprintf("%.2f",as.numeric(as.character(x)))))
names(newdf) <- names(df_results_full)
rownames(newdf) <- rownames(df_results_full)
knitr::kable(newdf, align=rep('c', ncol(newdf)))
```

|             |  BRL  | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL  |
|:------------|:-----:|:---------:|:--------:|:-------:|:-----:|:--------:|:-------:|
| TP          |  NA   |    NA     |    NA    |   NA    |  NA   |    NA    |   NA    |
| FP          |  NA   |    NA     |    NA    |   NA    |  NA   |    NA    |   NA    |
| FN          |  NA   |    NA     |    NA    |   NA    |  NA   |    NA    |   NA    |
| sensitivity |  NA   |    NA     |    NA    |   NA    |  NA   |    NA    |   NA    |
| FDP         |  NA   |    NA     |    NA    |   NA    |  NA   |    NA    |   NA    |
| minutes     | 17.00 |  1028.33  |  23.00   |  0.75   | 11.58 |  11.42   | 2771.33 |

``` r

colours <- grDevices::hcl.colors(length(methods), "Dark 3")
names(colours) <- methods
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "sensitivity", ylab = "FDP",
     main = "SHIW linked region by region (one point per region)")
for (method in methods) {
  d <- results_IREG[results_IREG$method == method, ]
  points(d$sensitivity, d$FDP, col = colours[method], pch = 19)
}
legend("bottomright", legend = methods, col = colours, pch = 19, bty = "n")
```

![plot of chunk
shiw-regional-blocking-plot](figure/shiw-regional-blocking-plot-1.png)

plot of chunk shiw-regional-blocking-plot

The 5 PIVs have little discriminating power on these data and the two
waves overlap on about 50% of their records (10.568 true links, 19.362
records in A, 20.016 records in B), so the task is hard for every
package. `BRL`, `multilink`, `reclin2` and `diyar` are not scalable (out
of memory errors).

We compare methods with blocking over the 20 regions: `reclin2` (EM on
the comparison patterns, threshold at 0.5 on the posterior) links a lot
of records reaching a high sensitivity at the price of a high FDP;
`fastLink`, `fedmatch` and `multilink` link about half of the true pairs
with a FDP around 0.7; `BRL` and `diyar` are the most conservative, with
a FDP around 0.3 and a sensitivity around 0.2; `FlexRL` performance is
between these conservative and permissive methods.

In terms of computational time however, `FlexRL` and `multilink` are
more than 10 times slower than the other methods. The results were
computed with the code above, on a machine with Apple M1 Max ARM chip,
8-core CPU, 32 GB of memory.
