# Scalable Record Linkage with FlexRL

``` r

library(FlexRL)
#> If you are happy with FlexRL, please cite us! If you are unhappy, please cite us anyway (and feel free to complain as well).
```

`FlexRL` is a package proposing flexible probabilistic **R**ecord
**L**inkage with diagnostic tools for downstream inference on linked
data. The record linkage algorithm
[`StEM()`](https://robachowyk.github.io/FlexRL/reference/StEM.md) uses a
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
         FlexRL    = c(list(data = prep_data, threads = 4),
                       args[intersect(names(args), names(formals(FlexRL::StEM)))]),
         BRL       = c(list(flds = PIVs, types = rep("bi", length(PIVs)), nIter = 1000),
                       args[intersect(names(args), c(names(formals(BRL::compareRecords)), names(formals(BRL::bipartiteGibbs))))]),
         fastLink  = c(list(varnames = PIVs, tol.em = 1e-06, threshold.match = 0.5, return.all = TRUE, n.cores = 4),
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
                            unique_key_2 = "unique_key_B", nthread = 4, 
                            multivar_settings = fedmatch::build_multivar_settings(
                              compare_type = rep("indicator", length(PIVs)),
                              wgts = rep(1 / length(PIVs), length(PIVs)))),
                       args[intersect(names(args), names(formals(fedmatch::merge_plus)))])
  )
}

# unix::rlimit_as(20e9) # may be needed on some OS

run_wrapper <- function(method, prep_data, arguments) {
  PIVs <- names(prep_data$PIVs_config)
  link_with <- get(paste0("link_with_", method), envir = asNamespace("FlexRL"))
  if (method == "diyar") {
    prep_data$encodedA <- prep_data$encodedA[, c(PIVs, "local_id", "source")]
    prep_data$encodedB <- prep_data$encodedB[, c(PIVs, "local_id", "source")]
  }
  if (Sys.info()["sysname"] == "Darwin") {
    tryCatch(link_with(prep_data$encodedA, prep_data$encodedB, arguments),
           error = function(e) {
             if (grepl("cannot allocate|memory|bad_alloc|too large|long vectors", conditionMessage(e), ignore.case = TRUE)) {
               message(method, ": memory error (", conditionMessage(e), "), results reported as NA")
               return(NULL)
             }
             stop(e)
           })
  } else if (Sys.info()["sysname"] %in% c("Linux","Windows")) {
      tryCatch(callr::r(function(link_with, A, B, arguments) {link_with(A, B, arguments)},
                      args = list(link_with = link_with, A = prep_data$encodedA,
                                  B = prep_data$encodedB, arguments = arguments)),
           error = function(e) {
             if (grepl("cannot allocate|memory|bad_alloc|too large|long vectors|crashed or was killed|signal 9|SIGKILL",
               conditionMessage(e), ignore.case = TRUE)) {
               message(method, ": memory error (", conditionMessage(e), "), results reported as NA")
               return(NULL)
             }
             stop(e)
           })
  } else {
    print("Unknown Operating System.")
  }
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

methods <- c("BRL", "multilink", "fastLink", "reclin2", "diyar", "fedmatch", "FlexRL")

run_all_wrappers <- function(prep_data, variables_type, methods, FDP_estimation = TRUE,
                             maxIter4CV = 2, n_repeats = 3, ...) {
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
    
    diag <- do.call(RL_diagnostics, c(list(fit, encodedA, encodedB, PIVs, variables_type, PIVs,
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

cat("Record Linkage\n")
#> Record Linkage
knitr::kable(df_results_NLTCS[1:5, ], digits = 2, align=rep('c', ncol(df_results_NLTCS)))
```

|             | Naive  |  BRL   | multilink | fastLink | reclin2 | diyar  | fedmatch | FlexRL |
|:------------|:------:|:------:|:---------:|:--------:|:-------:|:------:|:--------:|:------:|
| TP          | 320.00 | 211.00 |  221.00   |  265.00  | 320.00  | 126.00 |  266.00  | 224.00 |
| FP          | 163.00 | 15.00  |   26.00   |  78.00   | 163.00  |  5.00  |  223.00  | 28.00  |
| FN          | 12.00  | 121.00 |  111.00   |  67.00   |  12.00  | 206.00 |  66.00   | 108.00 |
| sensitivity |  0.96  |  0.64  |   0.67    |   0.80   |  0.96   |  0.38  |   0.80   |  0.67  |
| FDP         |  0.34  |  0.07  |   0.11    |   0.23   |  0.34   |  0.04  |   0.46   |  0.11  |

``` r

cat("FDP estimation\n")
#> FDP estimation
knitr::kable(df_results_NLTCS[6:8, ], digits = 2, align=rep('c', ncol(df_results_NLTCS)))
```

|                  | Naive | BRL  | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL |
|:-----------------|:-----:|:----:|:---------:|:--------:|:-------:|:-----:|:--------:|:------:|
| est. FDP score   |  NA   | 0.15 |    NaN    |   0.30   |  0.49   |  NaN  |   0.07   |  0.16  |
| est. FDP\* score |  NA   | 0.17 |    NaN    |   0.31   |   NA    |  NaN  |   0.07   |  0.18  |
| est. FDP\* synth |  NA   | 0.08 |   0.22    |   0.44   |   NA    | 0.03  |   0.13   |  0.16  |

``` r

cat("Linked divergence\n")
#> Linked divergence
knitr::kable(df_results_NLTCS[9:11, ], digits = 2, align=rep('c', ncol(df_results_NLTCS)))
```

|                 | Naive | BRL  | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL |
|:----------------|:-----:|:----:|:---------:|:--------:|:-------:|:-----:|:--------:|:------:|
| max SMD         |  NA   | 0.15 |   0.13    |   0.14   |  0.18   | 0.43  |   0.11   |  0.18  |
| min support IoU |  NA   | 0.88 |   0.88    |   0.88   |  0.88   | 0.88  |   0.88   |  0.88  |
| max MMD         |  NA   | 0.01 |   0.01    |   0.02   |  0.02   | 0.03  |   0.01   |  0.02  |

``` r

cat("Computational time and memory\n")
#> Computational time and memory
knitr::kable(df_results_NLTCS[12:13, ], digits = 2, align=rep('c', ncol(df_results_NLTCS)))
```

|             | Naive | BRL  | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL |
|:------------|:-----:|:----:|:---------:|:--------:|:-------:|:-----:|:--------:|:------:|
| time (secs) | 0.23  | 5.8  |    165    |   4.2    |  0.36   |   3   |   2.6    |  7.3   |
| memory (Mb) | 2.50  | 85.0 |    385    |  273.1   | 155.50  |  117  |  105.5   | 109.3  |

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

df_results_full <- data.frame(matrix(NA, nrow = 7, ncol = 0))
rownames(df_results_full) <- c("TP", "FP", "FN", "sensitivity", "FDP", "time (mins)", "memory (Mb)")
for (method in methods) {
  t0 <- Sys.time()           
  args <- list(n_iter = 100)
  MeM <- peakRAMNew(fit <- run_wrapper(method, prep_data_full, wrapper_args(method, prep_data_full, args)))
  df_results_full[1:5, method] <- evaluate_linkage(fit, true_pairs_full)
  df_results_full[6, method] <- if (is.null(fit)) NA else round(as.numeric(Sys.time() - t0, units = "mins"))
  df_results_full[7, method] <- if (is.null(fit)) NA else MeM$Peak_RAM_Used_MiB
}
```

``` r

cat("Record Linkage performance on large files")
#> Record Linkage performance on large files
knitr::kable(df_results_full, digits = 2, align=rep('c', ncol(df_results_full)))
```

|             | BRL | multilink | fastLink | reclin2 | diyar | fedmatch | FlexRL  |
|:------------|:---:|:---------:|:--------:|:-------:|:-----:|:--------:|:-------:|
| TP          | NA  |    NA     | 5001.00  |   NA    |  NA   | 5057.00  | 4118.00 |
| FP          | NA  |    NA     | 11235.00 |   NA    |  NA   | 14154.00 | 5552.00 |
| FN          | NA  |    NA     | 5567.00  |   NA    |  NA   | 5511.00  | 6450.00 |
| sensitivity | NA  |    NA     |   0.47   |   NA    |  NA   |   0.48   |  0.39   |
| FDP         | NA  |    NA     |   0.69   |   NA    |  NA   |   0.74   |  0.57   |
| time (mins) | NA  |    NA     |   8.00   |   NA    |  NA   |   5.00   |  6.00   |
| memory (Mb) | NA  |    NA     | 25918.30 |   NA    |  NA   | 1845.40  | 4033.60 |

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
    MeM <- peakRAMNew(fit <- run_wrapper(method, prep_data_r, wrapper_args(method, prep_data_r)))
    results_IREG[[length(results_IREG) + 1]] <- data.frame(
      IREG = r, n_A = nrow(prep_data_r$encodedA),
      n_B = nrow(prep_data_r$encodedB), n_links = length(true_pairs_r),
      method = method, round(t(evaluate_linkage(fit, true_pairs_r)),2),
      seconds = round(as.numeric(Sys.time() - t0, units = "secs")),
      memory = MeM$Peak_RAM_Used_MiB)
  }
}

results_IREG <- do.call(rbind, results_IREG)
```

``` r

cat("Record Linkage on regional subsets")
#> Record Linkage on regional subsets
knitr::kable(results_IREG, digits = 2, align=rep('c', ncol(results_IREG)))
```

| IREG | n_A  | n_B  | n_links |  method   |  TP  |   FP    |  FN  | sensitivity | FDP  | seconds | memory |
|:----:|:----:|:----:|:-------:|:---------:|:----:|:-------:|:----:|:-----------:|:----:|:-------:|:------:|
|  1   | 1360 | 1495 |   645   |    BRL    |  61  |   25    | 584  |    0.09     | 0.29 |   19    | 1886.3 |
|  1   | 1360 | 1495 |   645   | multilink | 153  |   399   | 492  |    0.24     | 0.72 |  1112   | 1651.5 |
|  1   | 1360 | 1495 |   645   | fastLink  | 321  |   812   | 324  |    0.50     | 0.72 |    4    | 182.4  |
|  1   | 1360 | 1495 |   645   |  reclin2  | 632  |  25520  |  13  |    0.98     | 0.98 |    1    | 667.5  |
|  1   | 1360 | 1495 |   645   |   diyar   |  47  |    5    | 598  |    0.07     | 0.10 |   12    | 874.4  |
|  1   | 1360 | 1495 |   645   | fedmatch  | 402  |   958   | 243  |    0.62     | 0.70 |   10    | 874.8  |
|  1   | 1360 | 1495 |   645   |  FlexRL   | 205  |   113   | 440  |    0.32     | 0.36 |   22    | 874.8  |
|  2   |  44  |  53  |   24    |    BRL    |  22  |   10    |  2   |    0.92     | 0.31 |    1    |  5.5   |
|  2   |  44  |  53  |   24    | multilink |  23  |   10    |  1   |    0.96     | 0.30 |    1    | 124.2  |
|  2   |  44  |  53  |   24    | fastLink  |  21  |   14    |  3   |    0.88     | 0.40 |    3    |  15.9  |
|  2   |  44  |  53  |   24    |  reclin2  |  23  |   50    |  1   |    0.96     | 0.68 |    1    |  2.9   |
|  2   |  44  |  53  |   24    |   diyar   |  10  |    0    |  14  |    0.42     | 0.00 |    1    |  17.5  |
|  2   |  44  |  53  |   24    | fedmatch  |  22  |   22    |  2   |    0.92     | 0.50 |    1    |  34.2  |
|  2   |  44  |  53  |   24    |  FlexRL   |  22  |    1    |  2   |    0.92     | 0.04 |    3    | 109.9  |
|  3   | 2103 | 2140 |   918   |    BRL    |  40  |   11    | 878  |    0.04     | 0.22 |   37    | 1120.9 |
|  3   | 2103 | 2140 |   918   | multilink | 121  |   585   | 797  |    0.13     | 0.83 |  3063   | 3385.0 |
|  3   | 2103 | 2140 |   918   | fastLink  | 450  |  1383   | 468  |    0.49     | 0.75 |    6    | 483.2  |
|  3   | 2103 | 2140 |   918   |  reclin2  | 911  |  54782  |  7   |    0.99     | 0.98 |    2    | 859.8  |
|  3   | 2103 | 2140 |   918   |   diyar   |  54  |    1    | 864  |    0.06     | 0.02 |   23    | 1405.2 |
|  3   | 2103 | 2140 |   918   | fedmatch  | 545  |  1558   | 373  |    0.59     | 0.74 |   17    | 768.7  |
|  3   | 2103 | 2140 |   918   |  FlexRL   | 272  |   148   | 646  |    0.30     | 0.35 |   33    | 768.6  |
|  4   | 437  | 541  |   229   |    BRL    |  81  |   29    | 148  |    0.35     | 0.26 |    3    | 276.7  |
|  4   | 437  | 541  |   229   | multilink |  96  |   83    | 133  |    0.42     | 0.46 |   75    | 768.3  |
|  4   | 437  | 541  |   229   | fastLink  | 141  |   233   |  88  |    0.62     | 0.62 |    6    | 768.4  |
|  4   | 437  | 541  |   229   |  reclin2  | 228  |  16275  |  1   |    1.00     | 0.99 |    1    | 153.9  |
|  4   | 437  | 541  |   229   |   diyar   |  10  |    0    | 219  |    0.04     | 0.00 |    3    | 768.8  |
|  4   | 437  | 541  |   229   | fedmatch  | 157  |   280   |  72  |    0.69     | 0.64 |    3    | 545.5  |
|  4   | 437  | 541  |   229   |  FlexRL   |  74  |   31    | 155  |    0.32     | 0.30 |    9    | 769.1  |
|  5   | 1181 | 1319 |   750   |    BRL    | 156  |   31    | 594  |    0.21     | 0.17 |   13    | 764.7  |
|  5   | 1181 | 1319 |   750   | multilink | 266  |   410   | 484  |    0.35     | 0.61 |   609   | 1268.9 |
|  5   | 1181 | 1319 |   750   | fastLink  | 365  |   710   | 385  |    0.49     | 0.66 |   26    | 1805.3 |
|  5   | 1181 | 1319 |   750   |  reclin2  | 598  | 341135  | 152  |    0.80     | 1.00 |    1    | 744.8  |
|  5   | 1181 | 1319 |   750   |   diyar   |  4   |    0    | 746  |    0.01     | 0.00 |   10    | 741.9  |
|  5   | 1181 | 1319 |   750   | fedmatch  | 349  |   832   | 401  |    0.47     | 0.70 |    9    | 213.3  |
|  5   | 1181 | 1319 |   750   |  FlexRL   | 195  |   100   | 555  |    0.26     | 0.34 |   19    | 709.5  |
|  6   | 459  | 513  |   238   |    BRL    |  98  |   35    | 140  |    0.41     | 0.26 |    3    | 277.0  |
|  6   | 459  | 513  |   238   | multilink | 122  |   84    | 116  |    0.51     | 0.41 |   71    | 746.1  |
|  6   | 459  | 513  |   238   | fastLink  | 153  |   202   |  85  |    0.64     | 0.57 |    5    | 505.8  |
|  6   | 459  | 513  |   238   |  reclin2  | 236  |  2752   |  2   |    0.99     | 0.92 |    1    | 152.1  |
|  6   | 459  | 513  |   238   |   diyar   |  6   |    0    | 232  |    0.03     | 0.00 |    3    | 650.0  |
|  6   | 459  | 513  |   238   | fedmatch  | 169  |   289   |  69  |    0.71     | 0.63 |    4    | 117.9  |
|  6   | 459  | 513  |   238   |  FlexRL   |  98  |   48    | 140  |    0.41     | 0.33 |    8    | 314.8  |
|  7   | 670  | 710  |   321   |    BRL    |  79  |   31    | 242  |    0.25     | 0.28 |    5    | 544.7  |
|  7   | 670  | 710  |   321   | multilink | 132  |   151   | 189  |    0.41     | 0.53 |   178   | 741.9  |
|  7   | 670  | 710  |   321   | fastLink  | 210  |   337   | 111  |    0.65     | 0.62 |    6    |  83.3  |
|  7   | 670  | 710  |   321   |  reclin2  | 316  |  6320   |  5   |    0.98     | 0.95 |    1    | 306.3  |
|  7   | 670  | 710  |   321   |   diyar   |  34  |    2    | 287  |    0.11     | 0.06 |    4    | 747.1  |
|  7   | 670  | 710  |   321   | fedmatch  | 227  |   443   |  94  |    0.71     | 0.66 |    5    | 141.8  |
|  7   | 670  | 710  |   321   |  FlexRL   | 127  |   54    | 194  |    0.40     | 0.30 |   11    | 409.2  |
|  8   | 1479 | 1524 |   898   |    BRL    | 140  |   23    | 758  |    0.16     | 0.14 |   18    | 742.7  |
|  8   | 1479 | 1524 |   898   | multilink | 290  |   421   | 608  |    0.32     | 0.59 |  1130   | 2010.0 |
|  8   | 1479 | 1524 |   898   | fastLink  | 452  |   865   | 446  |    0.50     | 0.66 |   37    | 2174.9 |
|  8   | 1479 | 1524 |   898   |  reclin2  | 721  | 339232  | 177  |    0.80     | 1.00 |    1    | 930.9  |
|  8   | 1479 | 1524 |   898   |   diyar   |  20  |    0    | 878  |    0.02     | 0.00 |   14    | 879.2  |
|  8   | 1479 | 1524 |   898   | fedmatch  | 510  |   969   | 388  |    0.57     | 0.66 |   12    | 217.9  |
|  8   | 1479 | 1524 |   898   |  FlexRL   | 299  |   104   | 599  |    0.33     | 0.26 |   24    | 721.5  |
|  9   | 1363 | 1500 |   711   |    BRL    |  44  |   10    | 667  |    0.06     | 0.19 |   18    | 897.3  |
|  9   | 1363 | 1500 |   711   | multilink | 107  |   371   | 604  |    0.15     | 0.78 |  1060   | 1891.1 |
|  9   | 1363 | 1500 |   711   | fastLink  | 343  |   857   | 368  |    0.48     | 0.71 |   38    | 2576.3 |
|  9   | 1363 | 1500 |   711   |  reclin2  |  8   | 188790  | 703  |    0.01     | 1.00 |    1    | 829.4  |
|  9   | 1363 | 1500 |   711   |   diyar   |  27  |    0    | 684  |    0.04     | 0.00 |   12    | 837.2  |
|  9   | 1363 | 1500 |   711   | fedmatch  | 338  |  1025   | 373  |    0.48     | 0.75 |   11    | 212.7  |
|  9   | 1363 | 1500 |   711   |  FlexRL   | 199  |   88    | 512  |    0.28     | 0.31 |   22    | 706.0  |
|  10  | 658  | 687  |   384   |    BRL    |  70  |   12    | 314  |    0.18     | 0.15 |    4    | 519.6  |
|  10  | 658  | 687  |   384   | multilink | 119  |   139   | 265  |    0.31     | 0.54 |   159   | 838.6  |
|  10  | 658  | 687  |   384   | fastLink  | 206  |   359   | 178  |    0.54     | 0.64 |    8    | 840.0  |
|  10  | 658  | 687  |   384   |  reclin2  | 321  | 102152  |  63  |    0.84     | 1.00 |    1    | 300.0  |
|  10  | 658  | 687  |   384   |   diyar   |  2   |    0    | 382  |    0.01     | 0.00 |    5    | 827.2  |
|  10  | 658  | 687  |   384   | fedmatch  | 238  |   420   | 146  |    0.62     | 0.64 |    5    | 126.7  |
|  10  | 658  | 687  |   384   |  FlexRL   | 113  |   58    | 271  |    0.29     | 0.34 |   11    | 374.1  |
|  11  | 762  | 872  |   486   |    BRL    | 119  |   30    | 367  |    0.24     | 0.20 |    6    | 767.1  |
|  11  | 762  | 872  |   486   | multilink | 182  |   213   | 304  |    0.37     | 0.54 |   220   | 839.7  |
|  11  | 762  | 872  |   486   | fastLink  | 249  |   422   | 237  |    0.51     | 0.63 |   10    | 843.6  |
|  11  | 762  | 872  |   486   |  reclin2  | 394  | 135933  |  92  |    0.81     | 1.00 |    1    | 438.6  |
|  11  | 762  | 872  |   486   |   diyar   |  4   |    0    | 482  |    0.01     | 0.00 |    6    | 841.5  |
|  11  | 762  | 872  |   486   | fedmatch  | 251  |   511   | 235  |    0.52     | 0.67 |    6    | 145.9  |
|  11  | 762  | 872  |   486   |  FlexRL   | 142  |   54    | 344  |    0.29     | 0.28 |   13    | 429.0  |
|  12  | 977  | 1142 |   455   |    BRL    |  38  |   14    | 417  |    0.08     | 0.27 |   11    | 804.0  |
|  12  | 977  | 1142 |   455   | multilink | 125  |   235   | 330  |    0.27     | 0.65 |   516   | 1085.7 |
|  12  | 977  | 1142 |   455   | fastLink  | 240  |   567   | 215  |    0.53     | 0.70 |   18    | 1323.5 |
|  12  | 977  | 1142 |   455   |  reclin2  | 375  |  3911   |  80  |    0.82     | 0.91 |    1    | 713.7  |
|  12  | 977  | 1142 |   455   |   diyar   |  6   |    0    | 449  |    0.01     | 0.00 |    8    | 716.1  |
|  12  | 977  | 1142 |   455   | fedmatch  | 254  |   722   | 201  |    0.56     | 0.74 |    8    | 172.2  |
|  12  | 977  | 1142 |   455   |  FlexRL   | 152  |   56    | 303  |    0.33     | 0.27 |   17    | 529.0  |
|  13  | 454  | 493  |   311   |    BRL    | 139  |   24    | 172  |    0.45     | 0.15 |    3    | 264.6  |
|  13  | 454  | 493  |   311   | multilink | 164  |   95    | 147  |    0.53     | 0.37 |   57    | 715.0  |
|  13  | 454  | 493  |   311   | fastLink  | 197  |   203   | 114  |    0.63     | 0.51 |    6    | 715.3  |
|  13  | 454  | 493  |   311   |  reclin2  | 310  |  17730  |  1   |    1.00     | 0.98 |    1    | 146.0  |
|  13  | 454  | 493  |   311   |   diyar   |  39  |    0    | 272  |    0.13     | 0.00 |    3    | 557.1  |
|  13  | 454  | 493  |   311   | fedmatch  | 198  |   255   | 113  |    0.64     | 0.56 |    3    | 105.6  |
|  13  | 454  | 493  |   311   |  FlexRL   | 137  |   33    | 174  |    0.44     | 0.19 |    8    | 279.1  |
|  14  | 280  | 299  |   183   |    BRL    | 101  |   22    |  82  |    0.55     | 0.18 |    1    | 107.1  |
|  14  | 280  | 299  |   183   | multilink | 113  |   44    |  70  |    0.62     | 0.28 |   19    | 716.2  |
|  14  | 280  | 299  |   183   | fastLink  | 131  |   106   |  52  |    0.72     | 0.45 |    4    | 367.1  |
|  14  | 280  | 299  |   183   |  reclin2  | 182  |  7282   |  1   |    0.99     | 0.98 |    1    |  55.6  |
|  14  | 280  | 299  |   183   |   diyar   |  2   |    0    | 181  |    0.01     | 0.00 |    2    | 325.4  |
|  14  | 280  | 299  |   183   | fedmatch  | 137  |   143   |  46  |    0.75     | 0.51 |    2    |  83.4  |
|  14  | 280  | 299  |   183   |  FlexRL   | 105  |   20    |  78  |    0.57     | 0.16 |    6    | 186.1  |
|  15  | 1939 | 2160 |  1069   |    BRL    |  25  |   22    | 1044 |    0.02     | 0.47 |   35    | 947.1  |
|  15  | 1939 | 2160 |  1069   | multilink |  69  |   648   | 1000 |    0.06     | 0.90 |  2229   | 3837.2 |
|  15  | 1939 | 2160 |  1069   | fastLink  | 276  |  1526   | 793  |    0.26     | 0.85 |   74    | 4640.2 |
|  15  | 1939 | 2160 |  1069   |  reclin2  | 1069 | 1572964 |  0   |    1.00     | 1.00 |    2    | 2191.0 |
|  15  | 1939 | 2160 |  1069   |   diyar   |  1   |    0    | 1068 |    0.00     | 0.00 |   21    | 1505.1 |
|  15  | 1939 | 2160 |  1069   | fedmatch  | 280  |  1659   | 789  |    0.26     | 0.86 |   16    | 291.4  |
|  15  | 1939 | 2160 |  1069   |  FlexRL   | 108  |   117   | 961  |    0.10     | 0.52 |   30    | 847.4  |
|  16  | 1130 | 1354 |   683   |    BRL    |  47  |   11    | 636  |    0.07     | 0.19 |   13    | 842.6  |
|  16  | 1130 | 1354 |   683   | multilink | 112  |   470   | 571  |    0.16     | 0.81 |   528   | 1382.1 |
|  16  | 1130 | 1354 |   683   | fastLink  | 266  |   774   | 417  |    0.39     | 0.74 |   26    | 1916.9 |
|  16  | 1130 | 1354 |   683   |  reclin2  | 681  | 404295  |  2   |    1.00     | 1.00 |    1    | 825.0  |
|  16  | 1130 | 1354 |   683   |   diyar   |  9   |    1    | 674  |    0.01     | 0.10 |   10    | 822.6  |
|  16  | 1130 | 1354 |   683   | fedmatch  | 290  |   840   | 393  |    0.42     | 0.74 |   10    | 209.6  |
|  16  | 1130 | 1354 |   683   |  FlexRL   | 154  |   64    | 529  |    0.23     | 0.29 |   19    | 689.2  |
|  17  | 264  | 276  |   193   |    BRL    | 116  |   45    |  77  |    0.60     | 0.28 |    1    |  95.2  |
|  17  | 264  | 276  |   193   | multilink | 122  |   78    |  71  |    0.63     | 0.39 |   13    | 829.0  |
|  17  | 264  | 276  |   193   | fastLink  | 130  |   100   |  63  |    0.67     | 0.43 |    4    | 280.5  |
|  17  | 264  | 276  |   193   |  reclin2  | 193  |  7266   |  0   |    1.00     | 0.97 |    1    |  49.5  |
|  17  | 264  | 276  |   193   |   diyar   |  16  |    1    | 177  |    0.08     | 0.06 |    2    | 312.2  |
|  17  | 264  | 276  |   193   | fedmatch  | 127  |   137   |  66  |    0.66     | 0.52 |    2    |  88.6  |
|  17  | 264  | 276  |   193   |  FlexRL   | 134  |   39    |  59  |    0.69     | 0.23 |    5    | 193.3  |
|  18  | 501  | 554  |   362   |    BRL    |  99  |   22    | 263  |    0.27     | 0.18 |    3    | 334.9  |
|  18  | 501  | 554  |   362   | multilink | 151  |   131   | 211  |    0.42     | 0.46 |   61    | 832.6  |
|  18  | 501  | 554  |   362   | fastLink  | 205  |   246   | 157  |    0.57     | 0.55 |    6    | 832.3  |
|  18  | 501  | 554  |   362   |  reclin2  | 311  | 118067  |  51  |    0.86     | 1.00 |    1    | 190.1  |
|  18  | 501  | 554  |   362   |   diyar   |  7   |    0    | 355  |    0.02     | 0.00 |    3    | 682.0  |
|  18  | 501  | 554  |   362   | fedmatch  | 194  |   307   | 168  |    0.54     | 0.61 |    4    | 120.6  |
|  18  | 501  | 554  |   362   |  FlexRL   | 121  |   49    | 241  |    0.33     | 0.29 |    9    | 341.8  |
|  19  | 1524 | 1625 |  1010   |    BRL    |  52  |    9    | 958  |    0.05     | 0.15 |   19    | 829.8  |
|  19  | 1524 | 1625 |  1010   | multilink | 138  |   574   | 872  |    0.14     | 0.81 |   999   | 2245.8 |
|  19  | 1524 | 1625 |  1010   | fastLink  | 354  |  1057   | 656  |    0.35     | 0.75 |   41    | 2855.0 |
|  19  | 1524 | 1625 |  1010   |  reclin2  | 1006 | 677113  |  4   |    1.00     | 1.00 |    2    | 1277.3 |
|  19  | 1524 | 1625 |  1010   |   diyar   |  5   |    0    | 1005 |    0.00     | 0.00 |   15    | 892.4  |
|  19  | 1524 | 1625 |  1010   | fedmatch  | 351  |  1173   | 659  |    0.35     | 0.77 |   12    | 259.5  |
|  19  | 1524 | 1625 |  1010   |  FlexRL   | 166  |   58    | 844  |    0.16     | 0.26 |   24    | 649.3  |
|  20  | 826  | 827  |   478   |    BRL    |  69  |   24    | 409  |    0.14     | 0.26 |    6    | 649.3  |
|  20  | 826  | 827  |   478   | multilink | 142  |   273   | 336  |    0.30     | 0.66 |   239   | 846.6  |
|  20  | 826  | 827  |   478   | fastLink  | 232  |   491   | 246  |    0.49     | 0.68 |   12    | 852.4  |
|  20  | 826  | 827  |   478   |  reclin2  | 478  |  71205  |  0   |    1.00     | 0.99 |    1    | 443.9  |
|  20  | 826  | 827  |   478   |   diyar   |  13  |    2    | 465  |    0.03     | 0.13 |    6    | 851.7  |
|  20  | 826  | 827  |   478   | fedmatch  | 245  |   580   | 233  |    0.51     | 0.70 |    6    | 163.4  |
|  20  | 826  | 827  |   478   |  FlexRL   | 140  |   63    | 338  |    0.29     | 0.31 |   13    | 508.2  |

Averaged over the regions and total computing time:

``` r

summary_IREG <- do.call(rbind, lapply(split(results_IREG, results_IREG$method), function(d) {
  data.frame(TP = mean(d$TP, na.rm=TRUE), FP = mean(d$FP, na.rm=TRUE), 
             FN = mean(d$FN, na.rm=TRUE), sensitivity = mean(d$sensitivity, na.rm=TRUE),
             FDP = mean(d$FDP, na.rm=TRUE), method = d$method[1],
             minutes = sum(d$seconds) / 60,
             memory = mean(d$memory))
}))
rownames(summary_IREG) <- summary_IREG$method
```

``` r

newdf <- as.data.frame(apply(t(summary_IREG)[c("TP","FP","FN","sensitivity","FDP","minutes","memory"),methods], 2, function(x)  sprintf("%.2f",as.numeric(as.character(x)))))
names(newdf) <- names(df_results_full)
rownames(newdf) <- rownames(df_results_full)
knitr::kable(newdf, align=rep('c', ncol(newdf)))
```

|             |  BRL   | multilink | fastLink |  reclin2  | diyar  | fedmatch | FlexRL |
|:------------|:------:|:---------:|:--------:|:---------:|:------:|:--------:|:------:|
| TP          | 80.00  |  137.00   |  247.00  |  450.00   | 16.00  |  264.00  | 148.00 |
| FP          | 22.00  |  270.70   |  563.20  | 204638.70 |  0.60  |  656.10  | 64.90  |
| FN          | 438.00 |  380.00   |  270.00  |   68.00   | 502.00 |  253.00  | 369.00 |
| sensitivity |  0.26  |   0.36    |   0.55   |   0.89    |  0.06  |   0.58   |  0.36  |
| FDP         |  0.23  |   0.58    |   0.63   |   0.97    |  0.02  |   0.67   |  0.29  |
| time (mins) |  3.65  |  205.65   |   5.67   |   0.38    |  2.72  |   2.43   |  5.10  |
| memory (Mb) | 634.00 |  1338.00  | 1203.00  |  564.00   | 763.00 |  245.00  | 520.00 |

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
