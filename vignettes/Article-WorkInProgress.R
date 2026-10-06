
# load data

df2016 <- read.csv("SHIW2016.csv", row.names = 1)
df2020 <- read.csv("SHIW2020.csv", row.names = 1)
df2016 <- df2016[!df2016$ID %in% df2016$ID[duplicated(df2016$ID)], ]
df2020 <- df2020[!df2020$ID %in% df2020$ID[duplicated(df2020$ID)], ]

PIVs_config <- list(
  ANASCI = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  SESSO  = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  STUDIO = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  STACIV = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  IREG   = list(dynamics = "stable", bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)))
PIVs <- names(PIVs_config)
PIVs_type <- list(ANASCI = TRUE, SESSO = FALSE, STUDIO = FALSE, STACIV = FALSE, IREG = FALSE)

prep_data <- prepare_data(df2016, df2020, "2016", "2020", PIVs_config, same_mistakes = TRUE, uniq_id = "ID")

RL_agreement(prep_data$encodedA, prep_data$encodedB, PIVs, prep_data$true_pairs)

truelinksdata <- merge(prep_data$encodedA, prep_data$encodedB, by = "ID")

# add a new dynamic piv "move"
prep_data$encodedA[,"date"] <- prep_data$encodedA[,"ETA"]
prep_data$encodedB[,"date"] <- rep(0, nrow(prep_data$encodedB))
prep_data$encodedA[,"change"] <- FALSE
prep_data$encodedB[,"change"] <- FALSE

time_difference <- abs(prep_data$encodedB[prep_data$true_pairs[[2]], "date"] - prep_data$encodedA[prep_data$true_pairs[[1]], "date"])

intercept <- rep(1, nrow(prep_data$true_pairs))
proba_same_H <- matrix(1, nrow(prep_data$true_pairs), 6)
cov <- cbind(intercept)
model_dynamics <- survival_model("exponential")
proba_same_H[, 6] <- model_dynamics$S(as.matrix(cov), log(0.03), time_difference)

xp <- exp(0.23 * (0:(30 - 1)))
prep_data$encodedA[,"MOVE"] <- sample(seq_len(30), nrow(prep_data$encodedA), replace = TRUE, prob = xp / sum(xp))
prep_data$encodedB[,"MOVE"] <- sample(seq_len(30), nrow(prep_data$encodedB), replace = TRUE, prob = xp / sum(xp))
prep_data$encodedB[prep_data$true_pairs[[2]],"MOVE"] <- prep_data$encodedA[prep_data$true_pairs[[1]],"MOVE"]

k = 6
for (i in seq_len(nrow(prep_data$true_pairs))) {
  is_not_changing <- stats::rbinom(1, 1, proba_same_H[i, k])
  if (!is_not_changing && !is.na(prep_data$encodedA[prep_data$true_pairs[i,1], "MOVE"])) {
    prep_data$encodedB[prep_data$true_pairs[i,2], "MOVE"] <- sample((1:30)[-c(prep_data$encodedA[prep_data$true_pairs[i,1], "MOVE"])], 1)
    prep_data$encodedB[prep_data$true_pairs[i,2], "change"] <- TRUE
  }
}
###

PIVs_config <- c(PIVs_config, list(MOVE = list(dynamics = "structured", bound_mistakes = c(NA, NA), fix_mistakes = c(0, 0), cond_hazard_cov = list(cov1 = c(), cov2 = c()))))
PIVs <- names(PIVs_config)
PIVs_type <- c(PIVs_type, MOVE = FALSE)

prep_data <- prepare_data(prep_data$encodedA, prep_data$encodedB, "2020", "2016", PIVs_config, same_mistakes = TRUE, uniq_id = "ID")
true_pairs <- do.call(paste, c(prep_data$true_pairs, list(sep = "_")))

fit <- StEM(data = prep_data, StEM_iter = 10, StEM_burnin = 5, gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10)

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

df_results_full <- data.frame(matrix(NA, nrow = 6, ncol = 0))
rownames(df_results_full) <- c("TP", "FP", "FN", "sensitivity", "FDP", "minutes")
for (method in methods) {
  t0 <- Sys.time()
  
  
  df_results_full[1:5, method] <- evaluate_linkage(fit, true_pairs_full)
  df_results_full[6, method] <- if (is.null(fit)) NA else round(as.numeric(Sys.time() - t0, units = "mins"))
}








time_difference <- truelinksdata$ETA.x
intercept <- rep(1, nrow(truelinksdata))
proba_same_H <- matrix(1, nrow(truelinksdata), length(PIVs)+1)
cov <- cbind(intercept)
model_dynamics <- survival_model("exponential")
proba_same_H[, length(PIVs)+1] <- model_dynamics$S(as.matrix(cov), log(0.03), time_difference)

plot(truelinksdata$ETA.x, proba_same_H[,6])

xp <- exp(0.23 * (0:(30 - 1)))
truelinksdata$move.x <- sample(seq_len(30), nrow(truelinksdata), replace = TRUE, prob = xp / sum(xp))
truelinksdata$move.y <- truelinksdata$move.x
truelinksdata$change <- FALSE

k = 6
for (i in seq_len(nrow(truelinksdata))) {
  is_not_changing <- stats::rbinom(1, 1, proba_same_H[i, k])
  if (!is_not_changing && !is.na(truelinksdata[i, "move.x"])) {
    truelinksdata[i, "move.y"] <- sample((1:30)[-c(truelinksdata[i, "move.x"])], 1)
    truelinksdata[i, "change"] <- TRUE
  }
}

plot(truelinksdata$ETA.x, truelinksdata$move.x == truelinksdata$move.y)
###

X <- cbind(intercept = rep(1, nrow(truelinksdata)))
times <- truelinksdata$ETA.y
Hequal <- truelinksdata$STUDIO.x == truelinksdata$STUDIO.y
test = survival_model("exponential")
# test$S(X, alpha, times)
# expo$S(X, alpha = log(0.3), times)
alphastar = stats::nlminb(test$init(ncol(X)), test$negloglik, X = X, times = times, Hequal = Hequal)$par
plot(times, test$S(X,alphastar, times))
plot(times, Hequal)

plot(truelinksdata[truelinksdata$STUDIO.x == truelinksdata$STUDIO.y,"ETA.x"])
plot(truelinksdata[truelinksdata$STUDIO.x != truelinksdata$STUDIO.y,"ETA.x"])

sort(truelinksdata$STUDIO.x == truelinksdata$STUDIO.y)


(truelinksdata, )

apply(df2016, 2, function(x) length(unique(x)))
apply(df2020, 2, function(x) length(unique(x)))

# ETA: age (realted to birth year)
# nonoc unemployment type
# sttp7 emplyment branch of activity

# in 2016
# ETA: age (realted to birth year) (cont)
# STUDIO: Educational qualification (cat)
# SETTP11: sector of activity
# QUALP10: Main employment, work status

# in 2020
# Y etapen expected age of retirement (cont)

