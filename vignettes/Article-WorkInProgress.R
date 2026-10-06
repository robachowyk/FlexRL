

pak::pak("robachowyk/FlexRL")
library(FlexRL)

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

# data pre-processing

PIVs_config <- c(PIVs_config, list(MOVE = list(dynamics = "structured", bound_mistakes = c(NA, NA), fix_mistakes = c(0, 0), cond_hazard_cov = list(cov1 = c(), cov2 = c()))))
PIVs <- names(PIVs_config)
PIVs_type <- c(PIVs_type, MOVE = FALSE)

prep_data <- prepare_data(prep_data$encodedA, prep_data$encodedB, "2020", "2016", PIVs_config, same_mistakes = TRUE, uniq_id = "ID")
true_pairs <- do.call(paste, c(prep_data$true_pairs, list(sep = "_")))

# run StEM

fit <- StEM(data = prep_data, StEM_iter = 10, StEM_burnin = 5, gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10)

# default results

delta_result = fit$Delta
delta_result = delta_result[delta_result$x>0.5, ]

linked_pairs    = do.call(paste, c(delta_result[,c("i","j")], list(sep = "_")))
true_positive   = length( intersect(linked_pairs, true_pairs) ) 
false_positive  = length( setdiff(linked_pairs, true_pairs) ) 
false_negative  = length( setdiff(true_pairs, linked_pairs) )
sensitivity     = true_positive / (true_positive + false_negative) 
fdp             = false_positive / (true_positive + false_positive)  

c(true_positive,false_positive,false_negative,sensitivity,fdp)

# convergence of the alpha!

apply(fit$alpha$MOVE, 2, mean)
log(0.03)

# linked data for inference

# 2020 data
prep_data$encodedA[delta_result$i, ]

# 2016 data
prep_data$encodedB[delta_result$j, ]

# combined 
combined <- data.frame( cbind( data.frame(prep_data$encodedA[delta_result$i, ]),
                                 data.frame(prep_data$encodedB[delta_result$j, ]) ) )

# in 2016
# ETA: age (realted to birth year) (cont)
# STUDIO: Educational qualification (cat)
# SETTP11: sector of activity
# QUALP10: Main employment, work status
# .y and .1
# in 2020
# outome etapen expected age of retirement (cont)
# .x and _

model <- lm(etapen.x ~ ETA.y + as.factor(STUDIO.y) + as.factor(SETTP11.y) + as.factor(QUALP10.y), data = truelinksdata)
model <- lm(etapen.x ~ ETA.y + as.factor(STUDIO.y) + as.factor(SETTP3.y) + as.factor(QUALP3.y), data = truelinksdata)
summary(model)$coef
summary(model)$coef[summary(model)$coef[,"Pr(>|t|)"] < 0.1,]
summary(model)$adj.r.squared

modelRL <- lm(etapen ~ ETA.1 + as.factor(STUDIO.1) + as.factor(SETTP11.1) + as.factor(QUALP10.1), data = combined)
modelRL <- lm(etapen ~ ETA.1 + as.factor(STUDIO.1) + as.factor(SETTP3.1) + as.factor(QUALP3.1), data = combined)
summary(modelRL)$coef
summary(modelRL)$coef[summary(modelRL)$coef[,"Pr(>|t|)"] < 0.1,]
summary(modelRL)$adj.r.squared

# run diagnostics

diag <- RL_diagnostics(fit = fit, encodedA = prep_data$encodedA, encodedB = prep_data$encodedB,
                       compare_vars = PIVs, vars_type_cont = PIVs_type, PIVs = PIVs, 
                       true_pairs = prep_data$true_pairs, FDP_estimation = TRUE, 
                       RL_method = "FlexRL", data = prep_data, StEM_iter = 10, StEM_burnin = 5, 
                       gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10,
                       maxIter4CV = 1, n_repeats = 2)

