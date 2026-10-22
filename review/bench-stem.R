# FlexRL StEM benchmark harness.
#
# Usage:  R_LIBS=<lib> Rscript review/bench-stem.R [scale]
#
# Produces, for the SAME fixed problem and iteration settings:
#   1. end-to-end StEM wall-clock,
#   2. an Rprof self-time breakdown of one smaller run.
# Run it before and after each PERF-PLAN.md phase; record numbers in §8.
# The package under test must be installed into the library given via R_LIBS.

scale <- suppressWarnings(as.numeric(commandArgs(TRUE)[1]))
if (is.na(scale)) scale <- 1

n_records <- round(c(10000, 15000) * scale)
n_links   <- round(7500 * scale)

set.seed(20261021)
PIVs_config <- list(
  V1 = list(dynamics = "stable",     bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  V2 = list(dynamics = "stable",     bound_mistakes = c(0.10, 0.10), fix_mistakes = c(NA, NA)),
  V3 = list(dynamics = "flexible",   bound_mistakes = c(NA, NA),     fix_mistakes = c(NA, NA)),
  V4 = list(dynamics = "structured", bound_mistakes = c(NA, NA),     fix_mistakes = c(0.03, 0.03),
            cond_hazard_cov = list(cov1 = c("Xe", "Xf"), cov2 = c()))
)
n_values <- c(5, 6, 7, 12)
p_mistake <- list(V1 = c(0.02, 0.02), V2 = c(0.02, 0.02), V3 = c(0.05, 0.05), V4 = c(0.02, 0.02))
p_missing <- list(V1 = c(0.005, 0.005), V2 = c(0.005, 0.005), V3 = c(0.005, 0.005), V4 = c(0.005, 0.005))
cond_hazard_params <- list(V1 = c(), V2 = c(), V3 = c(), V4 = log(c(0.7, 0.6, 0.5)))

cat(sprintf("scale=%g  nA=%d  nB=%d  n_links=%d\n", scale, n_records[1], n_records[2], n_links))

gen <- FlexRL::simulate_data(PIVs_config, n_values, n_records, n_links,
                             p_mistake, p_missing, cond_hazard_params, TRUE)
prep <- FlexRL::prepare_data(gen$data1, gen$data2, "fileA", "fileB", PIVs_config,
                             uniq_id = "entity_id")

args <- list(data = prep, StEM_iter = 3, StEM_burnin = 1,
             gibbs_iter = 4, gibbs_burnin = 2, n_post_sample = 10)

# --- 1. end-to-end ------------------------------------------------------------
set.seed(1)
t0 <- proc.time()
fit <- suppressWarnings(do.call(FlexRL::StEM, args))
t_end <- (proc.time() - t0)[["elapsed"]]
n_sweeps <- 3 * 4 + 10
cat(sprintf("StEM end-to-end: %.2f s  (%d sweeps, %.3f s/sweep)\n",
            t_end, n_sweeps, t_end / n_sweeps))

# --- 2. profiled run for a kernel breakdown -----------------------------------
args_p <- args
args_p$StEM_iter <- 2; args_p$StEM_burnin <- 1; args_p$n_post_sample <- 5
set.seed(1)
utils::Rprof(interval = 0.02)
fit2 <- suppressWarnings(do.call(FlexRL::StEM, args_p))
utils::Rprof(NULL)
p <- utils::summaryRprof()$by.self
total <- sum(p$self.time)
cat(sprintf("\nProfiled run total: %.1f s — top self-time:\n", total))
show <- min(15, nrow(p))
for (i in seq_len(show)) {
  cat(sprintf("  %-24s %7.2f s  %5.1f%%\n",
              paste(rownames(p)[i], collapse = " > "), p$self.time[i], 100 * p$self.time[i] / total))
}
pt <- utils::summaryRprof()$by.total
cat("\nTop total-time (call-tree attribution):\n")
show <- min(12, nrow(pt))
for (i in seq_len(show)) {
  cat(sprintf("  %-24s %7.2f s  %5.1f%%\n",
              paste(rownames(pt)[i], collapse = " > "), pt$total.time[i], 100 * pt$total.time[i] / total))
}
n_sweeps_p <- args_p$StEM_iter * args_p$gibbs_iter + args_p$n_post_sample
cat(sprintf("\nSUMMARY: scale=%g end-to-end %.2f s (%d sweeps, %.3f s/sweep) | profiled %.1f s (%d sweeps, %.3f s/sweep)\n",
            scale, t_end, n_sweeps, t_end / n_sweeps, total, n_sweeps_p, total / n_sweeps_p))
invisible(NULL)
