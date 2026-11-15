# Checks the analysis outputs against the reported findings.
# Run from the repository root after the analysis:  Rscript tests/check_results.R
stab <- read.csv("results/maize_stability.csv")
stopifnot(nrow(stab) == 17)
win <- read.csv("results/maize_winner_by_environment.csv")
stopifnot(nrow(win) == 12, length(unique(win[[2]])) == 8)

vc <- read.csv("results/rapeseed_variance_components.csv")
v <- setNames(vc$variance, vc$component)
rs <- read.delim("data/shafii.rapeseed.txt")
n_env <- length(unique(paste0(rs$loc, rs$year)))
n_rep <- length(unique(rs$rep))
H2 <- v["gen"] / (v["gen"] + v["gen:env"] / n_env + v["Residual"] / (n_env * n_rep))
stopifnot(abs(H2 - 0.08) < 0.02)

mz <- read.delim("data/kang.maize.txt")
ss <- anova(lm(yield ~ gen + env, data = mz))[["Sum Sq"]]
env_share <- ss[2] / sum(ss)
stopifnot(abs(env_share - 0.77) < 0.01)
cat("All checks passed: environment share =", round(env_share, 3), "; rapeseed H2 =", round(H2, 3), "\n")
