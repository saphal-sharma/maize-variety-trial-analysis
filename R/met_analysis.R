# Multi-environment variety trial (MET) analysis for plant breeding
#
# Part A  Maize hybrids, 17 hybrids x 12 environments (Kang & Gorman, 1989)
#         - environment means and two-way ANOVA
#         - stability: Finlay-Wilkinson regression and Wricke's ecovalence
#         - GGE biplot and "which-won-where" pattern
# Part B  Winter rapeseed, 6 varieties x 17 trials with 4 reps (Shafii & Price, 1998)
#         - linear mixed model, variance components, broad-sense heritability
#         - BLUPs for variety ranking
#
# Run from the repository root:  Rscript R/met_analysis.R
# Packages: dplyr, tidyr, ggplot2, lme4

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(lme4)
})
dir.create("results", showWarnings = FALSE)

# ------------------------------------------------------------------ Part A
mz <- read.delim("data/kang.maize.txt") |> select(gen, env, yield, year, loc)

env_means <- mz |> group_by(env) |> summarise(env_mean = mean(yield)) |> arrange(env_mean)
gen_means <- mz |> group_by(gen) |> summarise(mean_yield = mean(yield))

# Two-way ANOVA on genotype x environment means. Without replicates, the
# residual is the G x E interaction.
aov_tab <- anova(lm(yield ~ gen + env, data = mz))
ss <- aov_tab[["Sum Sq"]]
ss_share <- data.frame(source = c("Genotype", "Environment", "G x E (residual)"),
                       sum_sq = round(ss, 2), share_pct = round(ss / sum(ss) * 100, 1))

# Finlay-Wilkinson: regress each hybrid's yield on the environment mean.
# b > 1: responsive to good conditions; b < 1: relatively stable in poor conditions.
mz <- mz |> left_join(env_means, by = "env")
fw <- mz |> group_by(gen) |>
  summarise(fw_slope = coef(lm(yield ~ env_mean))[2],
            fw_r2 = summary(lm(yield ~ env_mean))$r.squared, .groups = "drop")

# Wricke's ecovalence: each hybrid's contribution to the G x E sum of squares.
grand <- mean(mz$yield)
mz <- mz |> left_join(gen_means, by = "gen") |>
  mutate(ge = yield - mean_yield - env_mean + grand)
eco <- mz |> group_by(gen) |> summarise(ecovalence = sum(ge^2), .groups = "drop") |>
  mutate(ecovalence_pct = 100 * ecovalence / sum(ecovalence))

stab <- gen_means |> left_join(fw, by = "gen") |> left_join(eco, by = "gen") |>
  mutate(rank_yield = rank(-mean_yield), rank_stability = rank(ecovalence)) |>
  arrange(desc(mean_yield)) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(stab, "results/maize_stability.csv", row.names = FALSE)

# GGE biplot: SVD of the environment-centred genotype x environment matrix
mat <- mz |> select(gen, env, yield) |> pivot_wider(names_from = env, values_from = yield)
G <- as.matrix(mat[, -1]); rownames(G) <- mat$gen
Gc <- sweep(G, 2, colMeans(G))
sv <- svd(Gc)
pc_var <- round(100 * sv$d^2 / sum(sv$d^2), 1)
gscore <- sv$u[, 1:2] %*% diag(sqrt(sv$d[1:2]))     # symmetric scaling
escore <- sv$v[, 1:2] %*% diag(sqrt(sv$d[1:2]))
gdf <- data.frame(name = rownames(G), PC1 = gscore[, 1], PC2 = gscore[, 2], type = "Hybrid")
edf <- data.frame(name = colnames(G), PC1 = escore[, 1], PC2 = escore[, 2], type = "Environment")
winners <- data.frame(env = colnames(G), best_hybrid = rownames(G)[apply(G, 2, which.max)],
                      best_yield = apply(G, 2, max))
write.csv(winners, "results/maize_winner_by_environment.csv", row.names = FALSE)

p1 <- ggplot() +
  geom_segment(data = edf, aes(0, 0, xend = PC1, yend = PC2), colour = "darkgreen",
               arrow = arrow(length = unit(0.15, "cm")), alpha = 0.7) +
  geom_text(data = edf, aes(PC1, PC2, label = name), colour = "darkgreen", size = 3, vjust = -0.5) +
  geom_point(data = gdf, aes(PC1, PC2), colour = "navy") +
  geom_text(data = gdf, aes(PC1, PC2, label = name), colour = "navy", size = 3, vjust = 1.4) +
  geom_hline(yintercept = 0, linewidth = 0.2) + geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(title = "GGE biplot: maize hybrids across 12 Louisiana environments",
       x = paste0("PC1 (", pc_var[1], "%)"), y = paste0("PC2 (", pc_var[2], "%)")) +
  theme_minimal()
ggsave("results/maize_gge_biplot.png", p1, width = 8, height = 6, dpi = 150)

p2 <- ggplot(stab, aes(mean_yield, fw_slope, label = gen)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = mean(stab$mean_yield), linetype = 2, colour = "grey50") +
  geom_point(aes(size = ecovalence_pct), colour = "steelblue", alpha = 0.7) +
  geom_text(size = 3, vjust = -1) +
  labs(title = "Yield vs responsiveness of maize hybrids",
       x = "Mean yield (t/ha)", y = "Finlay-Wilkinson slope (b)",
       size = "Share of G x E (%)") +
  theme_minimal()
ggsave("results/maize_yield_vs_stability.png", p2, width = 8, height = 5.5, dpi = 150)

# ------------------------------------------------------------------ Part B
rs <- read.delim("data/shafii.rapeseed.txt") |>
  mutate(env = paste0(loc, year), yield_t = yield / 1000)

# Genotype random: variance components for heritability
m_rand <- lmer(yield_t ~ 1 + (1 | gen) + (1 | env) + (1 | env:rep) + (1 | gen:env), data = rs)
vc <- as.data.frame(VarCorr(m_rand))[, c("grp", "vcov")]
v <- setNames(vc$vcov, vc$grp)
n_env <- n_distinct(rs$env)
n_rep <- n_distinct(rs$rep)
H2 <- v["gen"] / (v["gen"] + v["gen:env"] / n_env + v["Residual"] / (n_env * n_rep))

blup <- ranef(m_rand)$gen |> tibble::rownames_to_column("variety") |>
  rename(blup_t_ha = `(Intercept)`) |>
  mutate(predicted_yield_t_ha = fixef(m_rand)[1] + blup_t_ha) |>
  left_join(rs |> group_by(variety = gen) |> summarise(raw_mean_t_ha = mean(yield_t)), by = "variety") |>
  arrange(desc(predicted_yield_t_ha)) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(blup, "results/rapeseed_blups.csv", row.names = FALSE)

vc_out <- data.frame(component = names(v), variance = round(v, 4),
                     share_pct = round(100 * v / sum(v), 1))
write.csv(vc_out, "results/rapeseed_variance_components.csv", row.names = FALSE)

# ------------------------------------------------------------------ Print
cat("\n== Maize: share of variation\n"); print(ss_share)
cat("\nEnvironment means (t/ha), lowest to highest\n"); print(as.data.frame(env_means), digits = 3)
cat("\n== Maize hybrid stability (sorted by mean yield)\n"); print(as.data.frame(stab))
cat("\nGGE biplot variance explained: PC1", pc_var[1], "% PC2", pc_var[2], "%\n")
cat("\nWinning hybrid per environment\n"); print(table(winners$best_hybrid))
cat("\n== Rapeseed variance components\n"); print(vc_out)
cat(sprintf("\nBroad-sense heritability (entry-mean basis): %.2f  (%d environments, %d reps)\n",
            H2, n_env, n_rep))
cat("\nRapeseed BLUPs\n"); print(as.data.frame(blup))
