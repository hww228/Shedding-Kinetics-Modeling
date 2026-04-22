library(readxl)
library(dplyr)

raw_data_ct <- read_excel("./data/ct_data_1.xlsx") 
shedding.data.ct <- raw_data_ct %>%
  rename(subjects = patients, 
         ct_value = ctvalue) %>%
  mutate(subject = 1,
         days = as.integer(round(days))
  )
shedding.data.viral <- data.frame(
  subjects = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9), 
  days = c(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 7, 8, 9, 10, 11, 12, 14, 15, 16, 19, 20, 6, 8, 9, 12, 20, 21, 26, 6, 8, 9, 13, 16, 18, 21, 22, 6, 18, 19, 20, 23, 24, 26, 3, 4, 5, 8, 9, 12, 8, 9, 10),
  value = c(6.64, 7.16, 6.68, 4.67, 4.7, 3.67, 3.7, 3.73, 3.01, 2.42, 3.46, 2.8, 5.74, NA, 4.88, 4.77, 3.63, 3.25, 2.83, 3.84, 3.87, 3.11, 2.28, 6.1, 6.03, 6.79, 6.62, 4.68, 5.1, 5.79, 4.3, 4.89, 3.64, 4.71, 4.75, 4.09, 2.33, NA, NA, NA, 6.05, 7.23, 7.51, 5.57, 4.43, 5.15, 4.88, 4.25, 3.98, 3.53, 3.56, 4.68, 2.78, 2.99, 4.27, 2.57, 2.95, NA, 3.88, 3.78, 4.09, 4.23, NA, 4.58, NA, 2.57, 3.66, 2.59, NA, NA, 2.04, NA, 1.97, 4.68, 4.33, 3.61, 3.99, 2.71, NA, NA, NA, NA)
  )

#### woelf & wang
# ───────────────────────────────────────────────────────────────────────────────
# 3) Censoring rules
# ───────────────────────────────────────────────────────────────────────────────
censorlimit_viral <- 2   # viral LOQ on log10 scale
censorlimit_ct <- 40         # Ct right-censor threshold

# Viral: observed if >= LOQ; censored if NA or < LOQ
v <- shedding.data.viral$value
days_v <- shedding.data.viral$days
ind_obs_vl <- which(!is.na(v) & v >= censorlimit_viral)
ind_cen_vl <- which(is.na(v) | v < censorlimit_viral)

t_obs_vl <- days_v[ind_obs_vl]
c_obs_vl <- v[ind_obs_vl]
t_cen_vl <- days_v[ind_cen_vl]

# Ct: observed if < ct_lim; censored if NA or >= ct_lim
ct <- shedding.data.ct$ct_value
days_ct <- shedding.data.ct$days
ind_obs_ct <- which(!is.na(ct) & ct < censorlimit_ct)
ind_cen_ct <- which(is.na(ct) | ct >= censorlimit_ct)

t_obs_ct <- days_ct[ind_obs_ct]
ct_obs <- ct[ind_obs_ct]
t_cen_ct <- days_ct[ind_cen_ct]

# Quick checks (these prevent the earlier Stan constraint failure in expon1cen-style models)
cat("Viral: nobs =", length(t_obs_vl), " ncen =", length(t_cen_vl), "\n")
cat("Ct   : nobs =", length(t_obs_ct), " ncen =", length(t_cen_ct), "\n")
cat("min observed viral =", min(c_obs_vl), " (must be >= censlim_viral =", censorlimit_viral, ")\n")

# ───────────────────────────────────────────────────────────────────────────────
# 4) Build Stan data for expon1_combined.stan
# ───────────────────────────────────────────────────────────────────────────────
data_combined <- list(
  # viral (matches expon1cen.stan structure)
  nobs    = length(t_obs_vl),
  ncen    = length(t_cen_vl),
  t_obs   = as.vector(t_obs_vl),
  c_obs   = as.vector(c_obs_vl),
  t_cen   = as.vector(t_cen_vl),
  censlim = censorlimit_viral,
  
  # ct (uses ct_lim as you specified)
  nobs_ct = length(t_obs_ct),
  ncen_ct = length(t_cen_ct),
  t_obs_ct = as.vector(t_obs_ct),
  ct_obs   = as.vector(ct_obs),
  t_cen_ct = as.vector(t_cen_ct),
  ct_lim   = censorlimit_ct
)

# ───────────────────────────────────────────────────────────────────────────────
# 5) Fit combined Stan model
# ───────────────────────────────────────────────────────────────────────────────
library(rstan)
source("./models/plotting_helper.R")

fit_multi <- stan(
  file    = "./models/expon1_combined.stan",
  data    = data_combined,
  chains  = 4,
  warmup  = 1000,
  iter    = 2000,
  cores   = 4,
  refresh = 500
)

# Adjust pars to what exists in your Stan file.
# (Common: a0, c0, beta0, beta1, sig_obs, sig_ct)
print(fit_multi, pars = c("a0", "c0", "alpha0","alpha1", "sig_obs", "sig_ct"), probs = c(0.1, 0.5, 0.9))

# ───────────────────────────────────────────────────────────────────────────────
# 6) Extract posterior draws
# ───────────────────────────────────────────────────────────────────────────────
post <- rstan::extract(fit_multi, permuted = FALSE)

a0 <- as.vector(post[, , "a0"])
c0 <- as.vector(post[, , "c0"])
alpha0 <- as.vector(post[, , "alpha0"])
alpha1 <- as.vector(post[, , "alpha1"]) 

# If your combined Stan uses different names, change these:
sig_obs <- if ("sig_obs" %in% dimnames(post)[[3]]) as.vector(post[, , "sig_obs"]) else NULL
sig_ct  <- if ("sig_ct"  %in% dimnames(post)[[3]]) as.vector(post[, , "sig_ct"])  else NULL

# ───────────────────────────────────────────────────────────────────────────────
# 7) Helper functions for curves and summaries  (UPDATED for alpha0/alpha1 model)
# ───────────────────────────────────────────────────────────────────────────────

# Viral mean on log10 scale (same as Stan gene_obs)
mu_vl <- function(t, a0_draws, c0_draws) {
  c0_draws / log(10) - a0_draws * t / log(10)
}

# Ct forward mean: Ct = alpha0 + alpha1 * VL10  (NO DIVISION)
mu_ct_forward <- function(t, a0_draws, c0_draws, alpha0_draws, alpha1_draws) {
  alpha0_draws + alpha1_draws * mu_vl(t, a0_draws, c0_draws)
}

# Quantiles helper
curve_quantiles <- function(draws, probs = c(0.025, 0.5, 0.975)) {
  quantile(as.numeric(draws), probs = probs, na.rm = TRUE)
}

# Draw viral trajectory band
draw_vl_band <- function(xseq, col = "black") {
  y <- t(sapply(xseq, function(tt) curve_quantiles(mu_vl(tt, a0, c0))))
  lines(xseq, y[, 1], col = col, lty = 3, lwd = 2)
  lines(xseq, y[, 2], col = col, lty = 1, lwd = 2)
  lines(xseq, y[, 3], col = col, lty = 3, lwd = 2)
}

# Draw Ct trajectory band (UPDATED)
draw_ct_band <- function(xseq, col = "black") {
  y <- t(sapply(xseq, function(tt) curve_quantiles(mu_ct_forward(tt, a0, c0, alpha0, alpha1))))
  lines(xseq, y[, 1], col = col, lty = 3, lwd = 2)
  lines(xseq, y[, 2], col = col, lty = 1, lwd = 2)
  lines(xseq, y[, 3], col = col, lty = 3, lwd = 2)
}

# Ct -> viral load (log10) via inversion of forward model:
# Ct = alpha0 + alpha1*VL10  =>  VL10 = (Ct - alpha0)/alpha1
vl_from_ct_draws <- function(ct_value) {
  (ct_value - alpha0) / alpha1
}

# Posterior summary for Ct-derived viral point (names guaranteed)
vl_from_ct_summary <- function(ct_value) {
  q <- quantile(vl_from_ct_draws(ct_value),
                probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  out <- c(lo = q[1], med = q[2], hi = q[3])
  out
}

# Ct-derived LOQ line on viral scale (posterior median), using ct_lim from data
vl_loq_from_ct_med <- median(vl_from_ct_draws(data_combined$ct_lim), na.rm = TRUE)

# # ───────────────────────────────────────────────────────────────────────────────
# # 8) Plot 1: Ct fit (observed Ct < ct_lim)
# # ───────────────────────────────────────────────────────────────────────────────
# x_ct <- seq(0, 50, by = 0.2)
# 
# plot(
#   x    = data_combined$t_obs_ct,
#   y    = data_combined$ct_obs,
#   pch  = 16,
#   col  = "black",
#   xlab = "Days After Symptom Onset",
#   ylab = "Ct Value",
#   xlim = c(0, 50),
#   ylim = rev(range(c(data_combined$ct_obs, data_combined$ct_lim), na.rm = TRUE)),
#   main = "Ct Model Fit (Right-censored at Ct ≥ ct_lim)"
# )
# abline(h = data_combined$ct_lim, lty = 3, col = "black")
# draw_ct_band(x_ct, col = "black")
# text(2, data_combined$ct_lim - 0.5, "Limit of Quantification", adj = 0)
# 
# # ───────────────────────────────────────────────────────────────────────────────
# # 9) Plot 2: Viral load fit (observed viral >= LOQ)
# # ───────────────────────────────────────────────────────────────────────────────
# x_vl <- seq(0, 50, by = 0.2)
# 
# plot(
#   x    = data_combined$t_obs,
#   y    = data_combined$c_obs,
#   pch  = 16,
#   col  = "black",
#   xlab = "Days After Symptom Onset",
#   ylab = expression(paste("Log"[10], " Viral Load (copies/mL)")),
#   xlim = c(0, 50),
#   ylim = range(c(data_combined$c_obs, data_combined$censlim), na.rm = TRUE),
#   main = "Viral Load Model Fit (Left-censored at LOQ)"
# )
# abline(h = data_combined$censlim, lty = 3, col = "black")
# draw_vl_band(x_vl, col = "black")
# text(2, data_combined$censlim + 0.2, "Viral LOQ", adj = 0)
# 
# # ───────────────────────────────────────────────────────────────────────────────
# # 10) Plot 3: Ct-derived viral load points (with uncertainty)  (UPDATED)
# # ───────────────────────────────────────────────────────────────────────────────
# ct_summ <- t(vapply(data_combined$ct_obs,
#                     vl_from_ct_summary,
#                     FUN.VALUE = c(lo=0, med=0, hi=0)))
# 
# vl_ct_lo  <- ct_summ[, "lo"]
# vl_ct_med <- ct_summ[, "med"]
# vl_ct_hi  <- ct_summ[, "hi"]
# 
# ylim_ctvl <- range(c(vl_ct_lo, vl_ct_hi, vl_loq_from_ct_med), na.rm = TRUE)
# 
# plot(
#   x    = data_combined$t_obs_ct,
#   y    = vl_ct_med,
#   pch  = 16,
#   col  = "black",
#   xlab = "Days After Symptom Onset",
#   ylab = expression(log[10]*"(Viral Load from Ct)"),
#   xlim = c(0, 50),
#   ylim = ylim_ctvl,
#   main = "Ct-derived Viral Load (via posterior inversion)"
# )
# 
# # Draw uncertainty FIRST, then points so points are visible
# # segments(data_combined$t_obs_ct, vl_ct_lo, data_combined$t_obs_ct, vl_ct_hi, col = "black")
# points(data_combined$t_obs_ct, vl_ct_med, pch = 16, col = "black", cex = 1.1)
# 
# abline(h = vl_loq_from_ct_med, lty = 3, col = "black")
# text(2, vl_loq_from_ct_med + 0.2, "Ct-derived LOQ (median)", adj = 0)

# ───────────────────────────────────────────────────────────────────────────────
# 11) Plot 4: Combined viral plot (observed viral + Ct-derived viral)  (UPDATED)
# ───────────────────────────────────────────────────────────────────────────────
ylim_all <- range(c(data_combined$c_obs, vl_ct_lo, vl_ct_hi,
                    data_combined$censlim, vl_loq_from_ct_med), na.rm = TRUE)
ylim_modify <- c(-3,11)

par(mgp = c(2.2, 0.8, 0)) 
plot( 0, 0,
  type = "n",
  xlim = c(0, 50),
  ylim = ylim_modify,
  xlab = "Days After Symptom Onset",
  ylab = expression(log[10]*"(Viral Load)"),
  main = "Combined Viral Load Trajectory: Observed Viral + Ct-derived Viral",
  cex.lab = 1.4,   # axis label size
  cex.axis = 1.4,  # tick label size
  cex.main = 1.7   # title size
)

# LOQ lines
abline(h = data_combined$censlim, lty = 3, col = "red")
abline(h = vl_loq_from_ct_med, lty = 3, col = "blue")

# Observed viral
points(data_combined$t_obs, data_combined$c_obs, pch = 16, col = "red")

# Ct-derived: error bars first, then points
# segments(data_combined$t_obs_ct, vl_ct_lo, data_combined$t_obs_ct, vl_ct_hi, col = "black")
points(data_combined$t_obs_ct, vl_ct_med, pch = 16, col = "blue", cex = 1.1)

# Combined Viral trajectory band (a0,c0)
draw_vl_band(x_vl, col = "black")

# Original Viral trajectory band
graph.resp("red")

legend(
  "topright",
  legend = c("Observed Viral (woelfel2020virological)",
             "Ct-derived Viral (wang2020fecal)",
             "Viral LOQ",
             "Ct-derived LOQ",
             "Original Viral trajectory",
             "Combined Viral trajectory"),
  pch    = c(16, 16, NA, NA, NA, NA),
  lty    = c(NA, NA, 3, 3, 1, 1),
  col    = c("red", "blue", "red", "blue", "red","black"),
  lwd    = c(NA, NA, 2, 2, 2, 2),
  pt.cex = 1.3,
  cex = 1.3
)




