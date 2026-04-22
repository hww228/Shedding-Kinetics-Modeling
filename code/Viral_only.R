library(dplyr)
library(rstan)

# Optional: helper file (ok if you have it; otherwise comment out)


# ---- Viral data (your exact structure) ----
shedding.data.viral <- data.frame(
  subjects = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9),
  days = c(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 7, 8, 9, 10, 11, 12, 14, 15, 16, 19, 20, 6, 8, 9, 12, 20, 21, 26, 6, 8, 9, 13, 16, 18, 21, 22, 6, 18, 19, 20, 23, 24, 26, 3, 4, 5, 8, 9, 12, 8, 9, 10),
  value = c(6.64, 7.16, 6.68, 4.67, 4.7, 3.67, 3.7, 3.73, 3.01, 2.42, 3.46, 2.8, 5.74, NA, 4.88, 4.77, 3.63, 3.25, 2.83, 3.84, 3.87, 3.11, 2.28, 6.1, 6.03, 6.79, 6.62, 4.68, 5.1, 5.79, 4.3, 4.89, 3.64, 4.71, 4.75, 4.09, 2.33, NA, NA, NA, 6.05, 7.23, 7.51, 5.57, 4.43, 5.15, 4.88, 4.25, 3.98, 3.53, 3.56, 4.68, 2.78, 2.99, 4.27, 2.57, 2.95, NA, 3.88, 3.78, 4.09, 4.23, NA, 4.58, NA, 2.57, 3.66, 2.59, NA, NA, 2.04, NA, 1.97, 4.68, 4.33, 3.61, 3.99, 2.71, NA, NA, NA, NA)
)

# ---- Censoring rule for viral load ----
# In your dataset, NA looks like "below LOQ" (censored). We'll treat NA as censored.
censorlimit_viral <- 2

ind_obs <- which(!is.na(shedding.data.viral$value) & shedding.data.viral$value >= censorlimit_viral)
ind_cen <- which(is.na(shedding.data.viral$value) | shedding.data.viral$value < censorlimit_viral)

t_obs <- shedding.data.viral$days[ind_obs]
c_obs <- shedding.data.viral$value[ind_obs]
t_cen <- shedding.data.viral$days[ind_cen]

data_viral <- list(
  nobs    = length(t_obs),
  ncen    = length(t_cen),
  t_obs   = as.vector(t_obs),
  c_obs   = as.vector(c_obs),
  t_cen   = as.vector(t_cen),
  censlim = censorlimit_viral
)

print(data_viral$nobs)
print(data_viral$ncen)

source("./models/plotting_helper.R")
library(rstan)
fit_viral <- stan(
  file    = "./models/expon1cen.stan",
  data    = data_viral,
  chains  = 4,
  warmup  = 1000,
  iter    = 2000,
  cores   = 4,
  refresh = 500
)

print(fit_viral, pars = c("a0", "c0", "sig_obs"), probs = c(0.1, 0.5, 0.9))

# ---- Extract posterior draws ----
post_v <- rstan::extract(fit_viral, permuted = FALSE)
a0_v   <- as.vector(post_v[, , "a0"])
c0_v   <- as.vector(post_v[, , "c0"])

# ---- Posterior curve ----
vir_exponential <- function(tm, a0, c0) {
  c0 / log(10) - a0 * tm / log(10)
}

quant.resp <- function(t) {
  vals <- vir_exponential(t, a0_v, c0_v)
  quantile(vals, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

graph.resp <- function(color = "black") {
  x <- seq(from = 0, to = 50, by = 0.2)
  y <- array(NA, dim = c(length(x), 3))
  for (k in seq_along(x)) y[k, ] <- quant.resp(x[k])
  lines(x, y[, 1], col = color, lty = 3, lwd = 2)
  lines(x, y[, 2], col = color, lty = 1, lwd = 2)
  lines(x, y[, 3], col = color, lty = 3, lwd = 2)
}

# ---- Plot: Viral load fit ----
plot(
  x    = data_viral$t_obs,
  y    = data_viral$c_obs,
  pch  = 16,
  col  = "black",
  xlab = "Days After Symptom Onset",
  ylab = expression(paste("Log"[10], " Viral Load (copies/mL)")),
  xlim = c(0, 50),
  ylim = c(0, 10),
  main = "Viral Load Model Fit (Left-censored)"
)

lines(x = c(0, 50), y = rep(censorlimit_viral, 2), lty = 3, col = "black")
graph.resp("black")

legend(
  "topright",
  legend = c("Observed", "LOQ", "Median", "95% Credible Interval"),
  lty    = c(NA, 3, 1, 3),
  pch    = c(16, NA, NA, NA),
  lwd    = 2,
  pt.cex = 1
)

text(5, censorlimit_viral + 0.2, "Limit of Quantification", adj = 0)
