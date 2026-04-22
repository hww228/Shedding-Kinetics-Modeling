library(readxl)
library(dplyr)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# ───────────────────────────────────────────────────────────────────────────────
# 1) Load & prep Ct data
# ───────────────────────────────────────────────────────────────────────────────
raw_data_ct <- read_excel("./data/B_data.xlsx")  # change to A_data.xlsx if needed

shedding.data.ct <- raw_data_ct %>%
  rename(subjects = patients,
         ct_value = ctvalue) %>%
  mutate(
    subject = 1,
    days = as.integer(round(days))
  )

# ───────────────────────────────────────────────────────────────────────────────
# 2) Censoring rule (right-censor at Ct >= censorlimit_ct)
# ───────────────────────────────────────────────────────────────────────────────
censorlimit_ct <- 40

ind_obs <- which(!is.na(shedding.data.ct$ct_value) &
                   shedding.data.ct$ct_value < censorlimit_ct)

ind_cen <- which(is.na(shedding.data.ct$ct_value) |
                   shedding.data.ct$ct_value >= censorlimit_ct)

t_obs  <- shedding.data.ct$days[ind_obs]
ct_obs <- shedding.data.ct$ct_value[ind_obs]
t_cen  <- shedding.data.ct$days[ind_cen]

data_ct <- list(
  nobs   = length(t_obs),
  ncen   = length(t_cen),
  t_obs  = as.vector(t_obs),
  ct_obs = as.vector(ct_obs),
  t_cen  = as.vector(t_cen),
  ct_lim = censorlimit_ct     # Stan variable name stays ct_lim
)

cat("nobs =", data_ct$nobs, "\n")
cat("ncen =", data_ct$ncen, "\n")

# ───────────────────────────────────────────────────────────────────────────────
# 3) Fit Stan model (should sample cleanly)
# ───────────────────────────────────────────────────────────────────────────────
fit_ct <- stan(
  file    = "./models/expon1ctcen.stan",
  data    = data_ct,
  chains  = 4,
  warmup  = 1000,
  iter    = 2000,
  cores   = 4,
  refresh = 500
)

print(fit_ct, pars = c("ct0", "k", "sig_ct"), probs = c(0.1, 0.5, 0.9))

# ───────────────────────────────────────────────────────────────────────────────
# 4) Extract posterior draws
# ───────────────────────────────────────────────────────────────────────────────
post_ct <- rstan::extract(fit_ct, permuted = FALSE)

ct0   <- as.vector(post_ct[, , "ct0"])
k     <- as.vector(post_ct[, , "k"])
sig_ct <- as.vector(post_ct[, , "sig_ct"])

# ───────────────────────────────────────────────────────────────────────────────
# 5) Plot: observed Ct + fitted trajectory + LOQ line
# ───────────────────────────────────────────────────────────────────────────────
mu_ct_draws <- function(t) ct0 + k * t

curve_quantiles <- function(draws, probs = c(0.025, 0.5, 0.975)) {
  quantile(as.numeric(draws), probs = probs, na.rm = TRUE)
}

draw_ct_band <- function(xseq, col = "black") {
  y <- t(sapply(xseq, function(tt) curve_quantiles(mu_ct_draws(tt))))
  lines(xseq, y[, 1], col = col, lty = 3, lwd = 2)
  lines(xseq, y[, 2], col = col, lty = 1, lwd = 2)
  lines(xseq, y[, 3], col = col, lty = 3, lwd = 2)
}

xseq <- seq(0, 50, by = 0.2)

plot(
  x    = data_ct$t_obs,
  y    = data_ct$ct_obs,
  pch  = 16,
  col  = "black",
  xlab = "Days After Symptom Onset",
  ylab = "Ct Value",
  xlim = c(0, 50),
  ylim = rev(range(c(data_ct$ct_obs, censorlimit_ct), na.rm = TRUE)),
  main = "Ct Model Fit (Right-censored at Ct ≥ 40)"
)

abline(h = censorlimit_ct, lty = 3, col = "black")
draw_ct_band(xseq, col = "black")
text(2, censorlimit_ct - 0.5, "Limit of Quantification", adj = 0)

legend(
  "topright",
  legend = c("Observed Ct (<40)", "Ct LOQ (=40)", "Median", "95% Credible Interval"),
  pch    = c(16, NA, NA, NA),
  lty    = c(NA, 3, 1, 3),
  lwd    = c(NA, 2, 2, 2),
  bty    = "n"
)
