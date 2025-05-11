library(readxl)
library(dplyr)

raw_data_ct <- read_excel("./data/B_data.xlsx") #You can also change to A 
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

ind_ct <- which(!is.na(shedding.data.ct$ct_value))
ind_viral <- which(!is.na(shedding.data.viral$value))

data_multi_noncensored <- list(
  "nobs_ct"= length(ind_ct), #number of noncensored observations;
  "nobs_viral" = length(ind_viral),
  "t_obs_ct" = shedding.data.ct$days[ind_ct], #days after symptom onset;
  "t_obs_viral" = shedding.data.viral$days[ind_viral],
  "c_obs_ct" = shedding.data.ct$ct_value[ind_ct], #log10 concentration of SARS-CoV-2 RNA (copies per mL);
  "c_obs_viral" = shedding.data.viral$value[ind_viral]
)
rm(ind_ct) #remove indicators;
rm(ind_viral) #remove indicators;
data_multi_noncensored

library(rstan)
source("./models/plotting_helper.R")
fit_multi <- stan(
  file = "./models/expon1_combined.stan",  # Combined Stan model path
  data = data_multi_noncensored,
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 4,
  refresh = 500
)

print(fit_multi, pars = c("a0", "c0", "beta0", "beta1", "sigma_viral"), probs = c(0.1, 0.5, 0.9))

stan_post_multi <- extract(fit_multi, permuted = FALSE)
a0_combined <- as.vector(stan_post_multi[, , "a0"])
c0_combined <- as.vector(stan_post_multi[, , "c0"])
beta0_combined <- as.vector(stan_post_multi[, , "beta0"])
beta1_combined <- as.vector(stan_post_multi[, , "beta1"])


# ───────────────────────────────────────────────────────────────────────────────
# Ct_value-R plot
# ───────────────────────────────────────────────────────────────────────────────

# Convert Ct value to log10 concentration
RNA_c_combined <- as.numeric(mean(beta0_combined)) + as.numeric(mean(beta1_combined)) * data_multi_noncensored$c_obs_ct
censorlimit_ct <- 40 
censorlimit_RNA_combined <- as.numeric(mean(beta0_combined)) + as.numeric(mean(beta1_combined)) * censorlimit_ct

vir_exponential <- function(tm, a0, c0, beta0, beta1) {
  return((c0 / log(10) - a0 * tm / log(10) - beta0) / beta1)
}
quant.resp <- function(t) {
  vals <- vir_exponential(t, a0_combined, c0_combined, beta0_combined, beta1_combined)
  return(quantile(vals, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}
graph.resp <- function(color = "black") {
  x <- seq(from = 5, to = 50, by = 0.2)
  y <- array(NA, dim = c(length(x), 3))
  for (k in seq_along(x)) {
    y[k, ] <- quant.resp(x[k])
  }
  lines(x, y[, 1], col = color, lty = 3, lwd = 2)
  lines(x, y[, 2], col = color, lty = 1, lwd = 2)
  lines(x, y[, 3], col = color, lty = 3, lwd = 2)
}

# Base R plot
plot(data_multi_noncensored$t_obs_ct,
     data_multi_noncensored$c_obs_ct,
     pch = 16, col = "black",
     xlab = "Days After Symptom Onset",
     ylab = "Ct Value",
     ylim = rev(range(data_multi_noncensored$c_obs_ct, 45)),
     main = "Ct Value Model Fit")
lines(x = c(5, 50), y = rep(censorlimit_ct, 2), lty = 3, col = "black")
graph.resp("black")
legend("topright",
       legend = c("Quantifiable", "Negative", "Median", "95% Credible Interval"),
       lty = c(NA, NA, 1, 3),
       pch = c(16, 1, NA, NA),
       lwd = 2,
       pt.cex = 1)
text(12, 39.5, "Limit of Quantification")

# ───────────────────────────────────────────────────────────────────────────────
# Ct_value transfer to Viral Load R plot 
# ───────────────────────────────────────────────────────────────────────────────
ct.vir_exponential <- function(beta0, beta1, ct) {
  return(beta0 + beta1 * ct)
}

quant.resp.2 <- function(t) {
  vals_ct <- (c0_combined / log(10) - a0_combined * t / log(10) - beta0_combined) / beta1_combined
  vals_log10_viral <- beta0_combined + beta1_combined * vals_ct
  return(quantile(vals_log10_viral, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
}

graph.resp.2 <- function(color = "black") {
  x <- seq(from = 5, to = 50, by = 0.2)
  y <- array(NA, dim = c(length(x), 3))
  for (k in seq_along(x)) {
    y[k, ] <- quant.resp.2(x[k])
  }
  lines(x, y[, 1], col = color, lty = 3, lwd = 2)
  lines(x, y[, 2], col = color, lty = 1, lwd = 2)
  lines(x, y[, 3], col = color, lty = 3, lwd = 2)
}

plot(data_multi_noncensored$t_obs_ct,
     RNA_c_combined,
     pch = 16, col = "black",
     xlab = "Days After Symptom Onset",
     ylab = expression(log[10]*"(Viral Load from Ct)"),
     main = "Viral Load (Transformed from Ct)")
lines(x = c(5, 50), y = rep(censorlimit_RNA_combined, 2), lty = 3, col = "black")
graph.resp.2("black")

legend("topright",
       legend = c("Quantifiable", "Negative", "Median", "95% Credible Interval"),
       lty = c(NA, NA, 1, 3),
       pch = c(16, 1, NA, NA),
       lwd = 2,
       pt.cex = 1)
text(12, censorlimit_RNA_combined + 0.05, "Limit of Quantification")


# ───────────────────────────────────────────────────────────────────────────────
# viral load-R plot
# ───────────────────────────────────────────────────────────────────────────────
# customize the plot function;
censorlimit_viral <- 2
plot(data_multi_noncensored$t_obs_viral,
     data_multi_noncensored$c_obs_viral,
     pch = 16, col = "black",
     xlab = "Days After Symptom Onset",
     ylab = expression(paste("Log"[10], " Viral Load (copies/mL)")),
     ylim = c(0, 10),
     main = "Viral Load Model Fit")
lines(x = c(0, 30), y = rep(censorlimit_viral, 2), lty = 3, col = "black")
graph.resp.1("red", model = "exp")
legend("topright",
       legend = c("Quantifiable", "Negative", "Median", "95% Credible Interval"),
       lty = c(NA, NA, 1, 3),
       pch = c(16, 1, NA, NA),
       lwd = 2,
       pt.cex = 1.5)
text(5, 2.3, "Limit of Quantification")


# ───────────────────────────────────────────────────────────────────────────────
# Combined-R plot
# ───────────────────────────────────────────────────────────────────────────────

# Create unified plot
plot(-1, -1,
     xlim = c(0, 50),
     ylim = c(0, 10),
     xlab = "Days After Symptom Onset",
     ylab = expression(log[10]*"(Viral Load)"),
     main = "Combined Viral Load: Observed & Transformed",
     col = "white")

lines(x = c(0, 50), y = rep(censorlimit_viral, 2), lty = 3, col = "red")
lines(x = c(0, 50), y = rep(censorlimit_RNA_combined, 2), lty = 3, col = "black")

points(data_multi_noncensored$t_obs_viral,
       data_multi_noncensored$c_obs_viral,
       pch = 16, col = "red")
graph.resp.1("red", model = "exp")

points(data_multi_noncensored$t_obs_ct,
       RNA_c_combined,
       pch = 16, col = "black")
graph.resp.2("black") 

legend("topright",
       legend = c("Viral Load Obs", "Ct-derived Viral", "Median (Obs)", "Median (Ct-derived)",
                  "95% CI (Obs)", "95% CI (Ct-derived)", "LOQ Viral", "LOQ Ct-derived"),
       lty = c(NA, NA, 1, 1, 3, 3, 3, 3),
       pch = c(16, 16, NA, NA, NA, NA, NA, NA),
       col = c("red", "black", "red", "black", "red", "black", "red", "black"),
       lwd = 2,
       pt.cex = 1.2,
       cex = 0.8)

text(5, censorlimit_viral + 0.2, "LOQ (Viral Load)", col = "red")
text(5, censorlimit_RNA_combined + 0.2, "LOQ (Ct-derived)", col = "black")















