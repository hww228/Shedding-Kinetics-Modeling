#data for kim2020viral
shedding.data.kim <- data.frame(subjects = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2), 
                            day = c(5,6,7,8,9,10,11,12,14,15,16,5,6,7,8,9,10,11,12,14,15,16,16,17,19,21, 24, 16, 17, 19, 21, 24),
                            ct_value = c(NA, NA, NA, 36.27, NA, NA, 34.29, 33.45, NA, NA, NA, 37.66, 37.75, 35.54, 33.29, NA, 28.2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 30.45, NA, NA, NA))
censorlimit = 35
#Create an indicator for noncensored observations of Subject 1;
ind <- which(shedding.data.kim$subjects==1 & !is.na(shedding.data.kim$ct_value))
#Create input data for Stan model 
data_one_noncensored <- list(
  "nobs"=length(ind), #number of noncensored observations;
  "t" = shedding.data.kim$day[ind], 
  "ct_value" = shedding.data.kim$ct_value[ind] 
)
rm(ind) 
data_one_noncensored
library(rstan)
fit_one_noncensored <- stan(
  file = "./models/expon1ct.stan",  
  data = data_one_noncensored,   
  chains = 4,            
  warmup = 1000,          
  iter = 2000,        
  cores = 4,             
  refresh = 500            
)
print(fit_one_noncensored, pars=c("a0","c0","beta0","beta1","sigma"), probs=c(.1,.5,.9))
fit_summary <- summary(fit_one_noncensored, 
                       pars = c("a0","c0","beta0","beta1","sigma"), 
                       probs = c(.1,.5,.9))

fit_summary_matrix <- fit_summary$summary
mean_values <- fit_summary_matrix[, "mean"]
beta0_hat <- mean_values["beta0"]
beta1_hat <- mean_values["beta1"]

stan_post <- extract(fit_one_noncensored, permuted = FALSE)
a0 <- as.vector(stan_post[,,2])
c0 <- as.vector(stan_post[,,3])
beta0  <- as.vector(stan_post[,,4])
beta1  <- as.vector(stan_post[,,5])

RNA_c <- as.numeric(beta0_hat) + as.numeric(beta1_hat) * data_one_noncensored$ct_value

censorlimit_RNA = as.numeric(beta0_hat) + as.numeric(beta1_hat) * censorlimit
#graph

plot(-1,-1,xlim=c(0,30),ylim=c(-5,10),col="white",
     xlab="Days after Symptom Onset",ylab="Genome Copies per mL",main="Subject 1",
     yaxt="n",cex.axis=1.5,cex.lab=1.5,cex.main=2);
lines(x = c(0, 30), y = rep(censorlimit_RNA, 2), lty = 3, col = "black")
points(data_one_noncensored$t,RNA_c,col="black",cex=1.5,pch=16)
graph.resp.1(color="black",model="exp")
legend("topright", legend=c("Quantifiable","Negative","Median","95% Credible Interval"), lty=c(NA,NA,1,3), pch=c(16,1,NA,NA), lwd=2, pt.cex=2)
text(5,-4,"Limit of Quantification")
ticks.log(2,cex.axis=1)



par(mfrow = c(2, 1), mar = c(4, 5, 3, 2))  # two stacked plots

# Plot 1: Ct values
plot(data_one_noncensored$t, 
     data_one_noncensored$ct_value,
     pch = 16, col = "darkblue",
     xlab = "", ylab = "Ct Value",
     main = "Raw Ct Values", ylim = rev(range(data_one_noncensored$ct_value)))

# Plot 2: Estimated RNA concentrations
plot(data_one_noncensored$t, RNA_c,
     pch = 16, col = "black",
     xlab = "Days after Symptom Onset",
     ylab = "Estimated RNA Concentration",
     main = "Transformed Viral Load (from Ct)")
par(mfrow = c(1, 1)) 
