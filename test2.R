#data for wang2020fecal
library(readxl)
library(dplyr)
raw_data <- read_excel("./B data.xlsx") #You can also change to A: A data.xlsx → Corresponds to ORF1ab gene measurements; B data.xlsx → Corresponds to N gene measurements
shedding.data.wang <- raw_data %>%
  rename(subject = patients, 
         ct_value = ctvalue) %>%
  mutate(subject = 1,
    days = as.integer(round(days))
  )

censorlimit_2 = 40 
#Create an indicator for noncensored observations of Subject 1;
ind_2 <- which(shedding.data.wang$subject==1 & !is.na(shedding.data.wang$ct_value))
#Create input data for Stan model 
data_one_noncensored_2 <- list(
  "nobs"=length(ind_2), #number of noncensored observations;
  "t" = shedding.data.wang$days[ind_2], 
  "ct_value" = shedding.data.wang$ct_value[ind_2] 
)
rm(ind_2) 
data_one_noncensored_2
library(rstan)
fit_one_noncensored_2 <- stan(
  file = "./models/expon1ct.stan",  
  data = data_one_noncensored,   
  chains = 4,            
  warmup = 1000,          
  iter = 2000,        
  cores = 4,             
  refresh = 500            
)
print(fit_one_noncensored_2, pars=c("a0","c0","beta0","beta1","sigma"), probs=c(.1,.5,.9))
fit_summary_2 <- summary(fit_one_noncensored_2, 
                       pars = c("a0","c0","beta0","beta1","sigma"), 
                       probs = c(.1,.5,.9))

fit_summary_matrix_2 <- fit_summary_2$summary
mean_values_2 <- fit_summary_matrix_2[, "mean"]
beta0_hat_2 <- mean_values_2["beta0"]
beta1_hat_2 <- mean_values_2["beta1"]

stan_post_2 <- extract(fit_one_noncensored_2, permuted = FALSE)
a0_2 <- as.vector(stan_post_2[,,2])
c0_2 <- as.vector(stan_post_2[,,3])
beta0_2  <- as.vector(stan_post_2[,,4])
beta1_2  <- as.vector(stan_post_2[,,5])

RNA_c_2 <- as.numeric(beta0_hat_2) + as.numeric(beta1_hat_2) * data_one_noncensored_2$ct_value

censorlimit_RNA_2 = as.numeric(beta0_hat_2) + as.numeric(beta1_hat_2) * censorlimit_2
#graph

plot(-1,-1,xlim=c(0,30),ylim=c(-5,10),col="white",
     xlab="Days after Symptom Onset",ylab="Genome Copies per mL",main="All Subjects Treated As One in Wang",
     yaxt="n",cex.axis=1.5,cex.lab=1.5,cex.main=2);
# lines(x=c(0,30),y=c(censorlimit,censorlimit),lty=3,col="black");
lines(x = c(0, 30), y = rep(censorlimit_RNA_2, 2), lty = 3, col = "black")
points(data_one_noncensored_2$t,RNA_c_2,col="black",cex=1.5,pch=16)
graph.resp.1(color="black",model="exp")
legend("topright", legend=c("Quantifiable","Negative","Median","95% Credible Interval"), lty=c(NA,NA,1,3), pch=c(16,1,NA,NA), lwd=2, pt.cex=2)
text(4,-3.5,"Limit of Quantification")
ticks.log(2,cex.axis=1)


par(mfrow = c(2, 1)) 

# Plot 1: Ct values
plot(data_one_noncensored_2$t, 
     data_one_noncensored_2$ct_value,
     pch = 16, col = "darkblue",
     xlab = "", ylab = "Ct Value",
     main = "Raw Ct Values", ylim = rev(range(data_one_noncensored_2$ct_value)))

# Plot 2: Estimated RNA concentrations
plot(data_one_noncensored_2$t, RNA_c_2,
     pch = 16, col = "black",
     xlab = "Days after Symptom Onset",
     ylab = "Estimated RNA Concentration",
     main = "Transformed Viral Load (from Ct)")


