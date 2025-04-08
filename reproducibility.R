shedding.data_rep <- data.frame(subjects = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9), 
                            day = c(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 19, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 7, 8, 9, 10, 11, 12, 14, 15, 16, 19, 20, 6, 8, 9, 12, 20, 21, 26, 6, 8, 9, 13, 16, 18, 21, 22, 6, 18, 19, 20, 23, 24, 26, 3, 4, 5, 8, 9, 12, 8, 9, 10),
                            value = c(6.64, 7.16, 6.68, 4.67, 4.7, 3.67, 3.7, 3.73, 3.01, 2.42, 3.46, 2.8, 5.74, NA, 4.88, 4.77, 3.63, 3.25, 2.83, 3.84, 3.87, 3.11, 2.28, 6.1, 6.03, 6.79, 6.62, 4.68, 5.1, 5.79, 4.3, 4.89, 3.64, 4.71, 4.75, 4.09, 2.33, NA, NA, NA, 6.05, 7.23, 7.51, 5.57, 4.43, 5.15, 4.88, 4.25, 3.98, 3.53, 3.56, 4.68, 2.78, 2.99, 4.27, 2.57, 2.95, NA, 3.88, 3.78, 4.09, 4.23, NA, 4.58, NA, 2.57, 3.66, 2.59, NA, NA, 2.04, NA, 1.97, 4.68, 4.33, 3.61, 3.99, 2.71, NA, NA, NA, NA))
censorlimit = 2 ##log10 scale limit of quantification
#Create an indicator for noncensored observations of Subject 3;
ind <- which(shedding.data_rep$subjects==3 & !is.na(shedding.data_rep$value))
#Create input data for Stan model 
data_one_noncensored_rep <- list(
  "nobs"=length(ind), #number of noncensored observations;
  "t" = shedding.data_rep$day[ind], #days after symptom onset;
  "c" = shedding.data_rep$value[ind] #log10 concentration of SARS-CoV-2 RNA (copies per mL);
)
rm(ind) #remove the indicator;
data_one_noncensored_rep
fit_one_noncensored_rep <- stan(
  file = "./models/expon1.stan",  # Stan program
  data = data_one_noncensored_rep,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = 500             # no progress shown
)
print(fit_one_noncensored_rep, pars=c("a0","c0","sig_obs"), probs=c(.1,.5,.9))

#extract posterior samples
stan_post_rep <- extract(fit_one_noncensored_rep, permuted = FALSE)
a0 <- as.vector(stan_post_rep[,,2])
c0 <- as.vector(stan_post_rep[,,3])

#graph
plot(-1,-1,xlim=c(0,30),ylim=c(-5,10),col="white",
     xlab="Days after Symptom Onset",ylab="Genome Copies per mL",main="Subject 3",yaxt="n",cex.axis=1.5,cex.lab=1.5,cex.main=2);
lines(x=c(0,30),y=c(censorlimit,censorlimit),lty=3,col="black");
points(data_one_noncensored_rep$t,data_one_noncensored_rep$c,col="black",cex=1.5,pch=16)
graph.resp.1(color="black",model="exp")
legend("topright", legend=c("Quantifiable","Negative","Median","95% Credible Interval"), lty=c(NA,NA,1,3), pch=c(16,1,NA,NA), lwd=2, pt.cex=2)
text(5,2.5,"Limit of Quantification")
ticks.log(2,cex.axis=1)







