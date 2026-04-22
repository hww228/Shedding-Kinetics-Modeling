#data for wang2020fecal
library(readxl)
library(dplyr)
raw_data <- read_excel("./data/B_data.xlsx") #You can also change to A 
shedding.data.wang <- raw_data %>%
  rename(subject = patients, 
         ct_value = ctvalue) %>%
  mutate(subject = 1,
         days = as.integer(round(days))
  )

censorlimit_2 = 40 
#Create an indicator for noncensored observations of Subject 1;
ind_2 <- which(shedding.data.wang$subject==1 & !is.na(shedding.data.wang$ct_value))
ind_2_cen <- which(shedding.data.wang$subject==1 & !is.na(shedding.data.wang$ct_value) & shedding.data.wang$ct_value >=40)
#Create input data for Stan model 
data_one_censored_2 <- list(
  "nobs"=length(ind_2), #number of noncensored observations;
  "ncen" = length(ind_2_cen),
  "censlim"=censorlimit_2,
  "t_obs" = shedding.data.wang$days[ind_2], #days after symptom onset;
  "c_obs" = shedding.data.wang$ct_value[ind_2], #log10 concentration of SARS-CoV-2 RNA (copies per mL);
  "t_cen" = shedding.data.wang$days[ind_2_cen] #days after symptom onset;
)
rm(ind_2) 
rm(ind_2_cen)
data_one_censored_2
