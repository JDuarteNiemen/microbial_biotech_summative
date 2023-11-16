#Microbial biotech Summative

#Packages ----
library(tidyverse)
library(ggpubr)
library(janitor)
library(readxl)
library(patchwork)

#Import data ----
glucose <- read_excel("data/microbial_biotech_data_summative.xlsx", sheet="Sheet1")

fructose <- read_excel("data/microbial_biotech_data_summative.xlsx", sheet="Sheet2")

#Clean data ----
colnames(glucose)
colnames(fructose)

glucose <- glucose %>%
  rename("time"="Time (hr)",
         "glucose"="Glucose (mM)",
         "biomass"="Biomass (g/L)",
         "aca"="ACA (g/L)")

colnames(glucose)

fructose <- fructose %>%
  rename("time"="Time (hr)",
         "fructose"="Fructose (mM)",
         "biomass"="Biomass (g/L)",
         "aca"="ACA (g/L)")

colnames(fructose)

# Question 1. Calculating growth rates and doubling time ----
#calculate glucose growth rate
glucose <- glucose %>%
  mutate(u=log(biomass)/time)

#calculate fructose growth rate
fructose <- fructose %>%
  mutate(u=log(biomass)/time)

#plot glucose growth rate 
glucose_growth_rate_graph <- glucose %>%
  ggplot(aes(x=time,
             y=log(biomass),
             colour="red"))+
  geom_point()+
  stat_regline_equation()+
  theme_classic()

#plot fructose growth rate
fructose_growth_rate_graph <- fructose %>%
  ggplot(aes(x=time,
             y=log(biomass),
             colour="red"))+
  geom_point()+
  theme_classic()

#Set limits to straight line region (growth phase) and find gradient
linear_glucose_growth_rate_graph <- glucose_growth_rate_graph +
  xlim(3,18)+
  geom_smooth(method = "lm", se = FALSE, col = "red")+
  stat_regline_equation()

linear_glucose_growth_rate_graph

glucose_gradient <- 0.46

#Set limits to straight line region (growth phase) and find gradient
linear_fructose_growth_rate_graph <- fructose_growth_rate_graph +
  xlim(12,33)+
  geom_smooth(method="lm", se=FALSE, col="red")+
  stat_regline_equation()

linear_fructose_growth_rate_graph

fructose_gradient <-0.3

#Calculating doubling time for glucose
glucose_doubling_time <- log(2)/glucose_gradient

#Calculating doubling time for fructose
fructose_doubling_time <- log(2)/fructose_gradient


#Creating patchwork figure
glucose_growth_rate_graph + fructose_growth_rate_graph

#Question 2 ----
#setting values for log phase of each reaction
glucose
log_glucose_product_initial <- glucose[2,4]
log_glucose_product_final <- glucose[7,4]

log_glucose_substrate_initial <- glucose[2,2] 
log_glucose_substrate_final <- glucose[7,2]

fructose
log_fructose_product_initial <- fructose[5,4]
log_fructose_product_final <- fructose[12,4]

log_fructose_substrate_initial <- fructose[5,2]
log_fructose_substrate_final <- fructose[12,2]

#setting values for stationary phase
glucose
sta_glucose_product_initial <- glucose[8,4]
sta_glucose_product_final <- glucose[17,4]

sta_glucose_substrate_initial <- glucose[8,2]
sta_glucose_substrate_final <- glucose[17,2]

fructose
sta_fructose_product_initial <- fructose[13,4]
sta_fructose_product_final <- fructose[17,4]

sta_fructose_substrate_initial <- fructose[13,2]
sta_fructose_substrate_final <- fructose[17,2]


#calculating yields
#glucose log phase yield
glucose_log_yield <- log_glucose_product_final - log_glucose_product_initial / log_glucose_substrate_initial - log_glucose_substrate_final

#fructose log phase yield
fructose_log_yield <- log_fructose_product_final - log_fructose_product_initial / log_fructose_substrate_initial - log_fructose_substrate_final

#glucose stationary phase yield
glucose_sta_yield <- sta_glucose_product_final - sta_glucose_product_initial / sta_glucose_substrate_intial - sta_glucose_substrate_final

#fructose stationary phase yield
fructose_sta_yield <- sta_fructose_product_final - sta_fructose_product_initial / sta_fructose_substrate_initial - sta_glucose_substrate_final

cat("Glucose log phase yield:", glucose_log_yield)
cat("Fructose log phase yield:", fructose_log_yield)

cat("Glucose stationary phase yield:", glucose_sta_yield)
cat("Fructose stationary phase yield:", fructose_sta_yield)

#this method hasnt quite gone right, do it the same way as in the formative bit long but oh well
