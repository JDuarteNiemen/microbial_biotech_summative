# MICROBIAL BIOTECHNOLOGY FORMATIVE DATA HANDLING ----
##Packages ----
library(tidyverse)
library(ggpubr)
library(png)
### Import data ----
yeast <- read.csv("data/microbial_biotech_data.csv")
monod <- read.csv("data/microbial_biotech_data_2.csv")
q34 <- read.csv("data/estimate_manipulations_results.csv")
estimate_graph <- readPNG("data/estimates.png")
dilute <- read.csv("data/part_c_table.csv")
#cleaning data ----
colnames(yeast)
colnames(monod)

#rename coloumns
yeast <- yeast %>%
  rename("time"="Time..hr.",
         "fructose"="Fructose..mM.",
         "biomass"="Biomass..g.L.",
         "acetate"="Acetate..mM.")

monod <- monod %>%
  rename("glucose_conc"="Glucose.concentration..mM.",
         "growth_rate"="Growth.rate..h.1.")

#remove first row of NA
yeast <- yeast %>%
  drop_na()

#exploratory analysis ----
head(yeast)

#calculating growth rate, biomass/time ----

yeast <- yeast %>%
  mutate(u=log(biomass)/time)

yeast

#plotting growth rate

growth_rate_graph <- yeast %>%
  ggplot(aes(x=time,
             y=log(biomass),
             colour="red"))+
  geom_point()+
  geom_line()+
theme_classic()
#exponential phase = hour 2-16, stationary = 18-24

u_calc <- yeast %>%
  ggplot(aes(x=time,
             y=log(biomass),
             colour="red"))+
  xlim(2,16)+
  geom_point()+
  geom_line()+
  stat_regline_equation()+
  theme_classic()
#gradient = 0.35 = u
u_calc

#calculate doubling time ----
  
doubling_time <- log(2)/0.35 
  
#Calculating yields ----
  
#calculating yield coefficients for acetate product 
ysx_exp_coef <-  0.35/(150-114)
ysp_exp_coef <-  (yeast[9,4]-yeast[2,4])/(150-114)

#yield of acetate in exponential and stationary phase
yield_exp_ace <- (yeast[9,4]-yeast[2,4])/(yeast[2,2]-yeast[9,2])
yield_sta_ace <- (yeast[13,4]-yeast[10,4])/(yeast[10,2]-yeast[13,2])

#yield of biomass in exponential and stationary phase
yield_exp_bio <- (yeast[9,3]-yeast[2,3])/(yeast[2,2]-yeast[9,2])
yield_sta_bio <- (yeast[13,3]-yeast[10,3])/(yeast[10,2]-yeast[13,2])
#these are calculated using the cells in the yeast table specifying the initial and final volumes for either biomass, substrate(fructose) or product(acetate)
#following the equation product produced/substrate used.

#compiling different yields into a table for ease of access.
yields <- c("yield of acetate in exponential","yield of Acetate in Stationary","yield of Biomass in exponential","yield of Biomass in stationary")
values <-c(yield_exp_ace, yield_sta_ace, yield_exp_bio, yield_sta_bio)
yield_table <- data.frame("yield"=yields, "value(g/mM)"=values)


#Visualising data ----

monod #using 2nd dataset

monod_graph <- monod %>%
  ggplot(aes(x=glucose_conc,
             y=growth_rate,
             colour="red"))+
           geom_point()+
           theme_classic()
    
print(monod_graph)


nls_monod_graph <- monod_graph +
  xlab('Glucose Concentration (mM)')+
  ylab(bquote('Growth Rate' (Hr^-1)))+
 stat_smooth(se=FALSE)
#consider adding umax and 1/2umax and ks values using annotate()
           

print(nls_monod_graph)
   
#Calculating estimated Ks and umax ---- 

# Define the Monod model
monod_model <- function(glucose_conc, umax, ks) {
  umax * glucose_conc / (ks + glucose_conc)  
}  

#set values for umax and ks
initial_params <- list(umax = 0.7, ks = 0.5)


#performs non-linear regression
fit <- nls(growth_rate ~ monod_model(glucose_conc, umax, ks),
           data = monod,
           start = initial_params,
           control = nls.control(maxiter = 100))


# Extract the estimated parameters
estimated_umax <- coef(fit)["umax"]
estimated_ks <- coef(fit)["ks"]

# Print the results
cat("Estimated umax:", estimated_umax, "\n")
cat("Estimated ks:", estimated_ks, "\n")



#Lineweaver-Burke Plot ----

#caluclating reciprocals for glucose concentration (S) and growth rate (u)
lwbp <- monod %>%
  mutate(oneglucose_conc = (1/glucose_conc),
         onegrowth_rate = (1/growth_rate))




#plot lineweaver burke

lwbpg <- lwbp %>%
  ggplot(aes(x=oneglucose_conc,
             y=onegrowth_rate,
             colour="red"))+
  geom_line()+
  geom_point()+
  stat_regline_equation()+
  xlab('Glucose Concentration (mM)')+
  ylab(bquote('Growth Rate' (Hr^-1)))+
  theme_classic()


#umax is 1/yintercept
lwbumax <- 1/0.97

#Ks is gradient/1/yintercept
lwbks <- 3.1/lwbumax

#creating stupid new data frames ----
  
#create dataframe for question 3
monod_3 <- data.frame(monod)

#change value from 0.14 to 0.168
monod_3[1,2]=0.168

#create dataframe for question 4
monod_4 <- data.frame(monod)

monod_4[8,2]=0.84

#calculating new umax and ks values for question 3 ----

# Define the Monod model
monod_model <- function(glucose_conc, umax, ks) {
  umax * glucose_conc / (ks + glucose_conc)  
}  

#set values for umax and ks
initial_params <- list(umax = 0.7, ks = 0.5)


#performs non-linear regression
fit_3 <- nls(growth_rate ~ monod_model(glucose_conc, umax, ks),
           data = monod_3,
           start = initial_params,
           control = nls.control(maxiter = 100))


# Extract the estimated parameters
estimated_umax_3 <- coef(fit_3)["umax"]
estimated_ks_3 <- coef(fit_3)["ks"]

# Print the results
cat("Estimated Non-linear umax:", estimated_umax_3, "\n")
cat("Estimated Non-linear ks:", estimated_ks_3, "\n")

#_________________________________________________________

#Lineweaver burke method 
lwbp_3 <- monod_3 %>%
  mutate(oneglucose_conc = (1/glucose_conc),
         onegrowth_rate = (1/growth_rate))


#plot lineweaver burke
lwbpg_3 <- lwbp_3 %>%
  ggplot(aes(x=oneglucose_conc,
             y=onegrowth_rate,
             colour="red"))+
  geom_point()+
  stat_regline_equation()+
  geom_smooth(method="lm",
              se=FALSE)+
  theme_classic()

lwbpg_3

#umax is 1/yintercept
lwbumax_3 <- 1/1.2

#Ks is gradient/1/yintercept
lwbks_3 <- 3.1/lwbumax_3

cat("Estimated lineweaver-burke umax:", lwbumax_3, "\n")
cat("Estimated lineweaver-burke Ks:", lwbks_3, "\n")

#calculating new umax and ks values for question 4 ----
# Define the Monod model
monod_model <- function(glucose_conc, umax, ks) {
  umax * glucose_conc / (ks + glucose_conc)  
}  

#set values for umax and ks
initial_params <- list(umax = 0.7, ks = 0.5)


#performs non-linear regression
fit_4 <- nls(growth_rate ~ monod_model(glucose_conc, umax, ks),
           data = monod_4,
           start = initial_params,
           control = nls.control(maxiter = 100))


# Extract the estimated parameters
estimated_umax_4 <- coef(fit_4)["umax"]
estimated_ks_4 <- coef(fit_4)["ks"]

# Print the results
cat("Estimated Non-linear umax:", estimated_umax_4, "\n")
cat("Estimated Non-linear ks:", estimated_ks_4, "\n")

#_________________________________________________________

#Lineweaver burke method 
lwbp_4 <- monod_4 %>%
  mutate(oneglucose_conc = (1/glucose_conc),
         onegrowth_rate = (1/growth_rate))


#plot lineweaver burke
lwbpg_4 <- lwbp_4 %>%
  ggplot(aes(x=oneglucose_conc,
             y=onegrowth_rate,
             colour="red"))+
  geom_point()+
  stat_regline_equation()+
  geom_smooth(method="lm",
              se=FALSE)+
  theme_classic()

lwbpg_4

#umax is 1/yintercept
lwbumax_4 <- 1/0.92

#Ks is gradient/1/yintercept
lwbks_4 <- 3.1/lwbumax_4

cat("Estimated lineweaver-burke umax:", lwbumax_4, "\n")
cat("Estimated lineweaver-burke Ks:", lwbks_4, "\n")












#Part C ----
dilute

#changing coloumn names
dilute <- dilute %>%
  rename("dilution_rate"="Dilution.rate.D..h.1.",
         "biomass"="Biomass.x..g.L.",
         "fructose_conc"="Fructose.concentration..mM.")

colnames(dilute)

head(dilute)


dilute %>%
  ggplot(aes(x=fructose_conc,
             y=biomass,
             colour=dilution_rate))+
  geom_smooth()+
  theme_classic()


#changing fructose concetration to amount consumed


dilute

#calculating yield productfinal-productinitial/substrateinitial-substratefinal

biomass_yield_1 <- (dilute[1,2]-0)/(10-dilute[1,3])
biomass_yield_2 <- (dilute[2,2]-0)/(10-dilute[2,3])
biomass_yield_3 <- (dilute[3,2]-0)/(10-dilute[3,3])
biomass_yield_4 <- (dilute[4,2]-0)/(10-dilute[4,3])
biomass_yield_5 <- (dilute[5,2]-0)/(10-dilute[5,3])

biomass_yields <- c(biomass_yield_1, biomass_yield_2, biomass_yield_3, biomass_yield_4, biomass_yield_5)

optimal_dilution_rate <- data.frame(dilute$dilution_rate, biomass_yields)

optimal_dilution_rate

#optimal dilution rate is 0.2


