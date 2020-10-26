# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 4 - Modelling Additional Vaccine Effects

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
library(ggpubr) # subplots ggplot2

pop <- 10^6

parameters <- c(beta = 0.25, # infection rate
                gamma = 0.1, # natural recovery rate
                Cs = 0.3, # reduction in force of infection (lambda)
                Ci = 0.5,
                p = 0.6) # vaccination coverage

vacc_eff_S <- 1 - unname(parameters["Cs"]) # Vaccine effectivity in reducing susceptibility
vacc_eff_I <- 1 - unname(parameters["Ci"]) # Vaccine effectivity in reducing infectivity

ini_state <- c(S = unname((1-parameters["p"])) * (pop-1),
               I = 1,
               Iv = 0,
               R = 0,
               V = unname(parameters["p"])*(pop-1))

times <- 0:730 # 2 years

vacc_model <- function (time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I + Iv + R + V
    
    lambda <- (I + Ci*Iv) * beta/N
    
    dS <- -lambda*S
    dI <- lambda*S - gamma*I
    dIv <- Cs*lambda*V - gamma*Iv
    dR <- gamma * (I+Iv)
    dV <- -Cs*lambda*V
    
    return(list(c(dS, dI, dIv, dR, dV)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = vacc_model,
                            parms = parameters))

# What is the peak prevalence (number of infected people) with a vaccine coverage of 60%?
output$all_I <- output$I + output$Iv

peak_I <- round(max(output$all_I)) # number of infected people at its peak
print(paste("peak of infection =", peak_I))

graf <- ggplot(data = output,
               aes(x = time,
                   y = all_I)) +
  geom_line() +
  xlab("Time (days)") + ylab("Number of infected people") + 
  labs(title = paste("Combined leaky vaccine with coverage of", parameters["p"]*100,"%"))
# print(graf)

# For this vaccine, what is the minimum vaccination coverage required to interrupt transmission (Reff < 1)?
  # Reff = (1-p)*Ro + p*Cs*Ci*Ro
Ro <- unname(parameters["beta"]/parameters["gamma"])

pc <- (1 - 1/Ro) / (1 - unname(parameters["Cs"]*parameters["Ci"]))
print(paste("Minimum Vaccination Coverage =", round(pc*100,2),"%"))

parameters["p"] <- pc

new_state <- c(S = unname((1-parameters["p"])) * (pop-1),
               I = 1,
               Iv = 0,
               R = 0,
               V = unname(parameters["p"])*(pop-1))

output2 <- as.data.frame(ode(y = new_state,
                            times = times,
                            func = vacc_model,
                            parms = parameters))

graf2 <- ggplot(data = output2,
               aes(x = time,
                   y = I + Iv)) +
  geom_line() +
  xlab("Time (days)") + ylab("Number of infected people") + 
  labs(title = paste("Combined leaky vaccine with coverage of", round(parameters["p"]*100,2), "%"))

print(ggarrange(graf, graf2, nrow = 2))