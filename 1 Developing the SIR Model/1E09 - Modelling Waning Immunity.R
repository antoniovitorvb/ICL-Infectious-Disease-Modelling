# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 9 - Modelling Waning Immunity

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
library(ggpubr) # subplots ggplot2

population <- 10^6

ini_state <- c(S = population - 1,
               I = 1,
               R = 0)

parameters <- c(beta = 0.4*365, # Infection Rate
                gamma = 0.2*365, # Recovery rate
                sigma = 1/10)  # Duration of immunity = 10 years

times <- seq(from = 0, to = 100, by = 2/365) # 100 years in timesteps of 2 days

sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R # Total population size
    
    lambda <- beta * I/N
    
    dS <- -lambda*S + sigma*R
    dI <- lambda*S - gamma*I
    dR <- gamma*I - sigma*R
    
    return(list(c(dS,dI,dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")
output_long$prevalence <- output_long$value/sum(ini_state)

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = prevalence*100,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence [%]") +
  labs(title = paste("Prevalence of S, I and R over time\n
       Sigma = ", parameters["sigma"],"/ year")) +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green")) +
  xlim(c(0,5))


# Changing parameters
parameters <- c(beta = 0.4*365, # Infection Rate
                gamma = 0.2*365, # Recovery rate
                sigma = 1/0.5)  # Duration of immunity = half a year

times <- seq(from = 0, to = 5, by = 2/365) # 5 years in timesteps of 2 days

output2 <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

output2_long <- melt(as.data.frame(output2), id = "time")
output2_long$prevalence <- output2_long$value/sum(ini_state)

graf2 <- ggplot(data = output2_long,
               aes(x = time,
                   y = prevalence*100,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence [%]") +
  labs(title = paste("Prevalence of S, I and R over time\n
       Sigma = ", parameters["sigma"],"/ year")) +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green"))

print(ggarrange(graf, graf2, nrow = 2))