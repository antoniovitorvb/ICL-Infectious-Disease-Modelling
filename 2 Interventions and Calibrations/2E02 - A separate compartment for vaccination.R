# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 2 - A separate compartment for vaccination

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

pop <- 10^6

parameters <- c(beta = 0.6, # infection rate
                gamma = 1/10, # natural recovery rate
                p = 0.15) # vaccination coverage

ini_state <- c(S = unname((1-parameters["p"]) * (pop-1)),
               I = 1,
               R = 0, # naturally recovered
               V = unname(parameters["p"] * (pop-1))) # already vaccinated people 

Ro <- unname(parameters["beta"]/parameters["gamma"])
times <- 0:90

vacc_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R + V # Total population size
    
    lambda <- beta * I/N
    
    dS <- -lambda*S
    dI <- lambda*S - gamma*I
    dR <- gamma*I
    dV <- 0
    
    return(list(c(dS,dI,dR,dV)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = vacc_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id="time")
output_long$prevalence <- output_long$value/sum(ini_state)

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = prevalence*100,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Prevalence [%]") +
  labs(title = "Prevalence of S, I and R over time") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered", "Vaccinated"), values = c("blue", "red", "green", "black"))

print(graf)