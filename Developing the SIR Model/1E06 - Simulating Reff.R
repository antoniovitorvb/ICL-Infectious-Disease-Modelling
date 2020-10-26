# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 6 - Simulating Reff

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
library(ggpubr) # subplots ggplot2

# Modelling the epidemic
population <- 10^6

ini_state <- c(S = population - 1,
               I = 1,
               R = 0)
parameters <- c(beta = 1/2,
                gamma = 1/10)
obs_period <- 100
times <- 0:obs_period

sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R
    
    lambda <- beta * I/N
    
    dS <- -lambda*S
    dI <- lambda*S - gamma*I
    dR <- gamma*I
    
    return(list(c(dS,dI,dR)))
  })
}

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")

# Calculating the proportion in each compartment as a column in the long-format output
output_long$proportion <- output_long$value/sum(ini_state)

# Plots the proportion of people in the S, I and R compartments over time
graf <- ggplot(data = output_long,
               aes(x = time,
                   y = proportion*100, # %
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("proportion of the population [%]") +
  labs(title = paste("Proportion of Suscetible,Infected and Recovered over time when
                     \nBeta = ",parameters["beta"],", and Gamma = ", parameters["gamma"]),
       colour = "Compartment") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green")) +
  theme(legend.position = "bottom")

# Plots the Reff over time
output$reff <-  parameters["beta"] / parameters["gamma"] * output$S / sum(ini_state)

ReffG <- ggplot(data = output,
                aes(x = time,
                    y = reff)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Reff") +
  labs(title = "Effective reproduction number over time")

# subplotting each graphic with dim (2,1)
print(ggarrange(graf, ReffG, nrow = 2, ncol = 1))