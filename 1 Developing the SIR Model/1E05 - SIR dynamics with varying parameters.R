# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 5 - SIR dynamics with varying parameters

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

SIR <- function(population, beta, gamma, duration = 90){

population <- 10^6

ini_state <- c(S = population - 1,
               I = 1,
               R = 0)
parameters <- c(beta = beta,
                gamma = gamma)

times <- 0:duration

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

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))
output_long <- melt(as.data.frame(output), id = "time")

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = value,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Number of people") +
  labs(title = paste("Number of Susceptible, Infected and Recovered over time when
       \nBeta = ",round(parameters["beta"],3),", and Gamma = ", round(parameters["gamma"],3)),
       colour = "Compartment") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("black", "blue", "red"))

print(graf)

}