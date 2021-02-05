# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 3 - SIR Model with a constant Force of Infection

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

population <- 10^6

ini_state <- c(S = population - 1, I = 1, R = 0)
parameters <- c(lambda = 0.2, gamma = 1/10)
obs_period <- 60

times <- 0:obs_period

# This SIR Model doesn't consider reinfectations or deaths in a closed population

sir_model <- function(times, state, parameters) {
  
  with(as.list(c(state,parameters)),{
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

output_long <- melt(as.data.frame(output), id = "time") # 1st converts de data to a long format

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = value,
                   colour = variable,
                   group = variable)) + 
  geom_line() + 
  xlab("Time (days)") + 
  ylab("Number of people") + 
  labs(title = paste("Number of Susceptible, Infected and Recovered over time when 
                     \nlambda =",parameters["lambda"],"days^-1"),
       colour = "Legend")
print(graf)