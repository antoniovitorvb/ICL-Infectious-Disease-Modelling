# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 6 - Manual Calibration of a SIR Model part 2

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

data <-  data.frame(time = 1:14,
                    number_infected = c(3,8,26,76,225,298,258,233,189,128,68,29,14,4))

ini_state <- c(S = 762,
               I = 1,
               R = 0)

parameters <- c(beta = 1.7,
                gamma = 0.44)

times <- seq(from = 0, to = length(data$number_infected-1), by = 0.1)

sir_model <- function (time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    N <- S + I + R
    
    lambda <- beta*I/N
    
    dS <- -lambda*S
    dI <- lambda*S - gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

graf <- ggplot() +
  geom_line(data = output,
            aes(x = time, y = I)) + # Plot of the estimated parameters
  geom_point(data = data,
             aes(x = time, y = number_infected),
             colour = "red") + # Plot of data
  xlab("Time (days)") + ylab("Number of infected people") +
  labs(title = paste("Manual calibration for SIR Model with\n
                     beta =", round(parameters["beta"], 2), "and gamma =", round(parameters["gamma"], 2)))

print(graf)