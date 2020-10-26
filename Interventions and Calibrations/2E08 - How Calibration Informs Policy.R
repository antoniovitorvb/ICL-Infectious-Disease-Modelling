# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 8 - How Calibration Informs Policy

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

dataset <- read.csv("2E08_dataset.csv")

ini_state <- c(S = 499,
               I = 1,
               R = 0)

times <- 0:200 # 200 day outbreak

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

# DISTANCE FUNCTION 
SSQ <- function(parameters, dat) {  # takes as inputs the parameter values and dataset
  
  beta <- parameters[1]
  gamma <- parameters[2]
  
  # Simulate the model with initial conditions and timesteps defined above,
  # and parameter values from function call
  result <- as.data.frame(ode(y = ini_state, 
                              times = times, 
                              func = sir_model,
                              parms = c(beta = beta,      # ode() takes the values for
                                        # beta and gamma extracted from
                                        gamma = gamma)))  # the "parameters" input argument
  # of the SIR_SSQ() function
  
  # Calculate the sum of squares by comparing the model output with the matching
  # datapoints: This involves, for each timepoint with available data, 
  # calculating the difference between the number of infections
  # predicted by the model and the observed number of infections, squaring all these 
  # differences,and taking the sum of all squared differences
  ssq <- sum((result$I[result$time %in% dat$time] - dat$number_infected)^2)
  
  return(ssq)
}

output <- optim(par = c(0.1, 0.1), # choosen sensible starting points for parameters
                fn = SSQ,
                dat = dataset)

parameters <- c(beta = output$par[1],
                gamma = output$par[2])

result <- as.data.frame(ode(y = ini_state, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))

graf <- ggplot() +
  labs(title = paste("Infection Curve for\n", 
                     "beta =", round(parameters["beta"], 2), "and gamma =", round(parameters["gamma"], 2)),
       x = "Time (days)",
       y = "Number of Infected People",
       colour = "Compartment") +
  
  geom_line(data = result,
            aes(x = time, y = I,
                colour = "SSQ")) +
  geom_point(data = dataset,
             aes(x = time, y = number_infected,
                 colour = "Data"))
print(graf)