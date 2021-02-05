# INFECTIOUS DISEASE MODELLING
# Module 2 e-tivity 9 - Performing Maximum Likelihood Estimation

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

reported <- read.csv("2E09_reporteddata.csv")
fulldata <- read.csv("2E09_fulldata.csv")

ini_state <- c(S = 762,
               I = 1,
               R = 0)

times <- seq(from = 0, to = 14, by = 0.1)

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

log_pois <- function (parameters, dat) {
  beta <- parameters[1]
  gamma <- parameters[2]
  
  result <- as.data.frame(ode(y = ini_state,
                              times = times,
                              func = sir_model,
                              parms = c(beta = beta,
                                        gamma = gamma)))
  # Calculate log-likelihood using code block 4 from the previous etivity, accounting for the reporting rate of 60%:
  LL <- sum(dpois(x = dat$number_reported, lambda = 0.6 * result$I[result$time %in% dat$time], log = TRUE))
  
  return(LL)
}

# Calculating the Likelyhood
LH <- optim(par = c(1.7, 0.1), # starting values for beta and gamma
            fn = log_pois,
            dat = reported,
            control = list(fnscale = -1)) # tells optim() to look for the maximum number instead of the minimum (the default)

parameters <- c(beta = LH$par[1],
                gamma = LH$par[2])

output <- as.data.frame(ode(y = ini_state, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))

graf <- ggplot() +
  labs(title = paste("Model fit to epidemic curve with\n
                     beta =", round(parameters["beta"], 2), 
                     "and gamma =",round(parameters["gamma"],2)),
       x = "Time (days)",
       y = "Number of people") +
  scale_color_manual(values = c("black", "green", "red")) +
  
  geom_line(data = output,
            aes(x = time, y = I,
                colour = "Model")) +
  geom_point(data = fulldata,
             aes(x = time, y = number_infected,
                 colour = "Total cases")) +
  geom_point(data = reported,
             aes(x = time, y = number_reported,
                 colour = "Reported cases"))
print(graf)