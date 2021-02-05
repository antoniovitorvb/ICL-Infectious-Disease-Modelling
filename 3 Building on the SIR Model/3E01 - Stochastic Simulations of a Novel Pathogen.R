# INFECTIOUS DISEASE MODELLING
# Module 3 e-tivity 1 - Stochastic Simulations of a Novel Pathogen

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

ini_state <- c(S = 10^6,
               I = 1,
               R = 0)

parameters <- c(beta = 0.3,
                gamma = 0.4)

times <- 0:100

sir_model <- function(time, state, parameters){
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    
    lambda <- beta *I/N
    
    dS <- -lambda * S
    dI <- lambda*S - gamma*I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

graf <- ggplot() +
  geom_line(data = output,
            aes(x = time,
                y = I)) +
  labs(title = paste("Deterministic Model output for R0 =",
                     parameters["beta"]/parameters["gamma"]),
       x = "Time (days)",
       y = "Prevalende of Infection") +
  ylim(c(0,50))
print(graf)

# Simulating Stochastic algorithm:
library("GillespieSSA")

# Defining model and input
a <- c("beta*S*I/10^6", "gamma*I")
nu <- matrix(c(-1, 1, 0, 0, -1, 1),
             nrow = 3,
             ncol = 2)
tf <- 100

# simulating
sir_out <- ssa(ini_state, a, nu, parameters, tf=tf, simName = "SIR")
#View(sir_out)

while (sir_out$stats$nSteps == 1) {
  sir_out <- ssa(ini_state, a, nu, parameters, tf=tf, simName = "SIR")
  # Record number of simulations
  n_sims <- n_sims + 1
}

stoch_plot <- stoch_plot +
  geom_line(data = as.data.frame(sir_out$data),
            aes(x = t,
                y = I),
            col = sample(rainbow(20), 1)) +
  labs(title = paste("Stochastic model output for R0",
                     parameters["beta"]/parameters["gamma"]),
       subtitle = paste(n_sims, "simulations"),
       x = "Time (days)",
       y = "Prevalende of Infection") +
  xlim(c(0,20)) +
  ylim(c(0,20))
print(stoch_plot)