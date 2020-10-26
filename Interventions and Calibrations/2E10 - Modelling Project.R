# INFECTIOUS DISEASE MODELLING
# Module 2 e-tivity 10 - Modelling Project

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

# Model for the Source Conutry
pop <- 1000
p <- 0
ini_state <- c(S = (1-p) * (pop-1),
               I = 1,
               M = 0,
               R = p*(pop-1))

parameters <- c(beta = 0.31,
                gamma = 1/5*0.97,
                mu = 0.03*0.194/(1-0.03))

times <- seq(from = 0, to = 200, by = 1)

sir_model <- function(times, state, parameters){
  
  with(as.list(c(state,parameters)),{
    N <- S + I + M + R
    lambda <- beta*I/N
    
    dS <- -lambda*S
    dI <- lambda*S - (mu+gamma)*I
    dM <- mu*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dM, dR)))
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
  labs(title = paste("beta =", signif(parameters["beta"],3),
                     "gamma =", signif(parameters["gamma"],3),
                     "mu =", signif(parameters["mu"],3))) +
  geom_line()
print(graf)
print(max(output_long$prevalence[output_long$variable == "I"]))