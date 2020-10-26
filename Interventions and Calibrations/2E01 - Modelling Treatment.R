# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 1 - Modelling Treatment

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

pop <- 300000
parameters <- c(beta = 0.6,
                gamma = 1/5,
                h = 1/4,
                gT = 1/1.25)

ini_state <- c(S = pop - 1,
               I = 1,
               T = 0,
               R = 0)

times <- 0:60

treatment_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R + T
    
    lambda <- beta * (I+T) / N
    
    dS <- -lambda*S
    dI <- lambda*S - (gamma+h)*I
    dT <- h*I - gT*T
    dR <- gamma*I + gT*T
    
    return(list(c(dS,dI,dT,dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = treatment_model,
                            parms = parameters))
output_long <- melt(as.data.frame(output), id = "time")

#q1 How many people are infected at the peak of the epidemic?
output$IT <- output$I + output$T
q1 <- round(max(output$IT))

#q2 How rapidly does treatment need to be initiated in order to interrupt transmission, i.e. to bring R0 below 1?
  # Ro <- beta/(gamma+h) + h/(gamma+h) * beta/gT
h <- parameters["gT"] * (parameters["gamma"]-parameters["beta"]) / (parameters["beta"]-parameters["gT"])
h <- unname(h)  #remove name of parameters