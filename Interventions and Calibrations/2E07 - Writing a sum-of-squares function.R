# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 7 - Writing a sum-of-squares function

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

# FULL SOLUTION

data <-  data.frame(time = 1:14,
                    I = c(3,8,26,76,225,298,258,233,189,128,68,29,14,4))

ini_state <- c(S = 762,
               I = 1,
               R = 0)

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

# Sum-of-Squares function measures de distance between data and output values
SSQ <- function (parameters, dat) {
  with(as.list(c(parameters)), {
    result <- as.data.frame(ode(y = ini_state,
                                times = times,
                                func = sir_model,
                                parms = parameters))
    dat <- na.omit(dat)
    
    # delta = [model - data]^2
    deltas <- (result$I[result$time %in% dat$time] - dat$I)^2
    
    SSQ <- c(delta = sum(deltas))
    
    return(list(c(result,SSQ)))
  })  
}

parameters <- c(beta = 1.15,
                gamma = 0.02)
parameters2 <- c(beta = 1.7,
                 gamma = 0.45)

times <- seq(from = 0, to = 14, by = 0.1)

ssq1 <- as.data.frame(SSQ(parameters = parameters,
                          dat = data))
print(paste("ssq1 =", round(ssq1$delta[length(ssq1$time-1)],2), "where beta =", round(parameters["beta"],2), "and gamma =", round(parameters["gamma"],2)))

ssq2 <- as.data.frame(SSQ(parameters = parameters2,
                          dat = data))
print(paste("ssq2 =", round(ssq2$delta[length(ssq2$time-1)],2), "where beta =", round(parameters2["beta"],2), "and gamma =", round(parameters2["gamma"],2)))



graf <- ggplot() +
  labs(title = paste("Different plots for\n",
       "beta1 =", parameters["beta"], "and gamma1 =", parameters["gamma"],
                     "\n beta2 =", parameters2["beta"], "and gamma2 =", parameters2["gamma"]),
       x = "Time (days)",
       y = "Number of Infected People",
       colour = "Compartment") +
  geom_point(data = data,
             aes(x = time, y = I,
             colour = "red")) +
  geom_line(data = ssq1,
            aes(x = time, y = I,
            colour = "blue")) +
  geom_line(data = ssq2,
            aes(x = time, y = I,
            colour = "green")) +
  scale_color_identity(guide = "legend") +
  scale_colour_manual(labels = c("SSQ 1", "SSQ 2", "Data"),
                     values = c("blue", "green", "red"))
print(graf)