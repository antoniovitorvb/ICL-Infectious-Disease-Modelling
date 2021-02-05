# INFECTIOUS DISEASE MODELLING
# Module 3 e-tivity 2 - Developing an 2 age-structured group model

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

pop <- 10^6
N1 <- 0.2 * pop # 20% are children
N2 <- pop - N1 # the rest are adults

ini_state <- c(S1 = N1-1,
               I1 = 1,
               R1 = 0,
               S2 = N2,
               I2 = 0,
               R2 = 0)

parameters <- c(b = 0.05, # probability of infection of 5%
                c11 = 7,
                c12 = 6,
                c21 = 1,
                c22 = 10,
                gamma = 1/5) # Rate of recovery

times <- seq(from = 0, to = 90, by = 0.1)

sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N1 <- S1 + I1 + R1
    N2 <- S2 + I2 + R2
    
    # Defining the force of infection
    # acting on susceptible children:
    lambda1 <- b*c11*I1/N1 + b*c12*I2/N2
    # acting on susceptible adults:
    lambda2 <- b*c21*I1/N1 + b*c22*I2/N2
    
    # ODE for children:
    dS1 <- -lambda1 * S1
    dI1 <- lambda1*S1 - gamma*I1
    dR1 <- gamma*I1
    
    # ODE for adults:
    dS2 <- -lambda2* S2
    dI2 <- lambda2*S2 - gamma*I2
    dR2 <- gamma*I2
    
    return(list(c(dS1, dI1, dR1, dS2, dI2, dR2)))
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
  labs(colour = "Compartment",
       x = "Time (days)",
       y = "Number of people")
print(graf)

# What was the cumulative incidence of infection during this epidemic?
total_inc <- (output$S1[1] - output$S1[nrow(output)]) + (output$S2[1] - output$S2[nrow(output)])
print(paste("Total cumulative incidence =", round(total_inc,0)))

# What proportion of those infections occurred in children?
chil_prop <- (output$S1[1] - output$S1[nrow(output)]) / total_inc
print(paste("Proportion of infections in children = ", round(chil_prop*100,1),"%"))

# Which age group was most affected by the epidemic?
percentage <- c(children = (output$S1[1] - output$S1[nrow(output)]) / N1,
                adults = (output$S2[1] - output$S2[nrow(output)]) / N2)
print(paste("Most affected group was:", names(percentage[percentage == max(percentage)]),
            "with", round(percentage[percentage == max(percentage)]*100, 1),"%"))