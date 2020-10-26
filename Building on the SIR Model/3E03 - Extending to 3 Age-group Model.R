# INFECTIOUS DISEASE MODELLING
# Module 3 e-tivity 3 - Extending to 3 age-group model

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

pop <- 10^6
# Group 1 = children [20%]
# Group 2 = adults [65%]
# Group 3 = elders [15%]

ini_state <- c(S1 = 0.2*pop-1, S2 = 0.65*pop, S3 = 0.15*pop,
               I1 = 1, I2 = 0, I3 = 0,
               R1 = 0, R2 = 0, R3 = 0)

# Daily number of contacts (c11, c12, c13, c21...)
contacts <- c(7, 5, 1, 2, 9, 1, 1, 3, 2)
contact_matrix <- t(matrix(contacts, nrow = 3, ncol = 3)) # t() transposes the matrix for better writing contacts

parameters <- c(b = 0.05, # probability of infection of 5%
                contact_matrix = contact_matrix, # the age-specific average number of daily contacts
                gamma = 1/5) # rate of recovery

times <- seq(from = 0, to = 90, by = 0.1)

sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    age_groups <- nrow(contact_matrix)
    # assigning ini_state to each S, I and R group
    S <- state[1:age_groups]
    I <- state[(age_groups+1):(2*age_groups)]
    R <- state[(2*age_groups+1):(3*age_groups)]
    
    N <- S + I + R # N is a vector of length = age_groups
    
    # Defining the force of infection for each group
    lambda <- b * contact_matrix %*% as.matrix(I/N) # %*% is used to multiply matrices in R
    # lambda is also a vector of length = age_groups
    
    # ODEs:
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
# the output column names are adopted from the names we assigned in the ini_state vector
output_long <- melt(as.data.frame(output), id = "time")

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = value,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  labs(colour = "Compartment",
       x = "Time (days)",
       y = "Number of People")
print(graf)