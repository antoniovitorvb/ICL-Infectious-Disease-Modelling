# INFECTIOUS DISEASE MODELLING
# Module 3 e-tivity 4 - Coding a Vector-Borne Disease (VBD) Model

library(deSolve)  # solves the differential equations model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

Nh <- 10000
Nv <- 20000

ini_state <- c(Sh = Nh - 0.0028*Nh, # infection prevalence of 0.28% in humans and 0.057% in mosquitoes
               Ih = 0.0028*Nh,
               Rh = 0,
               Sv = Nv - 0.00057*Nv,
               Iv = 0.00057*Nv)

parameters <- c(a = 1,
                bv = 0.4,
                bh = 0.4,
                muv = 1/4,
                r = 0.167)

times <- 0:90

vbd_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    Nh <- Sh + Ih + Rh
    Nv <- Sv + Iv
    
    # ODEs for hosts
    dSh <- -a*bh*Sh*Iv/Nh
    dIh <- a*bh*Sh*Iv/Nh - r*Ih
    dRh <- r*Ih
    
    # ODEs for vectors
    dSv <- muv*Nv - a*bv*Sv*Ih/Nh - muv*Sv
    dIv <- a*bv*Sv*Ih/Nh - muv*Iv
    
    return(list(c(dSh, dIh, dRh, dSv, dIv)))
  })
}

output <- as.data.frame(ode(y = ini_state, 
                            times = times, 
                            func = vbd_model,
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
       y = "Number of Hosts/Vectors")
# print(graf)

# Based on the plot, how do assumptions of the mosquito biting activity affect the human infection prevalence?
a_values <- seq(0,1,by=0.1)                         # a vector of values for the biting rate, ranging from 0 to 1 per day
out_list <- vector(length(a_values), mode = "list")
out_R0 <- vector(length(a_values), mode = "list") # create empty lists to store the output for each value of a
for (i in seq_along(a_values)) {
  out_list[[i]] <- as.data.frame(ode(y = ini_state,
                                     times = times,
                                     func = vbd_model,
                                     parms = c(a = a_values[i],
                                               parameters[c("bv", "bh", "muv", "r")])))
  out_R0[[i]] <- a_values[i]^2 * parameters["bh"]*parameters["bv"]/(parameters["r"]*parameters["muv"]) * Nv/Nh
}
names(out_list) <- a_values # Rename list elements with te corresponding bitting rate
names(out_R0) <- a_values
# Extract the infected host column from the list and creating a dataframe by time
out_Ih <- cbind(time = out_list[[1]]$time,
                sapply(out_list, 
                       "[[", "Ih")) # fill in the name of your column for infected hosts
out_Ih_long <- melt(as.data.frame(out_Ih), id = "time")
out_R0 <- as.data.frame(out_R0)

Ih_a <- ggplot(data = out_Ih_long,
               aes(x = time,
                   y = value,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  labs(colour = "Bitting rate [1/day]",
       x = "Time (days)",
       y = "Human Infection Prevalence")
print(Ih_a)
max_R0 <- max(out_R0[out_R0 < 1])
print(paste("R0 =", round(max_R0, 2), "for a =", names(out_R0[length(out_R0[out_R0 < 1])])))
      
      
      