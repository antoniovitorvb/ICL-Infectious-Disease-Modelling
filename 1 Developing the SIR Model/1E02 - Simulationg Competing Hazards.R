# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 2 - Simulating Competing Hazards

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

ini_state <- c(I = 10^6, R = 0, M = 0)
parameters <- c(gamma = 0.1,  mu = 0.2)

obs_period <- 4*7
times <- 1:obs_period

cohort_model <- function(times, state, parameters){
  
  with(as.list(c(state,parameters)),{
    dI <- -(mu+gamma)*I
    dR <- gamma*I
    dM <- mu*I
    
    return(list(c(dI,dR,dM)))
  })
}

output <- as.data.frame(ode(y = ini_state, 
                            times = times, 
                            func = cohort_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")     

graf <- ggplot(data = output_long,   # specify object containing data to plot
       aes(x = time, 
           y = value, 
           colour = variable, 
           group = variable)) +       # assign columns to axes and groups
  geom_line() +                       # represent data as lines
  xlab("Time (days)")+                # add label for x axis
  ylab("Number of people") +          # add label for y axis
  labs(title = paste("Number infected and recovered over time when
  \ngamma =",parameters["gamma"],"days^-1
                     \nmu = ",parameters["mu"],"days^-1"),
       colour = "Legend")        # add legend title

# ggplot automatically adds the new compartment as a separate colour/group 
# because it is contained within the output dataframe
print(graf)