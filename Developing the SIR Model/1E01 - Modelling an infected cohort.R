# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 1 - Modelling an Infected Cohort

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

ini_infected <- 1000000
ini_recovered <- 0
gamma <- 1/10  # Recovery Rate
obs_period <- 28  # 4 weeks

ini_state <- c(I = ini_infected, R = ini_recovered)

time <- 0:obs_period  # OR seq(from = 0, to = obs_period, by = 1)

# ini_state
# gamma
# time

cohort_model <- function(time, state, gamma){
  
  with(as.list(c(state, gamma)), {
    dI <- -gamma*I
    dR <- gamma*I
    
    return(list(c(dI,dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = time,
                            func = cohort_model,
                            parms = gamma))

recover <- output[output$time == 28, c("time","R")]
recover

recover_per <- output[output$time==28,"R"]/(output[output$time==28,"R"]+output[output$time==28,"I"])
print(paste(as.character(signif(recover_per*100,4)),"%"))

# First turn the output dataset into a long format, 
# so that the number in each compartment at each timestep
# are all in the same column
output_long <- melt(as.data.frame(output), id = "time")

# Plot the number of people in each compartment over time
graf <- ggplot(data = output_long,          # specify object containing data to plot
       aes(x = time, y = value, colour = variable, group = variable)) +          # assign columns to axes and groups
  geom_line() +                     # represent data as lines
  xlab("Time (days)") +              # add label for x axis
  ylab("Number of people")  +       # add label for y axis

labs(title = paste("Number infected and recovered over time when gamma =",
                   gamma,"days^-1")) # add title

print(graf)