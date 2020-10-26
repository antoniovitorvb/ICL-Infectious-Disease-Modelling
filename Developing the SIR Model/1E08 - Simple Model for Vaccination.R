# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 8 - Simple Model for Vaccination

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots

population <- 10^6
p <- 0  # Vaccine coverage

ini_state <- c(S = (1 - p) * (population - 1),
               I = 1,
               R = p*(population -1))

parameters <- c(beta = 0.319, # Infection Rate
                gamma = 0.2) # Recovery rate

times <- 0:730

sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R # Total population size
    
    lambda <- beta * I/N
    
    dS <- -lambda*S
    dI <- lambda*S - gamma*I
    dR <- gamma*I
    
    return(list(c(dS,dI,dR)))
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
  geom_line() +
  xlab("Time (days)") +
  ylab("Prevalence [%]") +
  labs(title = "Prevalence of S, I and R over time") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green"))

print(graf)
print(max(output_long$prevalence[output_long$variable == "I"])*100)