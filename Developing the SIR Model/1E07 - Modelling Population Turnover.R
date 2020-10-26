# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 7 - Modelling Population Turnover

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
library(ggpubr) # subplots ggplot2

population <- 10^6

ini_state <- c(S = population - 1,
               I = 1,
               R = 0)
parameters <- c(beta = 0.4*365, # per year
                gamma = 0.2*365,
                mu = 1/70, # 70 years lifespan
                b = 1/70)

times <- seq(from = 0, to = 400, by = 3/365)   # from 0 to 400 YEARS in a 3 day interval

sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R
    
    lambda <- beta * I/N
    
    dS <- -lambda*S - mu*S + b*N
    dI <- lambda*S - gamma*I - mu*I
    dR <- gamma*I - mu*R
    
    return(list(c(dS,dI,dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

graf <- ggplot(data = output,
               aes(x = time,
                   y = I)) +
  geom_line() +
  xlab("Time (1 year)") +
  ylab("Number of Infected People") +
  labs(title = "Epidemic curve in the first year after introduction of an infected case") +
  xlim(c(0,1))

output_long <- melt(as.data.frame(output), id = "time")

graf2 <- ggplot(data = output_long,
               aes(x = time,
                   y = value,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Number of Infected People") +
  labs(title = "Epidemic curve over 4 generations") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "black", "green")) +
  theme(legend.position = "bottom")


# Changing lifespan for 4 weeks
parameters <- c(beta = 0.4,
                gamma = 0.2,
                mu = 1/(4*7), # 4 week lifespan
                b = 1/(4*7))
obs_period <- 365 # 1YEAR
times2 <- 0:obs_period

output2 <- as.data.frame(ode(y = ini_state,
                            times = times2,
                            func = sir_model,
                            parms = parameters))

output_long2 <- melt(as.data.frame(output2), id = "time")

output_long2$prevalence <- output_long2$value / sum(ini_state)

graf3 <- ggplot(data = output_long2,
               aes(x = time,
                   y = prevalence*100, # %
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Prevalence Proportion [%]") +
  labs(title = "Prevalence of infection, susceptibility and recovery over time",
       colour = "Compartment") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green")) +
  theme(legend.position = "bottom")

# Calculating the effective reproduction number
output2$reff <- parameters["beta"] / parameters["gamma"] *
  output2$S / (output2$S + output2$I + output2$R)

graf4 <- ggplot(data = output2,
                aes(x = time,
                    y = reff)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Reff") +
  labs(title = "Effective Reproduction Number over time")
  # ylim(c(0,2))

graf5 <- ggplot(data = output,
               aes(x = time,
                   y = I)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Number of Infected People") +
  labs(title = "Infection curve over 4 generations")

output_long$prevalence <- output_long$value/sum(ini_state)
graf6 <- ggplot(data = output_long,
                aes(x = time,
                    y = prevalence*100, # %
                    colour = variable,
                    group = variable)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence Proportion [%]") +
  labs(title = "Prevalence proportion over 4 generations",
       colour = "Compartment") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green")) +
  theme(legend.position = "bottom")

# subplotting each graphic with dim (2,2)
print(ggarrange(graf, graf2, graf4, graf5, graf3, graf6, nrow = 3, ncol = 2))