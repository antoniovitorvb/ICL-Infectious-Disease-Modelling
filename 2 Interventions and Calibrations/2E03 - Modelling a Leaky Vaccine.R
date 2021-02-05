# INFECTIOUS dISEASE MODELLING
# Module 2 e-tivity 3 - Modelling a Leaky Vaccine

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
library(ggpubr) # subplots ggplot2

pop <- 10^6

parameters <- c(beta = 1/4, # infection rate
                gamma = 1/10, # natural recovery rate
                Cs = 0.3, # reduction in the force of infection
                pc = 0.6) # Vaccination coverage

ini_state <- c(S = unname((1-parameters["pc"])) * (pop-1),
               I = 1,
               R = 0,
               V = unname(parameters["pc"]) * (pop-1))

times <- 0:730 # 2 years

vacc_model <- function (time, state, parameters) {
  with(as.list(c(state, parameters)),{
    N <- S + I + R + V
    
    lambda <- beta * I/N
    
    dS <- -lambda*S
    dI <- lambda*S + Cs*lambda*V -gamma*I
    dR <- gamma*I
    dV <- -Cs*lambda*V
    
    return(list(c(dS, dI, dR, dV)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = vacc_model,
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
  labs(title = paste("Prevalence of S, I and R over time with\n",
                     "leaky vaccine with coverage =", parameters["pc"]*100, "%")) +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered", "Vaccinated"),
                     values = c("blue", "red", "green", "black"))

# now changing pc to 86 % of coverage
parameters["pc"] = 0.86

new_state <- c(S = unname((1-parameters["pc"])) * (pop-1),
               I = 1,
               R = 0,
               V = unname(parameters["pc"]) * (pop-1))

output2 <- as.data.frame(ode(y = new_state,
                            times = times,
                            func = vacc_model,
                            parms = parameters))
output2_long <- melt(as.data.frame(output2), id = "time")
output2_long$prevalence <- output2_long$value/sum(new_state)

graf2 <- ggplot(data = output2_long,
               aes(x = time,
                   y = prevalence*100,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (days)") +
  ylab("Prevalence [%]") +
  labs(title = paste("Prevalence of S, I and R over time with\n",
                     "leaky vaccine with coverage =", parameters["pc"]*100, "%")) +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered", "Vaccinated"),
                     values = c("blue", "red", "green", "black"))

print(ggarrange(graf, graf2, nrow = 2))