# INFECTIOUS dISEASE MODELLING
# Module 1 e-tivity 10 - Combining the Mechanisms

library(deSolve)  # solves the model
library(reshape2) # changes hte shape of the model output
library(ggplot2)  # plots
# library(ggpubr) # subplots ggplot2

population <- 3*10^5
parameters <- c(beta = 1*365, # Infection rate
                gamma = 1/20*365, # Infection duration = 20 days
                mu = 1/3, # 3 year lifespan
                b = 1/3, # Population don't vary over time
                sigma = 0, # Life-long immunity
                p_vacc = 0) # No neonatal vaccination so far

ini_state <- c(S = (1-0.05)*population,
               I = 0.05 * population,
               R = 0)

times <- seq(from = 0, to = 6, by = 1/365)
sir_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)),{
    N <- S + I + R
    
    lambda <- beta * I/N
    
    dS <- -lambda*S - mu*S + (1-p_vacc)*b*N +sigma*R
    dI <- lambda*S - gamma*I - mu*I
    dR <- gamma*I - mu*R + p_vacc*b*N - sigma*R
    
    return(list(c(dS,dI,dR)))
  })
}

output <- as.data.frame(ode(y = ini_state,
                            times = times,
                            func = sir_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")
output_long$prevalence <- output_long$value / sum(ini_state)

graf <- ggplot(data = output_long,
               aes(x = time,
                   y = prevalence*100,
                   colour = variable,
                   group = variable)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence [%]") +
  labs(title = "Prevalence of S, I and R  animals over time") +
  scale_color_manual(labels = c("Susceptible", "Infected", "Recovered"), values = c("blue", "red", "green"))

print(graf)

#1 What is the endemic prevalence of the disease currently (the baseline prevalence), assuming permanent immunity?
  # Prevalence seems to have stabilized at 2 years, so the prevalence of the desease relates to the Infected cases
q1 <- output_long$prevalence[output_long$time ==2 & output_long$variable == "I"]

#2 What proportion of newborn animals would you need to vaccinate to reduce the prevalence by half, assuming life-long immunity?
  # Introducing neonatal vaccine after endemic equilibrium (time == 2)
new_state <- c(S = output$S[output$time == 2],
               I = output$I[output$time == 2],
               R = output$R[output$time == 2])
parameters["p_vacc"] <- 0.5 # 50% of newborns

output2 <- as.data.frame(ode(y = new_state,
                             times = times,
                             func = sir_model,
                             parms = parameters))

output2_long <- melt(as.data.frame(output2), id = "time")
output2_long$prevalence <- output2_long$value / sum(new_state)

graf2 <- ggplot(data = output2,
               aes(x = time,
                   y = I/sum(new_state)*100)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence [%]") +
  labs(title = paste("infecction prevalence when\n
       neonatal vaccine coverage =", parameters["p_vacc"]*100,"%"))
print(graf2)

q2 <- output2_long$prevalence[output2_long == tail(output2_long$time,1) & output2_long$variable == "I"]

#q3 Would it be possible to eliminate the disease from the population using neonatal vaccination under the assumption of lifelong immunity?

R0 <- parameters["beta"] / parameters["gamma"]
parameters["p_vacc"] <- 1 - 1/R0

output3 <- as.data.frame(ode(y = new_state,
                             times = times,
                             func = sir_model,
                             parms = parameters))

output3_long <- melt(as.data.frame(output3), id = "time")
output3_long$prevalence <- output3_long$value / sum(new_state)

graf3 <- ggplot(data = output3,
                aes(x = time,
                    y = I/sum(new_state)*100)) +
  geom_line() +
  xlab("Time (years)") +
  ylab("Prevalence [%]") +
  labs(title = paste("infecction prevalence when\n
       neonatal vaccine coverage =", parameters["p_vacc"]*100,"%"))
print(graf3)

q3 <- output3_long$prevalence[output3_long == tail(output3_long$time,1) & output3_long$variable == "I"]

#q4 If the average duration of immunity is only 1 year, how would this impact the proportional reduction in the prevalence with the vaccine coverage you obtained above compared to the baseline?

#q5 Would it be possible to eliminate the disease from the population using neonatal vaccination under these assumptions?

#q6 If an adjuvant (a vaccine promoter) was given along with the vaccine, that would extend the duration of immunity to 2.5 years on average, what vaccine coverage would be needed to reduce
#   the baseline prevalence by half? Would it be possible to eliminate the disease from the population under these assumptions using neonatal vaccination?

#q7 Based on your results, what overall recommendation would you give to the Minister?

#q8 Also provide some information to help the Minister interpret these results. Write down the assumptions in your modelling approach that you think might affect your results.
#   Are there any adaptations you could make to the model structure that would make it more realistic or that would allow you to answer more detailed questions?