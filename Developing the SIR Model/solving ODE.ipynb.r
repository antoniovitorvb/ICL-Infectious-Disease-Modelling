
# If you are working in RStudio, you will need to install if you haven't already:
# install.packages("ggplot2")
# install.packages("deSolve")

# Run the library() function to load these packages
library(ggplot2)
library(deSolve)


require(deSolve)

result <- ode(      y = state_vars  # initial state variables: state variables at first timepoint
              , times = times       # a vector of timepoints 
              ,  func = exp_pop_fn  # the differential equation itself, in the form of a function
              , parms = parms)      # parameters for the equation in func

# Now to see what each input actually is.

?ode

# y
state_vars  <- c(N = 1) # a vector with a named element

# times
times       <- seq(0, 40, by = 0.5) # the shorter the timesteps, the more accurate a solution

# func: a function 
exp_pop_fn <- function(time, state, parameters) { 
                                      
  N <- state['N'] # within this function, the state variables - in this case, just N
  
  dN <- parameters['alpha'] * N  # alpha is extracted from the 'parameters' argument
                                 # this argument gets fed in by the ode() function
                                 # specify when running ode()
  
  return(list(c(dN))) # tell the function what to return: 
                      # here, we need the calculated dN value, in a list.
  
  # if you have more than one state variable to calculate
  # tell the function to return derivatives in the same order
  # as you entered them (in the `state` argument)
                      
} 

# Remember that this function is an argument into another function; it doesn't do a lot by itself.
# The inputs come from running the ode() function, the output goes back into the ode() function.

# parms
parms <- c(alpha = log(2)) # alpha has been chosen so the population doubles each timestep

parms['alpha'] # you can see "parms" is a vector with one named element, alpha.
               # this argument 'parms' gets fed, by ode(), into the function that you specify to use as func
               # so it needs to contain whatever variables that function is expecting.

# For this example:
result <- ode(y = state_vars         # contains initial N
              , times = times        # the vector of timepoints 
              , func = exp_pop_fn    # the exponential equation, written as a function
              , parms = parms)       # parameters for the exponential equation: here, just alpha


head(as.data.frame(result)) # shows the first few lines of the results

result <- as.data.frame(result) # turn the output into a data.frame

# use ggplot to create a graph of times against the population size N
require(ggplot2) # if

expplot <- ggplot(data = result)
expplot <- expplot + geom_line(aes(x = time, y = N)
                               ,  colour = "blue")
expplot <- expplot + labs(x = "time (days)")
expplot # shows the graph

logistic_fn <- function(t, state, parameters) { # You'll need a K in the parameters argument

        N  <- state['N']  # still our only state variable
    
        dN <- parameters['alpha'] * N * (1 - (N / parameters['K'])) 
        # this line represents the differential equation

        return(list(c(dN)))

        }


parms['K'] <- 1000000 
# the vector 'parms' now includes an element named K, assigned the value 1000000 to represent carrying capacity

result_K <- ode(    y = state_vars   # still only contains initial N
              , times = times        # the vector of timepoints we already have
              ,  func = logistic_fn  # the logistic equation function
              , parms = parms)       # parms now includes K

result_K <- as.data.frame(result_K)

#check the output, and plot
tail(result_K) # to view the last 6 rows; note that N approaches K

logplot <- ggplot(data = result_K)
logplot <- logplot + geom_line(aes(x = time, y = N)
                               ,  colour = "blue")
logplot <- logplot + labs(x = "time (days)")
logplot


with(as.list(c(state_vars, parms)), {   # give with() a list made of state_vars and parms
    print(alpha)                        # to find anything referenced within the curly brackets,
    print(N)                            # R looks for names within the object given to with()
})  

# clearer code for the logistic function:

logistic_fn <- function(t, state, parameters) { # You'll need a K in the parameters argument

    with(as.list(c(state, parameters)), {           
        
        # you can now just refer to N directly, so no need to assign
        
        dN <- alpha * N * (1 - (N / K))  # this line represents the differential equation
        
        return(list(c(dN)))
    })

        
    }



exp(parms['alpha']) 
#  alpha was chosen so that at each whole timestep, the population doubles

t <- seq(1, 40, by = 1) # vector of times at which to calculate the population size
# these don't have to be the same as the timepoints as the ode() output was generated at

N_calc <- state_vars['N'] * 2^t # every day, the population doubles

# R automatically vectorises this expression, applying it to each element of 't' in turn, to create vector N_calc
# N_calc should be the same length as t. Make them into a dataframe
pop_df <- data.frame(times = t, N = N_calc)

require(ggplot2)
expplot  <- expplot + 
            geom_point(data = pop_df,    # specify different dataframe here
                        aes(y = N_calc, x = t)
                      , shape = 4
                      , size = 2
                      , colour = "red")  

expplot
