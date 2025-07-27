rm(list = ls())

library(dplyr)
library(readr)
library(MVSE)

# Load the climate variables
climate_vars <- read_csv('./HT_Departments.csv') %>% 
  mutate(month_date = floor_date(date, unit='month')) %>% 
  group_by(Departamento, month_date) %>% 
  summarise(meanT = mean(meanT),
            meanH = mean(meanH)) %>% 
  rename(H = meanH,
         T = meanT,
         date= month_date)

# Initialize an empty list to store results
model_results <- list()

# Loop through each state and run the model
for(state in unique(climate_vars$Departamento)) {
  
  # Filter the data for the current state
  state_data <- climate_vars %>% filter(Departamento == state)
  
  # Run the model for the current state
  denv_model <- mvse_model(model_name="DENV_COL", 
                           climate_data=state_data, 
                           model_category="denv_aegypti", 
                           warning=FALSE)
  
  # Fit the model
  denv_fit <- MVSE::sampling(denv_model, iter=10^6, warmup=0.3, 
                             verbose=FALSE, samples=10^3, seed=1)
  
  # Store the fit results along with the state name
  model_results[[state]] <- list(state = state, fit = denv_fit)
}

# Combine all results into a final list or dataframe
final_results <- do.call(rbind, lapply(model_results, function(x) data.frame(State = x$state, Fit = list(x$fit))))

# Now, `final_results` contains the fit outputs for each state with the state name as a variable

mcmc <- coda::as.mcmc(model_results$AMAZONAS$fit@sim$indexP[,2:1000])

coda::HPDinterval(mcmc, prob = 0.95)

s <- summary(mcmc)

m <- s$statistics[,'Mean']


#Test the loop
# Load necessary libraries
library(dplyr)
library(readr)
library(MVSE)
library(coda)

# Load the climate variables
climate_vars <- read_csv('./HT_Departments.csv') %>% 
  rename(H = meanH,
         T = meanT)

# Initialize an empty list to store results
model_results <- list()

# Initialize an empty dataframe to store the final results
final_results <- data.frame(State = character(),
                            Date = as.Date(character()),
                            Mean = numeric(),
                            HPD95_low = numeric(),
                            HPD95_hi = numeric(),
                            stringsAsFactors = FALSE)

# Loop through each state and run the model
for(state in unique(climate_vars$Departamento)) {
  
  # Filter the data for the current state
  state_data <- climate_vars %>% filter(Departamento == state)
  
  # Run the model for the current state
  denv_model <- mvse_model(model_name="DENV_COL", 
                           climate_data=state_data, 
                           model_category="denv_aegypti", 
                           warning=FALSE)
  
  # Fit the model
  denv_fit <- MVSE::sampling(denv_model, iter=10^6, warmup=0.3, 
                             verbose=FALSE, samples=10^3, seed=1)
  
  # Extract the MCMC samples for the second to 1000th parameter (modify as needed)
  dates <- denv_fit@sim$indexP[,1]
  mcmc <- coda::as.mcmc(denv_fit@sim$indexP[,2:1000])
  
  # Calculate the HPD intervals and mean
  hpd <- coda::HPDinterval(mcmc, prob = 0.95)
  s <- summary(mcmc)
  m <- s$statistics[,'Mean']
  
  # Add the results to the final dataframe
  temp_df <- data.frame(State = state,
                        Date = dates,  # Assuming you want to add the current date
                        Mean = m,
                        HPD95_low = hpd[,1],
                        HPD95_hi = hpd[,2])
  
  final_results <- rbind(final_results, temp_df)
}

# Now, `final_results` contains the mean, HPD intervals, and state names for all the states



#Test 2
# Load the climate variables
climate_vars <- read_csv('/DENV/2024/Climate_vars/HT_Departments.csv') %>% 
  rename(H = meanH,
         T = meanT)

# Initialize an empty list to store results
model_results <- list()

# Initialize an empty dataframe to store the final results
final_results <- data.frame(State = character(),
                            Date = as.Date(character()),
                            Mean = numeric(),
                            stringsAsFactors = FALSE)

# Loop through each state and run the model
for(state in unique(climate_vars$Departamento)) {
  
  # Filter the data for the current state
  state_data <- climate_vars %>% filter(Departamento == state)
  
  # Run the model for the current state
  denv_model <- mvse_model(model_name="DENV_COL", 
                           climate_data=state_data, 
                           model_category="denv_aegypti", 
                           warning=FALSE)
  
  # Fit the model
  denv_fit <- MVSE::sampling(denv_model, iter=10^6, warmup=0.3, 
                             verbose=FALSE, samples=10^3, seed=1)
  
  # Extract the MCMC samples for the second to 1000th parameter (modify as needed)
  dates <- denv_fit@sim$indexP[,1]
  mcmc <- coda::as.mcmc(denv_fit@sim$indexP[,2:ncol(denv_fit@sim$indexP)])
  
  # Calculate the HPD intervals and mean
  mean_P <- rowMeans(mcmc)
  
  # Add the results to the final dataframe
  temp_df <- data.frame(State = rep(state, length(dates)),
                        Date = dates,
                        Mean = mean_P)
  
  final_results <- rbind(final_results, temp_df)
}

write_csv(final_results, './P_index.csv')

# Now, `final_results` contains the mean, HPD intervals, and state names for all the states
