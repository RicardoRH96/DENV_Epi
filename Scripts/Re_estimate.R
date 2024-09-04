library(estimateR); library(tidyverse); library(lubridate); library(fitdistrplus); library(EpiNow2)
#Estimation of effective reproduction number of DENV
#Dengue parameters from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10416499/#Sec41
##Delay distribution for DENV reporting
### mu = 2.9, sigma = 0.83 lognormal
### Intrinsic incubation period lognormal with mean 5.9 d and sd of 1.6 d
### Human to Human via mosquito vector - gamma distribution with mean = 23 d and sd = 8.5 d.
### LOESS smoothing parameter sigma of 10 weeks

######################## Get delay distribution ###################
#There is no delay data for the last 3 months of 2024 (May to July)
setwd("~/Documents/DENV")
case_data <- read_csv('2024/Metadata/210_denv.csv')
range(case_data$INI_SIN)

delay_samples <- case_data %>% 
  filter(FEC_NOT - INI_SIN <= 365) %>% 
  mutate(reporting_delay = FEC_NOT - INI_SIN) %>% 
  mutate(reporting_delay=ifelse(reporting_delay==0,0.1,reporting_delay)) %>% 
  mutate(reporting_delay = case_when(reporting_delay < 0 ~ reporting_delay *-1,
                                     .default = reporting_delay)) %>%  ##Correct for older database where the onset of symptom date and the reporting date seems to be inverted
  filter(reporting_delay <= 100) #Cutoff of 100 days for reporting delay

delay_dists <- lapply(list("2007"=2007,"2008"=2008,"2009"=2009,"2010"=2010,"2011"=2011,"2012"=2012,"2013"=2013,"2014"=2014,
                           "2015"=2015,"2016"=2016,"2017"=2017,"2018"=2018,"2019"=2019,"2020"=2020,"2021"=2021,"2022"=2022,
                           "2023"=2023, '2007-2023'=c(2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023),
                           '2014-2023'=c(2014,2015,2016,2017,2018,2019,2020,2021,2022,2023)), function(year) {
  ds <- delay_samples %>% 
    filter(Year %in% year) %>% 
    pull(reporting_delay) %>%
    as.numeric()
  fitted <- fitdist(ds[sample.int(length(ds), 5000)], "lnorm") # need to subsample, too much data
  return(fitted$estimate)
})

d_delay_dists <- bind_rows(lapply(list("2007"=2007,"2008"=2008,"2009"=2009,"2010"=2010,"2011"=2011,"2012"=2012,"2013"=2013,"2014"=2014,
                                       "2015"=2015,"2016"=2016,"2017"=2017,"2018"=2018,"2019"=2019,"2020"=2020,"2021"=2021,"2022"=2022,
                                       "2023"=2023), function(x) {
  data.frame(reporting_delay = 1:200, p = dlnorm(1:200, meanlog = delay_dists[[as.character(x)]]["meanlog"], sdlog = delay_dists[[as.character(x)]]["sdlog"]))
}), .id = "Year")

write_csv(d_delay_dists, '2024/Metadata/DENV_ReportingDelay.csv')

#Now to the estimation part
distribution_onset_to_confirmation <- list(name = "lnorm", meanlog = delay_dists[['2023']]['meanlog'], sdlog = delay_dists[["2023"]]["sdlog"])
#Incubation period
mean_incubation = 5.9
sd_incubation = 1.6
distribution_incubation <- list(name = "lnorm", meanlog = convert_to_logmean(mean_incubation, sd_incubation), sdlog = convert_to_logsd(mean_incubation, sd_incubation))


#Generation Interval
mean_serial_interval <- 23
std_serial_interval <- 8.5
generation_time <- list(mean = mean_serial_interval, sd = std_serial_interval, dist="gamma", max = 15, fixed = T)

dengue_data <- case_data %>% 
  group_by(FEC_NOT) %>% 
  summarise(cases = n()) %>% 
  dplyr::rename(date = FEC_NOT) %>% 
  filter(date >= '2014-01-01') 

#2024 data 
denv_2024 <- read_csv('/Users/ricardorivero/Documents/DENV/2024/Metadata/2024_grouped.csv') %>% 
  mutate(month_date = floor_date(FEC_NOT, unit='months')) %>% 
  group_by(month_date) %>% 
  summarise(cases = sum(counts)) %>% 
  dplyr::rename(date = month_date)



#Fit estimateR
Re_estimates <- try(get_block_bootstrapped_estimate(
  incidence_data = dengue_data %>% mutate(confirm = ifelse(cases==0,0.1,cases)) %>% pull(cases),
  N_bootstrap_replicates = 100,
  smoothing_method = "LOESS",
  data_points_incl = 7*24,
  initial_Re_estimate_window = 7*8,
  delay = list(distribution_incubation,  distribution_onset_to_confirmation),
  estimation_window = 14,
  minimum_cumul_incidence = 200,
  mean_serial_interval = mean_serial_interval,
  std_serial_interval = std_serial_interval,
  ref_date = min(dengue_data$date),
  time_step = "day",
  combine_bootstrap_and_estimation_uncertainties = TRUE,
  output_Re_only = FALSE
), silent = TRUE)


# #Monthly fit
# Re_estimates_2024 <- get_block_bootstrapped_estimate(
#   incidence_data = denv_2024_adjusted %>% mutate(confirm = ifelse(cases==0,0.1,cases)) %>% pull(cases),
#   N_bootstrap_replicates = 100,
#   smoothing_method = "LOESS",
#   data_points_incl = 28,
#   initial_Re_estimate_window = 1,
#   delay = list(distribution_incubation,  distribution_onset_to_confirmation),
#   estimation_window = 1,
#   minimum_cumul_incidence = 200,
#   mean_serial_interval = mean_serial_interval,
#   std_serial_interval = std_serial_interval,
#   ref_date = min(denv_2024_adjusted$date),
#   time_step = "month",
#   combine_bootstrap_and_estimation_uncertainties = TRUE,
#   output_Re_only = FALSE
# )
# 
# 
# estimates_re_2024 <- Re_estimates_2024 %>% 
#   filter(date >= '2024-01-22' & date <= '2024-06-22')
# 
# all_estimates <- rbind(Re_estimates, estimates_re_2024)

#Dataframe wrangling
dengue_R_all <- Re_estimates %>% 
  dplyr::select(date, R = Re_estimate, lower = CI_down_Re_estimate, upper=CI_up_Re_estimate) %>% 
  mutate(Software = "EpiNow2")
  
#Plot
colour_palette <- viridis::viridis(7)
Re_plot <- dengue_R_all %>% 
ggplot(aes(x=date)) +
  geom_hline(yintercept=1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Software), alpha = 0.25) +
  geom_line(aes(y=R, color = Software), lwd = 1.1) +
  scale_colour_manual(values = colour_palette[c(1,3)], aesthetics = c("colour", "fill"),
                      guide='none') +
  scale_x_date(date_breaks = "1 years",
               date_labels = "%Y",
               limits = c(as.Date("2014-01-01"), as.Date("2024-01-01")),
               expand = c(0,0)) +
  coord_cartesian(ylim = c(0,2)) +
  theme_classic() +
  xlab("") +
  ylab("Effective reproductive number") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = c(.9, .9),
        legend.background = element_blank())
plot(Re_plot)


ggsave('2024/Re_estimate.jpeg', Re_plot, dpi = 1200)
