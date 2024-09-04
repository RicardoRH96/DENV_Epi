#DENV data filtering
library(tidyverse)

metadata <- read_tsv('D2/filtered_d2_metadata.tsv')

table(metadata$country) %>% View()

caribbean <- c('Cuba', 'Guadeloupe','Martinique','St. Barthélemy','Haiti','Puerto Rico')
central_am <- c('Nicaragua','Guatemala','Costa Rica','Honduras','El Salvador')
other_southAm <- c('French Guiana','Bolivia','Paraguay')
north_america <- c('Mexico','United States')

grouped_demes <- metadata %>% 
  mutate(country = case_when(country %in% caribbean ~ 'Caribbean',
                             country %in% central_am ~ 'Central America',
                             country %in% other_southAm ~ 'Rest of South America',
                             country %in% north_america ~ 'North America',
                             .default = country))

table(grouped_demes$genotype) %>% View()


#Check those with assigned lineages
regions_of_interest <- c('Central-America','South-America', 'North America', 'Caribbean')

lineages <- read_tsv('D2/Nate_seqs/metadata.tsv') %>% 
  janitor::clean_names() %>% 
  filter(region %in% regions_of_interest)

range(lineages$date)
