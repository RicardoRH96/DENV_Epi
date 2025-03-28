#Harmonize Colombia seqs data for the circular tree
library(tidyverse); library(readxl); library(janitor)

data <- read_tsv('~/Documents/DENV/2024/CircularPhylogeny/All_ColSeqs.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(lett_dep = str_split_i(str_split_i(virus_name, '/', 3), '-', 1),
         location_dept = str_split_i(location, ' / ', 3),
         location_dept = case_when(is.na(location_dept) ~ 'Colombia',
                                   location_dept == 'Valle' ~ 'Valle del Cauca',
                                   location_dept == 'Valle de Cauca' ~ 'Valle del Cauca',
                                   .default = location_dept)) %>% 
  select(virus_name, genotype, accession_id, collection_date, location, location_dept) %>% 
  mutate(yale_id = str_extract(virus_name, "Yale-[A-Za-z]{2}\\d{3,5}"))



yale_metadata1 <- readxl::read_xlsx('~/Documents/DENV/2024/Metadata/Yale_samples_metadata.xlsx') %>% 
  clean_names() %>% 
  rename(virus_name = sample_id) %>% 
  select(virus_name, collection_date, country, city)
yale_metadata2 <- read_csv('~/Documents/DENV/2024/D2/Nate_seqs/2024.09.23_DENV-2_Colombia_metadata.csv') %>% 
  janitor::clean_names() %>% 
  drop_na(city) %>% 
  rename(virus_name = yale_id) %>% 
  select(virus_name, collection_date, country, city)


yale <- rbind(yale_metadata1, yale_metadata2) %>% 
  mutate(virus_name = str_split_i(virus_name, '\\.', 1))
  
#Geocode the data
library(tidygeocoder)

geocoded <- data %>%
  mutate(city = sapply(str_split(location, " / "), function(x) x[4])) %>%
  left_join(yale %>% select(virus_name, city), by = c("yale_id" = "virus_name")) %>%
  mutate(city = coalesce(city.x, city.y),
         country = "Colombia") %>%
  select(-city.x, -city.y) %>%
  rename(state = location_dept) %>% 
  mutate(addr = paste(country, state,city, sep = ',')) %>% 
  geocode(addr, method = 'osm', lat = latitude , long = longitude)

geocoded <- geocoded %>% 
  mutate(name = paste(virus_name, accession_id, collection_date, sep = '|'))

write_tsv(geocoded,'~/Documents/DENV/2024/CircularPhylogeny/All_ColSeqs_locs_geocoded.tsv')

  '#8dd3c7'
