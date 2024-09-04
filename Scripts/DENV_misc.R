#Dengue lineages merge metadata
#Load libraries
library(tidyverse);library(ggstream); library(colorspace); library(ggtext); library(cowplot); library(fedmatch)

lineages <- read_csv('~/Documents/DENV/2024/Metadata/DENV2_Lineages.csv') %>% 
  select(name, assignment)

colnames(lineages)[1] <- 'accession_id'

theme_set(theme_classic(base_family = 'Helvetica',base_size = 14))


dataset <- read_tsv('~/Documents/DENV/2024/Metadata/denv_metadata_07_12_24.tsv') %>% 
  janitor::clean_names()


merged_data <- dataset %>%
  inner_join(lineages, by = 'accession_id') %>%
  select(-sampling_strategy, -vaccination_history, -last_vaccination_date, -passage, -specimen, -additional_host_information) %>%
  separate(location, into = c('continent', 'country', 'state', 'city', 'additional_loc'), sep = ' / ') %>%
  mutate(country = countrycode(country, origin = 'country.name', destination = 'country.name')) %>%
  filter(country == 'Colombia') %>%
  filter(nchar(collection_date) > 4) %>%
  filter(grepl("-", collection_date)) %>%
  mutate(collection_date = case_when(
    nchar(collection_date) == 7 ~ as.Date(paste0(collection_date, '-01')),
    TRUE ~ as.Date(collection_date)
  ))

# Summarise the data for muller plot
muller_data <- merged_data %>%
  mutate(epi_week = floor_date(collection_date, "month")) %>%
  filter(epi_week >= as.Date('2013-01-01')) %>%
  group_by(epi_week) %>%
  mutate(total_count = n()) %>%
  group_by(assignment, epi_week) %>%
  summarise(count = n(),
            percentage = (count / dplyr::first(total_count))) %>%
  ungroup()


muller_data %>% 
  ggplot(aes(epi_week, percentage, color = assignment, fill = assignment)) +
  geom_stream(
    geom = "polygon",
    color = "white",
    size = 0.8,
    type = 'proportional',
    bw = .5 # Controls smoothness
  ) +
  scale_color_manual(
    expand = c(0, 0),
    values = pnw_palette("Bay", 5),
    guide = "none"
  ) +
  scale_fill_manual(
    values = pnw_palette("Bay", 5),
    name = NULL
  ) +
  ylab('Genotype share') +
  xlab('Collection date') +
  theme(
    legend.position = 'top',
    plot.margin = unit(c(0, 0, 0, 0), "cm") # Remove margins
  )


#Georeferencing of the samples
install.packages('tidygeocoder')
library(tidygeocoder)

for_geocoding <- merged_data %>% 
  mutate(Yale_id = str_extract(virus_name, "Yale-[A-Za-z]{2}\\d{3,5}")) %>% 
  full_join(., filtered_yale, by='Yale_id') %>% 
  mutate(state.x = case_when(Yale_id != "NA" ~ state.y,
                             .default = state.x),
         city.x = case_when(Yale_id != "NA" ~ city.y,
                            .default = city.x))

write_csv(for_geocoding, '~/Documents/DENV/2024/Metadata/needs_geocoding.csv')

#Read and geocode
needs_geocoding <- read_csv('~/Documents/DENV/2024/Metadata/needs_geocoding.csv') %>% 
  dplyr::select(-state.y, -city.y) %>% 
  dplyr::rename(state = state.x,
                city = city.x,
                country = country.x)

geocoded_samples <- needs_geocoding %>% 
  mutate(addr = paste(country, state,city, sep = ',')) %>% 
  geocode(addr, method = 'osm', lat = latitude , long = longitude)

geocoded_samples$state <- ifelse(grepl("LET", geocoded_samples$virus_name), "Amazonas", geocoded_samples$state)

geocoded_samples <- geocoded_samples %>% 
  mutate(addr = paste(country, state,city, sep = ','),
         latitude = case_when(
                              addr == 'Colombia,NA,NA' ~ 4.59808,
                              state == 'Amazonas' ~ -4.2081,
                              .default = latitude),
         longitude = case_when(
                               addr == 'Colombia,NA,NA' ~ -74.076044,
                               state == 'Amazonas' ~ -69.9432,
                               .default = longitude))

#write geocoded metadata
write_csv(geocoded_samples, '~/Documents/DENV/2024/Metadata/D2_col_geocoded.csv')


#Multitype data

Metadata_msbd <- geocoded_samples %>% 
  select(accession_id, country,assignment,collection_date)

write_tsv(Metadata_msbd, '~/Documents/DENV/2024/D2/all_colombia/MSBD/D2_MSBD_METADATA.tsv')


#Yale samples for location metadata Verity
temp <- geocoded_samples %>% 
  filter(grepl('Yale', virus_name)) %>% 
  select(virus_name, accession_id, collection_date) %>% 
  mutate(Yale_id = str_extract(virus_name, "Yale-[A-Za-z]{2}\\d{3,5}"))

#Nate's metadata
Yale_metadata <- readxl::read_xlsx('~/Documents/DENV/2024/Metadata/Yale_samples_metadata.xlsx') %>% 
  janitor::clean_names() %>% 
  mutate(sample_id = gsub('.1','',sample_id))

ids <- temp$Yale_id

filtered_yale <- Yale_metadata %>% 
  filter(sample_id %in% ids) %>% 
  dplyr::rename(Yale_id = sample_id)

#Merge with other metadata
geocoded_v2 <- full_join(geocoded_samples, filtered_yale)


#Extract subdatasets for continuous phylogeography of lineages
D.2 <- geocoded_samples %>% 
  filter(assignment == '2III_D.2') %>% 
  dplyr::select(virus_name, accession_id, collection_date.x, country, state, city, latitude, longitude, assignment)

F.1.2 <- geocoded_samples %>% 
  filter(assignment == '2II_F.1.2') %>% 
  dplyr::select(virus_name, accession_id, collection_date.x, country, state, city, latitude, longitude, assignment)

F.1.2_continuous_phylo <- F.1.2 %>% 
  mutate(traits = paste(virus_name,accession_id,collection_date.x, sep = "/")) %>% 
  dplyr::select(traits, latitude, longitude) %>% 
  dplyr::rename(
                lat = latitude,
                long = longitude)

write_tsv(F.1.2_continuous_phylo, '2024/D2/all_colombia/Continuous_phylo/lineage_F12_latlong.tsv')


#Lambo stuff
dataset %>% filter(collection_date > as.Date('2010-01-01'), serotype == 'DENV2', country == 'Colombia')
