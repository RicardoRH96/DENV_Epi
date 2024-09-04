library(tidyverse); library(janitor); library(countrycode)

#Working directory
setwd('~/Documents/DENV/2024/')

metadata <- read_tsv('Metadata/denv_metadata_07_12_24.tsv') |> 
  janitor::clean_names() |> 
  select(-sampling_strategy, -vaccination_history, -last_vaccination_date,-passage,-specimen,-additional_host_information) |> 
  separate(location, into = c('continent','country','state','city','additional_loc'), sep=' / ') |> 
  mutate(country = countrycode::countryname(country, destination = 'country.name')) |> 
  mutate(un_subregion = countrycode::countryname(country, destination = 'un.regionsub.name'),
strain = paste(virus_name, accession_id, collection_date, sep = '/')) |> 
  rename(
  date = collection_date)


glimpse(metadata)

#structure of the header is: >hDenV1/Brazil/RJ-IAL-12/2023|EPI_ISL_17699800|2023-03-20 (>hDenV1/country/), some sequences have white-spaces in the header, remember to remove it.
#DENV1
metadata |> filter(serotype == 'DENV1') |> 
  distinct(strain, .keep_all = TRUE) |> 
  select(strain, date, country) |> 
  write_csv('D1/filtered_d1_metadata.csv')

#DENV2
metadata |> filter(serotype == 'DENV2') |> 
  distinct(strain, .keep_all = TRUE) |> 
  select(strain, date, country, genotype) |> 
  write_tsv('D2/filtered_d2_metadata.tsv')

#DENV3
metadata |> filter(serotype == 'DENV3') |> 
  distinct(strain, .keep_all = TRUE) |> 
  select(strain, date, country, genotype) |> 
  write_tsv('D3/filtered_d3_metadata.tsv')
#DENV4
metadata |> filter(serotype == 'DENV4') |> 
  distinct(strain, .keep_all = TRUE) |> 
  select(strain, date, country) |> 
  write_tsv('D4/filtered_d4_metadata.tsv')

