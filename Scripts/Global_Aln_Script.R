#Format metadata for the E gene global run
pacman::p_load(tidyverse, purrr, MoMAColors)

#set working directory
setwd('/DENV/E_gene/')
#Read metadata 
data <- read_tsv('E_aln_D2_global_metadata.tsv') %>% 
  janitor::clean_names() %>% 
  mutate(virus = paste(virus_name, accession_id, collection_date, sep = '|'))

#lineage assignation from Nextclade
lineageFiles <- dir('nextcladeAssignation', pattern = '.tsv', full.names = T)

lineages <- lineageFiles %>%
  map(read_tsv) %>% 
  reduce(rbind) %>% 
  distinct(seqName, .keep_all = T) %>% 
  select(seqName, clade, coverage)

#Join the metadata with the lineage assignation
fullData <- left_join(data, lineages, by = join_by(virus == seqName)) %>% 
  drop_na(clade) %>% 
  filter(coverage >= 0.6) %>% 
  rename(strain = virus)

#Get the countries and group everything except Colombia into continents

finalData <- fullData %>% 
  separate(location, into = c('continent','country'), sep = ' / ') %>% 
  mutate(location = case_when(country == 'Colombia' ~ 'Colombia',
                              .default = continent),
         location = case_when(location == 'South America'~ 'Rest of South America',
                              location == 'North America'~ 'North America & the Caribbean',
                              .default = location)) %>% 
  rename(date = collection_date)

write_tsv(finalData, 'metadata/global_metadata_d2.tsv')


#The lineages that are most abundant in colombia are F.1.1.2, D.2 and D.3
#Create datasets for those and subsample them
F1112 <- fullData %>% 
  separate(location, into = c('continent','country','state','additional_location'), sep = ' / ') %>% 
  filter(clade == '2II_F.1.1.2') %>% 
  rename(date = collection_date)

write_tsv(F1112, 'F112/F112_metadata.tsv')

D2_lineage <- fullData %>% 
  separate(location, into = c('continent','country','state','additional_location'), sep = ' / ') %>% 
  filter(clade == '2III_D.2') %>% 
  rename(date = collection_date)

write_tsv(D2_lineage, 'D2_lineage/D2_lineage.tsv')


#Extraction mutations



#Diversity of the global Dengue 2 alignment
palette <- moma.colors('Ohchi', 25)


lineage_abundance <- finalData %>% 
  group_by(clade) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  mutate(clade = case_when(count <= 20 ~ 'Minority lineages',
                           .default = clade)) %>% 
  summarise(count = sum(count), .by = 'clade') %>% 
  ggplot(aes(x=clade, y = count, fill = clade, color = clade)) +
  geom_bar(stat = 'identity') +
  scale_color_manual(values = palette, aesthetics = c('fill','color'), name = 'Lineage') +
  ylab('Sequence count') +
  xlab('Lineage') +
  theme_bw(base_family = 'Helvetica', base_size = 16)

ggsave('lineage_abundance_plot_SUPP.jpeg', lineage_abundance, width = 27, height = 10.5, units = 'in')
