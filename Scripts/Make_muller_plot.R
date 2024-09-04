#Load libraries
library(tidyverse);library(ggstream); library(colorspace); library(ggtext); library(cowplot); library(fedmatch)

theme_set(theme_classic(base_family = 'Helvetica',base_size = 14))


library(dplyr)
library(readr)

# # Load lineage data and select relevant columns
# lineages <- read_csv('~/Documents/DENV/2024/Metadata/D2_lineages.csv') %>%
#   select(name, assignment)
# 
# # Clean the 'name' column by removing 'Colombia' and extracting the first 23 characters
# lineages$name <- substr(gsub('Colombia', '', lineages$name), 1, 23)
# 
# # Load metadata for Dengue2
# metadata <- read_tsv("~/Documents/DENV/2024/D2/all_colombia/subsampled_d2_metadata.csv") %>%
#   filter(country == 'Colombia')
# 
# 
# 
# #load merged data
# merged <- inner_join(lineages, metadata, by = 'name')
# 
# merged %>% 
#   group_by(date, assignment) %>% 
#   summarise(number = n()) %>% 
#   View()


#Muller plot, genotype level (note: update when lineage assignation is available)

data <- read_tsv('~/Documents/DENV/2024/Metadata/colombia_serotype_genotype_08_10_2024.tsv') %>% 
  janitor::clean_names() %>% 
  filter(collection_date != "unknown") %>% 
  filter(nchar(collection_date) > 4) %>%
  filter(grepl("-", collection_date)) %>% 
  mutate(collection_date = case_when(nchar(collection_date)==7 ~ paste0(collection_date, '-01'),
                                     .default = collection_date),
         collection_date = as.Date(collection_date))


#Summarise the data for muller plot
muller_data <- data %>%
  filter(genotype != 'unassigned') %>%
  mutate(genotype_serotype = paste(serotype, genotype, sep = '-'),
         epi_week = floor_date(collection_date, "month")) %>%
  filter(epi_week >= '2014-01-01') %>%
  group_by(epi_week) %>%
  mutate(total_count = n()) %>%  # Calculate total count per epi_week
  group_by(genotype_serotype, epi_week) %>%
  summarise(count = n(),
            percentage = (count / first(total_count))) %>%  # Calculate percentage
  ungroup()
  
#palette
pal <- c(
  "#00496f", lighten("#00496f", .25, space = "HLS"),
  "#0f85a0", lighten("#0f85a0", .3, space = "HLS"),
  "#edd746", lighten("#edd746", .25, space = "HLS"),
  "#ed8b00", lighten("#ed8b00", .2, space = "HLS"),
  "#dd4124", lighten("#dd4124", .15, space = "HLS")
)

serotype_palette <- c('DENV1-V' = "#00496f",
                      'DENV2-II-Cosmopolitan' = pal[2],
                      'DENV2-III-Asian-American'="#edd746",
                      'DENV3-III'=pal[4],
                      'DENV4-II'='#dd4124')

muller_plot <- muller_data %>% 
  ggplot(aes(epi_week, percentage, color = genotype_serotype, fill = genotype_serotype)) +
  geom_stream(
    geom = "polygon",
    color = "white",
    size = 0.8,
    type = 'proportional',
    bw = .5 # Controls smoothness
  ) +
  scale_color_manual(
    expand = c(0, 0),
    values = serotype_palette,
    guide = "none"
  ) +
  scale_fill_manual(
    values = serotype_palette,
    name = NULL
  ) +
  ylab('Genotype share') +
  xlab('Collection date') +
  theme(
    legend.position = 'top',
    plot.margin = unit(c(0, 0, 0, 0), "cm") # Remove margins
  )


ggsave('~/Documents/DENV/2024/genotype_muller_plot.jpeg', dpi = 1200, height = 4.5, width = 12, units = 'in', device = 'jpeg')
