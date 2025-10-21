

# # Load lineage data and select relevant columns
# lineages <- read_csv('DENV/2024/Metadata/D2_lineages.csv') %>%
#   select(name, assignment)
# 
# # Clean the 'name' column by removing 'Colombia' and extracting the first 23 characters
# lineages$name <- substr(gsub('Colombia', '', lineages$name), 1, 23)
# 
# # Load metadata for Dengue2
# metadata <- read_tsv("DENV/2024/D2/all_colombia/subsampled_d2_metadata.csv") %>%
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

library(tidyverse); library(ggstream); library(colorspace); library(ggtext); library(cowplot); library(fedmatch)
theme_set(theme_classic(base_family = 'Helvetica', base_size = 14))
library(dplyr); library(readr); library(janitor); library(lubridate); library(stats); library(scales)

setwd("~/Documents/DENV/")

data <- read_tsv('./2024/Metadata/colombia_serotype_genotype_08_10_2024.tsv') %>%
  janitor::clean_names() %>%
  filter(collection_date != "unknown", nchar(collection_date) > 4, grepl("-", collection_date)) %>%
  mutate(collection_date = if_else(nchar(collection_date) == 7, paste0(collection_date, "-01"), collection_date),
         collection_date = as.Date(collection_date))

muller_data <- data %>%
  filter(genotype != "unassigned") %>%
  mutate(genotype_serotype = paste(serotype, genotype, sep = "_"),
         epi_week = floor_date(collection_date, "month")) %>%
  filter(epi_week >= as.Date("2014-01-01")) %>%
  group_by(epi_week) %>%
  mutate(total_count = n()) %>%
  group_by(genotype_serotype, epi_week) %>%
  summarise(count = n(), percentage = count / first(total_count), .groups = "drop")

# dynamic palette covering all strata; enforce your preferred hues
all_keys <- sort(unique(muller_data$genotype_serotype))
auto_cols <- setNames(qualitative_hcl(length(all_keys), palette = "Dark 3"), all_keys)
preferred <- c(
  "DENV1_V"  = "#00496f",
  "DENV2_II" = lighten("#00496f", .25, space = "HLS"),
  "DENV2_III"= "#edd746",
  "DENV3_III"= lighten("#ed8b00", .2, space = "HLS"),
  "DENV4_II" = "#dd4124"
)
serotype_palette <- replace(auto_cols, names(preferred), preferred)

muller_plot <- ggplot(muller_data, aes(epi_week, percentage, color = genotype_serotype, fill = genotype_serotype)) +
  geom_stream(geom = "polygon", color = "white", size = 0.8, type = "proportional", bw = .5) +
  scale_color_manual(values = serotype_palette, guide = "none") +
  scale_fill_manual(values = serotype_palette, name = NULL) +
  ylab("Genotype share") + xlab("Collection date") +
  theme(legend.position = "top", plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot(muller_plot)


key <- "DENV3_III"
bw_months <- 10
win_start <- as.Date("2022-01-01")
win_end   <- as.Date("2024-12-31")

# complete monthly grid of counts
all_keys <- sort(unique(muller_data$genotype_serotype))
month_grid <- seq(min(muller_data$epi_week), max(muller_data$epi_week), by = "1 month")

counts_wide <- muller_data %>%
  select(epi_week, genotype_serotype, count) %>%
  complete(epi_week = month_grid, genotype_serotype = all_keys, fill = list(count = 0)) %>%
  pivot_wider(names_from = genotype_serotype, values_from = count, values_fill = 0) %>%
  arrange(epi_week)

# Gaussian kernel smoothing on counts
smooth_series <- function(y, bw) {
  x <- seq_along(y)
  ks <- ksmooth(x, y, kernel = "normal", bandwidth = bw, x.points = x)
  pmax(ks$y, 0)
}
smoothed_counts <- counts_wide
smoothed_counts[,-1] <- lapply(counts_wide[,-1], smooth_series, bw = bw_months)
smoothed_counts$total <- rowSums(smoothed_counts[,-1, drop = FALSE])

# smoothed proportion and windowed maximum
d3_prop <- smoothed_counts[[key]] / smoothed_counts$total
d3_prop[!is.finite(d3_prop)] <- 0

in_win <- counts_wide$epi_week >= win_start & counts_wide$epi_week <= win_end
if (!any(in_win)) stop("No dates within the 2022–2024 window in the data.")

imax <- which.max(replace(d3_prop, !in_win, -Inf))
peak_date <- counts_wide$epi_week[imax]
peak_val  <- d3_prop[imax]

cat("Kernel-smoothed max for", key, "in 2022–2024 =",
    round(peak_val, 3), " (", scales::percent(peak_val), ") at ",
    format(peak_date, "%Y-%m-%d"), "\n", sep = "")

#ggsave('DENV/2024/genotype_muller_plot.jpeg', dpi = 1200, height = 4.5, width = 12, units = 'in', device = 'jpeg')
