#Create animated KDE map from dengue cases in Monteria
library(tidyterra)

libs <- c(
  "tidyverse", "terra", "tidyterra",
  "osmdata", "sf", "ggmap", "classInt",
  "gifski"
)

installed_libraries <- libs %in% rownames(
  installed.packages()
)

if(any(installed_libraries == F)){
  install.packages(
    libs[!installed_libraries]
  )
}

invisible(lapply(
  libs, library, character.only = T
))


#Get data

setwd("~/Documents/DENV/Database/modeling/")

data <- readr::read_csv("geocoded.csv")

#Select only Cordoba and transform dates to epiweeks and months
CORDOBA <- data %>% filter(ndep_resi == "CORDOBA") %>%
  drop_na(lon) %>%
  mutate(epi_week = floor_date(fec_not, unit = "week"),
         month = floor_date(fec_not, unit="month"))

monteria <- CORDOBA %>%
  tidyterra::filter(lon > -77 & lon < -74 & lat > 7 & lat < 9)

df <- st_as_sf(monteria, coords = c('lon','lat'), crs = 4326, agr = "constant" )

df <- df %>% mutate(epi_week = floor_date(fec_not, unit = "week"), 
                    month = floor_date(fec_not, unit="month"))

coordinates(monteria) <- ~ lon + lat

#grid monteria data
library(akima)
lon_grid <- seq(min(monteria$lon), max(monteria$lon), length.out = 100)
lat_grid <- seq(min(monteria$lat), max(monteria$lat), length.out = 100)
grid_points <- expand.grid(lon_grid, lat_grid)



#base R (takes a lot to compute)
kde_data <- ks::kde(monteria[, c("lon", "lat")], compute.cont = FALSE)
kde_df <- expand.grid(lon = kde_data$x$lon, lat = kde_data$x$lat)
kde_df$density <- as.vector(kde_data$H)

kde_plot <- ggmap(map) + #does not render
  geom_tile(data = kde_df, aes(x = lon, y = lat, fill = density), width = 0.01, height = 0.01) +
  scale_fill_gradient(low = "green", high = "red") +
  theme_bw()


cor <- st_as_sf(monteria, coords = c('lon','lat'), crs = 4326, agr = "constant" )
# 
# library(spatialEco)
# kde.output <- sp.kde(x=cor, bw=100,
#                                  standardize = TRUE, scale.factor = 1000)
# plot(kde.output)
# writeRaster(kde.output, "KDE_den.tif
#             ")
# 
# 
# projection(kde) <- CRS("+proj=longlat +datum=WGS84")
# 
# 
# 
# tmap::tm_shape(kde) + tmap::tm_raster("ud")

# v1 <- vect(CORDOBA)
# r1 <- rast(v1)
# nams <- names(v1)
# 
# 
# allrast <- lapply(nams, function(x) {
#   rasterize(v1, r1,
#             field = x,
#             touches = TRUE
#   )
# })
# allrast <- do.call("c", allrast)
# 
# 
# 
# coordinates(CORDOBA) <- ~ lon + lat
# raster <- terra::rast(CORDOBA)
location = "Monteria"
map = get_googlemap(center = location, zoom = 13, style = paste0("feature:all|element:labels|visibility:off"))

kde_plot <- ggmap(map) +
  geom_point(data = monteria,
             aes(x=lon, y=lat, group = interaction(month, lat), colour=factor(año)), alpha = 0.8)+
  stat_density2d(data = monteria, aes(x = lon, y = lat, fill = after_stat(density), alpha = after_stat(density)), 
                 geom = "tile", binwidth = c(0.01, 0.01), contour = FALSE) +
  scale_fill_gradient(low = "green", high = "red") +
  scale_alpha(range = c(0.1, 0.6)) +
  facet_wrap(~factor(año)) +
  theme_bw()

dev.off()

kde_plot

#Code for animated point map and leaflet interactive map (this section works)

library(leaflet)
library(leaflet.providers)
library(leaflet.extras2)

pal <- c('#555599', '#66BBBB', '#DD4444')
factpal <- colorFactor(pal, df$año)

m <- leaflet() %>%
  addTiles() %>%
  addProviderTiles("Stamen.TonerLite") %>%
  fitBounds(-76.946225,7.432629,-74.474301,8.769072) %>%
  addCircles(data = df, weight = 1, opacity = 1, fillColor = ~factpal(año), color = ~factpal(año), fillOpacity = 2)

#Animate the map
library(gganimate)
library(gifski)
library(ggplot2)

location = "Monteria"
map = get_googlemap(center = location, zoom = 13, style = paste0("feature:all|element:labels|visibility:off"))

data_map <- ggmap(map, extent = "panel") + geom_point(data = CORDOBA,
                                                      aes(x=lon, y=lat, group = interaction(month, lat), colour=factor(año)), alpha = 0.8) +
  stat_density_2d(data = CORDOBA,aes(x=lon, y=lat),
                  geom = "density_2d", n = 10)+
  scale_color_manual(values = c('#555599', '#66BBBB', '#DD4444'))+
  facet_wrap(vars(año))+
  labs(colour="Year")+
  xlab("Lat")+
  ylab("Lon")+
  theme(text = element_text(family = "Times New Roman", size = 14),
        legend.spacing.y = unit(1, "cm"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  guides(fill = guide_legend(byrow = TRUE))



density_map_dengue <- qmplot(x=lon, y=lat,
                             data = CORDOBA, geom = "blank",
                             legend = "topright") +
  stat_density2d(
    geom = "polygon", 
    alpha = .5,
    color = NA) +
  scale_fill_gradient2(low = "blue", 
                       mid = "green", 
                       high = "red")




map_with_data <- ggmap(map, extent = "panel") + 
  geom_point(data = CORDOBA,
             aes(x=lon, y=lat, group = interaction(month, lat), colour=factor(año)), alpha = 0.8) +
  scale_color_manual(values = c('#555599', '#66BBBB', '#DD4444'))+
  labs(colour="Year")+
  xlab("Lat")+
  ylab("Lon")+
  theme(text = element_text(family = "Times New Roman", size = 14),
        legend.spacing.y = unit(1, "cm"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  guides(fill = guide_legend(byrow = TRUE))


location2 = "Cordoba, Colombia"
map2 = get_map(location = location2, zoom = 9)
map_with_data2 = ggmap(map2) + geom_point(data = CORDOBA,
                                          aes(x=lon, y=lat, group = interaction(month, lat), colour=factor(año))) +
  scale_color_manual(values = c('#555599', '#66BBBB', '#DD4444'))

map_with_animation <- map_with_data +
  shadow_mark(color="red", alpha=0.3, size=0.2) +
  transition_time(CORDOBA$epi_week) +
  ggtitle('Epiweek: {frame_time}',
          subtitle = 'Frame {frame} of {nframes}')

num_months <- max(CORDOBA$month) - min(CORDOBA$month) + 1
animate(map_with_animation, nframes = num_weeks, fps = 30, res=600, height=10, width=10, units="in")
anim_save("denv_cases.gif")



