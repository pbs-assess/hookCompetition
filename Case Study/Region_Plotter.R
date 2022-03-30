library(sf)
library(rgeos)
library(tidyverse)
# COAST <- readOGR("~/OneDrive - The University Of British Columbia/Whale Project/Data for Joe/Coast shapefile",
# layer="coast_ALBERS")
# COAST <- rgeos::gSimplify(COAST, tol=1000)
COAST <- readRDS('Hires_BC_Coastline.rds')
load("C:/Users/WATSONJOE/OneDrive - DFO-MPO/Documents/hookcompetition/data/yelloweye_rockfish.rda")
# Convert to sf format
stations <- yelloweye_rockfish$set_counts
stations <-
  sf::st_as_sf(stations,
               coords=c('lon','lat'),
               crs=st_crs(4326))
plot(stations)
# Subset to keep only a single observation from each station (remove duplicates)
stations <-
  stations %>%
  group_by(station) %>%
  filter(row_number()==1 &
           standard == 'Y' &
           year > 1996)

stations <-
  sf::st_transform(stations,crs=COAST@proj4string)

COAST_sf <- st_as_sf(COAST)

stations %>%
  ggplot() +
  geom_sf(colour='yellow') +
  ylim(c(range(st_coordinates(stations)[,2]))) +
  xlim(c(range(st_coordinates(stations)[,1]))) +
  #theme(panel.background = element_blank()) +
  theme(panel.background = element_rect(fill='light blue')) +
  geom_sf(data=COAST_sf, fill='light green', alpha=0.8) +
  ggtitle('The locations of all 174 IPHC fishing stations') +
  xlab('Longitude') + ylab('Latitude')



