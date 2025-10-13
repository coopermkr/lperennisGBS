#''''''''''''''''''''''''''''''''''''''''''''''''
#' 
#' Mapping Sample Sites
#' @date 2025-06-18
#' @author Cooper Kimball-Rhines
#' 
#''''''''''''''''''''''''''''''''''''''''''''''''


# Load libraries
library(maptools)
library(tidyverse)
library(sf)
library(ggrepel)
library(cowplot)

# Load Northeast shapefile
northeast <- st_read("4.IBD/Northeastern_States/Northeast_State_Polygon.shp")


# Convert coords to lat long and filter out western MA
northeast <- st_transform(northeast, 4326)

northeast <- northeast |>
  filter(STATE_COD %in% c("NH", "NY", "MA", "VT"))

neMap <- ggplot(data = northeast) +
  geom_sf(color = "darkgrey")
neMap

# Load in site coords
sites <- read_csv("4.IBD/nesitecoords.csv")

palMap <- c("#E2A3C7", "#778da9", "#EC7D10", "#63A46C")


# Make sample site map
siteMap <- neMap +
  geom_point(data=sites, aes(x=Longitude, y=Latitude, color = State), 
             size=2.5, show.legend = FALSE) +
  geom_label_repel(data=sites, aes(x=Longitude, y=Latitude, label=Population), nudge_y = 0, nudge_x = -0.03) +
  theme_classic(base_size = 16) +
  xlab("Longitude") + ylab("Latitude") +
  scale_color_manual(values = palMap)

siteMap

inbredHier <- read_csv("3.popgen/inbreeding.csv") |>
  mutate(Region = c("Vermont", "New Hampshire", "New Hampshire", "New York",
                    "Florida", "Florida", "Mid", "Vermont", "New Hampshire",
                    "Mid", "New Hampshire", "Mid", "Massachusetts",
                    "Florida", "Mid", "Mid", "New York"),
  pop = factor(pop, levels = c("AL", "SA", "MO", "AB", "CN", "HK", "92-HK", "VT", "92-VT")))

# Just Northeast region Fis with confidence intervals
inbredRegions <- inbredHier |>
  filter(Region %in% c("Massachusetts", "New Hampshire", "New York", "Vermont")) |>
  ggplot(mapping = aes(x = pop, y = Fis, color = Region)) +
  geom_point(size = 2) +
  geom_errorbar(mapping = aes(ymin = ll, ymax = hl)) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(#title = "Population Inbreeding (95% CI)",
       x = "") +
  ylab(bquote(F[IS])) +
  guides(color = "none") +
  scale_color_manual(values = palMap)


# Inset inbreeding inside the map
inset <- ggdraw() +
  draw_plot(siteMap) +
  draw_plot(inbredRegions, x=0.12, y=0.18, height = 0.6, width = 0.43)

jpeg(filename = "4.IBD/inbreedingMap.jpg", height = 6, width = 8.5, units = "in", res = 300)
inset
dev.off()

