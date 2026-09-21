# 16 September 2026
# Global home range analysis

## Set up ----------------------------------------------------------------------

# load packages
library(adehabitatHR)
library(beepr)
library(lubridate)
library(purrr)
library(readxl)
library(scales)
library(sf)
library(tidyverse)

# load data
obs <- read_excel("data/PromObs_2008-2024_RB.xlsx")

# clean coordinates & timestamps
obs_clean <- obs %>%
  mutate(X_num = as.numeric(X), Y_num = as.numeric(Y)) %>%
  drop_na(X_num, Y_num) %>%
  
  # coordinate offset for AGD66 / AMG Zone 55
  mutate(X_AGD66 = X_num + 400000, Y_AGD66 = Y_num + 5600000) %>%
  
  # spatial & seasonal (Jul 10 to Jan 31) filters
  filter(X_num > 36612 & X_num < 38590, Y_num > 88032 & Y_num < 88929) %>%
  mutate(Date = as_date(Date), Year = year(Date), yday = yday(Date)) %>%
  filter(yday <= 31 | yday >= 191)


## Calculate global fixed bandwidth --------------------------------------------

# calculate median href across all years to eliminate small-sample-size inflation
obs_sf_init <- st_as_sf(obs_clean, coords = c("X_AGD66", "Y_AGD66"), crs = 20255)
obs_sp_init <- as(obs_sf_init, "Spatial")

k_ud_init <- kernelUD(obs_sp_init[, "Year"], h = "href", same4all = TRUE)
h_values <- map_dbl(k_ud_init, ~ .x@h$h)
global_h <- median(h_values)

cat("Global fixed bandwidth (median href):", round(global_h, 2), "meters\n")


## Define bootstrap function ---------------------------------------------------

get_boot_area <- function(year_df, h_val) {
  # resample observations with replacement
  boot_df <- year_df %>% slice_sample(prop = 1, replace = TRUE)
  
  # convert to SpatialPoints object
  sp_obj <- st_as_sf(boot_df, coords = c("X_AGD66", "Y_AGD66"), crs = 20255) %>% 
    as("Spatial")
  
  # run KDE with fixed bandwidth
  k_ud <- kernelUD(sp_obj, grid = 200, extent = 0.25, h = h_val)
  hr <- getverticeshr(k_ud, percent = 95, unin = "m", unout = "km2")
  
  return(hr$area)
}


## Spatial analysis (KDE) ------------------------------------------------------

set.seed(123)
n_boot <- 100

hr <- obs_clean %>%
  group_by(Year) %>%
  nest() %>%
  mutate(boot_areas = map(data, function(df){
    replicate(n_boot, get_boot_area(df, h_val = global_h), simplify = TRUE)
  })) %>%
  unnest(boot_areas) %>%
  group_by(Year) %>%
  summarise(
    Area_mean = mean(boot_areas),
    Area_median = median(boot_areas),
    Area_se = sd(boot_areas),
    Area_low = quantile(boot_areas, 0.025),
    Area_high = quantile(boot_areas, 0.975),
    .groups = "drop"
  )

print(hr)
beep(2)

# # save output
# write_csv(hr, "data/globalHR_to2024_hfix.csv")


## Plot results ----------------------------------------------------------------

# # load data
# hr <- read_csv("data/globalHR_to2024.csv")

hr %>% 
  mutate(Area = Area_mean*100, LCI = Area_low*100, UCI = Area_high*100) %>% 
  ggplot(aes(x = Year, y = Area)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), fill = "#7D9570", alpha = 0.2) +
  geom_line(color = "#7D9570", linewidth = 1) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(y = "Area (ha)", x = "Year") +
  theme_bw()

# ggsave("figures/areaVStime.png", width = 18.0, height = 10.0, units = c("cm"), dpi = 600)

# compare href & global h approaches
hr_og <- read_csv("data/globalHR_to2024.csv") %>%
  mutate(Model = "Dynamic bandwidth (href)")

hr_gh <- read_csv("data/globalHR_to2024_hfix.csv") %>%
  mutate(Model = "Global fixed bandwidth")

# combine & format for plotting
hr <- bind_rows(hr_og, hr_gh) %>% 
  mutate(Area = Area_mean * 100, LCI = Area_low * 100, UCI = Area_high * 100)

# generate comparison plot
ggplot(hr, aes(x = Year, y = Area, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  scale_color_manual(values = c("Dynamic bandwidth (href)" = "#D68D38", 
                                "Global fixed bandwidth" = "#7D9570")) +
  scale_fill_manual(values = c("Dynamic bandwidth (href)" = "#D68D38", 
                               "Global fixed bandwidth" = "#7D9570")) +
  labs(y = "Area (ha)", x = "Year", color = NULL, fill = NULL) +
  theme_bw()

# ggsave("figures/areaVStime_compare.png", plot = p_compare, width = 18.0, height = 10.0, units = "cm", dpi = 600)


## Multi-panel figure ----------------------------------------------------------

# generate baseline contours
obs_sf <- st_as_sf(obs_clean, coords = c("X_AGD66", "Y_AGD66"), crs = 20255)
obs_sp <- as(obs_sf, "Spatial")

k_ud_plot <- kernelUD(obs_sp[, "Year"], grid = 300, extent = 0.25, h = global_h, same4all = TRUE)
hr_95_plot <- getverticeshr(k_ud_plot, percent = 95, unin = "m", unout = "km2")

# merge bootstrap data
hr_polygons <- st_as_sf(hr_95_plot) %>%
  # st_set_crs(20255) %>% # if necessary
  st_transform(st_crs(obs_sf)) %>%
  transmute(Year = as.integer(as.character(id))) %>%
  left_join(hr_uncertainty, by = "Year") %>%
  mutate(
    title_label = paste0(
      Year, "\n", round(Area_mean, 3), " ",
      "[", round(Area_low, 3), " - ", round(Area_high, 3), "]"
    )
  )

obs_sf_labeled <- obs_sf %>%
  left_join(st_drop_geometry(hr_polygons), by = "Year") %>%
  filter(!is.na(title_label)) %>%
  st_as_sf() # %>%
  # st_set_crs(20255) # if necessary

# create plot
p_grid_boot <- ggplot() +
  geom_sf(data = hr_polygons, fill = "gray30", alpha = 0.3, color = "black") +
  geom_sf(data = obs_sf_labeled, size = 0.3, alpha = 0.4) +
  facet_wrap(~ title_label, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray92", color = NA),
    strip.text = element_text(face = "bold", size = 8, lineheight = 1.2),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ); p_grid_boot

ggsave("figures/globalHR_boot.pdf", plot = p_grid_boot, width = 10, height = 10.2)
ggsave("figures/globalHR_boot.png", width = 24.0, height = 24.5, units = "cm", dpi = 600)

