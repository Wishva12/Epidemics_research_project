library(spatstat)
library(sf)
library(ggplot2)
library(viridis)
library(dplyr)

data <- read.csv("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/All Data.csv")

clean_data <- data[complete.cases(data[, c("Dengue_A", "lon", "lat")]), ]
clean_data <- clean_data[clean_data$Dengue_A > 0, ]

district_data <- aggregate(cbind(Dengue_A, lon, lat) ~ RDHS, data = clean_data, FUN = mean)

sri_lanka_sf <- st_read("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/lka_admbnda_adm2_slsd_20220816.shp")

sri_lanka_bbox <- st_bbox(sri_lanka_sf)
win <- owin(xrange = c(sri_lanka_bbox$xmin, sri_lanka_bbox$xmax), 
            yrange = c(sri_lanka_bbox$ymin, sri_lanka_bbox$ymax))

pp <- ppp(x = district_data$lon, 
          y = district_data$lat, 
          window = win,
          marks = district_data$Dengue_A)

kde_result <- density(pp, sigma = 0.5, weights = marks(pp))  

plot(kde_result, main = "Kernel Density Estimation of Dengue Cases in Sri Lanka")
plot(pp, add = TRUE, pch = 20, col = "black", cex = 1.5)

sri_lanka_df <- st_coordinates(sri_lanka_sf)[,1:2] %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y)


expanded_points <- district_data[rep(row.names(district_data), 
                                     round(district_data$Dengue_A/10)), 
                                 c("lon", "lat", "RDHS", "Dengue_A")]

ggplot() +
  geom_sf(data = sri_lanka_sf, fill = "white", color = "black", size = 0.3) +
  stat_density_2d(data = expanded_points, 
                  aes(x = lon, y = lat, fill = after_stat(level)),
                  geom = "polygon", alpha = 0.6, contour = TRUE, bins = 10) +
  scale_fill_viridis_c(name = "Dengue\nDensity", option = "plasma") +
  geom_point(data = district_data, 
             aes(x = lon, y = lat, size = Dengue_A), 
             color = "red", alpha = 0.7) +
  scale_size_continuous(name = "Cases", range = c(1, 6)) +
  coord_sf(expand = FALSE) +
  labs(title = "Kernel Density Estimation of Dengue Cases",
       subtitle = "Sri Lankan Districts with Case Intensity",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        legend.title = element_text(size = 10))

ggplot() +
  geom_sf(data = sri_lanka_sf, fill = "white", color = "darkgray", size = 0.3) +
  stat_density_2d_filled(data = expanded_points, 
                         aes(x = lon, y = lat), 
                         alpha = 0.7, contour = TRUE, bins = 8) +
  scale_fill_viridis_d(name = "Density\nLevel", option = "inferno") +
  geom_point(data = district_data, 
             aes(x = lon, y = lat), 
             color = "white", size = 2, stroke = 1, shape = 21, fill = "red") +
  coord_sf(expand = FALSE) +
  labs(title = "Dengue Cases - Kernel Density Distribution",
       subtitle = "Higher density areas indicate clustering of cases",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

cat("KDE Analysis Summary:\n")
cat("Number of districts:", nrow(district_data), "\n")
cat("Total dengue cases:", sum(district_data$Dengue_A), "\n")
cat("Mean cases per district:", round(mean(district_data$Dengue_A), 2), "\n")
cat("Bandwidth (sigma) used:", 0.5, "\n")

density_at_points <- kde_result[list(x = district_data$lon, y = district_data$lat)]
high_density_threshold <- quantile(density_at_points, 0.75)
high_density_districts <- district_data[density_at_points > high_density_threshold, ]

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

plot(st_geometry(sri_lanka_sf), border = "gray", main = "Original Dengue Data")
points(district_data$lon, district_data$lat, 
       pch = 21, bg = heat.colors(10)[cut(district_data$Dengue_A, breaks = 10)], 
       cex = sqrt(district_data$Dengue_A)/5)
legend("topright", legend = c("Low", "Medium", "High"), 
       pt.bg = heat.colors(3), pch = 21, title = "Dengue Cases")

plot(kde_result, main = "KDE Surface")
plot(pp, add = TRUE, pch = 20, col = "white", cex = 0.8)

par(mfrow = c(1, 1))