library(spdep)
library(sp)
library(sf)
library(ggplot2)
library(maps)
library(mapdata)

data <- read.csv("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/All Data.csv")

clean_data <- data[complete.cases(data[, c("Dengue_A", "lon", "lat")]), ]
clean_data <- clean_data[clean_data$Dengue_A > 0, ]

district_data <- aggregate(cbind(Dengue_A, lon, lat) ~ RDHS, data = clean_data, FUN = mean)

coordinates <- cbind(district_data$lon, district_data$lat)
colnames(coordinates) <- c("lon", "lat")
spdf <- SpatialPointsDataFrame(coordinates, district_data)

k <- 3
coords_matrix <- coordinates(spdf)
nb <- knn2nb(knearneigh(coords_matrix, k = k))

lw <- nb2listw(nb, style = "B", zero.policy = TRUE)  

log_dengue <- log(district_data$Dengue_A)

global_g_result <- globalG.test(log_dengue, lw, zero.policy = TRUE)
print("Global Getis-Ord G test results:")
print(global_g_result)

global_g_stat <- globalG.test(log_dengue, lw, zero.policy = TRUE)$statistic
cat("Global G statistic:", global_g_stat, "\n")

nb_include_self <- include.self(nb)
lw_gi <- nb2listw(nb_include_self, style = "B", zero.policy = TRUE)

local_gi_result <- localG(log_dengue, lw_gi, zero.policy = TRUE)
print("Local Gi* test results summary:")
print(summary(local_gi_result))

sri_lanka_sf <- st_read("G:/Academic/4000 Level/Semester 1/Research/Coding from sir/lka_admbnda_adm2_slsd_20220816.shp")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

plot(st_geometry(sri_lanka_sf), border = "grey", main = "Sri Lankan Districts - Spatial Structure")
plot(nb, coordinates, add = TRUE, col = "red", lwd = 1.5)
points(coordinates, pch = 20, col = "darkblue", cex = 1.5)
text(coordinates, labels = district_data$RDHS, pos = 3, cex = 0.6, col = "black")

dengue_colors <- heat.colors(10, rev = TRUE)[cut(district_data$Dengue_A, breaks = 10)]
plot(st_geometry(sri_lanka_sf), border = "grey", main = "Dengue_A Cases Intensity")
points(coordinates, pch = 21, bg = dengue_colors, cex = 2)
plot(nb, coordinates, add = TRUE, col = "gray", lwd = 0.5)
legend("topright",
       legend = c("Low", "Medium", "High"),
       pt.bg = heat.colors(3, rev = TRUE),
       pch = 21, pt.cex = 1.5, title = "Dengue Cases", cex = 0.8)

log_colors <- terrain.colors(10)[cut(log_dengue, breaks = 10)]
plot(st_geometry(sri_lanka_sf), border = "grey", main = "log(Dengue_A) Distribution")
points(coordinates, pch = 21, bg = log_colors, cex = 2)
plot(nb, coordinates, add = TRUE, col = "gray", lwd = 0.5)

gi_colors <- ifelse(local_gi_result > 1.96, "red",          # Significant hot spots
                    ifelse(local_gi_result < -1.96, "blue",   # Significant cold spots
                           "lightgray"))                       # Not significant
plot(st_geometry(sri_lanka_sf), border = "grey", main = "Local Gi* Hot/Cold Spots")
points(coordinates, pch = 21, bg = gi_colors, cex = 2)
plot(nb, coordinates, add = TRUE, col = "gray", lwd = 0.5)
legend("topright",
       legend = c("Hot Spot", "Cold Spot", "Not Significant"),
       pt.bg = c("red", "blue", "lightgray"),
       pch = 21,
       pt.cex = 1.5,
       title = "Gi* Results",
       cex = 0.8)

par(mfrow = c(1, 1))


dengue_colors <- heat.colors(10, rev = TRUE)[cut(district_data$Dengue_A, breaks = 10)]
plot(st_geometry(sri_lanka_sf), border = "lightgray", lwd = 1.5,
     main = "Sri Lankan Districts: Dengue Cases and Spatial Connections",
     cex.main = 1.2, cex.lab = 1.1)
plot(nb, coordinates, add = TRUE, col = "red", lwd = 2)
points(coordinates, pch = 21, bg = dengue_colors, cex = 3)
text(coordinates, labels = district_data$RDHS, pos = 1, cex = 0.7, col = "black", font = 2)
legend("topright",
       legend = c("High Dengue", "Medium Dengue", "Low Dengue", "Spatial Links"),
       col = c("red", "orange", "yellow", "red"),
       pch = c(21, 21, 21, NA),
       lty = c(NA, NA, NA, 1),
       pt.bg = c("red", "orange", "yellow", NA),
       pt.cex = c(2, 2, 2, NA),
       lwd = c(NA, NA, NA, 2),
       title = "Legend",
       cex = 0.9)


print("Summary Statistics")
cat("Number of districts analyzed:", nrow(district_data), "\n")
cat("Mean Dengue_A cases per district:", round(mean(district_data$Dengue_A), 2), "\n")
cat("Median Dengue_A cases per district:", round(median(district_data$Dengue_A), 2), "\n")
cat("Standard deviation:", round(sd(district_data$Dengue_A), 2), "\n")
cat("Mean log(Dengue_A):", round(mean(log_dengue), 3), "\n")

significant_gi <- abs(local_gi_result) > 1.96  
hot_spots <- local_gi_result > 1.96
cold_spots <- local_gi_result < -1.96

cat("Number of significant hot spots:", sum(hot_spots), "\n")
cat("Number of significant cold spots:", sum(cold_spots), "\n")

plot(st_geometry(sri_lanka_sf), border = "grey", main = "Significant Hot/Cold Spots (Gi*)")
points(coordinates, pch = 21,
       bg = ifelse(hot_spots, "red", 
                   ifelse(cold_spots, "blue", "lightgray")),
       cex = ifelse(significant_gi, 2.5, 1.5))
plot(nb, coordinates, add = TRUE, col = "gray", lwd = 1)
text(coordinates[significant_gi, ],
     labels = district_data$RDHS[significant_gi],
     pos = 3, cex = 0.8, 
     col = ifelse(hot_spots[significant_gi], "darkred", "darkblue"), 
     font = 2)
legend("topright",
       legend = c("Hot Spot", "Cold Spot", "Not Significant"),
       pt.bg = c("red", "blue", "lightgray"),
       pch = 21,
       pt.cex = c(2.5, 2.5, 1.5),
       title = "Gi* Results")