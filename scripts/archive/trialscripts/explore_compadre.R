load("/Users/payal/Google Drive/Jobs/TB Discovery Project/Data/COMPADRE_v.4.0.1.RData")

# Explore data
ls()
head(compadre$metadata)
head(compadre$mat)

data <- compadre
# Plot by coordinates
plot(data$metadata$Lon, data$metadata$Lat, type='n') 
text(data$metadata$Lon, data$metadata$Lat, data$metadata$NumberPopulations)
  # not much use



coordinates(data$metadata) <- c(data$metadata$Lon, data$metadata$Lat)
