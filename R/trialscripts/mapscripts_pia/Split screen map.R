setwd("//uom-file3.unimelb.edu.au/6300/Users/plentini/Desktop/For Payal")
library(raster)

load("vic.Rdata")
load("region.Rdata")
load("PU.Rdata")
load("zones.Rdata")
load("links.Rdata")
load("studyveg.Rdata")
load("veg.Rdata")

col <- c("grey98", "darkseagreen4", "grey", "darkseagreen", "lightblue1", "skyblue2")
breaks <- c(0, 1, 2, 3, 4, 5, 6)

windows(10,8)
image(veg, col = col, axes = F, xlab = "", ylab = "", breaks = breaks)
plot(vic, add = 2, lwd = 2)
plot(region, add = 2, lwd = 2, col = rgb(0, 0, 0, 0.5))

windows(10,8)
plot(studyveg, legend = F, axes = F, box = F, col = c("grey98", "darkseagreen4"))
plot(region, add = T, lwd = 2)
plot(zones, add = T, lwd = 2)
plot(PU, add = T, lwd = 1)
plot(links, add = T, lwd = 1, col = "white", border = "red")



################################ PLOTTING THEM IN THE SAME WINDOW USING "split.screen" - the 'layout' function also works well for this! https://stackoverflow.com/questions/38810854/how-to-use-layout-function-in-r

windows(10,9)
split.screen(rbind(c(0, 0.7, 0, 0.5), c(0.3, 0.95, 0.4, 0.95)))
par(mar = c(0, 0, 0, 0))

screen(1)
col <- c("grey98", "darkseagreen4", "grey", "darkseagreen", "lightblue1", "skyblue2")
breaks <- c(0, 1, 2, 3, 4, 5, 6)
image(veg, col = col, axes = F, xlab = "", ylab = "", breaks = breaks)
plot(vic, add = 2, lwd = 2)
plot(region, add = 2, lwd = 2, col = rgb(0, 0, 0, 0.5))

screen(2)
par(mar = c(0, 0, 0, 0))
col <- c("grey98", "darkseagreen4", "grey", "darkseagreen", "lightblue1", "skyblue2")
breaks <- c(-1, 0, 1)
image(studyveg, axes = F, col = c("grey98", "darkseagreen4"), breaks = breaks, xlab = "", ylab = "")
plot(region, add = T, lwd = 2)
plot(zones, add = T, lwd = 2)
plot(PU, add = T, lwd = 1)
plot(links, add = T, lwd = 1, col = "white", border = "red")


close.screen(split.screen())

