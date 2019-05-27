library(raster)
library(fields)

setwd("//uom-file3.unimelb.edu.au/6300/Users/plentini/Desktop/For Payal")

load("area.RData")
load("action.RData")
load("cost.RData")
load("luse.RData")
load("BT.RData")

################# SET UP PLOT
windows(15,10)
par (mfrow = c (2, 2), mar=c(1, 1, 2.5, 4), xpd = T)
col = grey.colors(3000, start = 0, end = 1, gamma = 2.2, alpha = NULL)
plot(BT, box = F, axes = F, legend = F, zlim = c(0,1), col = col)
mtext("Brown Treecreeper", side = 3, at = 143.925, cex = 0.8, line = 0.5)
image.plot(zlim = c(0, 1), legend.only=TRUE, smallplot=c(0.90, 0.91, 0.1, 0.9), axis.args = c(cex.axis = 0.7), col = col)
plot(area, add =T)

plot(luse, box = F, axes = F, legend = F, col = c("darkgreen", "red", "gold1", "gray87"))
mtext("Land use", side = 3, at = 143.865, cex = 0.8, line = 0.5)
par(xpd = T)
legend(144.64 ,-36.85, c("Public Land", "Veg pasture", "Cropping", "Mod pasture"), pt.bg = c("darkgreen", "red", "gold1", "gray87"), col = "black", pch = 22, bty = "n", cex = 0.8)
plot(area, add =T)

colours <- c("white", "darkorchid3", "green3")
plot(action, box = F, axes = F, col = colours, legend = F)
par(xpd = T)
plot(area, add =T)
mtext("Action", side = 3, at = 143.865, cex = 0.8, line = 0.5)
legend(144.64 ,-36.85, c("fencing", "revegetation"), pt.bg = c("darkorchid3", "green3"), col = "black", pch = 22, bty = "n", cex = 0.8)

plot(cost, box = F, axes = F, legend = F, col = rev(terrain.colors(1000, alpha = 1)), zlim = c(0, 1000))
mtext("Cost", side = 3, at = 143.865, cex = 0.8, line = 0.5)
image.plot(zlim = c(0, 1000), legend.only=TRUE, smallplot=c(0.90, 0.91, 0.1, 0.9), axis.args = c(cex.axis = 0.7), col = rev(terrain.colors(1000, alpha = 1)))
plot(area, add =T)

dev.copy2pdf(file = "Action cost SDM.pdf")