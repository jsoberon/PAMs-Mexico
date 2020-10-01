# ------------------------------------------------------------------------------
# Project title: Visualizing Species Richness and Site Similarity from 
#                Presence-absence Matrices
# Authors: Jorge Soberon, Marlon E. Cobos, Claudia Nunez-Penichet
# Date: 1/10/2020 (dd/mm/yyyy)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Details:
#
# This script serves to replicate the analyses of for the paper ""Visualizing 
# Species Richness and Site Similarity from Presence-absence Matrices"
# 
# This script and the data used for analysis is stored in a GitHub repository:
# https://github.com/jsoberon/PAMs-Mexico
# 
# The data is a presence-absence matrix (PAM) for 1573 species of terrestrial 
# vertebrates of Mexico. Mexico was divided in a grid of 711 cells of 0.5 decimal
# degrees for purposes of representation. Species ranges used to construct the PAM were obtained from 
# the IUCN database.
#
# The analyses presented here are performed using base functions of R, and 
# other functions from the packages biosurvey, maps, and viridis (for colors).
# 
# biosurvey is a new R package available from GitHub (not yet in CRAN), 
# instructions for installing this package are provided below.
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Packages 
## installing remotes to install biosurvey using it
install.packages("remotes")

## installing biosurey (you may be asked to update other packages; not updating 
## those packages may prevent potential errors)
remotes::install_github("claununez/biosurvey")

## loading other packages
library(biosurvey)
library(maps)

install.packages("viridis")
library(viridis)

# function to produce old diagrams
source("https://raw.githubusercontent.com/jsoberon/PAMs-Mexico/master/Function_old_diagrams.R")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Data for analysis
## setting directory
setwd("R/Christen-Soberon") # change it to your working directory

## this is the presence absence matrix  
pam_url <- "https://github.com/jsoberon/PAMs-Mexico/blob/master/PAM_CS_MEX.RDATA?raw=true"
load(url(pam_url))

iucn[1:6, 1:8]
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Analysis 
## performing calculations to produce Christen-Soberon diagrams and their 
## geographic representations. Depending on your computer this process could be 
## time consuming (to make it faster you can use fewer 'randomization_iterations')
iucnCS <- prepare_PAM_CS(PAM = iucn, exclude_column = 1:2, significance_test = TRUE, 
                         randomization_iterations = 500, CL = 0.05, 
                         keep_randomizations = TRUE, parallel = TRUE)

## the list produced contains all elements needed to produce the following plots

# saving results as RData
save(iucn, iucnCS, file = "PAM_CS_MEX.RDATA")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Plotting the PAM and its geographic representation

## colors for richness
cols <- magma(length(unique(iucnCS$Richness_normalized)))
colfact <- as.factor(iucnCS$Richness_normalized)

## colors for dispersion field
cols1 <- viridis(length(unique(iucnCS$Dispersion_field_normalized)))
colfact1 <- as.factor(iucnCS$Dispersion_field_normalized)

## part of PAM
pampart <- t(as.matrix(iucn[15:1, 3:37]))

## plot configuration
x11()
layout(matrix(c(1, 1:3), nrow = 2, byrow = T))

## plotting part of a PAM
par(mar = c(0.5, 2.5, 3, 0.5))
image(t(as.matrix(iucn[, -(1:2)])), axes = FALSE)

axis(2, at = 0.5, labels = "Cells", tick = FALSE)
axis(3, at = 0.5, labels = "Species", tick = FALSE)
axis(3, at = 0.1, labels = "Presence-absence matrix\n", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species richness normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols[colfact], pch = 19)
legend_image <- as.raster(matrix(rev(cols), ncol = 1))
text(x = -115, y = 22, labels = "Richness")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
axis(3, at = -111, labels = "Geographic representation", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species dispersion field normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols1[colfact1], pch = 19)
legend_image <- as.raster(matrix(rev(cols1), ncol = 1))
text(x = -115, y = 22, labels = "Dispersion field")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Christen-Soberon diagrams vs previous diagram
x11()
par(mfrow = c(1, 2))
par(mar = c(4.5, 4.5, 3.5, 0.5))

## original diagram
rdp(iucn[, -(1:2)], view = 1, limits = 2)
title(main = "Christen diagram")

## simple Christen Soberon plot  
plot_PAM_CS(iucnCS, main = "Christen-Soberon diagram")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Christen-Soberon diagram options
colsg <- ifelse(iucnCS$S_significance_id == 0, "gray70", 
                ifelse(iucnCS$S_significance_id == 1, "gray35", "gray1"))

x11()
par(mfrow = c(2, 2), cex = 0.8)
par(mar = c(4.5, 4.5, 3.5, 0.5))

## simple Christen Soberon plot  
plot_PAM_CS(iucnCS, main = "Simple diagram")

## Christen Soberon plot showing randomized values   
plot_PAM_CS(iucnCS, main = "Diagram and random expectations",
            add_random_values = T)

## Christen Soberon plot showing significant   
plot_PAM_CS(iucnCS, main = "Diagram and significant values", add_significant = T, 
            col_significant_low = "gray35", col_significant_high = "gray1")

## geographic representation
map(regions = "Mexico")
title(main = "Geographic representation", line = 2)
points(iucn[, 1:2], col = colsg, pch = 19)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Christen-Soberon diagram, exploration in geography
## partitioning diagram area in blocks
n_cols <- 4
n_rows <- 4
data <- data.frame(iucnCS$Richness_normalized, iucnCS$Dispersion_field_normalized)
id <- paste(data[, 1], data[, 2])

xrange <- range(data[, 1])
xinter <- diff(xrange)/n_cols
yrange <- range(data[, 2])
yinter <- diff(yrange)/n_rows
xlb <- seq(xrange[1], xrange[2], xinter)
xlb[length(xlb)] <- xrange[2]
ylb <- seq(yrange[1], yrange[2], yinter)
ylb[length(ylb)] <- yrange[2]
blocks <- assign_blocks(data, 1, 2, n_cols, n_rows, xlb, ylb, block_type = "equal_area")
blocks <- blocks[match(id, paste(blocks[, 1], blocks[, 2])), ]

## sample of 4 blocks to plot
set.seed(1)
sam4 <- sample(unique(blocks$Block), 4)

## block colors 
colsbg1 <- ifelse(blocks$Block == sam4[1], "dodgerblue4", "gray55")
colsbg2 <- ifelse(blocks$Block == sam4[2], "dodgerblue4", "gray55")
colsbg3 <- ifelse(blocks$Block == sam4[3], "dodgerblue4", "gray55")
colsbg4 <- ifelse(blocks$Block == sam4[4], "dodgerblue4", "gray55")


## plotting
x11()
layout(matrix(1:18, nrow = 6), heights = c(1.2, 10, 10, 10, 10, 1), 
       widths = c(1, 10, 10))
par(cex = 0.7)

## Christen Soberon diagram in blocks 
### y labels
par(mar = rep(0, 4))
plot.new()
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new()

### title 1
plot.new(); text(0.5, 0.5, labels = "Blocks in the diangram", font = 2, cex = 1)

### CS plots
par(mar = c(2.5, 2, 0.5, 0.5))
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg1)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg2)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg3)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg4)

### x label
par(mar = rep(0, 4))
plot.new(); text(0.5, 0.5, labels = "Normalized richness")

### title 2
plot.new(); text(0.5, 0.5, labels = "Geographic representation", font = 2, cex = 1)

### mapping blocks
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg1, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg2, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg3, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg4, pch = 19)
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Exporting figures
## PAM 
jpeg(filename = "Figure_1.jpg", width = 166, height = 130, units = "mm", res = 600)

## plot configuration
layout(matrix(c(1, 1:3), nrow = 2, byrow = T))

## plotting part of a PAM
par(cex = 0.7)
par(mar = c(0.5, 2.5, 3, 0.5))
image(t(as.matrix(iucn[, -(1:2)])), axes = FALSE)

axis(2, at = 0.5, labels = "Cells", tick = FALSE)
axis(3, at = 0.5, labels = "Species", tick = FALSE)
axis(3, at = 0.1, labels = "Presence-absence matrix\n", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species richness normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols[colfact], pch = 19)
legend_image <- as.raster(matrix(rev(cols), ncol = 1))
text(x = -115, y = 21.5, labels = "Richness")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)
axis(3, at = -110.5, labels = "Geographic representation", font = 2, 
     tick = FALSE, cex = 1.2)

## plotting species dispersion field normalized in the geography
par(mar = c(0.5, 0.5, 0.5, 0.5))
map(regions = "Mexico")
points(iucn[, 1:2], col = cols1[colfact1], pch = 19)
legend_image <- as.raster(matrix(rev(cols1), ncol = 1))
text(x = -115, y = 21.5, labels = "Dispersion field")
text(x = -113, y = c(16, 20), labels = c("Low", "High"), cex = 0.8)
rasterImage(legend_image, -116, 16, -115, 20)

dev.off()



## Christen Soberon plots 
jpeg(filename = "Figure_2.jpg", width = 166, height = 90, units = "mm", res = 600)

par(mfrow = c(1, 2))
## original diagram
par(cex = 0.6)
rdp(mat = iucn[, -(1:2)], view = 1, limits = 2)
title(main = "Christen diagram")

## simple Christen Soberon plot  
par(cex = 0.6)
plot_PAM_CS(iucnCS, main = "Christen-Soberon diagram")

dev.off()



## Christen Soberon options 
jpeg(filename = "Figure_4.jpg", width = 166, height = 166, units = "mm", res = 600)

par(mfrow = c(2, 2), cex = 0.7)

## simple Christen Soberon plot  
plot_PAM_CS(iucnCS, main = "Simple diagram")

## Christen Soberon plot showing randomized values   
plot_PAM_CS(iucnCS, main = "Diagram and random expectations", 
            add_random_values = T)

## Christen Soberon plot showing significant   
plot_PAM_CS(iucnCS, main = "Diagram and significant values", add_significant = T, 
            col_significant_low = "gray35", col_significant_high = "gray1")

## geographic representation
map(regions = "Mexico")
title(main = "Geographic representation", line = 3.3)
points(iucn[, 1:2], col = colsg, pch = 19)

dev.off()



## Christen Soberon diagram, blocks explorations
jpeg(filename = "Figure_5.jpg", width = 140, height = 200, units = "mm", res = 600)
layout(matrix(1:18, nrow = 6), heights = c(1.2, 10, 10, 10, 10, 1), 
       widths = c(1, 10, 10))
par(cex = 0.7)

## Christen Soberon diagram in blocks 
### y labels
par(mar = rep(0, 4))
plot.new()
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new(); text(0.5, 0.5, labels = "Norm. dispersion field", srt = 90)
plot.new()

### title 1
plot.new(); text(0.5, 0.5, labels = "Blocks in the diangram", font = 2, cex = 1)

### CS plots
par(mar = c(2.5, 2, 0.5, 0.5))
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg1)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg2)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg3)
plot_PAM_CS(iucnCS, main = "", add_legend = FALSE, xlab = "", ylab = "", 
            col_all = colsbg4)

### x label
par(mar = rep(0, 4))
plot.new(); text(0.5, 0.5, labels = "Normalized richness")

### title 2
plot.new(); text(0.5, 0.5, labels = "Geographic representation", font = 2, cex = 1)

### mapping blocks
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg1, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg2, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg3, pch = 19)
map(regions = "Mexico"); points(iucn[, 1:2], col = colsbg4, pch = 19)

dev.off()
# ------------------------------------------------------------------------------
