# ------------------------------------------------------------------------------
# Unit 3: Plotting Mapping ESDA
# ------------------------------------------------------------------------------

## SET UP
# Clean the environment
rm(list = ls())

# Create a function to install missing packages and load all required packages
ipack <- function(pack) {
  create.pkg <- pack[!(pack %in% installed.packages()[, "Package"])]
  if (length(create.pkg))
    install.packages(create.pkg, dependencies = TRUE)
  sapply(pack, require, character.only = TRUE)
}

# Specifcy packages and install and load them
packages <- c("sp", "sf",  "spData", "spdep", "rgdal", "tmap", "tidyverse")
ipack(packages)


## DATA
# Choose a dataset from the spData package
data("house")
sf.house <- st_as_sf(house) ; rm(house)
nb.house <- LO_nb ; rm(LO_nb)
trmat.house <- trMat ; rm(trMat)

# ------------------------------------------------------------------------------
# Data on 25,357 single family homes sold in Lucas County, Ohio, 1993-1998 from 
# the county auditor, together with an nb neighbour object constructed as a 
# sphere of influence graph from projected coordinates.
# ------------------------------------------------------------------------------

## EXPLORATION
# Plot the houses according to their prices and year built
tmap_mode('plot')
m3 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, OH", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "price", title = "House prices", palette = "plasma", style = "order")
m4 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, OH", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "yrbuilt", title = "Houses, year built", palette = "plasma", style = "order")
tmap_arrange(m3, m4)

# ------------------------------------------------------------------------------
# It looks like there are clusters of similar house prices and age in the city 
# center as well as in some parts of the suburbs.
# ------------------------------------------------------------------------------


## GLOBAL SPATIAL AUTOCORRELATION
# Convert nb object to listw object
lw.house <- nb2listw(nb.house, style = "W")

# Test for global autocorrelation using Moran's I and Geary's C
options(scipen = 99999)
moran.test(sf.house$price, 
           lw.house, 
           randomisation = T, 
           zero.policy = T)  
geary.test(sf.house$price, 
           lw.house, 
           randomisation = T, 
           zero.policy = T)  

moran.test(sf.house$yrbuilt, 
           lw.house, 
           randomisation = T, 
           zero.policy = T)  
geary.test(sf.house$yrbuilt, 
           lw.house, 
           randomisation = T, 
           zero.policy = T)  

# ------------------------------------------------------------------------------
# Moran's I is expected to be 0 under the assumption that there is spatial 
# randomness. Moran's I for house prices using randomisation is 0.8 with a 
# tiny p-value, meaning that there is significant global positive autocorrelation 
# present in the data. Geary's C, which is more sensitive to local spatial 
# autocorrelation and inversely related to Moran's I, also indicates global 
# spatial autocorrelation. It's expected value is 1, whereas the range lies between
# 0 and an unspecified number larger than 1, where values below 1 indicate
# positive autorrelation. Hence, the value of 0.19 indicated strong positive 
# spatial autocorrelation, just like Moran's I.

# The same results adhere to the ages of the houses, while the positive 
# autocorrelation is slightly less strong (0.77 for Moran's I and 0.22 for Geary's C)

# Further, the same results are achieved using a Monte-Carlo simulation with the
# command moran.mc().
# ------------------------------------------------------------------------------

# Moran's scatterplot plotted manually for greater flexibility
par(mfrow = c(2,1), mar = c(6, 5, 4, 4))
sf.house$price <- scale(sf.house$price)
sf.house$l1.price <- lag.listw(lw.house, sf.house$price)
plot(x = sf.house$price, y = sf.house$l1.price, main = "Moran Scatterplot", 
     pch = 19, cex = 0.5, col = rgb(0,0,0,0.2),
     xlab = "House prices",
     ylab = "Spatially lagged house prices")
abline(h = 0, v = 0, lty = 3)
abline(lm(sf.house$l1.price ~ sf.house$price), lty = 3, lwd = 2, col = "green")
sf.house$yrbuilt <- scale(sf.house$yrbuilt)
sf.house$l1.yrbuilt <- lag.listw(lw.house, sf.house$yrbuilt)
plot(x = sf.house$yrbuilt, y = sf.house$l1.yrbuilt, main = "Moran Scatterplot", 
     pch = 19, cex = 0.5, col = rgb(0,0,0,0.2),
     xlab = "Houses year built",
     ylab = "Spatially lagged houses year built")
abline(h = 0, v = 0, lty = 3)
abline(lm(sf.house$l1.yrbuilt ~ sf.house$yrbuilt), lty = 3, lwd = 2, col = "green")

# ------------------------------------------------------------------------------
# I decided to manually calculate and plot the Moran scatterplot, in order to 
# leverage the additional flexibility to remove the outlier marks that 
# automatically come with the plot, such that the plot is better readable. 
# In contrast to using the command moran.plot(x, wlist), the manual plots are 
# centered. They show a clear positive correlation for both prices and yrbuilt.
# ------------------------------------------------------------------------------

par(mfrow = c(2,2), mar = c(2,2,2,1) + 0.1)
plot(sp.correlogram(neighbours = nb.house,
                    var = as.vector(sf.house$price),
                    order = 4,
                    zero.policy = T,
                    method = "I", #for moran's I
                    style = "W"),
     main = "Correlogram Moran's I - Price")
plot(sp.correlogram(neighbours = nb.house,
                    var = as.vector(sf.house$yrbuilt),
                    order = 4,
                    zero.policy = T,
                    method = "I", #for moran's I
                    style = "W"),
     main = "Correlogram Moran's I - Year built")
plot(sp.correlogram(neighbours = nb.house,
                    var = as.vector(sf.house$price),
                    order = 4,
                    zero.policy = T,
                    method = "C", #for moran's I
                    style = "W"),
     main = "Correlogram Geary's C - Price")
plot(sp.correlogram(neighbours = nb.house,
                    var = as.vector(sf.house$yrbuilt),
                    order = 4,
                    zero.policy = T,
                    method = "C", #for moran's I
                    style = "W"),
     main = "Correlogram Geary's C - Year built")

# ------------------------------------------------------------------------------
# The correlogram shows the spatial autocorrelations statistic with increasing 
# number of lags. For both price and year built, as well as according to both
# Moran's I and Geary's C, the data sets remains spatially autocorrelated at the
# fourth order, depsite gradually decreasing.
# ------------------------------------------------------------------------------

## LOCAL SPATIAL AUTOCORRELATION
# Test for local autocorrelation using the local moran statistic and adjusting 
# the p-values adjusted with the bonferroni method
price.lmi <- localmoran(as.vector(sf.house$price),
                       lw.house,
                       zero.policy = T,
                       alternative = "greater",
                       p.adjust.method = "bonferroni")
yrbuilt.lmi <- localmoran(as.vector(sf.house$yrbuilt),
                        lw.house,
                        zero.policy = T,
                        alternative = "greater",
                        p.adjust.method = "bonferroni")

# Add the values to the data.frame and create factor variables indicating significance
sf.house$price.lmi.level <- price.lmi[, 1]
sf.house$price.lmi.p_value <- price.lmi[, 5]
sf.house <- sf.house %>% 
  mutate(
    price.lmi.significance = 
      case_when(price.lmi.p_value < 0.01 ~ "p < 0.01",
                price.lmi.p_value < 0.02 & 
                  price.lmi.p_value > 0.01 ~ "p < 0.02",
                price.lmi.p_value < 0.05 & 
                  price.lmi.p_value > 0.02 ~ "p < 0.05",
                TRUE ~ "not significant")
    ) %>% 
  mutate(price.lmi.significance = as.factor(price.lmi.significance))
sf.house$yrbuilt.lmi.level <- yrbuilt.lmi[, 1]
sf.house$yrbuilt.lmi.p_value <- yrbuilt.lmi[, 5]
sf.house <- sf.house %>% 
  mutate(
    yrbuilt.lmi.significance = 
      case_when(yrbuilt.lmi.p_value < 0.01 ~ "p < 0.01",
                yrbuilt.lmi.p_value < 0.02 & 
                  yrbuilt.lmi.p_value > 0.01 ~ "p < 0.02",
                yrbuilt.lmi.p_value < 0.05 & 
                  yrbuilt.lmi.p_value > 0.02 ~ "p < 0.05",
                TRUE ~ "not significant")
  ) %>% 
  mutate(yrbuilt.lmi.significance = as.factor(yrbuilt.lmi.significance))

# Plot the significance of the local autocorrelations using tmap
tmap_mode("plot")
col1 <- viridisLite::magma(4, begin = 0.12, end = 0.83) ; print(col1)
m1 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, Year built", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "yrbuilt.lmi.significance", title = "Local Moran's I Significance", palette = c("not significant" = "#1B1043FF",
                                                                                                "p < 0.01" = "#FEAD77FF",
                                                                                                "p < 0.02" = "#DB476AFF",
                                                                                                "p < 0.05" = "#7B2382FF"), legend.show = F) +
  tm_add_legend(type = "fill", 
                col = c("#1B1043FF", "#FEAD77FF", "#DB476AFF", "#7B2382FF"),
                labels = levels(sf.house$yrbuilt.lmi.significance),
                title = "Local Moran's I")
m2 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, Price", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "price.lmi.significance", title = "Local Moran's I Significance", palette = c("not significant" = "#1B1043FF",
                                                                                              "p < 0.01" = "#FEAD77FF",
                                                                                              "p < 0.02" = "#DB476AFF",
                                                                                              "p < 0.05" = "#7B2382FF"), legend.show = F) +
  tm_add_legend(type = "fill", 
                col = c("#1B1043FF", "#FEAD77FF", "#DB476AFF", "#7B2382FF"),
                labels = levels(sf.house$price.lmi.significance),
                title = "Local Moran's I")
tmap_arrange(m1, m2)

# ------------------------------------------------------------------------------
# This plot shows in which areas there is significant spatial autocorrelation
# where black is insignificant and orange is highly significant.Given that my
# data set does not have any polygons, but only the coordinates of the houses,
# I chose to plot the insignificant values in a neutral color and highlight 
# significant ones, where the brighter the color, the more significant the value
# is.
# ------------------------------------------------------------------------------

# Plotting the local spatial autocorrelation in levels using tmap
m5 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, Price", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "price.lmi.level", title = "Local Moran's I", palette = "plasma", style = "order")
m6 <- tm_shape(sf.house) +
  tm_layout(title = "Lucas County, Year built", frame = F, title.position = c("center", "top")) +
  tm_dots(col = "yrbuilt.lmi.level", title = "Local Moran's I", palette = "plasma", style = "order")
tmap_arrange(m5, m6)

# ------------------------------------------------------------------------------
# As we would intuitively expect, there isn't really much negative autocorrelation
# anywhere in Lucas County, OH. In contrast, there are a few areas of serious
# autocorrelation which could be explained historically, meaning the entire city 
# centre was built around the same time and the entire area has experienced 
# significant deterioation since, while other positively autocorrelated areas in 
# the suburb are most likely large scale housing construction projects, where 
# all the houses look the same, were built at the same time, and have the same
# structure - your typically American suburbia. In the north and south, there 
# seems to be a good mix of available houses according to price and age.
# ------------------------------------------------------------------------------