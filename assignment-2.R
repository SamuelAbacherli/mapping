
# Specify and store the first-order contiguity matrix W and the corresponding 
# row-standardized weights matrix as matrix objects.
library(imager)
library(rgdal)
library(spdep)

img.austria <- load.image("../data/austria.png")
plot(img.austria, axes = F)

# Shapefile of Austria with MGI Austria Lambert projection
shp <- readOGR(dsn = "../data", layer ="NUTS_RG_03M_2013_4326_LEVL_2")  # reads in the shapefile
shp <- subset(shp, CNTR_CODE == "AT")
shp <- spTransform(shp, CRS("+init=epsg:31287"))
proj4string(shp)  # check the projection
plot(shp)

# Unstandardised first order queen contiguity matrix
nb.1q <- poly2nb(shp, row.names = shp$NUTS_NAME, queen = T)
W.matrix.1q <- nb2mat(nb.1q, style = "B")
colnames(W.matrix.1q) <- rownames(W.matrix.1q) ; View(W.matrix.1q)

# Check that first order rook is the same as queen
nb.1r <- poly2nb(shp, row.names = shp$NUTS_NAME, queen = F)
W.matrix.1r <- nb2mat(nb.1r, style = "B")
all(W.matrix.1r == W.matrix.1q)

# Standardised first order queen contiguity matrix
W.matrix.1q.std <- nb2mat(nb.1q, style = "W")
colnames(W.matrix.1q.std) <- rownames(W.matrix.1q.std) ; View(W.matrix.1q.std)

# Second order row standardized weights matrix as a matrix
W.list.2q.std = nblag(nb.1q, 2)
W.matrix.2q.std <- nb2mat(W.list.2q.std[[2]], style = "B") + 
  nb2mat(W.list.2q.std[[1]], style = "B")
t <- rowSums(W.matrix.2q.std)
W.matrix.2q.std <- W.matrix.2q.std / t
colnames(W.matrix.2q.std) <- rownames(W.matrix.2q.std) ; View(W.matrix.2q.std)
  
# First order queen matrix row standardized as a listw object
W.list.1q.std <- mat2listw(W.matrix.1q.std) ; W.list.1q.std

# Calculating the object characteristics by hand
nrow(W.matrix.1q)  # Number of regions
sum(rowSums(W.matrix.1q))  # Number of nonzero links
sum(rowSums(W.matrix.1q)) / (ncol(W.matrix.1q) * nrow(W.matrix.1q)) * 100  # Percentage nonzero weights
mean(rowSums(W.matrix.1q))  # Average number of links
nrow(W.matrix.1q)  # Number of zones n
nrow(W.matrix.1q) * ncol(W.matrix.1q)  # n * n (nn)
sum(W.matrix.1q.std)  # Global sum of weights S0

object.size(W.list.1q.std)
object.size(W.matrix.1q.std)

# The list consists of a style parameter, where "M" indicates a conversion from a matrix,
# a list of the neighbours of each area indicated through area indexes
# and a list of the corresponding weights

# Despite the listw format currently being larger than the matrix format, at a certain
# size the listw format is probably much more efficient in storing all information than the
# corresponding matrix would be. Furthermore, I assume that performing calculations
# on the weight matrix will be slower than on its corresponding listw format.
