# Basic use of rgplates using the GPlates Web Service
# Analytical Paleo Workshop
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
library("rgplates")
library(sf)
library(terra)

# the positions of erlangen
erlangen <- matrix(c(10.97, 49.58),  ncol=2)

# coordinate reconstruction to target age
reconstruct(erlangen, age=150)

# Plate reconstruction 
plates <- reconstruct("plates", age=150)
class(plates)

#  Just the plates
sfPlates <- as(plates, "sf")
plot(sfPlates$geometry, col="gray")
axis(1, at=seq(-180, 180, 30), label= paste(seq(-180, 180, 30), "°"))
axis(2, at=seq(-90, 90, 30), label= seq(-90, 90, 30))
box()



