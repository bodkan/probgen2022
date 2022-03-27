# This script contains coordinates of polygon population ranges and
# migration trajectories used in example 3. We store this separately
# to make the ex3.R script easier to read.

map <- world(xrange = c(-15, 60), yrange = c(20, 65), crs = 3035)

R1 <- region(
  "EHG range", map,
  polygon = list(c(26, 55), c(38, 53), c(48, 53), c(60, 53),
                 c(60, 60), c(48, 63), c(38, 63), c(26, 60))
)
R2 <- region(
  "Europe", map,
  polygon = list(
    c(-8, 35), c(-5, 36), c(10, 38), c(20, 35), c(25, 35),
    c(33, 45), c(20, 58), c(-5, 60), c(-15, 50)
  )
)
R3 <- region(
  "Anatolia", map,
  polygon = list(c(28, 35), c(40, 35), c(42, 40),
                 c(30, 43), c(27, 40), c(25, 38))
)
R4 <- join(R2, R3)
R5 <- region(
  "YAM range", map,
  polygon = list(c(26, 50), c(38, 49), c(48, 50),
                 c(48, 56), c(38, 59), c(26, 56))
)

ooa_trajectory <- list(c(40, 30), c(50, 30), c(60, 40), c(45, 55))