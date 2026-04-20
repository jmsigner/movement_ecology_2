#######################################################X
#------------- Movement Ecology with R 2 --------------X
#--------------Module 01 -- Introduction --------------X
#----------------Last updated 2026-04-19---------------X
#-------------------Code Walkthrough-------------------X
#######################################################X

# Load packages ----
library(tidyverse)
library(amt)
library(lubridate) 
library(sf)

# Load data ---- 

# We will use data from fishers that were tracked by Scott
# LaPoint and have been widely used for many methodological comparisons, examples
# etc. The data is saved in a `csv`-file in the data directory. 


# Create tracks ----

# The basic building block to work with the `amt` package are so called
# tracks. A track consists of a series of relocations. The function
# `make_track()` is used to create a track from a `data.frame` or `tibble`. It
# expects at least the data set, and coordinates (x and y). Optionally time
# stamps, additional columns and a coordinate reference system (crs) can be
# passed to the function using the EPSG code.

# Load data from the amt package

data(amt_fisher)
head(amt_fisher)

# Note, that a track is characterized by `x_`, `y_` and `t_`. We could add a crs
# using the crs argument with the EPSG code.

tr <- make_track(amt_fisher, x_, y_, t_, id = name, crs = 5070)
class(tr)

# ... Changing the CRS ----

# We can change the the CRS with the function `transform_coords()`. For example
# to change to geographic coordinates, we could just use:

tr |> transform_coords(4326)

# ... Visually inspect a track -----

tr |> inspect()

# ... Sampling rate ----

# The rate at which data are sampled for tracks can be different and irregular.
# To get an overview of the sampling rate the function
# `summarize_sampling_rate()` exists.

leroy <- filter(tr, id == "Leroy")

summarize_sampling_rate(leroy)

# This suggests that the median sampling rate is 15 min. We can now resample the
# track a 15 min interval.

leroy2 <- track_resample(leroy, rate = minutes(15), tolerance = seconds(60))
leroy2


# Fitting a single MCP home range -----

hr1 <- leroy2 |> hr_mcp()
plot(hr1)

# Dealing with many animals -----------
tr |> nest(data = -id)
tr |> nest(gps = -id)

# Adding home ranges
tr |> nest(data = -id) |> 
  mutate(hr = map(data, hr_mcp))

# Visualize
tr |> nest(data = -id) |> 
  mutate(hr = map(data, hr_mcp)) |> hr_to_sf(hr, id) |> 
  ggplot() + geom_sf(aes(fill = id))
