# Fitting SSFs with GAMs using `mgcv`

# Load packages ----
library(amt)
library(terra)
library(dplyr)
library(lubridate)
library(ggplot2)
library(mgcv)
library(gratia)
library(tidyr)
library(purrr)

# Load data ----
# Load the GPS data
gps <- read.csv("data/coyote_cougar.csv") %>% 
  mutate(t_ = ymd_hms(t_)) %>% 
  # Subset to the cougar F53
  filter(id == "F53")

# Load the habitat layers
hab <- rast("data/coyote_cougar_habitat.tif")
names(hab) <- c("elev", "tree", "biomass", "dist_to_road")

plot(hab)

# Prepare data ----
dat <- gps %>% 
  make_track(x_, y_, t_, all_cols = TRUE, crs = 32612) %>% 
  track_resample(rate = hours(4), tolerance = minutes(30)) %>% 
  steps_by_burst() %>% 
  random_steps(n_control = 30) %>% 
  extract_covariates(hab) %>% 
  # For mgcv, create "times" column and
  # explicitly create factor variables 
  mutate(times = 1,
         id = "F53",
         id_step_id_ = paste(id, step_id_, sep = "_"),
         id = factor(id),
         id_step_id_ = factor(id_step_id_)
  )

# Fit SSF ----
m <- gam(cbind(times, id_step_id_) ~ 
           # Habitat selection
           s(elev, bs = "cr", k = 50) +
           # Movement
           log(sl_) + sl_ + cos(ta_),
         family = cox.ph,
         data = dat,
         weights = case_)

summary(m)

draw(m)

# Predict log-RSS
x1 <- data.frame(elev = seq(2000, 3000, length.out = 100),
                 sl_ = 1000,
                 ta_ = 0)

x2 <- data.frame(elev = 2200,
                 sl_ = 1000,
                 ta_ = 0)

# Predict
g1 <- predict(m, newdata = x1, type = "link")
g2 <- predict(m, newdata = x2, type = "link")

# Log-rss
lr <- x1 %>% 
  mutate(log_rss = as.vector(g1) - as.vector(g2),
         rss = exp(log_rss))

# Plot
ggplot(lr, aes(x = elev, y = rss)) +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", 
             linetype = "dashed") +
  theme_bw()

# If you wanted to construct the confidence intervals with
# standard errors, you need the design matrix, coefficient
# vector, and variance-covariance matrix.

# Get the design matrix
predict(m, newdata = x1, type = "lpmatrix")
# Get the coefficient vector
coef(m)
# Get the variance-covariance matrix
vcov(m)

# Process data for multiple individuals ----
# Load the GPS data
multi <- read.csv("data/coyote_cougar.csv") %>% 
  mutate(t_ = ymd_hms(t_)) %>% 
  # Create a nested data.frame
  # See ?tidyr::nest (imported by `amt`)
  nest(data = x_:t_) %>% 
  # Create a new column with the track_xyt object
  mutate(trk = map(data, function(df) {
    df  %>% 
      make_track(x_, y_, t_, all_cols = TRUE, crs = 32612) %>% 
      track_resample(rate = hours(8), tolerance = minutes(30)) %>% 
      steps_by_burst() %>% 
      random_steps(n_control = 30) %>% 
      extract_covariates(hab)  %>% 
      return()
  })) %>% 
  select(id, trk) %>% 
  # Unnest back into a single data.frame
  unnest(cols = trk) %>% 
  mutate(id_step_id_ = factor(paste(id, step_id_, sep = "_")),
         id = factor(id),
         times = 1)

# Fit model with a random slope for elevation

m2 <- gam(cbind(times, id_step_id_) ~ 
            s(elev, bs = "cr", k = 20) +
            s(elev, id, bs = "re") +
            log(sl_) + sl_ + cos(ta_),
          data = multi,
          family = cox.ph,
          weights = case_)

summary(m2)

draw(m2)
