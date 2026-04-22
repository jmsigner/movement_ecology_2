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

# Plot from the 'gratia' package
draw(m)

# Predict log-RSS
x1 <- data.frame(elev = seq(2000, 3000, length.out = 100),
                 sl_ = 1000,
                 ta_ = 0)

x2 <- data.frame(elev = 2200,
                 sl_ = 1000,
                 ta_ = 0)

# Calculate the linear predictor
g1 <- as.vector(predict(m, newdata = x1, type = "link"))
g2 <- as.vector(predict(m, newdata = x2, type = "link"))

# log-RSS is just the difference in linear predictors
lr <- x1 %>% 
  mutate(log_rss = g1 - g2,
         # Exponentiate to get RSS
         rss = exp(log_rss))

# Plot
ggplot(lr, aes(x = elev, y = rss)) +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", 
             linetype = "dashed") +
  labs(x = "Elevation (m)", y = "RSS") +
  theme_bw()

# For models fitted with glm() or clogit(), amt::log_rss() 
# will calculate confidence intervals for you. We don't (yet)
# have support for models fitted with gam().

# If you wanted to construct the confidence intervals with
# standard errors, you need the design matrix, coefficient
# vector, and variance-covariance matrix.

# Get the design matrix
X1 <- predict(m, newdata = x1, type = "lpmatrix")
X2 <- predict(m, newdata = x2, type = "lpmatrix")

# Calculate a matrix of the difference between X1 and X2
D <- sweep(X1, 2, X2)

# Get the coefficient vector
B <- coef(m)
# Get the variance-covariance matrix
S <- vcov(m)

# Incidentally, we can also calculate log-RSS from the 
# difference matrix and the vector of betas
LR <- D %*% B
all.equal(as.vector(LR), lr$log_rss)

# Variance for the prediction
var_pred <- as.vector(diag(D %*% S %*% t(D)))
# SE for the prediction
SE <- sqrt(var_pred)

# Calculate 95% confidence interval
lr <- lr %>% 
  mutate(se = SE,
         log_lwr = log_rss - 1.96*se,
         log_upr = log_rss + 1.96*se) %>% 
  # Exponentiate to get RSS
  mutate(rss = exp(log_rss),
         lwr = exp(log_lwr),
         upr = exp(log_upr))

# Plot
ggplot(lr, aes(x = elev, y = rss)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "gray80") +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, color = "red", 
             linetype = "dashed") +
  geom_vline(xintercept = x2$elev, linetype = "dashed") +
  labs(x = "Elevation (m)", y = "RSS") +
  theme_bw()

# Example with multiple individuals ----
# Note there are only 3 individuals in the data
# Note that they are different species!
#   - C028 = coyote
#   - F53 & F64 = cougars

# Load the GPS data and format
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
# This takes ~ 1.5 minutes to fit.

system.time({
  m2 <- gam(cbind(times, id_step_id_) ~ 
              # Habitat selection
              s(elev, bs = "cr", k = 20) +
              # Random slopes
              s(elev, id, bs = "re") +
              # Movement
              log(sl_) + sl_ + cos(ta_),
            data = multi,
            family = cox.ph,
            weights = case_)
}) # 105 sec


summary(m2)

draw(m2)
