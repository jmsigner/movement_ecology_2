#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#---------------Module 05 -- Interpretation -----------X
#--------------- Last updated 2026-04-20 --------------X
#------------------ Code Walkthrough ------------------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(terra)

set.seed(1323)



# ... Calculating Relative Selection Strength (RSS) for Two Locations -----
#
# Calculating the relative use of location $s_1$ versus location $s_2$ is fairly
# straightforward when $s_1$ and $s_2$ share the same values for all but one
# covariate. For more complex scenarios, we have implemented a function in `amt`
# that will calculate the log-relative intensity [referred to as the
# log-Relative Selection Strength, or log-RSS, see Avagar et al 2017]. If you
# prefer to quantify relative-use (i.e., RSS), you can simply exponentiate the
# results.
 
# Here we demonstrate the use  of `log_rss()`: 

#   *  $s_1$: `elevation_end = 100`, `forest = 0`
#   *  $s_2$: `elevation_end = 100`, `forest = 1`

# We want to calculate the log-RSS for a step ending in $s_1$ versus a step
# ending in $s_2$, assuming $s_1$ and $s_2$ are equally accessible to the
# animal.

# The function `log_rss()` expects the two locations to be formatted as separate
# `data.frame` objects. Each `data.frame` must include all covariates used to
# fit the model. Note, Furthermore, `factor` variables should have the same
# `levels` as the original data.
 
# Let us start with creating the two positions.
s1 <- data.frame(
  elevation = 100,
  forest = TRUE)
 
# data.frame for s2; note the value for forest is different.
s2 <- data.frame(
  elevation = 100,
  forest = FALSE)

# Now that we have specified each location as a `data.frame`, we can pass them
# along to `log_rss()` for the calculation. The function will return an object
# of class `log_rss`, which is also more generally a `list`. The `list` element
# `"df"` contains a `data.frame` which contains the log-RSS calculation and
# could easily be used to make a plot when considering relative selection
# strength across a range of environmental characteristics.
lr1 <- log_rss(m1, x1 = s1, x2 = s2)

lr1$df

exp(lr1$df$log_rss)

# This is the same as: 
summary(m1)
exp(coef(m1)["forestTRUE"]) # See Module 5 for this

# `log_rss()` is designed to be able to consider several locations as `x1`,
# relative to a **single** location in `x2`. 

# Lets next consider the scenario where we want to compare locations that are
# outside of forests, but from a range of different elevations

# data.frame for s1
s1 <- data.frame(
  elevation = 10:150,
  forest = FALSE)
 
# data.frame for s2
s2 <- data.frame(
  elevation = 100,
  forest = FALSE)

lr1 <- log_rss(m1, x1 = s1, x2 = s2)

head(lr1$df)
exp(lr1$df)

ggplot(lr1$df, aes(elevation_x1, log_rss)) + 
  geom_line() +
  geom_hline(yintercept = 0, col = "red", lty = "dashed")

ggplot(lr1$df, aes(elevation_x1, exp(log_rss))) + 
  geom_line() +
  geom_hline(yintercept = 1, col = "red", lty = "dashed")

# Large-Sample Confidence Intervals
#
# We can use the standard errors from our fitted model to estimate these
# confidence intervals based on a normal approximation to the sampling
# distribution of $\hat{\beta}$.
lr1_ci_se <- log_rss(m1, s1, s2, ci = "se", ci_level = 0.95)
head(lr1_ci_se$df)

ggplot(lr1_ci_se$df, aes(elevation_x1, log_rss)) + 
  geom_hline(yintercept = 0, col = "red", lty = "dashed") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() 

# If want the RSS instead of log-RSS we have to exponentiate the results again.
# One way to do this is to use `mutate()` in combination with `across`.
lr1_ci_se$df |> mutate(across(log_rss:upr, exp)) |> 
  ggplot(aes(elevation_x1, log_rss)) + 
  geom_hline(yintercept = 1, col = "red", lty = "dashed") +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() 