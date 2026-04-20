#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#------------------Module 03 -- iSSF ------------------X
#--------------- Last updated 2026-04-17 --------------X
#------------------ Code Walkthrough ------------------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(terra)

set.seed(1323)


# Integrated step-selection analyses are implemented using the following steps:

# 1. Estimate a tentative selection-free movement kernel, using observed step-lengths and turn angles, giving.
# 2. Generate time-dependent available locations by simulating potential movements from the previously observed location.
# 3. Estimate $\beta$ using conditional logistic regression, with strata formed by combining time-dependent used and available locations.
# 4. Re-estimate the movement parameters.

# Simulate data for one animal ----

# ... Landscape ----

# Simulate covariates. We start with forest and simulate it as Gaussian random
# filed using the default settings. In order to create a binary landscape (0 and
# 1), we take 0.5 as a threshold which should lead on average to approx. 50 %
# forest.

r <- rast(xmin = -100, xmax = 100, ymin = -100, ymax = 100)
values(r) <- rnorm(ncell(r))

# smooth it 
r <- focal(r, w = matrix(1, 15, 15), fun = mean) 
r <- focal(r, w = matrix(1, 15, 15), fun = mean) 
r <- focal(r, w = matrix(1, 15, 15), fun = mean) 
forest <- as.numeric(r > 0)
plot(forest)


# We do the same for elevation, 
r <- rast(xmin = -100, xmax = 100, ymin = -100, ymax = 100)
values(r) <- rnorm(ncell(r))

# smooth it 
r <- focal(r, w = matrix(1, 15, 15), fun = mean) 
r <- focal(r, w = matrix(1, 15, 15), fun = mean) 
ele <- r
ele[] <- scales::rescale(ele[], c(0, 500))
plot(ele)

# Finally, we create a stack and give the layers a meaningful name.
covars <- c(forest, ele)
names(covars)
names(covars) <- c("forest", "elevation")
names(covars)
plot(covars)

# ... Movement 
# Next we model the movement of the animal using a exponential distribution for
# the step length and a uniform distribution for the turn angles.
curve(dexp(x, rate = 0.5), from = 0, to = 20)


# We will discuss simulations in much more detail in module 11.
dat1 <- redistribution_kernel(
  make_issf_model(coefs = c(forest_end = 0.5, elevation_end = 0.005), 
                  sl = make_exp_distr(rate = 0.5)), 
  start = make_start(c(0, 0)), 
  map = covars, n.control = 1e4) |> 
  simulate_path(n = 500)

# Inspect it visually
plot(ele)
points(dat1)

plot(forest)
points(dat1)
head(dat1)


# Preparing data for SSF ----
# ... Creating steps ----

# We have to convert the point representation of the track that we obtain from
# `amt::make_track()` to steps. Ideally, we would before ensure a regular
# sampling rate and go through the steps introduced earlier this week for data
# cleaning.
dat1 <- make_track(dat1, x_, y_, t_)
dat1 |> steps()

# Adding random steps. 
tmp <- dat1 |> steps() |> random_steps() 
tmp |> print(n = 15)

# Check the step-length and turn-angle distribution that was fitted to the data.
# We see later how we can plot this. As you notice, the step length distribution
# that was fitted to the data is a Gamma distribution and not an exponential.
# The reason for this is, that by default a Gamma distribution is used. We will
# see later how we can change this. But it is also no problem, since the
# exponential distribution is a special case of the Gamma distribution.
sl_distr(tmp)
ta_distr(tmp)

# We can vary the number of control steps
dat1 |> steps() |> random_steps(n_control = 20) 

# We can now extract the covariates
dat1 |> steps() |> 
  random_steps() |> 
  extract_covariates(covars) 

# Finally, we can some additional covariates
dat1 |> steps() |> 
  random_steps() |> 
  extract_covariates(covars) |> 
  mutate(log_sl_ = log(sl_), 
         cos_ta_ = cos(ta_))

# Everything at once
ssf.dat <- dat1 |> steps() |> 
  random_steps() |> 
  extract_covariates(covars) |> 
  mutate(log_sl_ = log(sl_), 
         cos_ta_ = cos(ta_))

ssf.dat

# - `x1_` and `y1_` are the coordinates associated with the starting location of the step.
# - `x2_` and `y2_` are the coordinates associated with the ending location of the step.
# - `sl_` and `ta_` are the step-length and turn angle associated with the step.
# - `t1_` and `t2_` are the time stamps associated with the start and end of the step.
# - `dt_` is the time duration associated with the step.
# - `step_id_` is a unique identifier associated with each step's observed and random locations.
# - `case_` is an indicator variable equal to `TRUE` for observed steps and `FALSE` for random steps.
# - `forest` and `elevation` are the covariate values at the **end** of the step. 


# ... Tentative Movement Parameters   ----
# Let's have a closer look at what is happening under the hood when applying the
# `random_steps()` function. This function conveniently fits tentative
# step-length and turn-angle distributions and then samples from these
# distributions to generate random steps. The default arguments lead to
# `random_steps()` fitting  a Gamma distribution for the step lengths and a von
# Mises distribution for the turn angles. It is possible to use other
# distributions.
#
# The tentative parameters in these statistical distributions (gamma, von Mises)
# are stored as attributes of the resulting object. We can view the parameters
# of the tentative step-length distribution using:
sl_distr(ssf.dat)

# Directly access the parameter of the distribution
sl_distr(ssf.dat)$params

# And plot them
curve(dgamma(x, shape = sl_distr_params(ssf.dat)$shape, 
             scale = sl_distr_params(ssf.dat)$scale), from = 0.1, to = 30
)
# Add the distribution that we used to simulate the data
curve(dexp(x, rate = 0.2), add = TRUE, col = "red")

# The function `random_steps()` has an argument called `sl_distr` that is set to
# gamma by default, but can be changed to any other suitable supported
# distribution.
sl_exp <- fit_distr(dat1 |> steps() |> pull(sl_), "exp")
dat1 |> steps() |> 
  random_steps(sl_distr = sl_exp) |> sl_distr()

curve(dexp(x, rate =  sl_exp$params$rate), from = 0, to = 30
)

# Add the distribution that we used to simulate the data
curve(dexp(x, rate = 0.2), add = TRUE, col = "red")

# Similarly, we can access the parameters of the tentative turn-angle
# distribution using:
ta_distr(ssf.dat)

# In the next module on issf 2, we will learn how we can update these
# distributions to adjust for the influence of habitat selection when estimating
# parameters of the movement kernel.
 
# Basic iSSF   ----
# We will fit a model with the same habitat covariates used for the simulations.
# Suppose we hypothesize that habitat selection at the step scale depends on
# forest (`forest`) and elevation (`elevation`). Furthermore, suppose we
# hypothesize that, in the absence of habitat selection, step lengths and turn
# angles follow a constant distribution, i.e., the selection-free movement
# kernel does not depend on covariates. We can fit that model as follows:
 
# Note the use of strata in the model formula. This indicates that each step
# (together with its random steps) forms a strata.
m1 <- ssf.dat |> 
  fit_issf(case_ ~ forest + elevation +
             strata(step_id_), model = TRUE)


# Interpreting Habitat-Selection Parameters ------
#
# We begin by exploring the coefficients in the fitted model using the
# `summary()` function:

summary(m1)


# Fitting an integrated SSF
m2 <- ssf.dat |> 
  fit_issf(case_ ~ forest + elevation + sl_ + log_sl_ + cos_ta_ +
             strata(step_id_), model = TRUE)

summary(m2)

