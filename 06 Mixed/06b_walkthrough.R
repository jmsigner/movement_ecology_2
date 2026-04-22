#######################################################X
#----Analysis of Animal Movement Data in R Workshop----X
#----------- Module 06 -- Multiple animals ------------X
#----------------Last updated 2026-04-20---------------X
#-------------------Code Walkthrough-------------------X
#-----Created by Johannes Signer and John Fieberg------X
#######################################################X

library(tidyverse)
library(amt)
library(broom)
library(patchwork)
library(here)
library(mvtnorm)
library(glmmTMB)
library(broom.mixed)
library(terra)
library(tictoc)
library(boot)
library(broom.mixed)
library(survival)

# Prepare data --------------
#' 
#' Prepare data for fitting SSFs to each of several individuals
#' or a mixed SSF
#' 
#' We will use consider a GPS telemetry dataset containing observations from 20 red deer (*Cervus elaphus*). 
#' 
#' We will consider data derived from the following sources:
#'  - landuse data from the ESA worldmap of the study area
#'  - elevation using a digital elevation model of the study area
#'  
#' ## Attribution
#' 
#' If you should consider these data in any future publications, please make sure to give proper attribution. For the ESA World Cover data: 
#' 
#' - Zanaga, D., Van De Kerchove, R., Daems, D., De Keersmaecker, W., Brockmann, C., Kirches, G., Wevers, J., Cartus, O., Santoro, M., Fritz, S. , Lesiv, M. , Herold, M., Tsendbazar, N.-E., Xu, P., Ramoino, F., & Arino, O. (2022). ESA WorldCover 10 m 2021 v200. 10.5281/zenodo.7254221.
#'
#' For the digital elevation model:
#' 
#' - https://doi.org/10.5270/ESA-c5d3d65
#' 
#' For the GPS data and any scripts: 
#' https://zenodo.org/records/18911728
#' And please talk to us (jsigner@uni-goettingen.de)
#'
dat <- read_csv("data/reddeer/gps/gps.csv")
head(dat)

#' Load environmental data
ele <- rast("data/reddeer/env/elevation.tif")
esa <- rast("data/reddeer/env/esa_worldmap.tif")

#' Create derived quantities from elevation (slope, aspect and terrain ruggedness index)
ter <- terrain(ele, v = c("slope", "TRI", "aspect"))
ter

#' And we want the distance to the next urban area
unique(esa)

#' Landcover clases: https://collections.sentinel-hub.com/worldcover/readme.html
#'
#' -  10: Tree cover
#' -  20: Shrubland
#' - 30: Grassland
#' -  40: Cropland
#' -  50: Built up
#' -  60: Bare, sparse vegetation
#' -  70: Snow and ice
#' -  80: permanent water
#' -  90: Wetland
#' - 100: Moss and lichen

#'  We are interested in class 50 (Built up).  Create a raster that sets all
#'  other class to NA.
built_up <- terra::subst(esa == 50, from = FALSE, to = NA)
dist_built_up <- distance(built_up)
names(dist_built_up) <- "dist_built_up"

#' Plot the rasters
plot(esa)
plot(built_up)
plot(dist_built_up)

#' We are also interested in distance to forest.  Use the same strategy for it.
forest <- terra::subst(esa == 10, from = FALSE, to = NA)
dist_forest <- distance(forest)
names(dist_forest) <- "dist_forest"
plot(dist_forest)

#' ## Prepare data for an iSSF. 
#' 
#' Rather than fit tentative step length and turn angle distributions to 
#' each individual separately, we will pool the data and fit a single gamma
#' distribution to the **pooled** step length data and a single von-Mises
#' distribution to the pooled turn angle data. This will make "fixed effects" 
#' parameters in mixedSSA's more interpretable = > they will correspond to
#' movement parameters for a "typical individual" (one with all random effects set to 0).  
#' 
#' 
#' We will make use of nested data frames in R, which will allow us to apply 
#' functions separately to data from each individual using the map functions in the Purr
#' package. If you are unfamiliar with `purrr` syntax, may want to view one or more of 
#' the tutorials, below, or make use of the [purrr cheat sheet](https://github.com/rstudio/cheatsheets/blob/master/purrr.pdf).
#'
#' - http://www.rebeccabarter.com/blog/2019-08-19_purrr/
#' - https://www.r-bloggers.com/2020/05/one-stop-tutorial-on-purrr-package-in-r/
#' - https://jennybc.github.io/purrr-tutorial/index.html
#'
#' We have to take the following steps:
#'
#' 1. Create a track
trk <- make_track(dat, x_, y_, t_, id = id, sex = sex, crs = 3035)

#' 2. Resample data: Lets resample to 200 min with a tolerance of 10 min
summarize_sampling_rate_many(trk, "id", time_unit = "min")

#' Nest the data set, so we can work with steps for each individual,
#' 
#' - Resample to 200 +/- 10 minutes
#' - Filter to only include animals with > 1000 steps
#' - Create a steps by burst version of the data
steps <- trk |> nest(data = -c(id, sex)) |> 
  mutate(
    data_resample = map(
      data, ~ track_resample(
        .x, rate = minutes(200), tolerance = minutes(10)) |> 
        filter_min_n_burst()), 
    # Calculate the number of relocations per animal
    n = map_int(data_resample, nrow)
  ) |> 
  filter(
    n > 1000
  ) |> 
  # Create steps by burst
  mutate(
    steps = map(data_resample, steps_by_burst)
  ) |> 
  # Select only relevant columns
  dplyr::select(id, sex, steps) |> 
  # Unnest
  unnest(cols = steps)

#' 3. Create random steps, by pairing each observed step with 10 random steps
#' generated from fitted gamma and von Mises distributions. 
rs <- steps |> 
  random_steps()

#' 4. Next, we need to create a unique id for each individual animal x step_id
#' combination. This will be our stratum indicator when fitting mixed SSFs.
rs <- rs |> mutate(id_step_id_ = paste(id, step_id_, sep=":"))

#' Now check that everything is correct (we should have 11 points per stratum)
rs %>% 
  as.data.frame() |>
  dplyr::count(id_step_id_) |>
  head(5)

#' Get rid of two steps.  
steps_to_rm <- filter(rs, is.na(x2_)) |> pull(id_step_id_) |> unique()
rs <- rs |> filter(!id_step_id_ %in% steps_to_rm)

#' 5. Annotate data with covariates.  We will annotate using values at both 
#' the start and end of the steps.  And, time of day at start and end of the steps.
rs <- rs |> 
  extract_covariates(esa, where = "both") |> 
  extract_covariates(dist_built_up, where = "both") |> 
  extract_covariates(ter, where = "both") |> 
  time_of_day(where = "both")

# Two step approach ------------
#' 
#' Demonstrate a two-step approach to inference when analyzing data
#' from multiple animals. 
#'
#' **Research questions**: Do red deer select/avoid human infrastructure at
#'  different times of days? Are there sex-specific differences in habitat
#'  selection? 
#'  
rs
rs <- rs |> mutate(
  female = as.numeric(sex == "f"), 
  night = as.numeric(tod_end_ == "night"), 
  y = as.numeric(case_)
)

#' Restrict inference to summer months
dat.summer <- filter(rs, month(t1_) %in% 5:8)

#' ## Fit models to individuals
#' 
#' Our models/analysis will allow for the following features:
#' 
#' - selection for built up areas depends on time of day
#' - movement and habitat-selection parameters will vary by individual since we fit
#' models separately to each individual
tic()
m.ind <- dat.summer |> nest(data = -c(id, sex)) |> 
  mutate(mod = map(data, ~ clogit(
    y ~ dist_built_up_end + dist_built_up_end:night + 
      log(sl_) + sl_ + cos(ta_) +
      strata(id_step_id_), data = .x) |> 
      tidy()
  ))
toc()


m.ind$mod[[1]]

#' ## Summarize results, t-based confidence intervals
#' 
#' Now, summarize fits to individual animals. First, we will unnest mod
#' so that we can see the individual coefficients.  We will then drop the data
#' column.
ind.coefs <- m.ind |>
  unnest(mod) |> 
  dplyr::select(-c(data))   
ind.coefs

#' Now, we can calculate mean coefficients across animals and their SE and sample
#' size (number of individuals), which can then be used for t-based confidence intervals
se <- function(x) sd(x) / sqrt(length(x))

two.step.summary <- ind.coefs |> 
  group_by(term) |>   
  summarize(mean = mean(estimate), se = se(estimate), n = n(), naive.var = var(estimate)) |>
  mutate(t.uci = mean + qt(0.975, df = n-1)*se,
         t.lci = mean + qt(0.025, df = n-1)*se)  
two.step.summary

#' Alternatively, we could get confidence intervals using a bootstrap rather than
#' have to rely on the central limit theorem (or a Normality assumption) for 
#' inference.
ind.coefs.wide <- ind.coefs |> 
  dplyr::select(term, estimate, id) |> 
  pivot_wider(names_from=term, values_from=estimate)

#' ## Bootstrap Confidence intervals
#' 
#' We can then resample individuals with replacement, calculate the mean coeffient,
#' rinse, repeat, to estimate sampling variability for the mean coefficient.
boot_mean <- function(data, indices) {
  d <- data[indices, ]  # Resample rows, select column
  apply(d, 2, FUN=mean, na.rm = TRUE)
}

bootresults <- boot(data = ind.coefs.wide[,-1], statistic = boot_mean,
                    R = 9999)


#' Let's use bias-corrected and accelerated intervals 
ci_list <- lapply(1:5, function(i) {
  boot.ci(bootresults, type = "bca", index = i)
})

bca_intervals <- sapply(ci_list, function(ci) {
  c(b.lci = ci$bca[4], b.uci = ci$bca[5])
})
bca_intervals

#' Now, let's merge these onto our two.step.summary data set
colnames(bca_intervals) <- colnames(ind.coefs.wide[,-1])
bca_intervals <- as.data.frame(t(bca_intervals)) %>%
  mutate(term = colnames(bca_intervals)) %>%
  arrange(term) 
two.step.summary <-left_join(two.step.summary, bca_intervals)

#' Plot results for habitat selection parameters
two.step.summary |> 
  filter(term %in% c("dist_built_up_end", "dist_built_up_end:night")) |> 
  ggplot(aes(x = term, y = mean)) + 
  
  # 1. Reference line
  geom_hline(yintercept = 0, lty = 2, color = "gray50") +
  
  # 2. Individual coefficients (Kept gray and outside aes() so they stay out of the legend)
  geom_jitter(data = subset(ind.coefs, term %in% c("dist_built_up_end", "dist_built_up_end:night")),
              aes(x = term, y = estimate), 
              color = "gray40", alpha = 0.4, width = 0.15, height = 0) +
  
  # 3. Red confidence intervals (Color moved INSIDE aes)
  geom_pointrange(aes(ymin = t.lci, ymax = t.uci, color = "t-distribution"), 
                  size = 1, alpha = 0.9, 
                  position = position_nudge(x = -0.2)) +
  
  # 4. Blue confidence intervals (Color moved INSIDE aes)
  geom_pointrange(aes(ymin = b.lci, ymax = b.uci, color = "Bootstrap"), 
                  size = 1, alpha = 0.9, 
                  position = position_nudge(x = 0.2)) +
  
  # 5. Define the colors for the legend
  scale_color_manual(name = "Estimation Method", 
                     values = c("t-distribution" = "red", "Bootstrap" = "blue")) +
  
  theme_light() +
  theme(legend.position = "bottom") +
  xlab("") + ylab("Coefficient Estimates")

# Mixed model --------------

#' **Purpose:** Demonstrate a mixedSSA (step-selection functions) fit to data
#' from multiple animals using a hierarchical model.
#'
#' **Research questions**: Do red deer select/avoid human infrastructure at
#'  different times of days? Are there sex-specific differences in habitat
#'  selection? 

#' Make sure id_step_id is a factor
dat.summer$id_step_id_ <- as.factor(dat.summer$id_step_id_)


#' ## glmmTMB: random slopes model
#' 
#' We will use the following process to help ensure we fit the model correctly:
#' 
#' 1. Set up model, but do not fit it
#' 2. Set random intercept variance to a large fixed value, set other variance components
#'  to 0
#' 3. Fit the model with Poisson likelihood 
#' 
#' There are many models we could fit.  In this first one, we will assume the coefficients
#' vary independently from one another since we have only 18 individuals. With more 
#' individuals, we could estimate covariances between the random effects.  
#' The "(0 + x | id)" is used to let R know that we want to model random coefficients
#' for x and that these coefficients should vary independently of the intercept 
#' or other random effects included in the model.

tic("mixed ssf")
mixedSSA.tmp <- glmmTMB(
  y ~ -1 + dist_built_up_end + dist_built_up_end:night + 
    log(sl_) + sl_ + cos(ta_) +
    (1 | id_step_id_) + # Note, this should be our first random effect!
    (0 + sl_|id) + (0+log(sl_) | id) + (0 + cos(ta_) | id) +
    (0 + dist_built_up_end | id) + (0 + dist_built_up_end:night | id),
  family=poisson(), data = dat.summer,
  doFit=FALSE)

#' Set variance of random intercept to 10^6 and fit the model
mixedSSA.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(mixedSSA.tmp$parameters$theta)
mixedSSA.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
MixedSSA <- glmmTMB:::fitTMB(mixedSSA.tmp)
summary(MixedSSA)
toc()


# Comparing the results
#' 
#' Models fit using glmmTMB
mixedssa1.means <-data.frame(term = rownames(summary(MixedSSA)$coef$cond),
                             mean = summary(MixedSSA)$coef$cond[,1],
                             se = summary(MixedSSA)$coef$cond[,2],
                             method = "MixedSSA1")

#' Get approximate Normal based CI
mixedssa1.means <- mixedssa1.means |> 
  mutate(lci = mean - 1.96 * se, uci = mean + 1.96 * se)

#' Two step
two.step.means.boot <- two.step.summary |> 
  dplyr::select(term, mean, se, b.lci, b.uci) |>
  mutate(method = "Two-step Bootstrap") |>
  rename(uci = b.uci, lci = b.lci)

two.step.t <- two.step.summary |> 
  dplyr::select(term, mean, se, t.lci, t.uci) |>
  mutate(method = "Two-step t-dist") |>
  rename(uci = t.uci, lci = t.lci)

#' Combine all and plot
all_means <- bind_rows(
  mixedssa1.means, two.step.means.boot, two.step.t)


#+ fig.alt = "Estimates of mean parameters and CI", out.width="100%", fig.height = 9, fig.width =12
ggplot(data = all_means, aes(x = method, y = mean)) +
  geom_point() +
  geom_errorbar( aes(x = method, ymin = lci, ymax = uci),
                 width=0.2, linewidth=1) + 
  xlab("") + ylab(expression(hat(beta))) +
  facet_wrap(~term, scales="free") +
  # Rotate the x-axis text so the large font doesn't overlap
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Angle text and align it
  )

#' Variance parameters
#' 
#' Note, we should expect naive estimates, formed by taking the variance of the individual
#' coefficients, will be biased high as this variability will be impacted by both
#' true variability and sampling error.  We can see that naive estimates are bigger
#' than those from glmmTMB
var.two.step<- data.frame(term = two.step.summary$term,
                          sd.hat = sqrt(two.step.summary$naive.var), 
                          method = "Naive two-step")

var.mixed1 <- data.frame(term =  tidy(MixedSSA, effects="ran_pars")$term[-1],
                         sd = tidy(MixedSSA, effects="ran_pars")$estimate[-1]) |>
  mutate(method = "Mixed SSA1") 


#' Combine into a single data frame for plotting
names(var.two.step)[2]<-"sd"
allvars <- rbind(var.two.step[,-1], var.mixed1[c(3, 4, 5, 2, 1),-1]) 
allvars$term = rep(var.two.step[,1], 2)

#' Plot the estimated variance parameters. Here we see that, as expected,
#' var(individual fits) > var(mixed model).
ggplot(data = allvars, aes(x = method, y = sd)) +
  geom_point() +
  xlab("") + ylab("Standard deviation of individual-specific parameters)") +
  facet_wrap(~term, scales="free")+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Angle text and align it
  )


#' ### Individual Coefficients
#' 
#' Here, we see that individual coefficients are shrunk back toward the mean
#' coefficient in the mixed effect approach.
#' 
#' Let's start with calculating confidence intervals for individual fits,
#' appending asympototic normal-based CI for parameters in models fit to individuals)
# 1. Prep individual coefs (Force id to character)
ind.coefs.se <- ind.coefs |>  
  mutate(
    uci = estimate + 1.96 * std.error,
    lci = estimate - 1.96 * std.error,
    method = "Individual Fits",
    id = as.character(id)) |> 
  dplyr::select(id, term, estimate, method, uci, lci)

# 2. Prep mixed coefs (Force id to character and use NA_real_)
mixed_coefs <- coef(MixedSSA)$cond$id[,-1]  
mixed_coefs$id <- as.character(rownames(mixed_coefs)) # <-- Enforce character type

mixed_coefs <- mixed_coefs %>% 
  pivot_longer(
    cols = c("dist_built_up_end", "log(sl_)", "sl_", "cos(ta_)", "dist_built_up_end:night"), 
    names_to = "term", 
    values_to = "estimate"
  ) %>%
  mutate(
    method = "MixedSSA1",
    uci = NA_real_, # <-- Use numeric NA
    lci = NA_real_
  ) %>%
  # Ensure column order matches exactly
  dplyr::select(id, term, estimate, method, uci, lci)

# 3. Combine safely
allests <- bind_rows(mixed_coefs, ind.coefs.se)

# 4. Plot
cbp1 <- c("#999999", "#E69F00")

# Define a dodge width so the points sit side-by-side
pd <- position_dodge(width = 0.5)

ggplot(data = allests, aes(x = id, y = estimate, col = method)) +
  
  # Add dodging to points
  geom_point(size = 3.5, position = pd) +
  
  # Apply dodging to error bars, letting ggplot silently drop the Mixed NAs
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2, size = 1, position = pd) +
  geom_hline(aes(yintercept = mean, col = method), data = mixedssa1.means, linetype = 2, lwd =1.2) +
  xlab("") + 
  ylab(expression(hat(beta))) + 
  facet_wrap(~term, scales = "free") +
  scale_colour_manual(values = cbp1) +
  theme_light() +
  # Update the theme to remove text and move the legend
  theme(
    axis.text.x = element_blank(),   # Removes the ID numbers (labels)
    axis.ticks.x = element_blank(),  # Removes the little tick marks
    legend.position = "bottom"       # Moves the legend below the plot
  )
