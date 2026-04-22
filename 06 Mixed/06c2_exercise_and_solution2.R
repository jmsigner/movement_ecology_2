library(tidyverse) 
library(survival) 
library(broom) 
library(broom.mixed) 
library(boot) 
library(infer)
library(multcomp)

#' In the models we fitted so fare we completely ignored differences between males and females. So, that seems like a nice next question to ask, do males and females differ?

#' Test if male and females differ in their selection for/against builtup areas during the day. 

#' Possible solutions ----------------------------------------------------------

# Using `ind.coefs` from the code walkthrough: 

# t-test: day
d1 <- ind.coefs |> filter(term == "dist_built_up_end")
t.test(estimate ~ sex, d1) # No difference

# bootstrap differences using the infer package
obs_diff <- d1 %>%
  specify(estimate ~ sex) %>%
  calculate(stat = "diff in means", order = c("m", "f"), na.rm. = TRUE)

obs_diff

boot_dist <- d1 %>%
  specify(estimate ~ sex) %>%
  generate(reps = 10000, type = "bootstrap") %>%
  calculate(stat = "diff in means", order = c("m", "f"), na.rm = TRUE)

get_confidence_interval(
  boot_dist,
  level = 0.95,
  type = "percentile"
)

boot_dist <- filter(boot_dist, !is.nan(stat))

get_p_value(
  boot_dist,
  obs_stat = obs_diff,
  direction = "two-sided"
)


visualize(boot_dist) +
  shade_p_value(obs_stat = obs_diff, direction = "two-sided")


# Use mixed ssfs with glmmTMB -------------------------

# We have to refit the model

tic()
mixedSSA.tmp2 <- glmmTMB(
  y ~ -1 + dist_built_up_end + dist_built_up_end:night +
    dist_built_up_end:female + dist_built_up_end:female:night +
    log(sl_) + sl_ + cos(ta_) +
    (1 | id_step_id_) + # Note, this should be our first random effect!
    (0 + sl_+ log(sl_) + cos(ta_) | id) +
    (0 + dist_built_up_end + dist_built_up_end:night | id),
  family=poisson(), data = dat.summer,
  doFit=FALSE)

#' Set variance of random intercept to 10^6 and fit the model
mixedSSA.tmp2$parameters$theta[1] <- log(1e3)
nvarparm<-length(mixedSSA.tmp2$parameters$theta)
mixedSSA.tmp2$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
MixedSSA2 <- glmmTMB:::fitTMB(mixedSSA.tmp2)
summary(MixedSSA2)
toc()

# Look at the difference between day an night
broom.mixed::tidy(MixedSSA2) |> 
  filter(term == "dist_built_up_end:female")

summary(MixedSSA2)

# We can use glht again to test for the differences
glht(
  MixedSSA2, 
  c(
    "dist_built_up_end = 0", # male day
    "dist_built_up_end + dist_built_up_end:female = 0", # female day
    "dist_built_up_end + dist_built_up_end:night = 0", # male night
    "dist_built_up_end + dist_built_up_end:female + dist_built_up_end:night + dist_built_up_end:night:female = 0" # female night
  )
) |> tidy(conf.int = TRUE) |> 
  ggplot(aes(contrast, estimate)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high))

# Now the differences

glht(
  MixedSSA2, 
  c(
    "dist_built_up_end:night = 0", # male
    "dist_built_up_end:night + dist_built_up_end:night:female = 0" # female
  )
) |> tidy(conf.int = TRUE) |> 
  mutate(sex = c("male", "female")) |> 
  ggplot(aes(sex, estimate)) + 
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, lty = 2)
