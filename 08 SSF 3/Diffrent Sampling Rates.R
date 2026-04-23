# Irregular sampling rates ---
library(amt)
library(tidyverse)

data("amt_fisher")
covars <- get_amt_fisher_covars()

# Check the sampling rates of all indiviuals
summarize_sampling_rate_many(amt_fisher, "name")

# Using the dynamic+model approach for Lupe.
f1 <- filter(amt_fisher, name == "Lupe")
summarize_sampling_rate(f1)

dt <- round(difftime(lead(f1$t_), f1$t_, units = "min"))
table(dt)
# We see that most samples have a sampling rate of 2min

f1 <- f1 |> track_resample(rate = minutes(2), tolerance = minutes(5))
table(f1$burst_)

# Now lets create random steps and ensure the time difference in minutes
steps_bursted <- steps_by_burst(f1, time_diff_units = "minuntes") |>
  mutate(dt_ = as.numeric(dt_ / 60) |> round()) |>
  filter(dt_ %in% 2:7)

# We see that for different time differences, we have different step length distributions
steps_bursted %>%
  dplyr::select(dt_, ta_, sl_) %>%
  subset(!is.na(dt_) & !is.na(sl_) & !is.na(ta_)) %>%
  pivot_longer(ta_:sl_, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Value, color = as.factor(dt_))) +
  geom_density() +
  facet_wrap(~ Metric, scales = "free") +
  scale_color_viridis_d(begin = 0.3, name = "Step-Duration (dt)") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "gray95", color = "white"))

# Fit for each sampling rate a gamma and a von Mises distribution
distributions <- tibble(dt_ = 2:7)
distributions$Gamma <- lapply(distributions$dt_, function(x) {
  resampled <- track_resample(f1, rate = minutes(x), tolerance = seconds(20))
  resampled <- steps_by_burst(resampled)
  gamma     <- fit_distr(resampled$sl_, dist_name = "gamma")
  return(gamma)
})

distributions$VonMises <- lapply(distributions$dt_, function(x) {
  resampled <- track_resample(f1, rate = minutes(x), tolerance = seconds(20))
  resampled <- steps_by_burst(resampled)
  vonmises  <- fit_distr(resampled$ta_, dist_name = "vonmises")
  return(vonmises)
})

steps_bursted_nested <- steps_bursted %>%
  nest(Steps = -dt_) %>%
  mutate(dt_ = as.numeric(dt_))

# We now have each time difference a movement kernel
steps_bursted_nested
distributions


# Next we will add the sampling-rate-specific movement kernel to the steps of this
# sampling rate. 
steps_bursted_nested <- left_join(steps_bursted_nested, distributions, by = "dt_")
steps_bursted_nested

# Now we can iterate over each sampling rate
steps_bursted_nested <- steps_bursted_nested |> mutate(
  rs = pmap(list(Steps, Gamma, VonMises), ~ {
    try(random_steps(..1, n_control = 25,
                     sl_distr = ..2, ta_distr = ..3))
  })
)

# Not for sampling rate of 5 and 7 there were to few steps
steps_bursted_nested

# Lets remove these sampling rates
steps_bursted_nested <- steps_bursted_nested |> filter(!map_lgl(rs, is, "try-error"))

rs <- steps_bursted_nested |> dplyr::select(dt_, rs) |> unnest(cols = rs)
class(rs) <- class(steps_bursted_nested$rs[[1]])
rs

rs <- extract_covariates(rs, covars$elevation)

# Fit a model with sampling rate of 2 minutes. This is what we would have done
# if we would have just resampled the data to 2 minutes. 
rs |> filter(dt_ == 2) |>
  fit_clogit(case_ ~ elevation + sl_ + strata(step_id_)) |> summary()

rs |>
  fit_clogit(case_ ~ elevation + sl_ + sl_:dt_ + strata(step_id_)) |> summary()
