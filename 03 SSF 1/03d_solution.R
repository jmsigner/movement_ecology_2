library(amt)

# 1. Load the `amt_fisher` dataset (it is shipeed with `amt` and you can get it with `data(amt_fisher)` or read it from the csv file `data/fisher.csv`). If you're working with your own data, format it as a `track_xyt` object using `amt::make_track()`.

data("amt_fisher")
amt_fisher_covar <- get_amt_fisher_covars()

fisher <- make_track(amt_fisher, x_, y_, t_, name = name)

# 2. Filter all data points for  `Leroy` and resample the data to 30 min. 

leroy <- filter(fisher, name == "Leroy")
summarize_sampling_rate(leroy)
leroy <- track_resample(leroy, rate = minutes(30), tolerance = minutes(2))

print(leroy, n = 100)
   
# 3. Use `Leroy` and create 15 random steps for each observed step. 

leroy <- leroy |> steps_by_burst() |> random_steps(n_control = 15)


# 4. Extract covariates at the end of each step for elevation. A raster with the elevation data is also shipped with `amt` and you can get the elevation by executing `data(amt_fisher_covar)`.  This returns a list with three covariates, we are only interested in elevation. To access elevation you can use: `amt_fisher$elevation`.
leroy <- leroy |> extract_covariates(amt_fisher_covar$elevation)

# 5. Fit a SSF an SSF, where you use elevation as the only covariate.
m1 <- leroy |> fit_ssf(case_ ~ elevation + strata(step_id_), model = TRUE)
summary(m1)

