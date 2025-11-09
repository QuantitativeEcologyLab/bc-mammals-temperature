This folder contains all custom functions created in this project. The outlier cleaning functions assume the telemetry data is in a data frame with a row for each individual that contains a nested data frame called `tel` with the telemetry data. The other columns indicate the animal ID (`animal`), the species (`species`), and the dataset name (`dataset_name`). The data should therefore look something like this:

```
> d
# A tibble: 434 × 4
   animal species     dataset_name       tel                   
   <chr>  <chr>       <chr>              <list>                
 1 BW003  Canis lupus Canis_lupus_boreal <tibble [15 × 16]>    
 2 BW006  Canis lupus Canis_lupus_boreal <tibble [61 × 16]>    
 3 BW008  Canis lupus Canis_lupus_boreal <tibble [11,846 × 16]>
...
```

The functions are:

- `check_animal.R`: runs a series of checks for detecting outlier locations. See appendix B in the `writing` folder.
- `detrend_speeds.R`: detrends the speeds to make uninformed speeds closer to 0 and reduce the mean-reversion tendency of `ctmm` speed estimates. it also reduces the degrees of freedom and widens the CIs by removing the dof for the mean
- `find_angle.R`: finds an angle in radians or degrees between three points, with coordinates given as two vectors of `x` and `y`
- `flag_outlier.R`: flags outliers by changing the value of the `outlier` column in the dataset based on some user-input cutoffs
- `get_legend.R`: a working version of `cowplot::get_legend()`
- `import_rda.R`: a function to import an Rda file without creating an object in the environment
- `labeller_perc.R`: a labelling function to convert axis breaks to percentages 
- `outlier_plots.R`: a funtion called by `check_animal()` that creates the diagnostic plots for checking telemetry data
- `plot_adj.R`: a function to plot `n_adj` points before and after any location that matches the selection criteria
- `remove_outlier_flags.R`: a function to clear the flags in the `outlier` column by resetting it to the values in the  `original_outliers` column. Since I only used this function to start over with an animal, it's not included in any of the scripts.
- `seq_range.R`: a wrapper for `seq(from = min(x), to = max(x), length.out = length(x))`, given a vector x
- `test_speeds.R`: a function to quickly fit a movement model for a given animal and check its mean speed. not used in the scripts.
