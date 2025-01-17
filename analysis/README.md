This folder contains all `R` scripts used for analyses in this project (except for custom functions, which are in the `functions` folder, and those used to create the manuscript which are in the `writing` folder). Below is a list of the scripts, in the order in which they should be run to replicate the analysis.

* `00-compile-goat-tracking-data.R` compiles the goat tracking data into a single `RDS` file.
* `00-crop-resource-rasters.R` creates the resource rasters necessary for the project using rasters from http://www.earthenv.org/landcover.
* `00-download-bc-dem.R` downloads Digital Elevation Models (DEMs) necessary for the project.
* `00-preview-boreal-data.R` previews the boreal data obtained from https://www.bcogris.ca/projects/boreal-caribou-telemetry-data/.

## Estimate animals' speeds and space use

* `01-compile-tracking-data.R` merges all telemetry data into a single tibble, which it saves as an `RDS` file.
* `01a-compile-goat-calibration-data.R` compiles the goat calibration data into a single `CSV` file.
* `02-remove-outliers.R` screens the telemetry data for outliers and removes any problematic locations that the error model could not account for (given the DOP values and estimated error. Since I only have calibration data for the goats, I assume the error to be 10 m for all other telemetries).
* `03-movement-models'parallel.R` fits variograms, continuous-time movement models, and utilization distributions using multiple `R` sessions.

## Obtain estimated historical weather data

* `04-dowload-era5-temperature-data.R` downloads estimated historical hourly weather data from the European Centre for Medium-Range Weather Forecasts

## Annotate telelemetry data with hourly weather data

* `05-add-speeds-and-temperature.R`

## Estimate effects of weather on amimal movment

* `06-temperature-hgams.R` fits two Hierarchical Generalized Additive Models (HGAMs) that estimate the effects of temperature on: (1) the probability of an animal moving, and (2) the animal's speed, given that it is moving.
* `06-temperature-rsfs.R` fits Hierarchical Resource Selection Functions while accounting for changes in resoure selection with temperature.

## Make predictions for each species' region

* `07-download-bc-climateNA-data.R` downloads climate projection data using the `climatenaR` package.
* `08-bind-climateNA-data.R` binds the climate projection data 
* `09-model-era5-temperature-data.R` models the ERA5 temperature data to simulate weather given the monthly average temperature.

## Create figures

The `analysis/figures` folder contains scripts used to create figures with data from previous scripts.
