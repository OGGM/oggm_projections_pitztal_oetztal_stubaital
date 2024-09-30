# OGGM Projections for Pitztal, Ötztal and Stubaital

[![DOI](https://zenodo.org/badge/791932856.svg)](https://zenodo.org/doi/10.5281/zenodo.11072662)

This repository contains OGGM projections for Pitztal, Ötztal and Stubaital. It is work in progress and will be updated regularly to accompany an upcoming publication. We don't expect any change to the results before submission.

We used [OGGM v1.6.1](https://zenodo.org/badge/latestdoi/43965645) calibrated it with high resolution, regionally aquired data (Hartl et al., in prep) for the 1980-2023 period. For the projections (2023-2100) we used all available CMIP6 climate projections, as described for the ['OGGM standard projections'](https://github.com/OGGM/oggm-standard-projections-csv-files). We aggregated these scenarios into global temperature increase levels +1.5°C, +2.0°C, +3.0°C and +4.0°C compared to preindustrial levels, similar to the methodology of [Rounce et al. (2023)](https://doi.org/10.1126/science.abo1324).

![Volume Evolution](_static/volume_evolution_and_elevation.png)
_Figure: Median (colored lines) and 5th and 95th percentile (shading) of the OGGM regional projections per future temperature scenario as percentage of 2017 glacier volume in the region of interest (a, b, c and d). The grey lines in each subplot show the median of the other three scenarios for reference. Temperature increase and number of climate model realizations n per scenario stated above the volume evolution plots. Insets (e, f, g and h) show the distribution of ice volume per 50 m elevation bands in different years for the four temperature scenarios, for the model run closest to the volume median of the scenario ensemble (colored lines). The shown models are +1.5°C INM-CM4-8 SSP1-2.6, +2.0°C NorESM2-MM SSP2-4.5, +3.0°C INM-CM5-0 SSP3-7.0 and +4.0°C NorESM2-MM SSP5-8.5._

## Data availability

The repository contains regional totals of glacier area and volume for each GCM realisation and temperature level. For glacier specific data, please contact us by opening an issue!

Further the script used for conducting the OGGM model runs is included.

