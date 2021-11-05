# hookCompetition
An R package for estimating relative abundance indices from longline data using a censored likelihood approach

## Installation

remotes::install_github('pbs-assess/hookCompetition', dep=T, build_vignettes = TRUE)

## Notes on how to use on non-IPHC data

To use the package on custom data, bypass the `read_Data_hookcomp()` step and ensure the following variable names are used inside the `sf` points object containing the catch counts:

- `N_it_species` (where species is the name of the species used (e.g. yelloweye)). This contains the catch counts of the species
- `effSkateIPHC` contains a measure of effort (e.g. number of hooks used etc.,)
- `prop_removed` contains the proportion of baits removed from the fishing event
- `year` the year of the fishing event
- `region_INLA` the region where the fishing event took place. An index (numeric variable) matching the corresponding row in `survey_boundaries`
- `station` containing the index of the fishing station where the fishing event took place. Useful if the same locations are visited each year (for multiple years)
- `twentyhooks` a bindary variable indicating years where a different number of hooks were used. Useful if 2 fishing protocols are used and differences want to be controlled for. Set to 0 if not relevant.
