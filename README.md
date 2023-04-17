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
- `twentyhooks` a binary variable indicating years where a different number of hooks were used. Useful if 2 fishing protocols are used and differences want to be controlled for. Set to 0 if not relevant.

## Simulation Study

The simulation study accompanying the paper can be found in the folder 'Simulation Study'

## Notes from Joe when Andy visited him in Sept 2022 and by Zoom in Oct 2022

Sean and Jillian -- some of this should be useful. I haven't had a chance to run anything more. I think for IPHC data it will use data as saved by my gfiphc package.
 
Joe pushed new `censored_index_fun_sdmTMB()` to hookCompetition, for operationally making indices.

Zoom chat 24/10/22 - IPHC data suggests fitting unstructured (no spatial) models with unique fixed effects for each year to get a coastwide index, plus spatiotemporal so can match spatial maps (see above function I think).

For generating new indices see Example.Rmd vignette, but for proper use (?) recommend using the functions that end in ...sdmTMB() to go into pacea, for example.

Use Censored_Longline_RCode (joe's GitHub) as submitted with manuscript, if have to re-make any figures. .rds files of results are saved. 

hookCompetition/Case Study is original code that then got tidied up into Censored_Longline_RCode/

Andy: ran hookCompetition/vignette/Example.Rmd at PSEC but was slow (think just due to connection): looks like code extracts new .rds files, but Andy ran from vignettes/ where it's saved them. They should match the ones sent by Joe in ../, maybe just move those to dummy folder so I have them, but I think they should be duplicated.
