#' Loads the IPHC catch count data for a species. It returns a list containing
#' two entries:
#' \code{Reduced_data_sp} is a \code{SpatialPointsDataFrame} containing all the
#' set-level catch counts of the species specified (e.g. stored as a variable
#' name `N_it_yelloweye`).
#' \code{survey_boundaries} are the survey boundaries projected into the same
#' CRS as the \code{Reduced_data_sp} object.
#'
#' @param species_vec is a character vector of the one or more species desired for modelling
#' @param at_PBS is a logical flag determining if the user has access to the gf data. If FALSE, yelloweye rockfish must be the specified species.
#' @param simple_species_vec is a character vector containing the desired shortened names of the species
#' @param data_wd is a character vector stating the working directory where the rds files will be saved
#' @param min_dist is a numeric stating the minimum distance (km) required for a temporary IPHC station to be from a permanent IPHC station to be included
#' @param min_dist_to_boundary is a numeric stating the minimum distance (km) required for an IPHC station to be from the survey boundary
#' @param years_all_vec is a numeric vector stating the years all hooks were enumerated
#' @param years_20_vec is a numeric vector stating the years only 20 hooks were enumerated
#' @param survey_boundaries is a sf polygons object containing user-specified survey boundaries for indices to be computed in. If NULL QCS, WCHG, WCVI, HS used. The data.frame *must* contain a factor variable `Region` with levels equal to the names of the regions
#' @export
#' @importFrom magrittr %>%
#' @import stats
read_Data_hookcomp <- function(species_vec='yelloweye rockfish', at_PBS=F,
                               simple_species_vec='yelloweye',
                               data_wd = './',
                               min_dist = 50,
                               min_dist_to_boundary=0,
                               years_all_vec = c(1995,1996,2003:2012, 2014:2019),
                               years_20_vec = c(1997:2002, 2013, 2020),
                               survey_boundaries=NULL
                               )
{
  # species_vec is a character vector of the one or more species desired for modelling
  # at_PBS is a logical flag determining if the user has access to the gf data
  # simple_species_vec is a character vector containing the desired shortened names of the species
  # data_wd is a character vector stating the working directory where the rds files will be saved
  # min_dist is a numeric stating the minimum distance required for a temporary IPHC station to be from a permanent IPHC station to be included
  # min_dist_to_boundary is a numeric stating the minimum distance required for an IPHC station to be from the survey boundary
  # years_all_vec is a numeric vector stating the years all hooks were enumerated
  # years_20_vec is a numeric vector stating the years only 20 hooks were enumerated
  # survey_boundaries is a SpatialPolygonsDataFrame object containing user-specified survey boundaries for indices to be computed in. If NULL QCS, WCHG, WCVI, HS used. The data.frame must contain a factor variable Region with levels equal to the names of the regions

  if(length(species_vec) != length(simple_species_vec))
  {
    stop('The length of species_vec needs to match the length of simple_species_vec')
  }
  if(!is.null(survey_boundaries))
  {
    if(!(c('Region') %in% names(survey_boundaries)))
    {
      stop('The data.frame *must* contain a factor variable `Region` with levels equal to the names of the regions')
    }
  }

  setwd(data_wd)

  if(at_PBS)
  {
    # This chunk will only work within PBS
    for(sp in species_vec)
    {
      if(sum(grepl(list.files(),pattern =paste0(gsub("/", "-", gsub(" ", "-", sp)),'.rds')))==0)
      {
        gfiphc::cache_pbs_data_iphc(sp)
        # Creates 'sp'.rds. if file is not already present in local directory
      }
    }

    if(sum(grepl(list.files(),pattern =paste0(gsub("/", "-", gsub(" ", "-", "hook with bait")),'.rds')))==0)
    {
      gfiphc::cache_pbs_data_iphc("hook with bait")
      # creates hook-with-bait.rds if file not already in local directory
    }

    if(sum(grepl(list.files(),pattern ="sets-other-years.rds"))==0)
    {
      saveRDS(gfiphc::get_iphc_sets_info(), "sets-other-years.rds", compress = TRUE)
      # extracts the set level data from GFBio for 2003 onwards (excluding 2013) if file not in directory
    }


    #saveRDS(gfiphc::get_iphc_skates_info(), "skates-other-years.rds", compress = TRUE)
    # extracts skate level data from GFBio for 2003 onwards (excluding 2013)
  }
  # If not within PBS, the code below requires the above rds files offline
  data_species_vec <- rep(NA, length(species_vec))
  for(i in 1:length(species_vec))
  {
    data_species_vec[i] <- paste0(gsub(" ", "-", species_vec[i]), ".rds")
    #sp_set_counts <- readRDS(paste0(gsub(" ", "-", sp), ".rds"))
  }

  sets_1996_2002 <- gfiphc::data1996to2002
  sets_1996_2002$obsHooksPerSet <- sets_1996_2002$hooksObserved
  sets_1996_2002$long = sets_1996_2002$lon
  sets_1996_2002$station <- as.character(sets_1996_2002$station)

  # read in 1995 sets?
  # Not yet available

  var_names <- c("year","station","lat","lon","E_it","usable","standard","N_it","N_it20","N_itAll")

  hook_with_bait <- suppressWarnings(tryCatch(
    {
      readRDS("hook-with-bait.rds")$set_counts
    },
    error=function(cond)
    {
      return((hookCompetition::hook_with_bait)$set_counts)
    },
    finally={

    }
  ))

  hook_with_bait$N_itAll_hook <- hook_with_bait$N_it
  hook_with_bait$N_it20_hook <- hook_with_bait$N_it20
  hook_with_bait$N_it_hook <- hook_with_bait$N_it
  hook_with_bait$N_it_hook[hook_with_bait$year %in% years_20_vec] <- hook_with_bait$N_it20_hook[hook_with_bait$year %in% years_20_vec]

  sets <- suppressWarnings(tryCatch(
    {
      readRDS("sets-other-years.rds")
    },
    error=function(cond)
    {
      return(hookCompetition::sets_other_years)
    },
    finally={

    }
  ))

  # skates <- suppressWarnings(tryCatch(
  #   {
  #     readRDS("skates-other-years.rds")$set_counts
  #   },
  #   error=function(cond)
  #   {
  #     return(hookCompetition::`skates-other-years`)
  #   },
  #   finally={
  #
  #   }
  # ))

  # loop through the datasets compute the neccessary data
  tmp <- vector('list',length(data_species_vec))

  for(i in 1:length(data_species_vec))
  {
    tmp[[i]] <- suppressWarnings(tryCatch(
      {
        readRDS(data_species_vec[i])$set_counts
      },
      error=function(cond)
      {
        if(species_vec[i] == 'yelloweye rockfish' & at_PBS == F)
        {
          return((hookCompetition::yelloweye_rockfish)$set_counts)
        }
        if(species_vec[i] != 'yelloweye rockfish' | at_PBS == T)
        {
          stop(paste0('the rds file for species ',sp,' is not available.'))
        }
      },
      finally={

      }
    ))

    tmp[[i]][,paste0('N_itAll_',simple_species_vec[i])] <- tmp[[i]]$N_it
    tmp[[i]][,paste0('N_it20_',simple_species_vec[i])] <- tmp[[i]]$N_it20
    tmp[[i]][,paste0('N_it_',simple_species_vec[i])] <- tmp[[i]]$N_it
    tmp[[i]][which(tmp[[i]]$year %in% years_20_vec),
             paste0('N_it_',simple_species_vec[i])] <-
      tmp[[i]]$N_it20[which(tmp[[i]]$year %in% years_20_vec)]

    ## 2012 is the year a bait change experiment occurred
    ## 2013 is the year that only the first 20 hooks of each skate were enumerated
    tmp[[i]][,paste0('N_itAll_',simple_species_vec[i])][tmp[[i]]$year == 2012,1] <- NA
    tmp[[i]][,paste0('N_it20_',simple_species_vec[i])][tmp[[i]]$year == 2012,1] <- NA
    tmp[[i]][,paste0('N_it_',simple_species_vec[i])][tmp[[i]]$year == 2012,1] <- NA

  }

  # update variable names with species names
  var_names_species <- vector('list',length(simple_species_vec))

  var_names_species[[1]] <- var_names
  var_names_species[[1]][grepl('N_it',var_names_species[[1]])] <-
    paste0(var_names_species[[1]][grepl('N_it',var_names_species[[1]])],'_',
           simple_species_vec[1])

  All_data_sp <- tmp[[1]][,var_names_species[[1]]]
  #All_data_sp_skate <- tmp[[1]][,var_names_species[[1]]]

  if(length(data_species_vec) > 1)
  {
    for(i in 2:length(data_species_vec))
    {
      # create new variable names
      var_names_species[[i]] <- var_names
      var_names_species[[i]][grepl('N_it',var_names_species[[i]])] <-
        paste0(var_names_species[[i]][grepl('N_it',var_names_species[[i]])],'_',
               simple_species_vec[i])

      # keep only the relevant variables specified by var_names vector
      All_data_sp <-
        All_data_sp %>%
        dplyr::inner_join(tmp[[i]][,var_names_species[[i]]])

      # All_data_sp_skate <-
      #   All_data_sp_skate %>%
      #   inner_join(tmp[[i]][,var_names_species[[i]]])

    }
  }
  #browser()

  All_data_sp <-
    All_data_sp %>%
    dplyr::inner_join(hook_with_bait[,c(1,2,3,4,5,8,11,12,13,14,15)]) %>%
    dplyr::left_join(sets[,c('year','setID','tripID','station','long','lat','obsHooksPerSet','deplHooksPerSet','effSkateIPHC','skatesCount','usable','iphcUsabilityDesc','standard','setInTrip')]) %>%
    dplyr::left_join(unique(sets_1996_2002[,c('year','station','obsHooksPerSet','usable')]),
              by=c('year','station'))

  # Merge the new usable.x and usable.y and obsHooksPerSet.x and obsHooksPerSet.y
  All_data_sp$obsHooksPerSet <- All_data_sp$obsHooksPerSet.y
  All_data_sp$obsHooksPerSet[is.na(All_data_sp$obsHooksPerSet)] <-
    All_data_sp$obsHooksPerSet.x[is.na(All_data_sp$obsHooksPerSet)]
  All_data_sp$usable <- All_data_sp$usable.y
  All_data_sp$usable[is.na(All_data_sp$usable)] <-
    All_data_sp$usable.x[is.na(All_data_sp$usable)]
  All_data_sp <-
    All_data_sp[,-c(which(names(All_data_sp) %in%
                            c('obsHooksPerSet.y','obsHooksPerSet.x',
                              'usable.y','usable.x')))]

  # All_data_sp <- sp::SpatialPointsDataFrame(
  #   coords = cbind(All_data_sp$lon,All_data_sp$lat),
  #   data = All_data_sp,
  #   proj4string = sp::CRS(SRS_string = 'EPSG:4326'))

  All_data_sp <- sf::st_as_sf(
    All_data_sp,
    coords = c("lon","lat"),
    crs = 4326)

  # Update the effective skate to match the number of hooks observed
  All_data_sp$effSkateIPHC[which(All_data_sp$year %in% years_all_vec)] <-
    All_data_sp$E_it[which(All_data_sp$year %in% years_all_vec)]
  All_data_sp$effSkateIPHC[which(All_data_sp$year %in% years_20_vec)] <-
    All_data_sp$E_it20[which(All_data_sp$year %in% years_20_vec)]

  if(is.null(survey_boundaries))
  {
    survey_boundaries <-
      sf::st_as_sf(
      hookCompetition::survey_boundaries
      )

    # survey_boundaries <- sp::SpatialPolygons(
    #   Srl=list(sp::Polygons(list(sp::Polygon(coords = survey_boundaries$HS)),ID=c('HS')),
    #            sp::Polygons(list(sp::Polygon(coords = survey_boundaries$QCS)),ID=c('QCS')),
    #            sp::Polygons(list(sp::Polygon(coords = survey_boundaries$WCHG)),ID=c('WCHG')),
    #            sp::Polygons(list(sp::Polygon(coords = survey_boundaries$WCVI)),ID=c('WCVI'))),
    #   proj4string = sp::CRS(SRS_string = 'EPSG:4326'))
    #
    # survey_df <- data.frame(Region = factor(names(survey_boundaries)))
    # rownames(survey_df) <- names(survey_boundaries)
    #
    # survey_boundaries <- sp::SpatialPolygonsDataFrame(
    #   survey_boundaries,
    #   data=survey_df)
  }

  # # Convert all objects to a 'Canadian CRS'
  # M_CRS <- CRS(SRS_string = 'EPSG:3005')
  #
  # survey_boundaries <- spTransform(survey_boundaries,M_CRS)
  # All_data_sp <- spTransform(All_data_sp,M_CRS)

  # Change to a new 'Canadian' CRS in units of km
#   KM_CRS <- sp::CRS('+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0
# +datum=NAD83 +units=km +no_defs')

  KM_CRS <- sf::st_crs('+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0
+datum=NAD83 +units=km +no_defs')

  # survey_boundaries <- sp::spTransform(survey_boundaries,KM_CRS)
  survey_boundaries <- sf::st_transform(survey_boundaries,KM_CRS)

  # All_data_sp <- sp::spTransform(All_data_sp,KM_CRS)
  All_data_sp <- sf::st_transform(All_data_sp,KM_CRS)

  #All_data_sp_skate <- spTransform(All_data_sp_skate,KM_CRS)

  # Where are all the stations that only record in a single year that lie within
  # min_dist km of the nearest permanent station?
  ind_multiyear <- All_data_sp$station %in% names(table(All_data_sp$station)[which(table(All_data_sp$station)>1)])

  ind_singleyear <- All_data_sp$station %in% names(table(All_data_sp$station)[which(table(All_data_sp$station)==1)])

  # ind_close <- apply(rgeos::gDistance(All_data_sp,
  #                              All_data_sp[which(ind_multiyear),],
  #                              byid = T), 2,
  #                    FUN = function(x){min(x)<min_dist})

  ind_close <- apply(sf::st_distance(All_data_sp,
                                      All_data_sp[which(ind_multiyear),],
                                      by_element = F), 1,
                     FUN = function(x){min(x)<min_dist})

  All_data_sp <- All_data_sp[which(ind_close),]

  # Create a variable for proportion of bait removed
  All_data_sp$prop_removed <-
    (All_data_sp$obsHooksPerSet-All_data_sp$N_it_hook)/All_data_sp$obsHooksPerSet

  # All_data_sp$region_INLA <-
  #   apply(rgeos::gDistance(survey_boundaries, All_data_sp, byid=T),
  #         MARGIN=1, FUN=function(x){
  #           which.min(x)})
  All_data_sp$region_INLA <-
    apply(sf::st_distance(survey_boundaries, All_data_sp, by_element = F),
          MARGIN=2, FUN=function(x){
            which.min(x)})

  # ind_close2 <- apply(rgeos::gDistance(All_data_sp,
  #                               survey_boundaries,
  #                               byid = T), 2,
  #                     FUN = function(x){min(x)>min_dist_to_boundary})
  ind_close2 <- apply(sf::st_distance(All_data_sp,
                                       survey_boundaries,
                                       by_element = F), 1,
                      FUN = function(x){min(x)>min_dist_to_boundary})

  All_data_sp$region_INLA[ind_close2] <- NA

  All_data_sp$twentyhooks=ifelse(
    is.na(as.numeric(as.matrix(All_data_sp[,paste0('N_itAll_',simple_species_vec[1])])[,1])),
    1,0)

  Reduced_data_sp <-
    All_data_sp[which(All_data_sp$station %in% (names(table(All_data_sp$station))[table(All_data_sp$station)>1])),]

  Reduced_data_sp <- Reduced_data_sp[which(!(Reduced_data_sp$year %in% c(2012))),]

  # NOTE THAT 2013 USED ONLY FIRST 20 HOOKS BUT THE OBSERVED NUMBER OF HOOKS IS NA.
  # MAP 1997'S EFFECTIVE SKATE VALUES TO 2013 VALUES TO OBTAIN THE NHOOKS.
  Reduced_data_sp[which(Reduced_data_sp$year==2013),]$obsHooksPerSet
  mean_effskate_perhook_1997 <-
    by(Reduced_data_sp[which(Reduced_data_sp$year==1997),]$effSkateIPHC,
       INDICES = as.factor(
         Reduced_data_sp[which(Reduced_data_sp$year==1997),]$obsHooksPerSet),
       FUN=mean
    )
  # map 2013 values of eff skate to the closest matching n hooks
  Reduced_data_sp[which(Reduced_data_sp$year==2013),]$obsHooksPerSet <-
    c(as.numeric(names(mean_effskate_perhook_1997))[
      apply(
        outer(Reduced_data_sp[which(Reduced_data_sp$year==2013),]$effSkateIPHC,mean_effskate_perhook_1997,'-')^2,
        1, which.min
      )
    ])
  # update prop_removed variable
  Reduced_data_sp[which(Reduced_data_sp$year==2013),]$prop_removed <-
    (Reduced_data_sp[which(Reduced_data_sp$year==2013),]$obsHooksPerSet -
       Reduced_data_sp[which(Reduced_data_sp$year==2013),]$N_it_hook) /
    Reduced_data_sp[which(Reduced_data_sp$year==2013),]$obsHooksPerSet

  return(list(
    Reduced_data_sp = Reduced_data_sp,
    #All_data_sp_skate=All_data_sp_skate,
    survey_boundaries=survey_boundaries
  ))

}
