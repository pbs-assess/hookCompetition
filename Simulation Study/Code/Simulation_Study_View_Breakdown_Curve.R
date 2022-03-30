# Simulation study demonstrating the censored hook competition method
# In this simulation study, the target species is the most aggressive.
library(INLA)
library(ggplot2)
library(tidyverse)
library(mgcv)
seed <- 25032022 # date
set.seed(seed)
nspecies <- 3
bite_funs <- cbind(rep('constant',3),
                   c('exp_decay','mixed', 'constant'))
soak_time <- 5
n_hooks <- 800
nstation <- 6
nyears <- 1000
saturation_level <- c(1,2,3)
mean_attract <- c('constant', 'linear target', 'linear aggressive')
cov_matrices <- list( neg=
                        diag(c(0.8, 0.2, 0.2)) %*%
                        matrix(c(1,0,0,0,1,-0.6,0,-0.6,1), 3,3,byrow = T) %*%
                        diag(c(0.7, 0.2, 0.2)),
  none=
  diag(c(0.8, 0.2, 0.2)) %*%
  matrix(c(1,0,0,0,1,0,0,0,1), 3,3,byrow = T) %*%
  diag(c(0.7, 0.2, 0.2)),
  positive=
  diag(c(0.8, 0.2, 0.2)) %*%
    matrix(c(1,0,0,0,1,0.6,0,0.6,1), 3,3,byrow = T) %*%
    diag(c(0.7, 0.2, 0.2)))
n_sim <- 1
hook_sat_level <- 0.85#0.5 # true proportion at which saturation effects begin

# create saturation function
sat_fun <- function(sat_effect, sat_level, hook_sat_level=0.85)
{
  val <- rep(1, length(sat_level))

  val[which(sat_level>hook_sat_level)] <-
    (1 - sat_effect*((sat_level[which(sat_level>hook_sat_level)]-hook_sat_level)/(1-hook_sat_level)))
  return(val)
}

# try the hook competition adjusted method
comp_factor_fun <- function(prop_hook, n_hook)
{
  prop <- prop_hook
  # if all hooks saturated - map to 1 hook
  prop[which(prop == 0)] <- 1 / n_hook[which(prop == 0)]
  return(-log(prop)/(1-prop))
}

# sample the unadjusted bite times for each of the 'attracted' fish for each fishing event
bite_samp <- function(bite_fun, n)
{
  if(bite_fun=='exp_decay')
  {
    val <- rexp(n)
  }
  if(bite_fun=='constant')
  {
    val <- runif(n, min = 0, max=soak_time)
  }
  if(bite_fun=='mixed')
  {
    val <- c(rexp(n),
             runif(n, min = 0, max=soak_time))[
               rbinom(n=n, size = 1, prob=0.5)+1
             ]
  }
  return(val)
}

Results <- data.frame(
  #nsim = rep(1:n_sim, each=nstation*2*3*4*6),
  sat_level = rep(rep(c('low','high','very high'), times=n_sim),each=3*3*6*nstation*nyears*nspecies),
  mean_attract = rep(rep(mean_attract, times=n_sim*3), each=3*6*nstation*nyears*nspecies),
  correlation = rep(rep(c('negative','none','positive'),times=n_sim*3*3), each=6*nstation*nyears*nspecies),
  #sat_effects = rep(rep(c('no saturation','saturation'),times=n_sim*2*2), each=2*4*nstation),
  #bite_fun = rep(rep(c('constant','mixed'), times=n_sim*2*2*2), each=4*nstation),
  #model=rep(rep(c('naive','adjust','censored_upper85','censored_upper95','censored_upper100','censored'), times=n_sim*3*3*3), each=nstation*nyears*nspecies),

  Prop_Sat=rep(0, times=n_sim*3*3*3*6*nstation*nyears*nspecies),
  N_Bite=rep(0, times=n_sim*3*3*3*6*nstation*nyears*nspecies),
  Species=rep(0, times=n_sim*3*3*3*6*nstation*nyears*nspecies))

results_mapper <- function(i,j,k)
{
  return(
    which(Results$sat_level == c('low','high','very high')[i] &
            Results$mean_attract == mean_attract[j] &
            Results$correlation == c('negative','none','positive')[k])
         )
}

for(i in 1:length(saturation_level))
  {
    for(j in 1:length(mean_attract))
    {
      for(k in 1:length(cov_matrices))
      {
        # NOTE THAT WE NOW SCALE THE MEAN OF THE TARGET SPECIES BY THE EXPONENTIAL OF THE
        # LOG-NORMAL VARIANCE TO ENSURE THE SAME MEAN FROM THE PREVIOUS EXPERIMENT
        # DOESN'T AFFECT RELATIVE ABUNDANCE!!
        if(mean_attract[j] == 'constant')
        {
          mean_bite_gen =
            cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i],
                  rep(400,6),
                  (rep(100, 6)/exp(cov_matrices[[k]][3,3]/2)))
          mean_bite =
            cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i],
                  rep(400,6),
                  (rep(100, 6)))
        }
        if(mean_attract[j] == 'linear target')
        {
          mean_bite_gen =
            cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i],
                  rep(400,6),
                  (c(120, 140, 160, 180, 200, 220)-100)/exp(cov_matrices[[k]][3,3]/2))
          mean_bite =
            cbind(c(200, 200, 200, 200, 200, 200)*saturation_level[i],
                  rep(400,6),
                  (c(120, 140, 160, 180, 200, 220)-100))
        }
        if(mean_attract[j] == 'linear aggressive')
        {
          mean_bite_gen =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  rep(100, 6))
          mean_bite =
            cbind(c(120, 140, 160, 180, 200, 280)*saturation_level[i],
                  rep(400,6),
                  rep(100, 6))
        }
          # for(l in 1:dim(bite_funs)[2])
          # {
          saturation_effect <- rep(0,3)#c(0, 0.2, 0.8)
          # sample the number of each species that WOULD bite at each station for each year if hooks were available
          nbite <- data.frame(bites = rep(0, times=nspecies*nstation*nyears),
                              attracted = rep(0, times=nspecies*nstation*nyears),
                              species=rep(1:3, each=nstation*nyears),
                              station=rep(1:nstation, times=nspecies*nyears),
                              year=rep(rep(1:nyears, each=nstation), times=nspecies))
          nbite$attracted <- rpois(dim(nbite)[1],
                                   lambda = as.numeric(mean_bite_gen[cbind(nbite$station,nbite$species)])*
                                     exp(as.numeric(rmvn(n=nstation*nyears, mu = rep(0,nspecies), V=cov_matrices[[k]])[cbind(rep(1:(nstation*nyears),times=nspecies),nbite$species)])))
                                     #exp(rnorm(dim(nbite)[1], mean = 0, sd=sd_log_bite[nbite$species])))

          for(i2 in 1:nstation)
          {
            for(j2 in 1:nyears)
            {
              bite_time_1 <- bite_samp(bite_funs[1],sum(nbite$attracted[nbite$species==1 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_1)>soak_time)
              {
                bite_time_1[bite_time_1>soak_time] <-
                  bite_samp(bite_funs[1],sum(bite_time_1>soak_time))
              }
              # repeat for species 2 from uniform distribution
              bite_time_2 <- bite_samp(bite_funs[2],sum(nbite$attracted[nbite$species==2 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_2)>soak_time)
              {
                bite_time_2[bite_time_2>soak_time] <-
                  bite_samp(bite_funs[2],sum(bite_time_2>soak_time))
              }

              # repeat for species 3 from uniform distribution
              bite_time_3 <- bite_samp(bite_funs[3],sum(nbite$attracted[nbite$species==3 &
                                                                          nbite$station==i2 &
                                                                          nbite$year==j2]))
              # truncate them to 0-5 interval
              while(max(bite_time_3)>soak_time)
              {
                bite_time_3[bite_time_3>soak_time] <-
                  bite_samp(bite_funs[3],sum(bite_time_3>soak_time))
              }

              # Now we sample the first n_hooks*n_hook_sat_level unadjusted
              if((length(bite_time_1) + length(bite_time_2) + length(bite_time_3)) <= n_hooks*hook_sat_level)
              {
                nbite$bites[nbite$species==1 &
                              nbite$station==i2 &
                              nbite$year==j2] <- length(bite_time_1)
                nbite$bites[nbite$species==2 &
                              nbite$station==i2 &
                              nbite$year==j2] <- length(bite_time_2)
                nbite$bites[nbite$species==3 &
                              nbite$station==i2 &
                              nbite$year==j2] <- length(bite_time_3)
              }
              if((length(bite_time_1) + length(bite_time_2) + length(bite_time_3)) > n_hooks*hook_sat_level)
              {
                species_ind <- c(rep(1, length(bite_time_1)),rep(2,length(bite_time_2)),rep(3,length(bite_time_3)))
                # Now we sample the first n_hooks*n_hook_sat_level unadjusted
                all_times <- c(bite_time_1,bite_time_2,bite_time_3)
                time_ind <- sort.int(all_times, index.return = T, decreasing = F)$ix
                # for the remaining hooks we sample/thin according to the sat_fun
                current_sat_level <- hook_sat_level
                counter <- round(n_hooks*hook_sat_level)
                for(k2 in (round(n_hooks*hook_sat_level)+1):(length(all_times)))
                {
                  flag <- T
                  if(species_ind[time_ind[k2]]==1)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[1],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                    flag <- F
                  }
                  if(species_ind[time_ind[k2]]==2 & flag)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[2],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                    flag <- F
                  }
                  if(species_ind[time_ind[k2]]==3 & flag)
                  {
                    time_ind[k2] <- ifelse(rbinom(n=1,size=1,
                                                  prob=sat_fun(sat_effect = saturation_effect[3],
                                                               sat_level = current_sat_level,
                                                               hook_sat_level = hook_sat_level))==1,
                                           time_ind[k2], NA)
                  }
                  if(!is.na(time_ind[k2]))
                  {
                    counter <- counter + 1
                    current_sat_level <- counter/n_hooks
                  }
                  if(counter==n_hooks)
                  {
                    time_ind <- time_ind[1:k2]
                    break
                  }
                }
                time_ind <- time_ind[!is.na(time_ind)]
                nbite$bites[nbite$species==1 &
                              nbite$station==i2 &
                              nbite$year==j2] <- sum(species_ind[time_ind]==1)
                nbite$bites[nbite$species==2 &
                              nbite$station==i2 &
                              nbite$year==j2] <- sum(species_ind[time_ind]==2)
                nbite$bites[nbite$species==3 &
                              nbite$station==i2 &
                              nbite$year==j2] <- sum(species_ind[time_ind]==3)
              }

            }
          }

          nbite <-
            nbite %>%
            group_by(station, year) %>%
            mutate(prop_sat=sum(bites/n_hooks),
                   composition=bites/(sum(bites)))

          # ggplot(nbite, aes(x=station, y=bites, group=factor(species), colour=factor(species))) +
          #   geom_point() + geom_smooth() + geom_hline(yintercept = mean_bite[1,2:3]*exp(0.5*sd_log_bite[c(2,3)]^2))
          #
           # ggplot(nbite, aes(x=prop_sat, y=composition, group=factor(species), colour=factor(species))) +
           #   geom_point() + geom_smooth()
          #
          # ggplot(nbite[nbite$species!=2,], aes(x=prop_sat, y=composition, group=factor(species), colour=factor(species))) +
          #   geom_point() + geom_smooth()
          #
          # ggplot(nbite[nbite$species!=2,], aes(x=prop_sat, y=composition, group=factor(species), colour=factor(species))) +
          #   geom_point() + geom_smooth() + facet_wrap(~factor(station))

          Results[results_mapper(i,j,k),'Prop_Sat'] <- nbite$prop_sat
          Results[results_mapper(i,j,k),'N_Bite'] <- nbite$bites
          Results[results_mapper(i,j,k),'Composition'] <- nbite$composition
          Results[results_mapper(i,j,k),'Species'] <- nbite$species

      }
    }
}


Results %>% filter(Species!=2 & mean_attract=='constant') %>%
  mutate(sat_level=factor(sat_level, levels=c('low', 'high', 'very high'), ordered = T),
         Species=ifelse(Species==3, 'Target','Aggressive')) %>%
  ggplot(aes(x=Prop_Sat, y=N_Bite, group=factor(Species), colour=factor(Species))) +
  geom_smooth() + geom_vline(xintercept = hook_sat_level) + facet_grid(mean_attract+sat_level~correlation) +
  coord_cartesian(ylim=c(0,NA), xlim=c(NA,NA)) +
  ggtitle('The trend between ')

Results %>% filter(Species!=2 & mean_attract=='constant') %>%
  mutate(sat_level=factor(sat_level, levels=c('low', 'high', 'very high'), ordered = T)) %>%
  ggplot(aes(x=Prop_Sat, y=Composition, group=factor(Species), colour=factor(Species))) +
  geom_smooth() + geom_vline(xintercept = hook_sat_level) + facet_grid(mean_attract+sat_level~correlation)
