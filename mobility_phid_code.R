


# ====================== MOBILITY MODELS IN R: A SHORT TUTORIAL =============================

# This tutorial script focuses on using population estimates and geographical data sources to 
# predict relative population fluxes between pairs of administrative districts, and exploring ways of summarising these.

# The script will focus on implementing and exploring two types of general human mobility models:
# Gravity model: in this case, a naive (unparameterised) gravity model
# Radiation model: whose relative flux derivation is designed to be parameter-free

# We will estimate both sets of metrics across all second-level administrative units in the United Kingdom
# Which will hopefully give a good opportunity to see and critique some of the predictions that the models make
# (as a reminder that it is always worth keeping in mind the limitations of these kinds of models)

# The pipelines involves several stages:
# 1. Access administrative unit boundaries and populations for the country of interest
# 2. Develop a matrix of pairwise distances between each pair of administrative units (here, great circle distances)
# 3. Make and visualise predictions from naive gravity and radiation models
# 4. Visualise and summarise.



# ----------- Housekeeping and project setup -------------------

setwd("C:/Users/roryj/Documents/PhD/202007_lshtm_dengue/teaching/phid_mobility/")
library(raster); library(rgdal); library(maptools); library(sf); library(sp)
library(dplyr); library(exactextractr); library(ggplot2); library(geosphere)

# mapping theme
maptheme <- theme_classic() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        legend.title = element_text(size=10), 
        strip.background = element_blank(),
        strip.text = element_text(size=16))

# uk borders for mapping
uk_shp <- sf::st_read("./data/uk_shp/gadm36_GBR_0.shp")


# ------------ 1. Access input data for the United Kingdom ------------------

# Gridded population at 100m resolution for 2020
# Sourced from WorldPop
pop <- raster::raster("./data/uk_pop/gbr_ppp_2020_constrained.tif")

# UK admin level 2 shapefile, sourced from Office of National Statistics
shp <- sf::st_read("data/uk_shp/LPA_MAY_2021_UK_BFC_V2.shp")
shp$area_km2 <- as.vector(sf::st_area(shp)/10^6)
shp <- sf::st_transform(shp, crs=4326)

# Extract population for each admin unit...
exp <- exactextractr::exact_extract(pop, shp, fun="sum")

# and combine into a big project dataframe
dd <- shp %>%
  sf::st_drop_geometry() %>%
  dplyr::select(LPA21NM, LONG, LAT) %>%
  dplyr::rename("admin" = 1, "long"=2, "lat"=3) %>%
  dplyr::mutate(
    id = 1:nrow(shp),
    population = exp
  )

# Let's check everything has extracted properly by visualising population density over space - this looks fine.
# Although one thing we can see very quickly (as usual with adminstrative unit data) is that admin unit sizes are extremely variable
# This could have consequences for interpreting mobility flows between locations - this is always worth bearing in mind during spatial analyses.
shp %>%
  dplyr::mutate(pop = exp, popdens = (pop / area_km2) ) %>%
  ggplot() + 
  geom_sf(aes(fill=log(popdens)), col=NA) + 
  maptheme + 
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white")




# -------------- 2. Now we need to calculate pairwise distances between all administrative units -------------------

# First, let's build a dataframe of all combinations of source (area_i) and destination (area_j) locations
# Here, we are using the administrative unit centroid as our "location" from which to calculate distance
# Although remember that the differences in admin unit sizes mean that this may be more appropriate for some places (e.g. London boroughs)
# than others (e.g. the Yorkshire Dales). 
# One could also experiment with calculating distance using different metrics - for example, travel time between most populated locations 
# within each administrative unit.
# But now, for simplicity, we'll stick with lat-lon centroids.

# Build the dataframe (with n_districts^2 number of rows)
dd_i <- dd %>% rename("area_i"=admin, "id_i" = id, "x_i"=long, "y_i"=lat, "population_i"=population)
dd_j <- dd %>% rename("area_j"=admin, "id_j" = id, "x_j"=long, "y_j"=lat, "population_j"=population)
dists <- expand.grid(dd_i$area_i, dd_j$area_j) %>%
  dplyr::rename("area_i"=1, "area_j"=2) %>%
  left_join(dd_i) %>%
  left_join(dd_j)

# Use the haversine distance formula (accounting for Earth's curvature) to calculate the distance between the centroids in metres
# And then we'll divide that by 1000 to give the distance in kilometres
haversineCalc <- function(x){ geosphere::distHaversine(p1 = c(dists$x_i[x], dists$y_i[x]), p2 = c(dists$x_j[x], dists$y_j[x])) }
dists$dist_ij <- sapply(1:nrow(dists), haversineCalc) / 1000


# Some visualisation: 

# Let's take a look at the distribution of distances from any given location
dists %>% 
  dplyr::filter(area_i %in% sample(unique(area_i), 12, replace=FALSE)) %>%
  ggplot() + 
  geom_histogram(aes(x=dist_ij), col="black", fill="grey50") + 
  facet_wrap(~area_i) +
  theme_bw() + xlab("Distance (km)")

# And let's visualise how this looks geographically - we can see that the great circle distance is probably an extreme simplification
# in some contexts, particularly where there are coastlines or large geographical features involved...

dists %>%
  dplyr::filter(area_i == "Waltham Forest LPA") %>%
  ggplot() + 
  geom_sf(data=uk_shp, fill="grey85", col=NA) + 
  geom_point(aes(x_i, y_i)) + 
  geom_segment(aes(x=x_i, y=y_i, xend=x_j, yend=y_j, group=area_j, col=dist_ij), alpha=0.8, size=0.7) + 
  scale_color_viridis_c(option="F", direction=-1, end=0.9) + 
  maptheme

dists %>%
  dplyr::filter(area_i == "Yorkshire Dales National Park LPA") %>%
  ggplot() + 
  geom_sf(data=uk_shp, fill="grey85", col=NA) + 
  geom_point(aes(x_i, y_i)) + 
  geom_segment(aes(x=x_i, y=y_i, xend=x_j, yend=y_j, group=area_j, col=dist_ij), alpha=0.8, size=0.7) + 
  scale_color_viridis_c(option="F", direction=-1, end=0.9) + 
  maptheme




# ------------- 3. Use mobility models to predict pairwise relative fluxes: gravity model ------------------

# Here, we will use a naive gravity model that assumes all of the scaling parameters (alpha, beta and lambda) are 1
# i.e. Tij_grav ~ (pop_i * pop_j)/dist_ij
# This is a simplified version of a gravity model that consequently does not incorporate empirical information about
# the functional form of the distance decay in probability of movement, or about the weighting of source and destination attractiveness
# This kind of naive approach has been used successfully (/usefully) as a benchmark to test against the benefit of richer mobility data
# For example, Wesolowski et al, 2010: https://www.pnas.org/content/112/38/11887
# Although it tends to underperform against either mobile phone data or a parameterised gravity model, in general incorporating
# some mobility information into either statistical or dynamical models provides explanatory/predictive benefit over not including anything

# scaling coefficients are all set to 1 (i.e. naive)
g_alpha <- 1
g_beta <- 1
g_lambda <- 1

# gravity model formula, predict using population in thousands
gravityMod <- function(pop_i, pop_j, dist_ij, alpha, beta, lambda){
  t_ij = (pop_i^alpha * pop_j^beta) / (dist_ij^lambda)
  return(t_ij)
}
dists$Tij_grav <- gravityMod(pop_i=dists$population_i/1000, pop_j=dists$population_j/1000, dist_ij=dists$dist_ij, alpha=g_alpha, beta=g_beta, lambda=g_lambda)

# mobility models estimate fluxes *between* locations; self to self evaluates as Inf (because distance is 0) so set these examples to NA
dists$Tij_grav[ which(dists$area_i == dists$area_j) ] <- NA



# Now we can do some simple visualisation: here is the matrix of pairwise predicted fluxes
dists %>%
  ggplot() +
  geom_raster(aes(x = area_i, y=area_j, fill=Tij_grav)) + 
  scale_fill_gradientn(colors=viridisLite::turbo(200)) +
  theme(axis.text = element_blank())

# Here, we can visualise gravity from particular locations: can see that most of the predicted flux from the Yorkshire Dales
# Under gravity assumptions is towards nearby cities: Newcastle, Leeds, Sheffield
shp %>%
  dplyr::left_join(
    dists %>% dplyr::filter(area_i == "Yorkshire Dales National Park LPA") %>% dplyr::select(area_j, Tij_grav),
    by=c("LPA21NM"="area_j")
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=Tij_grav), color=NA) +
  maptheme + 
  scale_fill_viridis_c()

# Whereas the predicted flux from Waltham Forest is dominated by other London boroughs (although note the massively different scales)
shp %>%
  dplyr::left_join(
    dists %>% dplyr::filter(area_i == "Waltham Forest LPA") %>% dplyr::select(area_j, Tij_grav),
    by=c("LPA21NM"="area_j")
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=Tij_grav), color=NA) +
  maptheme + 
  scale_fill_viridis_c()




# ------------- 4. Use mobility models to predict pairwise relative fluxes: radiation model ------------------

# The radiation model was formalised by Simini et al 2012: https://www.nature.com/articles/nature10856
# It is based on the same assumptions as the gravity model but also accounts for the competitive effects of other population
# centres within an intervening radius.
# This is formulated based on an economic model of commuter flows that assumes that the number of jobs is proportional to the population
# of that location, and that commuters would choose the closest job that offers benefits better than those of their home location
# In practice, this means that predictions are not necessarily dominated by the most populated areas - other, more nearby places
# that offer opportunities, may be favoured by commuters.

# Because the radiation model requires estimating Sij (the total population within radius dist_ij around a focal location)
# Calculating this is rather more computationally intensive

# Firstly, create a dataframe containing areas, populations and distances, and use the data.table package for fast lookup
# By setting a lookup key...
distances <- dists[ , c("area_i", "area_j",  "dist_ij", "population_j")] 
distances <- data.table::data.table(distances)
distances$lookup <- distances$area_i
distances <- data.table::setkey(distances, lookup)

# Then we calculate S_ij (the population in radius dist_ij excluding source and destination populations).
# This is the memory intensive bit as it scales exponentially
# with the number of districts in your country (here, we need to do ~153,000 subsets and calculations)

# We run this function against each pair of locations
calc_Sij <- function(x){
  
  focal <- distances[ x ,]
  radius <- distances[ .(focal$lookup) ] %>% 
    dplyr::filter(dist_ij <= focal$dist_ij) %>% # get all admin units whose distance <= dist_ij
    dplyr::filter(!area_j %in% c(focal$area_i, focal$area_j)) # then exclude source and destination
  
  return(
    data.frame(
      area_i = focal$area_i,
      area_j = focal$area_j,
      Sij = sum(radius$population_j) # calculate total population
    )
  )
}

# Run the calculation - this will take quite a few minutes.
# radius_pop <- do.call(
#   rbind.data.frame,
#   lapply(1:nrow(distances), calc_Sij)
#)
#write.csv(radius_pop, "./output/radius_population.csv", row.names=FALSE)

# Here's one I prepared earlier
radius_pop <- read.csv("./output/radius_population.csv")

# Combine the S_ij population totals with the 
dists <- dists %>% dplyr::left_join(radius_pop)

# Then we predict radiation flux using the model as formulated in Simini et al 2012
# r_Ti = proportion persons in the source location that commute, here assume all persons (i.e. set as constant across all locations)
# pop_i = population of i
# pop_j = population of j
# s_ij = total population of radius between i and j

radiationMod <- function(r_Ti, pop_i, pop_j, s_ij){
  T_i <- r_Ti * pop_i
  T_ij <- T_i * ( (pop_i * pop_j) / ((pop_i + s_ij) * (pop_i + pop_j + s_ij)) )
  return(T_ij)
}
dists$Tij_rad <- radiationMod(r_Ti = 1, pop_i=dists$population_i, pop_j = dists$population_j, s_ij = dists$Sij)
dists$Tij_rad[ which(dists$area_i == dists$area_j) ] = NA # remove cases where area_i == area_i


# Now let's do the same visualisation of the fluxes matrix
dists %>%
  ggplot() +
  geom_raster(aes(x = area_i, y=area_j, fill=Tij_rad)) + 
  scale_fill_gradientn(colors=viridisLite::turbo(200)) +
  theme(axis.text = element_blank())

# And visualise what commuter flux looks like between locations: most of the outflow from Yorkshire Dales goes to much closer places
# than the gravity model would predict
shp %>%
  dplyr::left_join(
    dists %>% dplyr::filter(area_i == "Yorkshire Dales National Park LPA") %>% dplyr::select(area_j, Tij_rad),
    by=c("LPA21NM"="area_j")
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=Tij_rad), color=NA) +
  maptheme + 
  scale_fill_viridis_c()

# Whereas the predicted commuter outflow from Waltham Forest is *even more* dominated by other London boroughs
# Which is probably an accurate reflection of the reality that commuter flows tend to come intoward to London, rather than outwards
shp %>%
  dplyr::left_join(
    dists %>% dplyr::filter(area_i == "Waltham Forest LPA") %>% dplyr::select(area_j, Tij_rad),
    by=c("LPA21NM"="area_j")
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=Tij_rad), color=NA) +
  maptheme + 
  scale_fill_viridis_c()



# ---------------- 4. Finally, let's calculate some summary metrics ----------------

# We'll summarise by destination to calculate mean influx based on gravity and radiation estimates
fluxus <- dists %>%
  dplyr::select(area_j, population_j, Tij_grav, Tij_rad) %>%
  dplyr::group_by(area_j) %>%
  dplyr::summarise(population = head(population_j, 1), 
                   `Gravity influx` = mean(Tij_grav, na.rm=TRUE),
                   `Radiation influx` = mean(Tij_rad, na.rm=TRUE))

# Then we can visualise these over space - how do the two models differ in their predictions of commuter inflows?
# The major difference is that the radiation model predicts a much more even distribution of commuter influx across
# cities in the UK regardless of their size (and probably reflects regional economic hubs more accurately)
# In contrast, the gravity model's predictions are much more geographically restricted and focused on particularly densely populated cities
# That are also more "central" overall = particularly Birmingham, London, Leeds, Manchester, but very little in (for example)
# Edinburgh, Glasgow, Newcastle..

# Another interesting outcome of this analysis is noting how the size and shape of districts might influence predictions differently
# between models; for example just visually, the gravity model looks rather more sensitive to this (predicting, for example, larger
# inflows into areas that have high populations purely by virtue of their size)

shp %>%
  full_join(
    fluxus %>% reshape2::melt(id.vars = 1:2),
    by=c("LPA21NM"="area_j")
  ) %>%
  ggplot() + 
  geom_sf(aes(fill=value), col=NA) + 
  maptheme + 
  scale_fill_gradientn(colors=viridisLite::turbo(200), na.value="white") +
  facet_wrap(~variable)


# That's it! Feel free to ping me any questions: Rory.Gibb@lshtm.ac.uk



# ================================== ENDS =======================================

