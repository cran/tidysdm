## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# set options to limit threads used by imported libraries
options(cores=2)
options(mc.cores=2)
# xgboost uses data.table
data.table::setDTthreads(2)

## ----libraries----------------------------------------------------------------
library(tidysdm)

## ----load_presences-----------------------------------------------------------
 data(lacerta)
 lacerta

## ----download_presences, eval = FALSE-----------------------------------------
#  # download presences
#  library(rgbif)
#  occ_download_get(key = "0068808-230530130749713", path = tempdir())
#  # read file
#  library(readr)
#  distrib <- read_delim(file.path(tempdir(), "0068808-230530130749713.zip"))
#  # keep the necessary columns and rename them
#  lacerta <- distrib %>% select(gbifID, decimalLatitude, decimalLongitude) %>%
#    rename(ID = gbifID, latitude = decimalLatitude, longitude = decimalLongitude)

## ----cast_to_sf---------------------------------------------------------------
library(sf)
lacerta <- st_as_sf(lacerta, coords = c("longitude", "latitude"))
st_crs(lacerta) <- 4326

## ----land_mask, eval=FALSE----------------------------------------------------
#  library(pastclim)
#  download_dataset(dataset = "WorldClim_2.1_10m")
#  land_mask <-
#    get_land_mask(time_ce = 1985, dataset = "WorldClim_2.1_10m")
#  
#  # Iberia peninsula extension
#  iberia_poly <-
#    terra::vect(
#      "POLYGON((-9.8 43.3,-7.8 44.1,-2.0 43.7,3.6 42.5,3.8 41.5,1.3 40.8,0.3 39.5,
#       0.9 38.6,-0.4 37.5,-1.6 36.7,-2.3 36.3,-4.1 36.4,-4.5 36.4,-5.0 36.1,
#      -5.6 36.0,-6.3 36.0,-7.1 36.9,-9.5 36.6,-9.4 38.0,-10.6 38.9,-9.5 40.8,
#      -9.8 43.3))"
#    )
#  
#  crs(iberia_poly) <- "lonlat"
#  # crop the extent
#  land_mask <- crop(land_mask, iberia_poly)
#  # and mask to the polygon
#  land_mask <- mask(land_mask, iberia_poly)

## ----echo=FALSE, results='hide'-----------------------------------------------
library(pastclim)
set_data_path(on_CRAN = TRUE)
# Iberia peninsula extension
iberia_poly <- terra::vect("POLYGON((-9.8 43.3,-7.8 44.1,-2.0 43.7,3.6 42.5,3.8 41.5,1.3 40.8,0.3 39.5,0.9 38.6,-0.4 37.5,-1.6 36.7,-2.3 36.3,-4.1 36.4,-4.5 36.4,-5.0 36.1,-5.6 36.0,-6.3 36.0,-7.1 36.9,-9.5 36.6,-9.4 38.0,-10.6 38.9,-9.5 40.8,-9.8 43.3))")

crs(iberia_poly) <- "lonlat"
gdal(warn = 3)
land_mask <- rast(system.file("extdata/lacerta_land_mask.nc",
  package = "tidysdm"
))

## ----fig.width=6, fig.height=4------------------------------------------------
library(tidyterra)
library(ggplot2)
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_1985)) +
  geom_sf(data = lacerta) + scale_fill_gradient(na.value = "transparent")


## ----thin_by_cell-------------------------------------------------------------
set.seed(1234567)
lacerta <- thin_by_cell(lacerta, raster = land_mask)
nrow(lacerta)

## ----plot_thin_by_cell, fig.width=6, fig.height=4-----------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_1985)) +
  geom_sf(data = lacerta) + scale_fill_gradient(na.value = "transparent")

## ----thin_by_dist-------------------------------------------------------------
set.seed(1234567)
lacerta_thin <- thin_by_dist(lacerta, dist_min = km2m(20))
nrow(lacerta_thin)

## ----plot_thin_by_dist, fig.width=6, fig.height=4-----------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_1985)) +
  geom_sf(data = lacerta_thin) + scale_fill_gradient(na.value = "transparent")

## ----download_background, eval=FALSE------------------------------------------
#  # download presences
#  library(rgbif)
#  # download file
#  occ_download_get(key = "0121761-240321170329656", path = tempdir())
#  # read file
#  library(readr)
#  backg_distrib <- readr::read_delim(file.path(tempdir(), "0121761-240321170329656.zip"))
#  
#  # keep the necessary columns
#  lacertidae_background <- backg_distrib %>% select(gbifID, decimalLatitude, decimalLongitude) %>%
#    rename(ID = gbifID, latitude = decimalLatitude, longitude = decimalLongitude)
#  
#  lacertidae_background <- st_as_sf(lacertidae_background, coords = c("longitude", "latitude"))
#  st_crs(lacertidae_background) <- 4326

## ----echo=FALSE---------------------------------------------------------------
data("lacertidae_background")
lacertidae_background <- st_as_sf(lacertidae_background, coords = c("longitude", "latitude"))
st_crs(lacertidae_background) <- 4326

## ----background_to_raster-----------------------------------------------------
lacertidae_background_raster <- rasterize(lacertidae_background, land_mask, fun = "count")

plot(lacertidae_background_raster)


## ----sample_background--------------------------------------------------------
set.seed(1234567)
lacerta_thin <- sample_background(data = lacerta_thin, raster = lacertidae_background_raster,
                  n = 3 * nrow(lacerta_thin),
                  method = "bias",
                  class_label = "background",
                  return_pres = TRUE)


## ----plot_sample_pseudoabs, fig.width=6, fig.height=4-------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_1985)) +
  geom_sf(data = lacerta_thin, aes(col = class)) + scale_fill_gradient(na.value = "transparent")

## ----load_climate, eval=FALSE-------------------------------------------------
#  climate_vars <- get_vars_for_dataset("WorldClim_2.1_10m")

## ----echo=FALSE, results='hide'-----------------------------------------------
gdal(warn = 3)
climate_present <- rast(system.file("extdata/lacerta_climate_present_10m.nc",
  package = "tidysdm"
))
climate_vars <- names(climate_present)

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset("WorldClim_2.1_10m")

## ----climate_worldclim, eval=FALSE--------------------------------------------
#  climate_present <- pastclim::region_slice(
#    time_ce = 1985,
#    bio_variables = climate_vars,
#    data = "WorldClim_2.1_10m",
#    crop = iberia_poly
#  )

## -----------------------------------------------------------------------------
lacerta_thin <- lacerta_thin %>%
  bind_cols(terra::extract(climate_present, lacerta_thin, ID = FALSE))

## ----fig.height=11, fig.width=7-----------------------------------------------
lacerta_thin %>% plot_pres_vs_bg(class)

## -----------------------------------------------------------------------------
lacerta_thin %>% dist_pres_vs_bg(class)

## ----climate_variables--------------------------------------------------------

suggested_vars <- c("bio06", "bio05", "bio13", "bio14", "bio15")

## ----fig.width=7, fig.height=8------------------------------------------------
pairs(climate_present[[suggested_vars]])


## ----choose_var_cor_keep------------------------------------------------------
climate_present <- climate_present[[suggested_vars]]

vars_uncor <- filter_collinear(climate_present, cutoff = 0.7, method = "cor_caret")
vars_uncor

## -----------------------------------------------------------------------------
lacerta_thin <- lacerta_thin %>% select(all_of(c(vars_uncor, "class")))
climate_present <- climate_present[[vars_uncor]]
names(climate_present) # added to highlight which variables are retained in the end

## ----recipe-------------------------------------------------------------------
lacerta_rec <- recipe(lacerta_thin, formula = class ~ .)
lacerta_rec

## -----------------------------------------------------------------------------
lacerta_thin %>% check_sdm_presence(class)

## ----workflow_set-------------------------------------------------------------
lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf(),
      # boosted tree model (gbm) specs with tuning
      gbm = sdm_spec_boost_tree(),
      # maxent specs with tuning
      maxent = sdm_spec_maxent()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

## ----training_cv, fig.width=6, fig.height=4-----------------------------------
library(tidysdm)
set.seed(100)
#lacerta_cv <- spatial_block_cv(lacerta_thin, v = 5)

lacerta_cv <- spatial_block_cv(data = lacerta_thin, v = 3, n = 5)
autoplot(lacerta_cv)

## ----tune_grid----------------------------------------------------------------
set.seed(1234567)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
    resamples = lacerta_cv, grid = 3,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## ----autoplot_models, fig.width=7, fig.height=4-------------------------------
autoplot(lacerta_models)

## -----------------------------------------------------------------------------
lacerta_ensemble <- simple_ensemble() %>%
  add_member(lacerta_models, metric = "boyce_cont")
lacerta_ensemble


## ----autoplot_ens, fig.width=7, fig.height=4----------------------------------
autoplot(lacerta_ensemble)

## -----------------------------------------------------------------------------
lacerta_ensemble %>% collect_metrics()

## ----plot_present, fig.width=6, fig.height=4----------------------------------
prediction_present <- predict_raster(lacerta_ensemble, climate_present)
ggplot() +
  geom_spatraster(data = prediction_present, aes(fill = mean)) +
  scale_fill_terrain_c() +
  # plot presences used in the model
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

## ----plot_present_best, fig.width=6, fig.height=4-----------------------------

prediction_present_boyce <- predict_raster(lacerta_ensemble, climate_present,
  metric_thresh = c("boyce_cont", 0.7),
  fun = "median"
)
ggplot() +
  geom_spatraster(data = prediction_present_boyce, aes(fill = median)) +
  scale_fill_terrain_c() +
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

## -----------------------------------------------------------------------------
lacerta_ensemble <- calib_class_thresh(lacerta_ensemble,
  class_thresh = "tss_max", 
  metric_thresh = c("boyce_cont", 0.7)
)

## ----fig.width=6, fig.height=4------------------------------------------------
prediction_present_binary <- predict_raster(lacerta_ensemble,
  climate_present,
  type = "class",
  class_thresh = c("tss_max"), 
  metric_thresh = c("boyce_cont", 0.7)
)
ggplot() +
  geom_spatraster(data = prediction_present_binary, aes(fill = binary_mean)) +
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset("WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m")

## ----eval=FALSE---------------------------------------------------------------
#  get_time_ce_steps("WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m")

## ----echo=FALSE---------------------------------------------------------------
c(2030, 2050, 2070, 2090)

## ----eval=FALSE---------------------------------------------------------------
#  get_vars_for_dataset("WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m")

## ----echo=FALSE---------------------------------------------------------------
climate_vars[-length(climate_vars)]

## ----eval=FALSE---------------------------------------------------------------
#  climate_future <- pastclim::region_slice(
#    time_ce = 2090,
#    bio_variables = vars_uncor,
#    data = "WorldClim_2.1_HadGEM3-GC31-LL_ssp245_10m",
#    crop = iberia_poly
#  )

## ----echo=FALSE, results='asis'-----------------------------------------------
gdal(warn = 3)
climate_future <- rast(system.file("extdata/lacerta_climate_future_10m.nc",
  package = "tidysdm"
))

## ----plot_future, fig.width=6, fig.height=4-----------------------------------
prediction_future <- predict_raster(lacerta_ensemble, climate_future)

ggplot() +
  geom_spatraster(data = prediction_future, aes(fill = mean)) +
  scale_fill_terrain_c()

## -----------------------------------------------------------------------------
climate_future_clamped <- clamp_predictors(climate_future, 
                                           training = lacerta_thin,
                                           .col= class)
prediction_future_clamped <- predict_raster(lacerta_ensemble, 
                                            raster = climate_future_clamped)

ggplot() +
  geom_spatraster(data = prediction_future_clamped, aes(fill = mean)) +
  scale_fill_terrain_c()

## -----------------------------------------------------------------------------
lacerta_mess_future <- extrapol_mess(x = climate_future, 
                                      training = lacerta_thin, 
                                      .col = "class")

ggplot() + geom_spatraster(data = lacerta_mess_future) + 
  scale_fill_viridis_b(na.value = "transparent")

## -----------------------------------------------------------------------------
# subset mess 
lacerta_mess_future_subset <- lacerta_mess_future
lacerta_mess_future_subset[lacerta_mess_future_subset >= 0] <- NA
lacerta_mess_future_subset[lacerta_mess_future_subset < 0] <- 1

# convert into polygon
lacerta_mess_future_subset <- as.polygons(lacerta_mess_future_subset)

# plot as a mask 
ggplot() + geom_spatraster(data = prediction_future) + 
  scale_fill_viridis_b(na.value = "transparent") + geom_sf(data = lacerta_mess_future_subset, fill= "lightgray", alpha = 0.5, linewidth = 0.5)


## -----------------------------------------------------------------------------
bio05_prof <- lacerta_rec %>%
  step_profile(-bio05, profile = vars(bio05)) %>%
  prep(training = lacerta_thin)

bio05_data <- bake(bio05_prof, new_data = NULL)

bio05_data <- bio05_data %>%
  mutate(
    pred = predict(lacerta_ensemble, bio05_data)$mean
  )

ggplot(bio05_data, aes(x = bio05, y = pred)) +
  geom_point(alpha = .5, cex = 1)

## -----------------------------------------------------------------------------
# empty object to store the simple ensembles that we will create
ensemble_list <- list()
set.seed(123) # make sure you set the seed OUTSIDE the loop
for (i_repeat in 1:3) {
  # thin the data
  lacerta_thin_rep <- thin_by_cell(lacerta, raster = climate_present)
  lacerta_thin_rep <- thin_by_dist(lacerta_thin_rep, dist_min = 20000)
  # sample pseudo-absences
  lacerta_thin_rep <- sample_pseudoabs(lacerta_thin_rep,
    n = 3 * nrow(lacerta_thin_rep),
    raster = climate_present,
    method = c("dist_min", 50000)
  )
  # get climate
  lacerta_thin_rep <- lacerta_thin_rep %>%
    bind_cols(terra::extract(climate_present, lacerta_thin_rep, ID = FALSE))
  # create folds
  lacerta_thin_rep_cv <- spatial_block_cv(lacerta_thin_rep, v = 5)
  # create a recipe
  lacerta_thin_rep_rec <- recipe(lacerta_thin_rep, formula = class ~ .)
  # create a workflow_set
  lacerta_thin_rep_models <-
    # create the workflow_set
    workflow_set(
      preproc = list(default = lacerta_thin_rep_rec),
      models = list(
        # the standard glm specs
        glm = sdm_spec_glm(),
        # maxent specs with tuning
        maxent = sdm_spec_maxent()
      ),
      # make all combinations of preproc and models,
      cross = TRUE
    ) %>%
    # tweak controls to store information needed later to create the ensemble
    option_add(control = control_ensemble_grid())

  # train the model
  lacerta_thin_rep_models <-
    lacerta_thin_rep_models %>%
    workflow_map("tune_grid",
      resamples = lacerta_thin_rep_cv, grid = 10,
      metrics = sdm_metric_set(), verbose = TRUE
    )
  # make an simple ensemble and add it to the list
  ensemble_list[[i_repeat]] <- simple_ensemble() %>%
    add_member(lacerta_thin_rep_models, metric = "boyce_cont")
}

## -----------------------------------------------------------------------------
lacerta_rep_ens <- repeat_ensemble() %>% add_repeat(ensemble_list)
lacerta_rep_ens

## ----fig.width=6, fig.height=4------------------------------------------------
lacerta_rep_ens <- predict_raster(lacerta_rep_ens, climate_present,
  fun = c("mean", "median")
)
ggplot() +
  geom_spatraster(data = lacerta_rep_ens, aes(fill = median)) +
  scale_fill_terrain_c()

