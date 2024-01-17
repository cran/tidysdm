## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# set options to limit threads used by imported libraries
options(cores=2)
options(mc.cores=2)
# xgboost uses data.table
data.table::setDTthreads(1)

## -----------------------------------------------------------------------------
library(tidysdm)

## ----load_presences-----------------------------------------------------------
data(horses)
horses

## -----------------------------------------------------------------------------
library(sf)
horses <- st_as_sf(horses, coords = c("longitude", "latitude"))
st_crs(horses) <- 4326

## ----echo=FALSE, results='hide'-----------------------------------------------
library(pastclim)
set_data_path(on_CRAN = TRUE)

## -----------------------------------------------------------------------------
library(pastclim)
land_mask <- pastclim::get_land_mask(time_bp = 0, dataset = "Example")
europe_poly <- vect(region_outline$Europe)
crs(europe_poly) <- "lonlat"
land_mask <- crop(land_mask, europe_poly)
land_mask <- mask(land_mask, europe_poly)

## ----cast_to_sf, fig.width=6, fig.height=4------------------------------------
library(tidyterra)
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  scale_fill_terrain_c() +
  geom_sf(data = horses, aes(col = time_bp))

## ----thin_by_dist-------------------------------------------------------------
set.seed(123)
horses <- thin_by_dist_time(horses,
  dist_min = km2m(100),
  interval_min = y2d(2000),
  time_col = "time_bp",
  lubridate_fun = pastclim::ybp2date
)
nrow(horses)

## ----plot_thinned, fig.width=6, fig.height=4----------------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  scale_fill_terrain_c() +
  geom_sf(data = horses, aes(col = time_bp))

## ----load_climate-------------------------------------------------------------
library(pastclim)
climate_vars <- c("bio01", "bio10", "bio12")
climate_full <- pastclim::region_series(
  bio_variables = climate_vars,
  data = "Example",
  crop = region_outline$Europe
)

## ----thin_by_cell-------------------------------------------------------------
set.seed(123)
horses <- thin_by_cell_time(horses,
  raster = climate_full,
  time_col = "time_bp",
  lubridate_fun = pastclim::ybp2date
)
nrow(horses)

## ----fig.width=6, fig.height=4------------------------------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  scale_fill_terrain_c() +
  geom_sf(data = horses, aes(col = time_bp))

## -----------------------------------------------------------------------------
set.seed(123)
horses <- sample_pseudoabs_time(horses,
  n_per_presence = 3,
  raster = climate_full,
  time_col = "time_bp",
  lubridate_fun = pastclim::ybp2date,
  method = c("dist_min", km2m(70))
)

## ----fig.width=6, fig.height=4------------------------------------------------
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  scale_fill_terrain_c() +
  geom_sf(data = horses, aes(col = class))

## ----climate_for_locations----------------------------------------------------
horses_df <- horses %>%
  dplyr::bind_cols(sf::st_coordinates(horses)) %>%
  mutate(time_bp = date2ybp(time_step)) %>%
  as.data.frame() %>%
  select(-geometry)
# get climate
horses_df <- location_slice_from_region_series(horses_df,
  region_series = climate_full
)

# add the climate reconstructions to the sf object, and remove the time_step
# as we don't need it for modelling
horses <- horses %>%
  bind_cols(horses_df[, climate_vars]) %>%
  select(-time_step)

## ----recipe-------------------------------------------------------------------
horses_rec <- recipe(horses, formula = class ~ .)
horses_rec

## -----------------------------------------------------------------------------
horses_rec$var_info

## ----workflow_set-------------------------------------------------------------
horses_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = horses_rec),
    models = list(
      # the standard glm specs  (no params to tune)
      glm = sdm_spec_glm(),
      # the standard sdm specs (no params to tune)
      gam = sdm_spec_gam(),
      # rf specs with tuning
      rf = sdm_spec_rf(),
      # boosted tree model (gbm) specs with tuning
      gbm = sdm_spec_boost_tree()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # set formula for gams
  update_workflow_model("default_gam",
    spec = sdm_spec_gam(),
    formula = gam_formula(horses_rec)
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

## ----training_cv, fig.width=6, fig.height=4-----------------------------------
library(tidysdm)
set.seed(1005)
horses_cv <- spatial_block_cv(horses, v = 5)
autoplot(horses_cv)

## ----tune_grid----------------------------------------------------------------
set.seed(123)
horses_models <-
  horses_models %>%
  workflow_map("tune_grid",
    resamples = horses_cv, grid = 5,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## ----fig.width=6, fig.height=4------------------------------------------------
autoplot(horses_models)

## -----------------------------------------------------------------------------
horses_ensemble <- simple_ensemble() %>%
  add_member(horses_models, metric = "boyce_cont")

## ----fig.width=6, fig.height=4------------------------------------------------
autoplot(horses_ensemble)

## ----get_lgm------------------------------------------------------------------
climate_lgm <- pastclim::region_slice(
  time_bp = -20000,
  bio_variables = climate_vars,
  data = "Example",
  crop = region_outline$Europe
)

## ----plot_lgm, fig.width=6, fig.height=4--------------------------------------
prediction_lgm <- predict_raster(horses_ensemble, climate_lgm)
ggplot() +
  geom_spatraster(data = prediction_lgm, aes(fill = mean)) +
  scale_fill_terrain_c()

