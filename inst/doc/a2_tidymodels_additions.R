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
lacerta_ensemble

## -----------------------------------------------------------------------------
explainer_lacerta_ens <- explain_tidysdm(lacerta_ensemble)

## ----vip, fig.width=6, fig.height=4-------------------------------------------
library(DALEX)
vip_ensemble <- model_parts(explainer = explainer_lacerta_ens)
plot(vip_ensemble)

## ----pdp, fig.width=6, fig.height=4-------------------------------------------
pdp_bio05 <- model_profile(explainer_lacerta_ens, N = 500, variables = "bio05")
plot(pdp_bio05)

## -----------------------------------------------------------------------------
explainer_list <- explain_tidysdm(tidysdm::lacerta_ensemble, by_workflow = TRUE)

## ----profile, fig.width=6, fig.height=4---------------------------------------
profile_list <- lapply(explainer_list, model_profile,
  N = 500,
  variables = "bio05"
)
plot(profile_list)

## -----------------------------------------------------------------------------
library(tidysdm)
library(sf)
lacerta_thin <- readRDS(system.file("extdata/lacerta_climate_sf.RDS",
  package = "tidysdm"
))

## ----initial_split, fig.width=6, fig.height=4---------------------------------
set.seed(1005)
lacerta_initial <- spatial_initial_split(lacerta_thin,
  prop = 1 / 5, spatial_block_cv
)
autoplot(lacerta_initial)

## -----------------------------------------------------------------------------
check_splits_balance(lacerta_initial, class)

## ----training_cv, fig.width=6, fig.height=4-----------------------------------
set.seed(1005)
lacerta_training <- training(lacerta_initial)
lacerta_cv <- spatial_block_cv(lacerta_training,
  v = 5,
  cellsize = grid_cellsize(lacerta_thin),
  offset = grid_offset(lacerta_thin)
)
autoplot(lacerta_cv)

## -----------------------------------------------------------------------------
check_splits_balance(lacerta_cv, class)

## ----recipe-------------------------------------------------------------------
lacerta_rec_all <- recipe(lacerta_thin, formula = class ~ .)
lacerta_rec_uncor <- lacerta_rec_all %>%
  step_rm(all_of(c(
    "bio01", "bio02", "bio03", "bio04", "bio07", "bio08",
    "bio09", "bio10", "bio11", "bio12", "bio14", "bio16",
    "bio17", "bio18", "bio19", "altitude"
  )))

lacerta_rec_uncor

## ----workflow_set-------------------------------------------------------------
lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(
      uncor = lacerta_rec_uncor, # recipe for the glm
      all = lacerta_rec_all, # recipe for the random forest
      all = lacerta_rec_uncor # recipe for svm
    ),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf(),
      # svm specs with tuning
      svm = parsnip::svm_poly(
        cost = tune(),
        degree = tune()
      ) %>%
        parsnip::set_engine("kernlab") %>%
        parsnip::set_mode("classification")
    ),
    # make all combinations of preproc and models,
    cross = FALSE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  # note that we use the bayes version as we will use a Bayes search (see later)
  option_add(control = stacks::control_stack_bayes())

## -----------------------------------------------------------------------------
rf_param <- lacerta_models %>%
  # extract the rf workflow
  extract_workflow("all_rf") %>%
  # extract its parameters dials (used to tune)
  extract_parameter_set_dials() %>%
  # give it the predictors to finalize mtry
  finalize(x = st_drop_geometry(lacerta_thin) %>% select(-class))

# now update the workflowset with the new parameter info
lacerta_models <- lacerta_models %>%
  option_add(param_info = rf_param, id = "all_rf")

## ----tune_grid----------------------------------------------------------------
set.seed(1234567)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_bayes",
    resamples = lacerta_cv, initial = 8,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## -----------------------------------------------------------------------------
autoplot(lacerta_models)

## ----build_stack, fig.width=6, fig.height=4-----------------------------------
library(stacks)
set.seed(1005)
lacerta_stack <-
  # initialize the stack
  stacks() %>%
  # add candidate members
  add_candidates(lacerta_models) %>%
  # determine how to combine their predictions
  blend_predictions() %>%
  # fit the candidates with non-zero weights (i.e.non-zero stacking coefficients)
  fit_members()

autoplot(lacerta_stack, type = "weights")

## ----predict_test-------------------------------------------------------------
lacerta_testing <- testing(lacerta_initial)

lacerta_test_pred <-
  lacerta_testing %>%
  bind_cols(predict(lacerta_stack, ., type = "prob"))

## ----assess_test--------------------------------------------------------------
sdm_metric_set()(data = lacerta_test_pred, truth = class, .pred_presence)

## ----eval=FALSE---------------------------------------------------------------
#  download_dataset("WorldClim_2.1_10m")
#  climate_vars <- lacerta_rec_all$var_info %>%
#    filter(role == "predictor") %>%
#    pull(variable)
#  
#  climate_present <- pastclim::region_slice(
#    time_ce = 1985,
#    bio_variables = climate_vars,
#    data = "WorldClim_2.1_10m",
#    crop = iberia_poly
#  )

## ----echo=FALSE, results='hide'-----------------------------------------------
terra::gdal(warn = 3)
climate_present <- terra::rast(
  system.file("extdata/lacerta_climate_present_10m.nc", package = "tidysdm")
)
climate_vars <- lacerta_rec_all$var_info %>%
  filter(role == "predictor") %>%
  pull(variable)
if (!all(climate_vars %in% names(climate_present))) {
  stop("mismatched variables in the raster")
}

## ----plot_present, fig.width=6, fig.height=4----------------------------------
prediction_present <- predict_raster(lacerta_stack, climate_present,
  type = "prob"
)
library(tidyterra)
ggplot() +
  geom_spatraster(data = prediction_present, aes(fill = .pred_presence)) +
  scale_fill_terrain_c() +
  # plot presences used in the model
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

