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

## -----------------------------------------------------------------------------
library(tidysdm)
lacerta_thin <- readRDS(system.file("extdata/lacerta_climate_sf.RDS",
  package = "tidysdm"
))
lacerta_thin$bio05[37] <- NA

## -----------------------------------------------------------------------------
lacerta_rec <- recipe(lacerta_thin, formula = class ~ .) %>%
  step_rm(all_of(c(
    "bio01", "bio02", "bio03", "bio04", "bio07", "bio08",
    "bio09", "bio10", "bio11", "bio12", "bio14", "bio16",
    "bio17", "bio18", "bio19", "altitude"
  )))

lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

## ----error=TRUE---------------------------------------------------------------
set.seed(100)
lacerta_cv <- spatial_block_cv(lacerta_thin, v = 5)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
    resamples = lacerta_cv, grid = 3,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## -----------------------------------------------------------------------------
lacerta_prep <- lacerta_rec %>% prep(lacerta_thin)
lacerta_prep

## -----------------------------------------------------------------------------
lacerta_thin <- readRDS(system.file("extdata/lacerta_climate_sf.RDS",
  package = "tidysdm"
))
suggested_vars <- c("bio05", "bio06", "bio13", "bio14", "bio15")
lacerta_rec_sel <- recipe(lacerta_thin, formula = class ~ .) %>%
  step_select(all_of(suggested_vars))

## ----error=TRUE---------------------------------------------------------------
lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec_sel),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

set.seed(100)
lacerta_cv <- spatial_block_cv(lacerta_thin, v = 5)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
    resamples = lacerta_cv, grid = 3,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## -----------------------------------------------------------------------------
lacerta_prep_sel <- lacerta_rec_sel %>% prep(lacerta_thin)
lacerta_prep_sel

## -----------------------------------------------------------------------------
lacerta_thin <- readRDS(system.file("extdata/lacerta_climate_sf.RDS",
  package = "tidysdm"
))

lacerta_rec <- recipe(lacerta_thin, formula = class ~ .) %>%
  step_rm(all_of(c(
    "bio01", "bio02", "bio03", "bio04", "bio07", "bio08",
    "bio09", "bio10", "bio11", "bio12", "bio14", "bio16",
    "bio17", "bio18", "bio19", "altitude"
  )))

lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # the standard gam specs
      gam = sdm_spec_gam()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # set formula for gams
  update_workflow_model("default_gam",
    spec = sdm_spec_gam(),
    formula = gam_formula(lacerta_rec)
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

## -----------------------------------------------------------------------------
set.seed(100)
lacerta_cv <- spatial_block_cv(lacerta_thin, v = 5)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
    resamples = lacerta_cv, grid = 3,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## -----------------------------------------------------------------------------
lacerta_thin <- readRDS(system.file("extdata/lacerta_climate_sf.RDS",
  package = "tidysdm"
))
set.seed(123)
lacerta_thin <- lacerta_thin[sample(
  1:nrow(lacerta_thin),
  nrow(lacerta_thin) / 5
), ]

lacerta_rec <- recipe(lacerta_thin, formula = class ~ .) %>%
  step_rm(all_of(c(
    "bio01", "bio02", "bio03", "bio04", "bio07", "bio08",
    "bio09", "bio10", "bio11", "bio12", "bio14", "bio16",
    "bio17", "bio18", "bio19", "altitude"
  )))

lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # the standard gam specs
      gam = sdm_spec_gam(),
      # rf specs with tuning
      rf = sdm_spec_rf()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # set formula for gams
  update_workflow_model("default_gam",
    spec = sdm_spec_gam(),
    formula = gam_formula(lacerta_rec)
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

## -----------------------------------------------------------------------------
set.seed(100)
lacerta_cv <- spatial_block_cv(lacerta_thin, v = 3)
lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
    resamples = lacerta_cv, grid = 3,
    metrics = sdm_metric_set(), verbose = TRUE
  )

## -----------------------------------------------------------------------------
gam_results <- extract_workflow_set_result(lacerta_models, id = "default_gam")
gam_results

## -----------------------------------------------------------------------------
gam_results$.notes[2]

## -----------------------------------------------------------------------------
problem_split <- gam_results$splits[2][[1]]
summary(training(problem_split))

## -----------------------------------------------------------------------------
gam_workflow <- extract_workflow(lacerta_models, id = "default_gam")
faulty_gam <- fit(gam_workflow, training(problem_split))

