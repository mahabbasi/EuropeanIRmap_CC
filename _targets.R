library(targets)
library(tarchetypes)

tar_option_set(packages = c("sf", "data.table", "here"),
               memory = "transient", garbage_collection = TRUE)

tar_source("src/used_libraries.R")
tar_source("src/pre_processing_functions.R")
tar_source("src/make_prediction_functions.R")
tar_source("src/post_processing_fuctions.R")
tar_source("src/CAST_funs_editted.R")
#----------------------- define the paths ----------------------------------------
root_dir <- here::here()
data_dir <- file.path(root_dir, 'data')
res_dir <- file.path(root_dir, 'results')

#Create new directories if they are not existed
watergap_extractdir <- file.path(data_dir, 'watergap_netrextract')

lapply(list(watergap_extractdir), function(path) {
              if (!dir.exists(path)) {
                dir.create(path, recursive = TRUE)
              }
  })


# ---------------------- preparing high-resolution predicots ------------
# a tibble of high-resolution predictors combination for mapping
params <- tribble(
  ~name, ~n_window,  ~Fun, ~pred_name, ~out_name,
  "mean3", 3  ,"mean",  'mean_p3', 'q_mean_p3_eu', 
  "mean12", 12, "mean",  'mean_p12','q_mean_p12_eu',
  "min3", 3  ,"min",  'min_p3', 'q_min_p3_eu', 
  "minn12", 12, "min",  'min_p12','q_min_p12_eu'
) 

prepare_hr_predictors <- tar_plan(
  
  ## gfdl-esm4 GCM --------
  tar_target(
    name = 'hr_flow_gfdl_hist',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/GFDL-ESM4/hist"),
                                out_dir = file.path(watergap_extractdir, 'GFDL-ESM4/historical'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'gfdl-esm4',
                                phase = 'hist')
  ),
  tar_target(
    name = 'hr_flow_gfdl_ssp126',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/GFDL-ESM4/future/ssp126"),
                                out_dir = file.path(watergap_extractdir, 'GFDL-ESM4/future/ssp126'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'gfdl-esm4',
                                phase = 'ssp126')
  ),
  tar_target(
    name = 'hr_flow_gfdl_ssp585',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/GFDL-ESM4/future/ssp585"),
                                out_dir = file.path(watergap_extractdir, 'GFDL-ESM4/future/ssp585'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'gfdl-esm4',
                                phase = 'ssp585')
  ),
  
  ## ipsl-cm6a-lr GCM --------
  tar_target(
    name = 'hr_flow_ipsl_hist',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/IPSL-CM6A-LR/hist"),
                                out_dir = file.path(watergap_extractdir, 'IPSL-CM6A-LR/historical'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ipsl-cm6a-lr',
                                phase = 'hist')
  ),
  tar_target(
    name = 'hr_flow_ipsl_ssp126',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/IPSL-CM6A-LR/future/ssp126"),
                                out_dir = file.path(watergap_extractdir, 'IPSL-CM6A-LR/future/ssp126'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ipsl-cm6a-lr',
                                phase = 'ssp126')
  ),
  tar_target(
    name = 'hr_flow_ipsl_ssp585',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/IPSL-CM6A-LR/future/ssp585"),
                                out_dir = file.path(watergap_extractdir, 'IPSL-CM6A-LR/future/ssp585'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ipsl-cm6a-lr',
                                phase = 'ssp585')
  ),
  
  ## MPI-ESM1-2-HR GCM -----
  tar_target(
    name = 'hr_flow_mpi_hist',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MPI-ESM1-2-HR/hist"),
                                out_dir = file.path(watergap_extractdir, 'MPI-ESM1-2-HR/historical'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mpi-esm1-2-hr',
                                phase = 'hist')
  ),
  tar_target(
    name = 'hr_flow_mpi_ssp126',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MPI-ESM1-2-HR/future/ssp126"),
                                out_dir = file.path(watergap_extractdir, 'MPI-ESM1-2-HR/future/ssp126'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mpi-esm1-2-hr',
                                phase = 'ssp126')
  ),
  tar_target(
    name = 'hr_flow_mpi_ssp585',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MPI-ESM1-2-HR/future/ssp585"),
                                out_dir = file.path(watergap_extractdir, 'MPI-ESM1-2-HR/future/ssp585'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mpi-esm1-2-hr',
                                phase = 'ssp585')
  ),
  # MRI-ESM2-0 GCM -----------
  tar_target(
    name = 'hr_flow_mri_hist',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MRI-ESM2-0/hist"),
                                out_dir = file.path(watergap_extractdir, 'MRI-ESM2-0/historical'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mri-esm2-0',
                                phase = 'hist')
  ),
  tar_target(
    name = 'hr_flow_mri_ssp126',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MRI-ESM2-0/future/ssp126"),
                                out_dir = file.path(watergap_extractdir, 'MRI-ESM2-0/future/ssp126'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mri-esm2-0',
                                phase = 'ssp126')
  ),
  tar_target(
    name = 'hr_flow_mri_ssp585',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/MRI-ESM2-0/future/ssp585"),
                                out_dir = file.path(watergap_extractdir, 'MRI-ESM2-0/future/ssp585'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'mri-esm2-0',
                                phase = 'ssp585')
  ),
  
  # UKESM1-0-LL GCM -----------
  tar_target(
    name = 'hr_flow_ukesm_hist',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/UKESM1-0-LL/hist"),
                                out_dir = file.path(watergap_extractdir, 'UKESM1-0-LL/historical'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ukesm1-0-ll',
                                phase = 'hist')
  ),
  tar_target(
    name = 'hr_flow_ukesm_ssp126',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/UKESM1-0-LL/future/ssp126"),
                                out_dir = file.path(watergap_extractdir, 'UKESM1-0-LL/future/ssp126'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ukesm1-0-ll',
                                phase = 'ssp126')
  ),
  tar_target(
    name = 'hr_flow_ukesm_ssp585',
    command = extract_dsflow_hr(ncs_dir = file.path(data_dir, "raw",
                                                    "Climate_forcing/UKESM1-0-LL/future/ssp585"),
                                out_dir = file.path(watergap_extractdir, 'UKESM1-0-LL/future/ssp585'),
                                prpts_path = file.path(root_dir, 
                                                       'data/others/dryvernet_pourpoints_eu/dryvernet_pourpoints_eu.csv'),
                                model_name = 'ukesm1-0-ll',
                                phase = 'ssp585')
  ),
  
  # compute high-resolution predictors -------------------
  tar_map(
    values = params,
    names = name,
    
    # GFDL-ESM4 
    tar_target(gfdl_hist,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/GFDL-ESM4/historical'),
                                     out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/historical'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=1985, end_year=2014,
                                     start_date='1985-01-01', end_date='2014-12-01')),
    tar_target(gfdl_ssp126,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/GFDL-ESM4/future/ssp126'),
                                     out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/future/ssp126'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    tar_target(gfdl_ssp585,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/GFDL-ESM4/future/ssp585'),
                                     out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/future/ssp585'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    # IPSL-CM6A-LR 
    tar_target(ipsl_hist,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/IPSL-CM6A-LR/historical'),
                                     out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/historical'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=1985, end_year=2014,
                                     start_date='1985-01-01', end_date='2014-12-01')),
    tar_target(ipsl_ssp126,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/IPSL-CM6A-LR/future/ssp126'),
                                     out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/future/ssp126'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    tar_target(ipsl_ssp585,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/IPSL-CM6A-LR/future/ssp585'),
                                     out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/future/ssp585'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    # MPI-ESM1-2-HR 
    tar_target(mpi_hist,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MPI-ESM1-2-HR/historical'),
                                     out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/historical'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=1985, end_year=2014,
                                     start_date='1985-01-01', end_date='2014-12-01')),
    tar_target(mpi_ssp126,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MPI-ESM1-2-HR/future/ssp126'),
                                     out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/future/ssp126'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    tar_target(mpi_ssp585,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MPI-ESM1-2-HR/future/ssp585'),
                                     out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/future/ssp585'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    # MRI-ESM2-0 
    tar_target(mri_hist,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MRI-ESM2-0/historical'),
                                     out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/historical'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=1985, end_year=2014,
                                     start_date='1985-01-01', end_date='2014-12-01')),
    tar_target(mri_ssp126,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MRI-ESM2-0/future/ssp126'),
                                     out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/future/ssp126'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    tar_target(mri_ssp585,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/MRI-ESM2-0/future/ssp585'),
                                     out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/future/ssp585'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    # UKESM1-0-LL
    tar_target(ukesm_hist,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/UKESM1-0-LL/historical'),
                                     out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/historical'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=1985, end_year=2014,
                                     start_date='1985-01-01', end_date='2014-12-01')),
    tar_target(ukesm_ssp126,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/UKESM1-0-LL/future/ssp126'),
                                     out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/future/ssp126'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01')),
    tar_target(ukesm_ssp585,
               compute_hr_predictors(in_path = file.path(data_dir, 'watergap_netrextract/UKESM1-0-LL/future/ssp585'),
                                     out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/future/ssp585'),
                                     n_window=n_window, Fun=Fun, pred_name=pred_name,out_name=out_name,
                                     start_year=2041, end_year=2100,
                                     start_date='2041-01-01', end_date='2100-12-01'))
    
  ),
  
  # high-resolution inter-annual predictors (sd and cv) -----------
  # GFDL-ESM4 GCM
  tar_target(
    name = 'gfdl_hist_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/GFDL-ESM4/historical'),
                                                out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/historical'),
                                                start_year=1985, end_year=2014)
  ),
  tar_target(
    name = 'gfdl_ssp126_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/GFDL-ESM4/future/ssp126'),
                                                out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/future/ssp126'),
                                                start_year=2041, end_year=2100)
  ),
  tar_target(
    name = 'gfdl_ssp585_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/GFDL-ESM4/future/ssp585'),
                                                out_path=file.path(data_dir, 'predictors/HR/GFDL-ESM4/future/ssp585'),
                                                start_year=2041, end_year=2100)
  ),
  # IPSL-CM6A-LR GCM
  tar_target(
    name = 'ipsl_hist_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/IPSL-CM6A-LR/historical'),
                                                out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/historical'),
                                                start_year=1985, end_year=2014)
  ),
  tar_target(
    name = 'ipsl_ssp126_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/IPSL-CM6A-LR/future/ssp126'),
                                                out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/future/ssp126'),
                                                start_year=2041, end_year=2100)
  ),
  tar_target(
    name = 'ipsl_ssp585_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/IPSL-CM6A-LR/future/ssp585'),
                                                out_path=file.path(data_dir, 'predictors/HR/IPSL-CM6A-LR/future/ssp585'),
                                                start_year=2041, end_year=2100)
  ),
  # MPI-ESM1-2-HR GCM
  tar_target(
    name = 'mpi_hist_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MPI-ESM1-2-HR/historical'),
                                                out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/historical'),
                                                start_year=1985, end_year=2014)
  ),
  tar_target(
    name = 'mpi_ssp126_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MPI-ESM1-2-HR/future/ssp126'),
                                                out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/future/ssp126'),
                                                start_year=2041, end_year=2100)
  ),
  tar_target(
    name = 'mpi_ssp585_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MPI-ESM1-2-HR/future/ssp585'),
                                                out_path=file.path(data_dir, 'predictors/HR/MPI-ESM1-2-HR/future/ssp585'),
                                                start_year=2041, end_year=2100)
  ),
  # MRI-ESM2-0 GCM
  tar_target(
    name = 'mri_hist_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MRI-ESM2-0/historical'),
                                                out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/historical'),
                                                start_year=1985, end_year=2014)
  ),
  tar_target(
    name = 'mri_ssp126_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MRI-ESM2-0/future/ssp126'),
                                                out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/future/ssp126'),
                                                start_year=2041, end_year=2100)
  ),
  tar_target(
    name = 'mri_ssp585_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/MRI-ESM2-0/future/ssp585'),
                                                out_path=file.path(data_dir, 'predictors/HR/MRI-ESM2-0/future/ssp585'),
                                                start_year=2041, end_year=2100)
  ),
  # UKESM1-0-LL GCM
  tar_target(
    name = 'ukesm_hist_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/UKESM1-0-LL/historical'),
                                                out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/historical'),
                                                start_year=1985, end_year=2014)
  ),
  tar_target(
    name = 'ukesm_ssp126_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/UKESM1-0-LL/future/ssp126'),
                                                out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/future/ssp126'),
                                                start_year=2041, end_year=2100)
  ),
  tar_target(
    name = 'ukesm_ssp585_interannual',
    command = compute_hr_interannual_predictors(in_path = file.path(data_dir, 
                                                                    'watergap_netrextract/UKESM1-0-LL/future/ssp585'),
                                                out_path=file.path(data_dir, 'predictors/HR/UKESM1-0-LL/future/ssp585'),
                                                start_year=2041, end_year=2100)
  )
)

# ---------------------- preparing low-resolution predicots ---------------------------
prepare_lr_predictors <- tar_plan(
  
  # GFDL-ESM4 --------
  tar_target(
    name = 'gfdl_hist_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/gfdl_csv/wetdays",
                                                             'gfdl_esm4_r1i1p1f1_w5e5_historical_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/historical'),
                                          period_phase = 'hist', var_name='wetdays')),
  tar_target(
    name = 'gfdl_ssp126_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/gfdl_csv/wetdays",
                                                             'gfdl_esm4_r1i1p1f1_w5e5_ssp126_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/future/ssp126'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'gfdl_ssp585_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/gfdl_csv/wetdays",
                                                             'gfdl_esm4_r1i1p1f1_w5e5_ssp585_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/future/ssp585'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'gfdl_hist_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/gfdl_csv/qrd",
                                                             'watergap2_2e_gfdl_esm4_w5e5_historical_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/historical'),
                                          period_phase = 'hist', var_name='qrdifoverql')),
  tar_target(
    name = 'gfdl_ssp126_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/gfdl_csv/qrd",
                                                             'watergap2_2e_gfdl_esm4_w5e5_ssp126_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/future/ssp126'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  tar_target(
    name = 'gfdl_ssp585_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/gfdl_csv/qrd",
                                                             'watergap2_2e_gfdl_esm4_w5e5_ssp585_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/GFDL-ESM4/future/ssp585'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  
  # IPSL-CM6A ------
  tar_target(
    name = 'ipsl_hist_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/ipsl_csv/wetdays",
                                                             'ipsl_cm6a_lr_r1i1p1f1_w5e5_historical_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/historical'),
                                          period_phase = 'hist', var_name='wetdays')),
  tar_target(
    name = 'ipsl_ssp126_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ipsl_csv/wetdays",
                                                             'ipsl_cm6a_lr_r1i1p1f1_w5e5_ssp126_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/future/ssp126'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'ipsl_ssp585_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ipsl_csv/wetdays",
                                                             'ipsl_cm6a_lr_r1i1p1f1_w5e5_ssp585_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/future/ssp585'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'ipsl_hist_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ipsl_csv/qrd",
                                                             'watergap2_2e_ipsl_cm6a_lr_w5e5_historical_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/historical'),
                                          period_phase = 'hist', var_name='qrdifoverql')),
  tar_target(
    name = 'ipsl_ssp126_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ipsl_csv/qrd",
                                                             'watergap2_2e_ipsl_cm6a_lr_w5e5_ssp126_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/future/ssp126'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  tar_target(
    name = 'ipsl_ssp585_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ipsl_csv/qrd",
                                                             'watergap2_2e_ipsl_cm6a_lr_w5e5_ssp585_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/IPSL-CM6A/future/ssp585'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  
  # MPI-ESM1-2-HR ------
  tar_target(
    name = 'mpi_hist_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/mpi_csv/wetdays",
                                                             'mpi_esm1_2_hr_r1i1p1f1_w5e5_historical_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/historical'),
                                          period_phase = 'hist', var_name='wetdays')),
  tar_target(
    name = 'mpi_ssp126_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mpi_csv/wetdays",
                                                             'mpi_esm1_2_hr_r1i1p1f1_w5e5_ssp126_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/future/ssp126'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'mpi_ssp585_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mpi_csv/wetdays",
                                                             'mpi_esm1_2_hr_r1i1p1f1_w5e5_ssp585_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/future/ssp585'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'mpi_hist_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mpi_csv/qrd",
                                                             'watergap2_2e_mpi_esm1_2_hr_w5e5_historical_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/historical'),
                                          period_phase = 'hist', var_name='qrdifoverql')),
  tar_target(
    name = 'mpi_ssp126_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mpi_csv/qrd",
                                                             'watergap2_2e_mpi_esm1_2_hr_w5e5_ssp126_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/future/ssp126'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  tar_target(
    name = 'mpi_ssp585_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mpi_csv/qrd",
                                                             'watergap2_2e_mpi_esm1_2_hr_w5e5_ssp585_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MPI-ESM1-2-HR/future/ssp585'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  # MRI-ESM2-0 ------
  tar_target(
    name = 'mri_hist_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/mri_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_historical_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/historical'),
                                          period_phase = 'hist', var_name='wetdays')),
  tar_target(
    name = 'mri_ssp126_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mri_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_ssp126_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/future/ssp126'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'mri_ssp585_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mri_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_ssp585_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/future/ssp585'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'mri_hist_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mri_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_historical_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/historical'),
                                          period_phase = 'hist', var_name='qrdifoverql')),
  tar_target(
    name = 'mri_ssp126_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mri_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_ssp126_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/future/ssp126'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  tar_target(
    name = 'mri_ssp585_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/mri_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_ssp585_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/MRI-ESM2-0/future/ssp585'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  
  # UKESM1-0-LL ------
  tar_target(
    name = 'ukesm_hist_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/ukesm_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_historical_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/historical'),
                                          period_phase = 'hist', var_name='wetdays')),
  tar_target(
    name = 'ukesm_ssp126_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ukesm_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_ssp126_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/future/ssp126'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'ukesm_ssp585_wetdays',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ukesm_csv/wetdays",
                                                             'mri_esm2_0_r1i1p1f1_w5e5_ssp585_wetdays.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/future/ssp585'),
                                          period_phase = 'future', var_name='wetdays')),
  tar_target(
    name = 'ukesm_hist_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ukesm_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_historical_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/historical'),
                                          period_phase = 'hist', var_name='qrdifoverql')),
  tar_target(
    name = 'ukesm_ssp126_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ukesm_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_ssp126_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/future/ssp126'),
                                          period_phase = 'future', var_name='qrdifoverql')),
  tar_target(
    name = 'ukesm_ssp585_qrd',
    command = select_modify_lr_predictors(in_dir = file.path(data_dir,"predictors/LR/ukesm_csv/qrd",
                                                             'watergap2_2e_mri_esm2_0_w5e5_ssp585_2015soc_from_histsoc_qrdifoverql.csv'),
                                          reachids_dir = file.path(data_dir, "others/european_reaches_DRYVER_RIVID.csv"),
                                          out_dir = file.path(data_dir, 'predictors/LR/UKESM1-0-LL/future/ssp585'),
                                          period_phase = 'future', var_name='qrdifoverql')
    )
  
)


# ---------------------- run RF model for the five GCMs ------------
run_rf_gcms <- tar_plan(
  # GFDL-ESM4 -------
  tar_target(
    name = 'gfdl_hist_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/GFDL-ESM4/historical"),
                                    path_HR = file.path(data_dir, "predictors/HR/GFDL-ESM4/historical"),
                                    start_year = 1985, end_year = 2014,
                                    outdir = file.path(res_dir, 'predictions/GFDL-ESM4/historical'),
                                    period_phase = 'hist')),
  tar_target(
    name = 'gfdl_ssp126_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/GFDL-ESM4/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/GFDL-ESM4/future/ssp126"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/GFDL-ESM4/future/ssp126'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'gfdl_ssp126_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/GFDL-ESM4/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/GFDL-ESM4/future/ssp126"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/GFDL-ESM4/future/ssp126'),
                                    period_phase = 'farfuture')),
  tar_target(
    name = 'gfdl_ssp585_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/GFDL-ESM4/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/GFDL-ESM4/future/ssp585"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/GFDL-ESM4/future/ssp585'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'gfdl_ssp585_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/GFDL-ESM4/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/GFDL-ESM4/future/ssp585"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/GFDL-ESM4/future/ssp585'),
                                    period_phase = 'farfuture')),
  # IPSL-CM6A-LR -------
  tar_target(
    name = 'ipsl_hist_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/IPSL-CM6A-LR/historical"),
                                    path_HR = file.path(data_dir, "predictors/HR/IPSL-CM6A-LR/historical"),
                                    start_year = 1985, end_year = 2014,
                                    outdir = file.path(res_dir, 'predictions/IPSL-CM6A-LR/historical'),
                                    period_phase = 'hist')),
  tar_target(
    name = 'ipsl_ssp126_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/IPSL-CM6A-LR/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/IPSL-CM6A-LR/future/ssp126"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/IPSL-CM6A-LR/future/ssp126'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'ipsl_ssp126_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/IPSL-CM6A-LR/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/IPSL-CM6A-LR/future/ssp126"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/IPSL-CM6A-LR/future/ssp126'),
                                    period_phase = 'farfuture')),
  tar_target(
    name = 'ipsl_ssp585_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/IPSL-CM6A-LR/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/IPSL-CM6A-LR/future/ssp585"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/IPSL-CM6A-LR/future/ssp585'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'ipsl_ssp585_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/IPSL-CM6A-LR/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/IPSL-CM6A-LR/future/ssp585"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/IPSL-CM6A-LR/future/ssp585'),
                                    period_phase = 'farfuture')),
  # MPI-ESM1-2-HR -------
  tar_target(
    name = 'mpi_hist_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MPI-ESM1-2-HR/historical"),
                                    path_HR = file.path(data_dir, "predictors/HR/MPI-ESM1-2-HR/historical"),
                                    start_year = 1985, end_year = 2014,
                                    outdir = file.path(res_dir, 'predictions/MPI-ESM1-2-HR/historical'),
                                    period_phase = 'hist')),
  tar_target(
    name = 'mpi_ssp126_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MPI-ESM1-2-HR/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/MPI-ESM1-2-HR/future/ssp126"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/MPI-ESM1-2-HR/future/ssp126'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'mpi_ssp126_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MPI-ESM1-2-HR/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/MPI-ESM1-2-HR/future/ssp126"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/MPI-ESM1-2-HR/future/ssp126'),
                                    period_phase = 'farfuture')),
  tar_target(
    name = 'mpi_ssp585_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MPI-ESM1-2-HR/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/MPI-ESM1-2-HR/future/ssp585"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/MPI-ESM1-2-HR/future/ssp585'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'mpi_ssp585_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MPI-ESM1-2-HR/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/MPI-ESM1-2-HR/future/ssp585"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/MPI-ESM1-2-HR/future/ssp585'),
                                    period_phase = 'farfuture')),
  # MRI-ESM2-0 -------
  tar_target(
    name = 'mri_hist_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MRI-ESM2-0/historical"),
                                    path_HR = file.path(data_dir, "predictors/HR/MRI-ESM2-0/historical"),
                                    start_year = 1985, end_year = 2014,
                                    outdir = file.path(res_dir, 'predictions/MRI-ESM2-0/historical'),
                                    period_phase = 'hist')),
  tar_target(
    name = 'mri_ssp126_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MRI-ESM2-0/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/MRI-ESM2-0/future/ssp126"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/MRI-ESM2-0/future/ssp126'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'mri_ssp126_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MRI-ESM2-0/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/MRI-ESM2-0/future/ssp126"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/MRI-ESM2-0/future/ssp126'),
                                    period_phase = 'farfuture')),
  tar_target(
    name = 'mri_ssp585_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MRI-ESM2-0/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/MRI-ESM2-0/future/ssp585"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/MRI-ESM2-0/future/ssp585'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'mri_ssp585_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/MRI-ESM2-0/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/MRI-ESM2-0/future/ssp585"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/MRI-ESM2-0/future/ssp585'),
                                    period_phase = 'farfuture')),
  # UKESM1-0-LL -------
  tar_target(
    name = 'ukesm_hist_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/UKESM1-0-LL/historical"),
                                    path_HR = file.path(data_dir, "predictors/HR/UKESM1-0-LL/historical"),
                                    start_year = 1985, end_year = 2014,
                                    outdir = file.path(res_dir, 'predictions/UKESM1-0-LL/historical'),
                                    period_phase = 'hist')),
  tar_target(
    name = 'ukesm_ssp126_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/UKESM1-0-LL/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/UKESM1-0-LL/future/ssp126"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/UKESM1-0-LL/future/ssp126'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'ukesm_ssp126_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/UKESM1-0-LL/future/ssp126"),
                                    path_HR = file.path(data_dir, "predictors/HR/UKESM1-0-LL/future/ssp126"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/UKESM1-0-LL/future/ssp126'),
                                    period_phase = 'farfuture')),
  tar_target(
    name = 'ukesm_ssp585_nearfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/UKESM1-0-LL/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/UKESM1-0-LL/future/ssp585"),
                                    start_year = 2041, end_year = 2070,
                                    outdir = file.path(res_dir, 'predictions/UKESM1-0-LL/future/ssp585'),
                                    period_phase = 'nearfuture')),
  tar_target(
    name = 'ukesm_ssp585_farfuture_run',
    command = runmodels_over_period(path_model1 = file.path(data_dir,'rf_models/rftuned_step1.qs'),
                                    path_model2 = file.path(data_dir,'rf_models/rftuned_step2.qs'),
                                    path_static = file.path(data_dir,'predictors/statics/static_preds_net_eu.fst'),
                                    path_LR = file.path(data_dir, "predictors/LR/UKESM1-0-LL/future/ssp585"),
                                    path_HR = file.path(data_dir, "predictors/HR/UKESM1-0-LL/future/ssp585"),
                                    start_year = 2071, end_year = 2100,
                                    outdir = file.path(res_dir, 'predictions/UKESM1-0-LL/future/ssp585'),
                                    period_phase = 'farfuture')
    )
)

# -------------------------- post-processsing ----------------------
create_final_products <- tar_plan(
  tar_target(
    name = 'create_shp_fig2ab',
    command = compute_qtot_fig2ab(in_dir=file.path(data_dir, 'predictors/LR/LR_qtot'),
                                  out_dir=file.path(res_dir, 'figures/figure2'),
                                  eu_countries_shp_dir=file.path(data_dir, 'shp/eu_countries_shp'))
  ),
  tar_target(
    name = 'create_shp_fig2cd',
    command = compute_hr_pc_fig2cd(in_dir=file.path(data_dir, 'predictors/HR'),
                                  out_dir=file.path(res_dir, 'figures/figure2'),
                                  eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'create_shp_fig2ef',
    command = compute_lowflow_pc_fig2ef(in_dir=file.path(data_dir, 'predictors/HR'),
                                        out_dir=file.path(res_dir, 'figures/figure2'),
                                        eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'create_shp_fig3',
    command = compute_np_changes_fig3(in_dir=file.path(res_dir, 'predictions/GCMs'),
                                      out_dir=file.path(res_dir, 'figures/figure3'),
                                      eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'create_shp_figS2',
    command = compute_np_changes_fig3(in_dir=file.path(res_dir, 'predictions/GCMs'),
                                      out_dir=file.path(res_dir, 'figures/figureS2'),
                                      eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'),
                                      rcp_name='RCP2.6')
  ),
  tar_target(
    name = 'create_plot_fig4',
    command = ggcdf_np_climzones_fig4(in_dir=file.path(res_dir,'predictions/GCMs'),
                                      out_dir=file.path(res_dir, 'figures/figure4'),
                                      in_shp_climcode_path=file.path(data_dir,
                                                                     'shp/climatecode_shp/eu_nets_withclimcode_final.shp'),
                                      climate_zone_png_dir=file.path(data_dir, 'others'))
  ),
  tar_target(
    name = 'create_plot_fig5',
    command = ggline_seasonality_fig5(in_dir=file.path(res_dir,'predictions/GCMs'),
                                      out_dir=file.path(res_dir, 'figures/figure5'))
  ),
  tar_target(
    name = 'create_plot_fig6',
    command = ggbar_changes_fig6(out_dir=file.path(res_dir, 'figures/figure6'))
  ),
  tar_target(
    name = 'create_shp_fig7',
    command = compute_interannual_vari_fig7(in_dir=file.path(res_dir,'predictions/GCMs'),
                                            out_dir=file.path(res_dir, 'figures/figure7'),
                                            eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'create_shp_figS5',
    command = compute_interannual_vari_fig7(in_dir=file.path(res_dir,'predictions/GCMs'),
                                            out_dir=file.path(res_dir, 'figures/figureS5'),
                                            eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'),
                                            rcp_name='RCP2.6')
  ),
  tar_target(
    name = 'create_plot_figS3',
    command = ggline_seasonality_figS3(in_dir=file.path(res_dir,'predictions/GCMs'),
                                       out_dir=file.path(res_dir, 'figures/figureS3'),
                                       in_shp_climcode_path=file.path(data_dir,
                                                                      'shp/climatecode_shp/eu_nets_withclimcode_final.shp'))
  ),
  tar_target(
    name = 'create_shp_figS4',
    command = create_shp_figS4(in_dir=file.path(res_dir,'predictions/GCMs/future/ssp585/far_future'),
                               out_dir=file.path(res_dir, 'figures/figureS4'),
                               eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'generate_table4',
    command = compute_status_shifts_table4(in_dir=file.path(res_dir,'predictions/GCMs'),
                                           in_shp_climcode_path=file.path(data_dir,
                                                                          'shp/climatecode_shp/eu_nets_withclimcode_final.shp'),
                                           threshold=354)
  ),
  tar_target(
    name = 'generate_table3',
    command = compute_changes_table3(in_dir=file.path(res_dir, 'predictions/GCMs'),
                                     eu_net_shp_dir=file.path(data_dir, 'shp/eu_nets/dryver_net_eu_final.shp'))
  ),
  tar_target(
    name = 'generate_tableS3',
    command = compute_status_shifts_table4(in_dir=file.path(res_dir,'predictions/GCMs'),
                                           in_shp_climcode_path=file.path(data_dir,
                                                                          'shp/climatecode_shp/eu_nets_withclimcode_final.shp'),
                                           threshold=346)
  ),
  tar_target(
    name = 'reanalysis_aoa',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reanalysis'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reanalysis'))
  ),
  tar_target(
    name = 'reference_aoa_gdal',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reference/gfdl'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reference/gfdl'))
  ),
  tar_target(
    name = 'reference_aoa_ipsl',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reference/ipsl'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reference/ipsl'))
  ),
  tar_target(
    name = 'reference_aoa_mpi',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reference/mpi'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reference/mpi'))
  ),
  tar_target(
    name = 'reference_aoa_mri',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reference/mri'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reference/mri'))
  ),
  tar_target(
    name = 'reference_aoa_ukesm',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/reference/ukesm'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/reference/ukesm'))
  ),
  tar_target(
    name = 'future_aoa_gdal',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/2080s/gfdl'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/2080s/gfdl'))
  ),
  tar_target(
    name = 'future_aoa_ipsl',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/2080s/ipsl'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/2080s/ipsl'))
  ),
  tar_target(
    name = 'future_aoa_mpi',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/2080s/mpi'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/2080s/mpi'))
  ),
  tar_target(
    name = 'future_aoa_mri',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/2080s/mri'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/2080s/mri'))
  ),
  tar_target(
    name = 'future_aoa_ukesm',
    command = analyze_aoa(in_dir=file.path(data_dir, 'predictors/min_five_GCMs/2080s/ukesm'),
                          in_model=file.path(data_dir, 'rf_models/rftuned_step1.qs'),
                          in_data_obs = file.path(data_dir, 'others/model_data.qs'),
                          out_dir=file.path(res_dir, 'aoa/2080s/ukesm'))
  ),
  tar_target(
    name = 'diff_aoa_ref_and_reanalysis',
    command = compute_aoa_dif(in_dir_reanalysis=file.path(res_dir, 'aoa/reanalysis'),
                              in_dir_phase=file.path(res_dir, 'aoa/reference'),
                              out_dir=file.path(res_dir, 'aoa/reference/differences'))
  ),
  tar_target(
    name = 'diff_aoa_2080s_and_reanalysis',
    command = compute_aoa_dif(in_dir_reanalysis=file.path(res_dir, 'aoa/reanalysis'),
                              in_dir_phase=file.path(res_dir, 'aoa/2080s'),
                              out_dir=file.path(res_dir, 'aoa/2080s/differences'))
  ),
  tar_target(
    name = 'create_figure9_ref',
    command = create_fig9_and_S7(in_dir=file.path(res_dir, 'aoa/reference/differences'),
                                 out_dir=file.path(res_dir, 'figures/figure9'))
  ),
  tar_target(
    name = 'create_figure9_2080s',
    command = create_fig9_and_S7(in_dir=file.path(res_dir, 'aoa/2080s/differences'),
                                 out_dir=file.path(res_dir, 'figures/figure9'))
  ),
  tar_target(
    name = 'create_figureS7_ref',
    command = create_fig9_and_S7(in_dir=file.path(res_dir, 'aoa/reference/differences'),
                                 out_dir=file.path(res_dir, 'figures/figureS7'),
                                 figure9=FALSE,
                                 phase='reference')
  ),
  tar_target(
    name = 'create_figureS7_2080s',
    command = create_fig9_and_S7(in_dir=file.path(res_dir, 'aoa/2080s/differences'),
                                 out_dir=file.path(res_dir, 'figures/figureS7'),
                                 figure9=FALSE,
                                 phase='2080s')
  )
)

# --------------------------- run the pipeline ---------------------
list(
  # prepare_hr_predictors, # if the monthly HR streamflow is available, active this
  # prepare_lr_predictors, # if the LR predictors are available, active this
  # run_rf_gcms,           # if the two previous plans have been executed, active this
  create_final_products    # the link for the required data is provided in data avialability statement of the paper
)



















