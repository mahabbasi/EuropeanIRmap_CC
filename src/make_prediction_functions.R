
# ----------------------- Functions for applying the sequantial model to reaches -------------
# Costum functions to import data and the tuned models, and lastely execute the sequential models

#' combine all predictors
#' 
#' @description This function imports and combines all the 23 predictors form different 
#' directories.
#' 
#'  @param path_static (character) path to the static predictors in {fst} format
#'  @param path_LR (character) path to the low-resolution predictors in {fst} format
#'  @param path_HR (character) path to the high-resolution predictors in {fst} format
#'  @param num_mon (numeric) an integer value indicates the number of calendar month.
#'  i.e. 1,2,...,12
#'  @param num_year (numeric) an integer value indicates the number of year, e.g. 1985
#'  @param period_phase (character) a character which indicates the phase of the period.
#'  options: `hist`, `nearfuture`, and `farfuture`
#'    
#' @return a (1,23) data.frame, which represents the values for a specific month/year.
#' 
#'  @export
#'   
import_data <- function(path_static, path_LR,
                        path_HR, num_mon,
                        num_year, period_phase = 'hist'){
  
  # Static predictors -----
  static_pred <- fst::read_fst(path = path_static,
                               columns = c("glacier_fraction", "land_cover", "slope", "drainage_area",
                                           "pot_nat_vegetation", "karst_fraction", "karst_status",
                                           'lka_pc_use', 'ppd_pk_cav', 'ppd_pk_uav', 'ire_pc_cse', 
                                           'ire_pc_use', 'dor_pc_pva',
                                           paste0("ai_", num_mon))) %>%
    as.data.table() %>% 
    rename(., ai = paste0("ai_", num_mon))
  
  upa <- static_pred %>% pull(drainage_area)
  
  # Low resolution predictors ----
  path_qrdif_ql_ratio <-  file.path(path_LR, "gwr_to_runoff_ratio_eu.fst")
  path_wetdays <-  file.path(path_LR, "wetdays_net_eu.fst") 
  
  qrdif_ql_ratio <- fst::read_fst(path = path_qrdif_ql_ratio,
                                  columns = paste0("gwr_", num_mon, "_", num_year)) %>%
    as.data.table() %>% 
    replace(is.na(.), 0) %>% 
    `colnames<-`("gwr_to_runoff_ratio")
  
  wet_day <- fst::read_fst(path = path_wetdays,
                           columns = paste0("wetdays_", num_mon, "_", num_year)) %>%
    as.data.table() %>% 
    `colnames<-`("wet_days")
  # High resolution predictors ------
  path_Q <- file.path(path_HR, "watergap_dis_net_eu.fst")
  path_mean_p3 <-  file.path(path_HR, "q_mean_p3_eu.fst")
  path_mean_p12 <-  file.path(path_HR, "q_mean_p12_eu.fst")
  path_min_p3 <-  file.path(path_HR, "q_min_p3_eu.fst")
  path_min_p12 <-  file.path(path_HR, "q_min_p12_eu.fst")
  
  Q <- fst::read_fst(path = path_Q,
                     columns = c('DRYVER_RIVID', paste0("Q_", num_mon, "_", num_year))) %>%
    as.data.table()
  Q <- unique(Q, by='DRYVER_RIVID')
  Q <- Q[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("Q")
  Q <- Q/upa
  mean_p3m <- fst::read_fst(path = path_mean_p3,
                            columns = c('DRYVER_RIVID', paste0("mean_p3_", num_mon, "_", num_year))) %>% 
    as.data.table()
  mean_p3m <- mean_p3m[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("mean_p3m")
  mean_p3m <- mean_p3m/upa
  mean_p12m <- fst::read_fst(path = path_mean_p12,
                             columns = c('DRYVER_RIVID', paste0("mean_p12_", num_mon, "_", num_year))) %>% 
    as.data.table()
  mean_p12m <- mean_p12m[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("mean_p12m")
  mean_p12m <- mean_p12m/upa
  min_p3m <- fst::read_fst(path = path_min_p3,
                           columns = c('DRYVER_RIVID', paste0("min_p3_", num_mon, "_", num_year))) %>% 
    as.data.table() 
  min_p3m <- min_p3m[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("min_p3m")
  min_p3m <- min_p3m/upa
  min_p12m <- fst::read_fst(path = path_min_p12,
                            columns = c('DRYVER_RIVID', paste0("min_p12_", num_mon, "_", num_year))) %>% 
    as.data.table()
  min_p12m <- min_p12m[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("min_p12m")
  min_p12m <- min_p12m/upa
  if (period_phase == 'nearfuture'){
    path_cv <-  file.path(path_HR, "q_iav_cv_eu_nearfuture.fst")
    path_sd <-  file.path(path_HR, "q_iav_sd_eu_nearfuture.fst")  
    
  }else if (period_phase == 'farfuture') {
    path_cv <-  file.path(path_HR, "q_iav_cv_eu_farfuture.fst")
    path_sd <-  file.path(path_HR, "q_iav_sd_eu_farfuture.fst")  
  }else{
    path_cv <-  file.path(path_HR, "q_iav_cv_eu.fst")
    path_sd <-  file.path(path_HR, "q_iav_sd_eu.fst")  
  }
  
  q_cv <- fst::read_fst(path = path_cv,
                        columns = c('DRYVER_RIVID', paste0("cv_", num_mon))) %>% 
    as.data.table()
  q_cv <- q_cv[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("cv") %>% 
    replace(is.na(.), 5)
  q_sd <- fst::read_fst(path = path_sd,
                        columns = c('DRYVER_RIVID', paste0("sd_", num_mon))) %>% 
    as.data.table()
  q_sd <- q_sd[order(DRYVER_RIVID)] %>% .[,2] %>% 
    `colnames<-`("sd")
  outdt <- cbind(Q, mean_p3m, mean_p12m, min_p3m, min_p12m, q_sd, q_cv,
                 wet_day, qrdif_ql_ratio, static_pred)
}

#' Load the RF sequential models
#' 
#' @param path_model1 (character) path to the saved step one random forest model in {qs} format
#' @param path_model2 (character) path to the saved step two random forest model in {qs} format
#' 
#' @return a list that includes both RF models
#' 
#' @export
#' 
load_models <- function(path_model1, path_model2){
  
  model_step1 <- qs::qread(file = path_model1)
  model_step2 <- qs::qread(file = path_model2)
  out <- list(model_step1 = model_step1,
              model_step2 = model_step2)
}

#' Run the models for reaches
#' 
#' @param model_step1 (qs) a {qs} file includes the framework of the step one RF model
#' @param model_step2 (qs) a {qs} file includes the framework of the step two RF model
#' @param data_dt (data.frame) a data.frame (1,23) resulted from `import_data`
#' @param threshold (numeric) a value that classify the perennial and non-perennial condition
#' based on the output probability of the RF model. `default =` 0.5.
#' 
#' @return the function gives a list of 1: the results of step one RF model; 2: the results of
#' the step two RF model; and 3: the IDs of the reaches predicted non-perennial by step two.
#' 
#' @export
#' 
execute_models_nets <- function(model_step1, model_step2,
                                data_dt, threshold = 0.5){
  
  res_step1 <- model_step1$predict_newdata(data_dt)
  res_step1 <- res_step1$set_threshold(threshold) %>% 
    as.data.table()
  rows_id_step2 <- res_step1[response != 0, row_ids]
  
  res_step2 <- model_step2$predict_newdata(data_dt %>% 
                                             slice(rows_id_step2)) %>% 
    as.data.table()
  out <- list(res_step1 = res_step1,
              res_step2 = res_step2,
              rows_id_step2 = rows_id_step2)
  return(out)
}

#' Executing the predictive models over the given period
#' 
#' @description This function combines all the previous functions in this script to 
#' implement the predictive models over the reaches for the given period
#' 
#' @param path_model1 (character) path to the saved step one random forest model in {qs} format
#' @param path_model2 (character) path to the saved step two random forest model in {qs} format
#' @param path_static (character) path to the static predictors in {fst} format
#' @param path_LR (character) path to the low-resolution predictors in {fst} format
#' @param path_HR (character) path to the high-resolution predictors in {fst} format
#' @param start_year (numeric) an integer number that presents the start year of the interested period
#' @param end_year (numeric) an integer number that presents the end year of the interested period
#' @param outdir (character) path to the directory that the outputs are saved
#' @param period_phase (character) a character which indicates the phase of the period.
#'  options: `hist`, `nearfuture`, and `farfuture`
#'  
#'  @return The function returns `NULL`. But, two files in {fst} format are saved in the
#'  `outdir`. The files contain the streamflow intermittence class and the streamflow 
#'  intermittence probability for the reaches.
#'  
#'  @export
#'  
runmodels_over_period <- function(path_model1,path_model2,path_static,path_LR,path_HR,
                                  start_year=1981, end_year=2019, outdir, period_phase = 'hist'){
  
  seq_models <- load_models(path_model1 = path_model1,
                            path_model2 = path_model2)
  
  
  # load the DRYvER_id of the reaches over Europe
  id <- fst::read_fst(path = path_static,
                      columns = "DRYVER_RIVID") %>% 
    as.data.table()
  
  # create an empty matrix
  years <- start_year:end_year
  res_nets_mat <- matrix(NA, nrow = dim(id)[1], ncol = (length(years)*12+1))
  res_nets_mat[,1] <- id$DRYVER_RIVID
  
  res_nets_mat_prob <- matrix(NA, nrow = dim(id)[1], ncol = (length(years)*12+1))
  res_nets_mat_prob[,1] <- id$DRYVER_RIVID
  # run sequential model for the whole historical period
  for(i in seq_along(years)){
    cat("The number of no-flow days for reaches in", as.character(years[i]),
        "is undergoing...\n")
    for(j in 1:12){
      
      data_dt <- import_data(path_static = path_static,
                             path_LR = path_LR,
                             path_HR = path_HR,
                             num_mon = j,
                             num_year = years[i],
                             period_phase = period_phase)
      
      results_net <- execute_models_nets(model_step1 = seq_models$model_step1,
                                         model_step2 = seq_models$model_step2,
                                         data_dt = data_dt)
      res_nets_mat[,((i-1)*12+j+1)] <- 
        as.numeric(levels(results_net$res_step1$response))[results_net$res_step1$response]
      res_nets_mat[results_net$rows_id_step2, ((i-1)*12+j+1)] <- 
        as.numeric(levels(results_net$res_step2$response))[results_net$res_step2$response]
      #only the probability of class 0-perennial for step1
      res_nets_mat_prob[,((i-1)*12+j+1)] <- round(results_net$res_step1$prob.0, digits = 2)
    }
    
  }
  
  date_colnames <- seq.Date(as.Date(paste0(start_year,'-01-01')),
                            as.Date(as.Date(paste0(end_year,'-12-01'))),
                            'month') %>% 
    as.character()
  
  res_nets_mat_dt <- res_nets_mat %>%
    as.data.table() %>%
    `colnames<-`(c('DRYVER_RIVID', date_colnames)) %>%
    as.data.table()
  
  res_nets_mat_prob <- res_nets_mat_prob %>%
    as.data.table() %>%
    `colnames<-`(c('DRYVER_RIVID', date_colnames)) %>%
    as.data.table()
  
  if (!dir.exists(outdir)){
    dir.create(outdir, recursive = TRUE)
  }
  # write the data tables
  fst::write_fst(res_nets_mat_dt, 
                 path = file.path(outdir,paste0("results_nets_class_",period_phase, '.fst')))
  
  
  fst::write_fst(res_nets_mat_prob, 
                 path = file.path(outdir,paste0("results_nets_prob_",period_phase, '.fst')))
  
  return(NULL)
}

