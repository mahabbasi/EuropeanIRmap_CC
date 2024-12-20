### ------------------ Functions to do the analysis and create the final shapefiles --------------------
## -------------- Compute the LR total runoff in the future for the figure 2a;2b-------------
#' @title Compute the changes in the LR total runoff figure 2a;2b
#' 
#' @description This function computes the median changes of the LR total runoff HR the five
#' GCMs in the 2080s under RCP8.5 for figure 2a;2b.
#' 
#' @param `in_dir` the path to the LR total runoff of the five GCMs for the reference period and 
#' the future period.
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_countries_shp_dir` the path to the shapefile of the European countries required to mask
#' the final output.
#' 
#' @export
#' 
compute_qtot_fig2ab <- function(in_dir, out_dir, eu_countries_shp_dir) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  file_lists <- list.files(path = in_dir, pattern = '.nc', full.names = TRUE)
  file_lists_hist <- grep('historical', file_lists, value=TRUE)
  file_lists_rcp8.5 <- grep('ssp585', file_lists, value=TRUE)
  
  # read shapefile of the European countries
  shap_countries <- terra::vect(file.path(eu_countries_shp_dir,
                                          'europe_countires_final.shp'))
  # compute the median of the LR qtot
  perchange_rasters <- lapply(seq_along(file_lists_hist), function(i){
    
    print(i)
    # read the raster for future and historical between the specific periods
    inp_hist <- terra::rast(file_lists_hist[i])[[1621:1980]]
    inp_rcp8.5 <- terra::rast(file_lists_rcp8.5[i])[[673:1032]]
    # compute the average of grid cell over the specific periods
    mean_hist <- terra::mean(inp_hist)
    mean_rcp8.5 <- terra::mean(inp_rcp8.5)
    # compute the percentage change
    change_d <- (mean_rcp8.5 - mean_hist)/mean_hist * 100
    
  }) %>% 
    do.call(c,.)
  
  median_changes <- terra::median(perchange_rasters)
  median_changes_masked <- median_changes %>%
    terra::crop(., shap_countries) %>%
    terra::mask(., shap_countries)
  
  terra::writeRaster(median_changes_masked, file.path(out_dir,
                                                      'median_qtot_pc_fig2a.tif'), overwrite=TRUE)
  
  ## Signal to noise ratio ------------
  r_median <- terra::app(perchange_rasters, median, na.rm = TRUE)
  r_range <- app(perchange_rasters, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  snr <- r_median/r_range
  
  snr_masked <- snr %>%
    terra::crop(., shap_countries) %>%
    terra::mask(., shap_countries)
  
  terra::writeRaster(snr_masked, file.path(out_dir, 
                                           'SNR_qtot_pc_fig2b.tif'), overwrite=TRUE)
  
  return(NULL)
  
}
## -------------- Compute the HR streamflow of WaterGAP in the future 2c;2d-------------
#' @title Compute the changes HR streamflow figure 2c;2d
#' 
#' @description This function computes the median changes of the HR streamflow of rhe five
#' GCMs in the 2080s under RCP8.5 for figure 2c;2d.
#' 
#' @param `in_dir` the path to the HR streamflow of the five GCMs for the reference period and 
#' the future period.
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' 
#' @export
#' 
compute_hr_pc_fig2cd <- function(in_dir, out_dir, eu_net_shp_dir) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  
  GCM_names <- c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
  
  m <- lapply(seq_along(GCM_names), function(i){
    print(i)
    input_hist_dir <- file.path(in_dir, GCM_names[i], 'historical/watergap_dis_net_eu.fst')
    input_rcp8.5_dir <- file.path(in_dir, GCM_names[i], 'future/ssp585/watergap_dis_net_eu.fst')
    
    # read the fst files of the HR streamflow
    hist_fst <- fst::read_fst(input_hist_dir) %>% as.data.table()
    rcp8.5_fst <- fst::read_fst(input_rcp8.5_dir) %>% as.data.table()
    
    hist_fst[,hist_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = 2:ncol(hist_fst)]
    rcp8.5_fst[,future_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = 2:ncol(rcp8.5_fst)]
    com_dt <- data.table(hist_fst[,.(DRYVER_RIVID, hist_mean)], rcp8.5_fst[,.(future_mean)])
    com_dt[, per_change := (future_mean - hist_mean)/hist_mean *100]
    
    if (i ==1){
      com_dt[,.(DRYVER_RIVID, per_change)]
    } else{
      com_dt[,.(per_change)]
    }
  }) %>% 
    do.call('cbind',.)
  
  
  m[, median := apply(.SD, 1, median, na.rm=TRUE), .SDcols=2:ncol(m)]
  m[, max := apply(.SD, 1, max, na.rm=TRUE), .SDcols=2:ncol(m)]
  m[, min := apply(.SD, 1, min, na.rm=TRUE), .SDcols=2:ncol(m)]
  snr_streamflow <- m[,SNR:=median/(max-min)]
  setnames(m, 'DRYVER_RIVID', 'DRYVER_RIV')
  hist_joined_st <- left_join(reach_shp, m[,.(DRYVER_RIV, median, min, max, SNR)], by='DRYVER_RIV')
  sf::write_sf(hist_joined_st,
               dsn = file.path(out_dir, 'pc_SNR_hr_streamflow_rcp8.5_far_fig2cd.shp'))
  
  return(NULL)
}

## -------------- Compute the HR low flow of WaterGAP in the future figure 2e;2f -------------
#' @title Compute the changes of the HR low flow Q90 figure 2e;2f
#' 
#' @description This function computes the median changes of the HR streamflow of rhe five
#' GCMs in the 2080s under RCP8.5 for figure 2c;2d.
#' 
#' @param `in_dir` the path to the HR streamflow of the five GCMs for the reference period and 
#' the future period.
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' 
#' @export
#' 
compute_lowflow_pc_fig2ef <- function(in_dir, out_dir, eu_net_shp_dir) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  
  GCM_names <- c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
  
  m <- lapply(seq_along(GCM_names), function(i){
    print(i)
    input_hist_dir <- file.path(in_dir, GCM_names[i], 'historical/watergap_dis_net_eu.fst')
    input_rcp8.5_dir <- file.path(in_dir, GCM_names[i], 'future/ssp585/watergap_dis_net_eu.fst')
    
    # read the fst files of the HR streamflow
    hist_fst <- fst::read_fst(input_hist_dir) %>% as.data.table()
    rcp8.5_fst <- fst::read_fst(input_rcp8.5_dir) %>% as.data.table()
    
    quantile_10_hist <- hist_fst[,apply(.SD, 1, function(row) quantile(row, probs = 0.1)),
                                 .SDcols=-'DRYVER_RIVID']
    rcp8.5_fst[,'DRYVER_RIVID':=hist_fst[,DRYVER_RIVID]]
    quantile_10_rcp8.5 <- rcp8.5_fst[,apply(.SD, 1, function(row) quantile(row, probs = 0.1)),
                                           .SDcols=-'DRYVER_RIVID']
    
    out_dif <- data.table(DRYVER_RIVID=hist_fst[,DRYVER_RIVID],
                          pc_p10=(quantile_10_rcp8.5-quantile_10_hist)/quantile_10_hist * 100)
    
    if (i ==1){
      out_dif[,.(DRYVER_RIVID, pc_p10)]
    } else{
      out_dif[,.(pc_p10)]
    }
    
  }) %>% 
    do.call('cbind',.)
  
  m[, median := apply(.SD, 1, median, na.rm=TRUE), .SDcols=2:ncol(m)]
  m[, max := apply(.SD, 1, max, na.rm=TRUE), .SDcols=2:ncol(m)]
  m[, min := apply(.SD, 1, min, na.rm=TRUE), .SDcols=2:ncol(m)]
  setnames(m, 'DRYVER_RIVID', 'DRYVER_RIV')
  hist_joined_st <- dplyr::left_join(reach_shp, m[,.(DRYVER_RIV, median)], by='DRYVER_RIV')
  sf::write_sf(hist_joined_st,
               dsn = file.path(out_dir, 'change_p10_streamflow_rcp8.5_far_fig2e.shp'))
  
  ## signal-to-noise ratio --------
  m[,SNR:=median/(max-min)]
  hist_joined_st <- dplyr::left_join(reach_shp, m[,.(DRYVER_RIV, SNR)], by='DRYVER_RIV')
  
  hist_joined_st <- hist_joined_st %>% 
    dplyr::mutate(., SNR =ifelse(is.na(SNR), 1000, SNR)) %>% 
    dplyr::mutate(., SNR =ifelse(SNR < -10, -1000, SNR))
  
  sf::write_sf(hist_joined_st,
               dsn = file.path(out_dir, 'SNR_pc_p10_streamflow_rcp8.5_far_fig2f.shp'))
  
  
  return(NULL)
}

## ----------- Compute the changes of the median non-perennial reach-months figure 3------------
#' @title Compute the changes of median non-perennial reaches-month compared to 1985-2014
#' 
#' @description This function computes the median non-perennial reaches-month of the five
#' GCMs in the reference period for figure 3a; the median of the changes for the 2050s and 
#' the 2080s under RCP8.5 (figures 3b, 3c) and the corresponding Signal-to-Noise Ratio (figures 3e, 3f).
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and two future periods under RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' @param rcp_name (character) this presents the name of the RCP, and could be either RCP8.5 or RCP2.6.
#' 
#' @export
#' 
compute_np_changes_fig3 <- function(in_dir, out_dir, eu_net_shp_dir, rcp_name='RCP8.5'){
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  # create the inputs paths 
  if (rcp_name == 'RCP8.5') {
    in_hist_path <- file.path(in_dir, 'historical') 
    in_nearfuture_path <- file.path(in_dir, 'future/ssp585/near_future') 
    in_farfuture_path <- file.path(in_dir, 'future/ssp585/far_future')
  } else {
    in_hist_path <- file.path(in_dir, 'historical') 
    in_nearfuture_path <- file.path(in_dir, 'future/ssp126/near_future') 
    in_farfuture_path <- file.path(in_dir, 'future/ssp126/far_future')
  }

  
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  # define a costumized function to compute the number of intermittent months
  compute_number_inter_mon <- function(in_list){
    lapply(seq_along(in_list), function(i){
      in_dt <- fst::read_fst(in_list[[i]]) %>%
        as.data.table()
      col_counts_ratio <- rowSums(in_dt[, .SD[, -1, with = FALSE] > 0])
    }) %>% 
      do.call('cbind',.) %>% 
      `colnames<-`(c('gfdl', 'ipsl', 'mpi',
                     'mri', 'ukesm')) %>% 
      as.data.table()
  }
  
  # analyze the median of GCMs for reference period
  hist_path_list <- list.files(path = in_hist_path, pattern = '.fst', full.names = TRUE)
  hist_number_inter_mon_dt <- compute_number_inter_mon(hist_path_list)
  hist_number_inter_mon_dt[,median_row := apply(.SD, 1, median),
                           .SDcols=names(hist_number_inter_mon_dt)]
  hist_number_inter_mon_dt[,median_row := round((median_row/360) * 100, digits = 2)]
  
  in_reach_ids <- fst::read_fst(hist_path_list[1]) %>% 
    as.data.table() %>% .[,.(DRYVER_RIVID)]
  
  hist_number_inter_mon_dt[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  hist_joined_st <- left_join(reach_shp, hist_number_inter_mon_dt, by='DRYVER_RIV')
  sf::write_sf(hist_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_reference_fig3a.shp'))
  # analyze the changes of percentage point for the 2050s 
  nearfuture_path_list <- list.files(path = in_nearfuture_path,
                                     pattern = '.fst', full.names = TRUE)
  
  nearfuture_number_inter_mon_dt <- compute_number_inter_mon(nearfuture_path_list)
  nearfuture_inter_pp <- (nearfuture_number_inter_mon_dt - hist_number_inter_mon_dt[,1:5])/360 * 100
  nearfuture_inter_pp[,median_row := apply(.SD, 1, median), .SDcols=names(nearfuture_inter_pp)]
  nearfuture_inter_pp[, range_row := do.call(pmax, .SD) - do.call(pmin, .SD),
                      .SDcols = names(nearfuture_inter_pp)[1:5]]
  
  nearfuture_inter_pp[,S2N:=median_row/range_row]
  
  nearfuture_inter_pp[, (names(nearfuture_inter_pp)) := lapply(.SD, function(x) ifelse(is.nan(x), 0, x)),
                      .SDcols = names(nearfuture_inter_pp)]
  nearfuture_inter_pp[, (names(nearfuture_inter_pp)) := lapply(.SD, function(x) ifelse(is.infinite(x), 0, x)),
                      .SDcols = names(nearfuture_inter_pp)]
  nearfuture_inter_pp <- round(nearfuture_inter_pp, digits = 2)
  
  nearfuture_inter_pp[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  nearfuture_joined_medianpp_st <- left_join(reach_shp, nearfuture_inter_pp[,.SD, .SDcols=-'S2N'], by='DRYVER_RIV')
  nearfuture_joined_s2n_st <- left_join(reach_shp, nearfuture_inter_pp, by='DRYVER_RIV')
  
  if (rcp_name == 'RCP8.5') {
    sf::write_sf(nearfuture_joined_medianpp_st,
                 dsn = file.path(out_dir, 'median_inter_2050s_pp_fig3b.shp'))
    
    sf::write_sf(nearfuture_joined_s2n_st,
                 dsn = file.path(out_dir, 'SNR_inter_2050s_pp_fig3d.shp'))
  } else {
    sf::write_sf(nearfuture_joined_medianpp_st,
                 dsn = file.path(out_dir, 'median_inter_2050s_pp_figS2b.shp'))
    
    sf::write_sf(nearfuture_joined_s2n_st,
                 dsn = file.path(out_dir, 'SNR_inter_2050s_pp_figS2d.shp'))
  }

  
  # analyze the changes of percentage point for the 2080s 
  
  farfuture_path_list <- list.files(path = in_farfuture_path,
                                    pattern = '.fst', full.names = TRUE)
  
  farfuture_number_inter_mon_dt <- compute_number_inter_mon(farfuture_path_list)
  farfuture_inter_pp <- (farfuture_number_inter_mon_dt - hist_number_inter_mon_dt[,1:5])/360 * 100
  farfuture_inter_pp[,median_row := apply(.SD, 1, median), .SDcols=names(farfuture_inter_pp)]
  farfuture_inter_pp[, range_row := do.call(pmax, .SD) - do.call(pmin, .SD),
                     .SDcols = names(farfuture_inter_pp)[1:5]]
  
  farfuture_inter_pp[,S2N:=median_row/range_row]
  
  farfuture_inter_pp[, (names(farfuture_inter_pp)) := lapply(.SD, function(x) ifelse(is.nan(x), 0, x)),
                     .SDcols = names(farfuture_inter_pp)]
  farfuture_inter_pp[, (names(farfuture_inter_pp)) := lapply(.SD, function(x) ifelse(is.infinite(x), 0, x)),
                     .SDcols = names(farfuture_inter_pp)]
  farfuture_inter_pp <- round(farfuture_inter_pp, digits = 2)
  
  farfuture_inter_pp[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  farfuture_joined_medianpp_st <- left_join(reach_shp, farfuture_inter_pp[,.SD, .SDcols=-'S2N'], by='DRYVER_RIV')
  farfuture_joined_s2n_st <- left_join(reach_shp, farfuture_inter_pp, by='DRYVER_RIV')
  
  if (rcp_name == 'RCP8.5') {
    sf::write_sf(farfuture_joined_medianpp_st,
                 dsn = file.path(out_dir, 'median_inter_2080s_pp_fig3c.shp'))
    
    sf::write_sf(farfuture_joined_s2n_st,
                 dsn = file.path(out_dir, 'SNR_inter_2080s_pp_fig3e.shp'))
  } else {
    sf::write_sf(farfuture_joined_medianpp_st,
                 dsn = file.path(out_dir, 'median_inter_2080s_pp_figS2c.shp'))
    
    sf::write_sf(farfuture_joined_s2n_st,
                 dsn = file.path(out_dir, 'SNR_inter_2080s_pp_figS2e.shp'))
  }

  
  return(NULL)
}

## ---------- Plot the CDF of the non-perennial figure 4 -----------------------
#' @title Plot the CDF of the non-perennial for different Climate zones
#' 
#' @description This function generates a plot that depicts the the monthly percentage
#' of non-perennial reaches for the five GCMs calculated over the three 360-months periods in
#' Europe and the six European climate zones.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and two future periods (the 2050s, the 2080s) under RCP2.6 and RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `in_shp_climcode_path` the path to the shapefile of the climate zone of the 
#' European network.
#' @param `climate_zone_png_dir` the path to the png file of the six climate zones over Europe.
#' 
#' @export
#' 
ggcdf_np_climzones_fig4 <- function(in_dir, in_shp_climcode_path, out_dir,
                                         climate_zone_png_dir){
  
  if(!file.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  in_hist_dt <- list.files(file.path(in_dir, 'historical'),
                           pattern = '.fst', full.names = TRUE)
  in_ssp126_near_dt <- list.files(file.path(in_dir, 'future/ssp126/near_future'),
                                  pattern = '.fst', full.names = TRUE)
  in_ssp126_far_dt <- list.files(file.path(in_dir, 'future/ssp126/far_future'),
                                 pattern = '.fst', full.names = TRUE)
  in_ssp585_near_dt <- list.files(file.path(in_dir, 'future/ssp585/near_future'),
                                  pattern = '.fst', full.names = TRUE)
  in_ssp585_far_dt <- list.files(file.path(in_dir, 'future/ssp585/far_future'),
                                 pattern = '.fst', full.names = TRUE)
  
  in_shp_climcode <- read_sf(in_shp_climcode_path)
  climcode_dt <- in_shp_climcode %>%
    st_drop_geometry() %>%
    as.data.table() %>% 
    .[,.(DRYVER_RIV,clim_code)]
  
  # combine all the scenarios for the individual GCM 
  in_list_gfdl <- list(in_hist_dt[1], in_ssp126_near_dt[1], in_ssp126_far_dt[1],
                       in_ssp585_near_dt[1], in_ssp585_far_dt[1])
  in_list_ipsl <- list(in_hist_dt[2], in_ssp126_near_dt[2], in_ssp126_far_dt[2],
                       in_ssp585_near_dt[2], in_ssp585_far_dt[2])
  in_list_mpi <- list(in_hist_dt[3], in_ssp126_near_dt[3], in_ssp126_far_dt[3],
                      in_ssp585_near_dt[3], in_ssp585_far_dt[3])
  
  in_list_mri <- list(in_hist_dt[4], in_ssp126_near_dt[4], in_ssp126_far_dt[4],
                      in_ssp585_near_dt[4], in_ssp585_far_dt[4])
  
  in_list_ukesm <- list(in_hist_dt[5], in_ssp126_near_dt[5], in_ssp126_far_dt[5],
                        in_ssp585_near_dt[5], in_ssp585_far_dt[5])
  
  inner_function <- function(in_list, model_name='gfdl', climcode_dt=climcode_dt,
                             eu_overall=TRUE, climate_zone=5){
    
    if (eu_overall) {
      # compute streamflow intermittence colwise 
      intermittent_colwise <- lapply(seq_along(in_list), function(i){
        in_dt <- fst::read_fst(in_list[[i]]) %>%
          as.data.table()
        row_counts <- colSums(in_dt[, lapply(.SD, function(x) x > 0)])/dim(in_dt)[1] * 100
      }) %>% 
        do.call('cbind',.) %>% 
        `colnames<-`(paste0(c('hist', 'ssp126_near', 'ssp126_far',
                              'ssp585_near', 'ssp585_far'),'_', model_name)) %>% 
        as.data.table() %>% .[-1]
      
    } else {
      cat('Climate zone', climate_zone, 'is undergoing....\n')
      intermittent_colwise <- lapply(seq_along(in_list), function(i){
        in_dt <- fst::read_fst(in_list[[i]]) %>%
          as.data.table()
        setnames(in_dt, 'DRYVER_RIVID', 'DRYVER_RIV')
        in_dt_withclim <- left_join(in_dt, climcode_dt, by='DRYVER_RIV')
        class_num <- in_dt_withclim[clim_code==climate_zone,.N]
        row_counts <- colSums(in_dt_withclim[clim_code==climate_zone,
                                             lapply(.SD, function(x) x > 0)])/class_num * 100
      }) %>% 
        do.call('cbind',.) %>% 
        `colnames<-`(paste0(c('hist', 'ssp126_near', 'ssp126_far',
                              'ssp585_near', 'ssp585_far'),'_', model_name)) %>% 
        as.data.table() %>% .[-c(1, 362)]
      
    }
  }

  ## plot the CDF of the five GCMs and three different scenarios 
  ggcdf_plot <- function(inter_gfdl,inter_ipsl,inter_mpi,inter_mri,inter_ukesm,
                         region_name= 'Europe', xlab=NULL, ylab=NULL,xloc=0.7, yloc=0.2,
                         col_region='black', legend = NULL){
    # Combining data tables and converting to long format
    combined_data <- cbind(inter_gfdl, inter_ipsl,inter_mpi, inter_mri, inter_ukesm)
    #plot the cdf of intermittent reaches over time
    theme_set(theme_bw(16))
    data_long <- melt(combined_data, measure.vars = names(combined_data),
                      variable.name = "scenario", value.name = "value")
    
    color_values <- c("hist_gfdl" = "black","hist_ipsl" = "black","hist_mpi" = "black",
                      "hist_mri" = "black","hist_ukesm" = "black",
                      "ssp126_near_gfdl" = "#92c5de","ssp126_near_ipsl" = "#92c5de","ssp126_near_mpi" = "#92c5de",
                      "ssp126_near_mri" = "#92c5de","ssp126_near_ukesm" = "#92c5de",
                      "ssp126_far_gfdl" = "#0571b0","ssp126_far_ipsl" = "#0571b0","ssp126_far_mpi" = "#0571b0",
                      "ssp126_far_mri" = "#0571b0","ssp126_far_ukesm" = "#0571b0",
                      "ssp585_near_gfdl" = "#f4a582","ssp585_near_ipsl" = "#f4a582","ssp585_near_mpi" = "#f4a582",
                      "ssp585_near_mri" = "#f4a582","ssp585_near_ukesm" = "#f4a582",
                      "ssp585_far_gfdl" = "#ca0020","ssp585_far_ipsl" = "#ca0020","ssp585_far_mpi" = "#ca0020",
                      "ssp585_far_mri" = "#ca0020","ssp585_far_ukesm" = "#ca0020")
    
    linetype_values <- c("hist_gfdl" = "solid","hist_ipsl" = "solid","hist_mpi" = "solid",
                         "hist_mri" = "solid","hist_ukesm" = "solid",
                         "ssp126_near_gfdl" = "solid","ssp126_near_ipsl" = "solid","ssp126_near_mpi" = "solid",
                         "ssp126_near_mri" = "solid","ssp126_near_ukesm" = "solid",
                         "ssp126_far_gfdl" = "dashed","ssp126_far_ipsl" = "dashed","ssp126_far_mpi" = "dashed",
                         "ssp126_far_mri" = "dashed","ssp126_far_ukesm" = "dashed",
                         "ssp585_near_gfdl" = "dashed","ssp585_near_ipsl" = "dashed","ssp585_near_mpi" = "dashed",
                         "ssp585_near_mri" = "dashed","ssp585_near_ukesm" = "dashed",
                         "ssp585_far_gfdl" = "solid","ssp585_far_ipsl" = "solid","ssp585_far_mpi" = "solid",
                         "ssp585_far_mri" = "solid","ssp585_far_ukesm" = "solid")
    
    # Creating CDF plot
    p <- ggplot(data_long, aes(x = value, color = scenario)) +
      stat_ecdf(geom = "step", pad = FALSE, linewidth=0.5) +
      labs(x = xlab,
           y = ylab) +
      scale_color_manual(values = color_values) +
      # scale_linetype_manual(values = linetype_values) +
      annotation_custom(grid::textGrob(region_name, xloc, yloc, gp = grid::gpar(fontsize=14, fontface="bold", col=col_region)))+
      theme(panel.grid.minor = element_blank(),axis.text=element_text(colour="black"),
            panel.grid.major = element_line(colour = "black", linewidth=0.1, linetype = 'dashed' ))+
      guides(linetype='none', color='none')
    
    return(p)
  }
  
  #for all the European reaches
  inter_gfdl_a <- inner_function(in_list = in_list_gfdl, model_name='gfdl')
  inter_ipsl_a <- inner_function(in_list = in_list_ipsl, model_name='ipsl')
  inter_mpi_a <- inner_function(in_list = in_list_mpi, model_name='mpi')
  inter_mri_a <- inner_function(in_list = in_list_mri, model_name='mri')
  inter_ukesm_a <- inner_function(in_list = in_list_ukesm, model_name='ukesm')

  p2a <- ggcdf_plot(inter_gfdl_a, inter_ipsl_a, inter_mpi_a, inter_mri_a,
                    inter_ukesm_a, region_name= 'Europe', ylab = '',
                    xloc = 0.2, yloc = 0.9)
  #for all the reaches within climate_zone=5
  inter_gfdl_b <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 5)
  inter_ipsl_b <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 5)
  inter_mpi_b <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 5)
  inter_mri_b <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 5)
  inter_ukesm_b <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 5)
  p2b <- ggcdf_plot(inter_gfdl_b, inter_ipsl_b, inter_mpi_b, inter_mri_b, inter_ukesm_b,
                    region_name= 'mediterranean/\nsemi-arid', col_region='#A80000')
  
  #for all the reaches within climate_zone=14
  inter_gfdl_c <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 14)
  inter_ipsl_c <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 14)
  inter_mpi_c <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 14)
  inter_mri_c <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 14)
  inter_ukesm_c <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 14)
  p2c <- ggcdf_plot(inter_gfdl_c, inter_ipsl_c, inter_mpi_c, inter_mri_c, inter_ukesm_c,
                    region_name= 'humid subtropical', ylab = '', col_region='#FFA77F')
  
  #for all the reaches within climate_zone=15
  inter_gfdl_d <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 15)
  inter_ipsl_d <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 15)
  inter_mpi_d <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 15)
  inter_mri_d <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 15)
  inter_ukesm_d <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 15)
  p2d <- ggcdf_plot(inter_gfdl_d, inter_ipsl_d, inter_mpi_d, inter_mri_d, inter_ukesm_d,
                    region_name= 'temperate oceanic', col_region='#D1FF73')
  
  #for all the reaches within climate_zone=25
  inter_gfdl_e <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 25)
  inter_ipsl_e <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 25)
  inter_mpi_e <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 25)
  inter_mri_e <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 25)
  inter_ukesm_e <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 25)
  p2e <- ggcdf_plot(inter_gfdl_e, inter_ipsl_e, inter_mpi_e, inter_mri_e, inter_ukesm_e,
                    region_name= 'humid continental', ylab = '', col_region='#70A800')
  
  #for all the reaches within climate_zone=27
  inter_gfdl_f <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 27)
  inter_ipsl_f <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 27)
  inter_mpi_f <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 27)
  inter_mri_f <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 27)
  inter_ukesm_f <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 27)
  p2f <- ggcdf_plot(inter_gfdl_f, inter_ipsl_f, inter_mpi_f, inter_mri_f, inter_ukesm_f,
                    region_name= 'sub-arctic', xlab = "", col_region='#73DFFF')
  
  #for all the reaches within climate_zone=30
  inter_gfdl_g <- inner_function(in_list = in_list_gfdl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='gfdl', climate_zone = 30)
  inter_ipsl_g <- inner_function(in_list = in_list_ipsl, climcode_dt = climcode_dt,
                                 eu_overall = FALSE,
                                 model_name='ipsl', climate_zone = 30)
  inter_mpi_g <- inner_function(in_list = in_list_mpi, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mpi', climate_zone = 30)
  inter_mri_g <- inner_function(in_list = in_list_mri, climcode_dt = climcode_dt,
                                eu_overall = FALSE,
                                model_name='mri', climate_zone = 30)
  inter_ukesm_g <- inner_function(in_list = in_list_ukesm, climcode_dt = climcode_dt,
                                  eu_overall = FALSE,
                                  model_name='ukesm', climate_zone = 30)
  
  p2g <- ggcdf_plot(inter_gfdl_g, inter_ipsl_g, inter_mpi_g, inter_mri_g, inter_ukesm_g,
                    region_name= 'polar/alpine', xlab = "", ylab = '', col_region='#005CE6')
  
  # read the png file of the climate zone classes
  map_image <- png::readPNG(file.path(climate_zone_png_dir, "climate_zones_classes.png"))
  
  # combine the subplots and the imported PNG file into one figure 
  map_grob <- rasterGrob(map_image, width = unit(1, "npc"), height = unit(1, "npc"))
  combined_plot <- plot_grid(
    p2a,
    ggdraw(map_grob),p2b,p2c, p2d, p2e, p2f, p2g,
    labels = c("a", "b", 'c', 'd', 'e', 'f', 'g', 'h'),
    label_size = 14, label_fontface = 'plain',
    ncol = 2,
    # rel_widths = c(0.1, 0.1),  # Relative widths of columns
    # rel_heights = c(0.01, 0.01),  # Relative heights of rows (if using nrow)
    align = 'v'  # Align plots horizontally and vertically
    # axis = 'tb'  # Align axes to the top and bottom
  )
  # Add common x and y labels with additional margins
  final_plot <- ggdraw() +
    draw_plot(combined_plot, x = 0, y = 0, width = 1, height = 1) +
    draw_label("Monthly proportion of non-perennial reaches [%]", x = 0.5, y = 0, hjust = 0.5, 
               vjust = 0, size = 16) +
    draw_label("F (x)", x = -0.05, y = 0.5, hjust = 0, vjust = 0.5, angle = 90, size = 16)
  
  # Adjust margins to avoid overlap
  combined_plot_with_margins <- ggdraw(final_plot) +
    theme(plot.margin = margin(20, 20, 40, 40))  # Adjust margins as needed
  # saveRDS(combined_plot_with_margins, file = 'combined_plot_with_margins.rds')
  ggsave(combined_plot_with_margins, filename = file.path(out_dir, 'CDF_non_perennial_reach_all_fig4.png'),
         units = 'in', dpi = 200, width = 8, height = 8)
  
  return(NULL)
}

## --------- Seasonality of streamflow intermittence figure 5--------------
#' @title Plot the seasonality of streamflow intermittence
#' 
#' @description This function generates a plot that depicts Percentage of European reaches 
#' that are non-perennial in each calendar month of the reference and future periods,
#' distinguishing RCP2.6 and RCP8.5.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and two future periods (the 2050s, the 2080s) under RCP2.6 and RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `in_shp_climcode_path` the path to the shapefile of the climate zone of the 
#' European network.
#' @param `climate_zone_png_dir` the path to the png file of the six climate zones over Europe.
#' 
#' @export
#' 
ggline_seasonality_fig5 <- function(in_dir, out_dir){
  
  if(!file.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  # inner function that computes the statistics of intermittency for each reach individually
  compute_statistics_inter_colwise <- function(in_path, col_names){
    
    path_list <- list.files(path = in_path, pattern = '.fst', full.names = TRUE)
    date_names <- fst::read_fst(path_list[1]) %>% 
      names() %>% .[-1]
    
    # compute streamflow intermittence colwise 
    intermittent_colwise <- lapply(seq_along(path_list), function(i){
      in_dt <- fst::read_fst(path_list[[i]]) %>%
        as.data.table()
      row_counts <- colSums(in_dt[, lapply(.SD, function(x) x > 0)])/dim(in_dt)[1] * 100
    }) %>% 
      do.call('rbind',.) %>% .[,-1] %>% 
      `colnames<-`(date_names)  %>% 
      as.data.table()
    intermittent_colwise[,models:=c('gfdl', 'ipsl', 'mpi',
                                    'mri', 'ukesm')]
    melted_dt <- melt(intermittent_colwise, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    out_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    setnames(out_dt, names(out_dt),
             c('month', col_names))
    return(out_dt)
  }
  
  in_hist_path <- file.path(in_dir,'historical')
  in_ssp126_nearfuture_path <- file.path(in_dir,'future/ssp126/near_future')
  in_ssp126_farfuture_path <- file.path(in_dir,'future/ssp126/far_future')
  in_ssp585_nearfuture_path <- file.path(in_dir,'future/ssp585/near_future')
  in_ssp585_farfuture_path <- file.path(in_dir,'future/ssp585/far_future')
  
  hist_dt <- compute_statistics_inter_colwise(in_path = in_hist_path,
                                              col_names = c('y1', 'y1hi', 'y1lo'))
  ssp126_nearfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp126_nearfuture_path,
                                                           col_names = c('y2', 'y2hi', 'y2lo'))
  ssp126_farfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp126_farfuture_path,
                                                          col_names = c('y3', 'y3hi', 'y3lo'))
  ssp585_nearfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp585_nearfuture_path,
                                                           col_names = c('y4', 'y4hi', 'y4lo'))
  ssp585_farfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp585_farfuture_path,
                                                          col_names = c('y5', 'y5hi', 'y5lo'))
  
  combined_dt <- cbind(hist_dt, ssp126_nearfuture_dt[,-1],
                       ssp126_farfuture_dt[,-1], ssp585_nearfuture_dt[,-1],
                       ssp585_farfuture_dt[,-1])
  
  # Long format data frame:
  longDF = tidyr::pivot_longer(combined_dt , cols=!month , names_to="line" , values_to="y" )
  longDF$fill = NA
  longDF$fill[grep( "1" , longDF$line  )] = "y1"
  longDF$fill[grep( "2" , longDF$line  )] = "y2"
  longDF$fill[grep( "3" , longDF$line  )] = "y3"
  longDF$fill[grep( "4" , longDF$line  )] = "y4"
  longDF$fill[grep( "5" , longDF$line  )] = "y5"
  
  longDF1 <- longDF %>% 
    mutate(line = gsub("\\d", "", line)) %>% 
    tidyr::pivot_wider(id_cols = c(month, fill), names_from = line, values_from = y)
  longDF1 <- longDF1 %>% 
    mutate(fill=case_when(
      fill == 'y1'~'Reference',
      fill == 'y2'~'RCP2.6 2050s',
      fill == 'y3'~'RCP2.6 2080s',
      fill == 'y4'~'RCP8.5 2050s',
      fill == 'y5'~'RCP8.5 2080s'
    ))
  
  # abbrevation of the months for x axis
  abb <- c("J","F","M","A","M","J","J","A","S","O","N","D") 
  # Define the order of legend labels
  desired_order <- c("Reference", "RCP2.6 2050s", "RCP2.6 2080s",
                     'RCP8.5 2050s', 'RCP8.5 2080s')
  
  shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(longDF1$month))), by = 1),
                        max = seq(from = 1.5, to = max(as.numeric(as.factor(longDF1$month))) + 0.5, by = 1),
                        col = c(0,1))
  
  outp <- ggplot(longDF1) +
    geom_rect(data = shading,
              aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf,
                  fill = factor(col), alpha = 0.1)) +
    scale_fill_manual(values = c("white", "gray85"))+
    geom_errorbar(aes(x=month,ymin=ylo, ymax=yhi, color = factor(fill, desired_order)), width=1,
                  position = position_dodge(0.8))+
    geom_line(aes(x=month, y=y, color = factor(fill, desired_order)), linewidth=1) +
    scale_color_manual(values = c('Reference' = "black", 'RCP2.6 2050s' = "#92c5de",
                                  'RCP2.6 2080s'='#0571b0', 'RCP8.5 2050s'='#f4a582',
                                  'RCP8.5 2080s'='#ca0020')) + 
    theme_bw(18)+
    scale_x_continuous(breaks=1:12, labels=abb) +
    theme(legend.position = c(0.35, 0.95),
          legend.justification = c('right', 'top'))+
    guides(color = guide_legend(ncol = 1, title=''),
           fill ='none', alpha='none')+
    labs(y="Fraction of non-perennial reaches [%]" , x="Calendar month")
  
  ggsave(outp, filename = file.path(out_dir, 'all_scenarios_inone_errorbar.png'),
         dpi = 400, units = 'in', height = 8, width = 8)
  
  return(NULL)
  
}

## ---------------- The percentage of non-perennial reach-months figure 6-----------
#' @title Plot the percentage of reach-months in the four intermittence classes
#' 
#' @description This function generates a plot that depicts Percentage of reach-months for Europe
#' and four climate zones in the reference period and the 2080s under RCP8.5.
#' 
#' @param `out_dir` the path to the location that the final plot is stored.
#' 
#' @export
#' 
ggbar_changes_fig6 <- function(out_dir){
  
  
  if(!file.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  
  long_dt <- data.frame(
    class = rep(c('Perennial', '1-5 no-flow days', '6-15 no-flow days',
                  '16-29 no-flow days', '30-31 no-flow days'), times = 14),
    period = rep(c("Reference","Reference","Reference","Reference","Reference",
                   "Future", "Future", "Future", "Future", "Future"), times=7),
    climate_zones = rep(c("Europe", "mediterranean / semi-arid", 'humid subtropical',
                          'temperate oceanic', 'humid continental', 'sub-arctic',
                          'polar/alpine'), each = 10),
    median_value = c(96.4, 0.09,0.1, 0.27, 3.1, 95.1, 0.1, 0.12, 0.29, 4.3, 84.4, 0.19, 0.24, 0.45, 14.3,
                     80.5, 0.20, 0.26, 0.40, 18.7, 95.5, 0.21, 0.19, 0.48, 3.6, 93.7, 0.22, 0.28, 0.47, 5.3,
                     98.7, 0.12, 0.05, 0.28, 0.86, 98.2, 0.13,0.07, 0.29, 1.3, 98.1, 0.06, 0.11, 0.30, 1.5,
                     97.4, 0.07,0.13, 0.35, 2.1, 99.9, 0.01, 0.01, 0.03, 0.03, 99.9, 0.01, 0.01, 0.03, 0.03,
                     100, 0, 0, 0, 0, 100, 0, 0, 0, 0),
    lower_bound = c(96.25, 0.087, 0.097, 0.25, 3.01, 94.75, 0.093, 0.11, 0.26, 4.03, 84.4, 0.19,0.23,0.43,13.8,
                    77.7,0.17,0.21,0.39,18.5,95.4,0.17,0.12,0.41,3.4,92.5,0.21,0.25,0.42,5.2,98.6,0.1,0.04,0.25,0.7,
                    97.9,0.12,0.06,0.26,1.2,97.9,0.05,0.1,0.27,1.4,97,0.06,0.1,0.32,1.7,100,0.01,0,0.02,0.02,99.8,0,
                    0.01,0.02,0.01,100,0,0,0,0,100,0,0,0,0),
    upper_bound = c(96.45,0.094,0.1,0.28,3.2,95.38,0.11,0.125,0.30,4.72,85.3,0.2,0.25,0.52,14.7,80.6,0.22,0.27,0.50,21.4,
                    95.7,0.24,0.21,0.49,3.8,93.8,0.22,0.32,0.48,6.5,98.9,0.13,0.05,0.30,0.9,98.3,0.14,0.08,0.30,1.5,
                    98.2,0.06,0.12,0.32,1.6,97.7,0.08,0.15,0.38,2.4,99.9,0.01,0.02,0.03,0.04,99.9,0.02,0.05,0.09,0.08,
                    100, 0, 0, 0, 0, 100, 0, 0, 0, 0)
  )
  long_dt <- long_dt %>% filter(class != 'Perennial')
  long_dt$class <- factor(long_dt$class, levels = c('1-5 no-flow days', '6-15 no-flow days',
                                                    '16-29 no-flow days', '30-31 no-flow days'))
  
  long_dt_editted <- long_dt[1:40,]
  
  # Create a combined label for ClimateZone and Period
  long_dt_editted$Zone_Period <- paste(long_dt_editted$climate_zones, long_dt_editted$period, sep = "_")
  long_dt_editted$Zone_Period <- factor(long_dt_editted$Zone_Period,
                                        levels = c("Europe_Reference", "Europe_Future",
                                                   "mediterranean / semi-arid_Reference",
                                                   "mediterranean / semi-arid_Future",
                                                   "humid subtropical_Reference",
                                                   "humid subtropical_Future",
                                                   "temperate oceanic_Reference",
                                                   "temperate oceanic_Future",
                                                   "humid continental_Reference",
                                                   "humid continental_Future"))
  
  # Calculate the cumulative position for each class segment
  long_dt_editted <- long_dt_editted %>%
    group_by(Zone_Period) %>%
    arrange(desc(class)) %>%
    mutate(CumulativeValue = cumsum(median_value),
           ymin = CumulativeValue - median_value,  # Base of each segment
           ymax = CumulativeValue,
           label_y = (ymin + ymax) / 2)    
  
  # Define colors for each class
  color_mapping <- c(
    '1-5 no-flow days' = '#5ab4ac',
    '6-15 no-flow days' = '#c7eae5',
    '16-29 no-flow days' = '#f6e8c3',
    '30-31 no-flow days' = '#d8b365'
  )
  
  # Create the stacked bar plot with reference and future periods next to each other
  outp <- ggplot(long_dt_editted, aes(x = Zone_Period, y = median_value, fill = class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65) +
    # Adding error bars
    geom_errorbar(aes(ymin = ymin + lower_bound, ymax = ymin + upper_bound), 
                  width = 0.3, 
                  position = position_dodge(width = 0.8), 
                  color = "black") +
    scale_fill_manual(values = color_mapping) +
    theme_bw(base_size = 16) +
    labs(y = "Fraction of reach-months in the four intermittence classes [%]",
         x = "") +
    theme(legend.title = element_text(colour = "black",
                                      size =  16),
          legend.text =  element_text(colour = "black",
                                      size =  14),
          legend.position = c(0.70, 0.75),
          axis.title = element_text(colour = "black",
                                    size =  16),
          axis.text = element_text(colour = "black",
                                   size =  14)) +
    guides(shape = 'none', alpha='none',
           fill = guide_legend(nrow=2, title="Intermittence classes")) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = c("Europe_Reference" = "Reference",
                                "Europe_Future" = "Future",
                                "mediterranean / semi-arid_Reference" = "Reference",
                                "mediterranean / semi-arid_Future" = "Future",
                                "humid subtropical_Reference" = "Reference",
                                "humid subtropical_Future" = "Future",
                                "temperate oceanic_Reference" = "Reference",
                                "temperate oceanic_Future" = "Future",
                                "humid continental_Reference" = "Reference",
                                "humid continental_Future" = "Future"))
  
  ggsave(outp, filename = paste0(out_dir, "classes_fraction_stackedbar_Figure6_editted.png"),
         units = "in", height = 8, width = 9, dpi = 500)
  
  return(NULL)
}

##---------------- Compute the inter-annual variability of the non-perennial reaches figure 7 -----------
#' @title Compute the inter-annual variability of the non-perennial reaches
#' 
#' @description This function computes 
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and the 2080s under RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' @param rcp_name (character) either RCP8.5 or RCP2.6
#' 
#' @export
#' 
compute_interannual_vari_fig7 <- function(in_dir, eu_net_shp_dir, out_dir, rcp_name='RCP8.5'){
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (rcp_name == 'RCP8.5') {
    in_rcp8.5_path <-  file.path(in_dir, 'future/ssp585/far_future')
  } else {
    in_rcp8.5_path <-  file.path(in_dir, 'future/ssp126/far_future')
  }
  
  in_hist_path <- file.path(in_dir, 'historical')
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  
  # define a costumized function to compute the quantile (90th-10th) of intermittent months
  compute_quantile_inter_mon <- function(in_list, p90=FALSE){
    
    lapply(seq_along(in_list), function(i, p90){
      print(i)
      in_dt <- fst::read_fst(in_list[[i]]) %>%
        as.data.table()
      long_data <- melt(in_dt,id.vars = 'DRYVER_RIVID')
      long_data[,variable:=fasttime::fastDate(variable)]
      long_data[,c('year'):=year(variable)]
      long_data[,value:=ifelse(value>0, 1, 0)]
      absolute_values <- long_data[, .(sum_values = sum(value)), by = .(DRYVER_RIVID, year)]
      wide_data <- dcast(absolute_values, DRYVER_RIVID ~ year, value.var = "sum_values")
      quantile_90 <- wide_data[,apply(.SD, 1, function(row) quantile(row, probs = 0.9)),
                               .SDcols=-'DRYVER_RIVID']
      quantile_10 <- wide_data[,apply(.SD, 1, function(row) quantile(row, probs = 0.1)),
                               .SDcols=-'DRYVER_RIVID']
      
      
      if (p90) {
        quantile_90
      }
      else{
        differ_quantiles <- quantile_90 - quantile_10
      }
    }, p90=p90) %>% 
      do.call('cbind',.) %>% 
      `colnames<-`(c('gfdl', 'ipsl', 'mpi',
                     'mri', 'ukesm')) %>% 
      as.data.table()
  }
  
  #the quantile of interannual variability for historical period for 5 GCMs
  hist_path_list <- list.files(path = in_hist_path, pattern = '.fst', full.names = TRUE)
  in_reach_ids <- fst::read_fst(hist_path_list[1]) %>% 
    as.data.table() %>% .[,.(DRYVER_RIVID)]
  hist_quantiles <- compute_quantile_inter_mon(in_list = hist_path_list)
  hist_p90 <- compute_quantile_inter_mon(in_list = hist_path_list, p90=TRUE)
  
  # save the values for the reference period 
  hist_med_p90 <- data.table::copy(hist_p90)
  hist_med_quantiles <- data.table::copy(hist_quantiles)
  
  hist_med_p90[,median_p90 := apply(.SD, 1, median), .SDcols=names(hist_med_p90)]
  hist_med_quantiles[,median_row := apply(.SD, 1, median),
                     .SDcols=names(hist_med_quantiles)]
  
  hist_med_quantiles[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  
  med_joined_st <- left_join(reach_shp,
                             hist_med_quantiles[,.(DRYVER_RIV, median_row)],
                             by='DRYVER_RIV')
  sf::write_sf(med_joined_st,
               dsn = file.path(out_dir, 'median_inter_var_ref_fig7b.shp'))
  
  hist_med_p90[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  
  med_joined_st <- left_join(reach_shp,
                             hist_med_p90[,.(DRYVER_RIV, median_p90)],
                             by='DRYVER_RIV')
  sf::write_sf(med_joined_st,
               dsn = file.path(out_dir, 'Idry_p90_ref_fig7a.shp'))
  
  #the quantile of interannual variability for rcp8.5 period for 5 GCMs
  rcp8.5_path_list <- list.files(path = in_rcp8.5_path, pattern = '.fst', full.names = TRUE)
  rcp8.5_quantiles <- compute_quantile_inter_mon(in_list = rcp8.5_path_list)
  rcp8.5_p90 <- compute_quantile_inter_mon(in_list = rcp8.5_path_list, p90=TRUE)

  #difference quantiles for rcp8.5
  differ_quantiles <- rcp8.5_quantiles - hist_quantiles
  differ_quantiles[,median_row := apply(.SD, 1, median), .SDcols=names(differ_quantiles)]
  differ_quantiles[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  
  med_joined_st <- left_join(reach_shp,
                             differ_quantiles[,.(DRYVER_RIV, median_row)],
                             by='DRYVER_RIV')
  if (rcp_name == 'RCP8.5') {
    
    sf::write_sf(med_joined_st,
                 dsn = file.path(out_dir, paste0('Iiv_p90p10_', rcp_name, '_fig7d.shp')))
  } else {
    sf::write_sf(med_joined_st,
                 dsn = file.path(out_dir, paste0('Iiv_p90p10_', rcp_name, '_figS5b.shp')))
  }
  
  differ_p90 <- rcp8.5_p90 - hist_p90
  differ_p90[,median_p90 := apply(.SD, 1, median), .SDcols=names(differ_p90)]
  differ_p90[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  p90_joined_st <- left_join(reach_shp,
                             differ_p90[,.(DRYVER_RIV,  median_p90)],
                             by='DRYVER_RIV')
  
  if (rcp_name == 'RCP8.5') {
    sf::write_sf(p90_joined_st,
                 dsn = file.path(out_dir, paste0('Idry_p90_', rcp_name, '_fig7c.shp')))
  } else {
    sf::write_sf(p90_joined_st,
                 dsn = file.path(out_dir, paste0('Idry_p90_', rcp_name, '_figS5a.shp')))
  }

  return(NULL)
}
## ------------ Compute the percentage of non-perennial reaches within climate zones figure S3-----------------
#' @title Plot the seasonality of streamflow intermittence witht the six climate zones
#' 
#' @description This function generates a plot that depicts Percentage of European reaches 
#' that are non-perennial in each calendar month of the reference and future periods,
#' distinguishing RCP2.6 and RCP8.5 for the climate zones.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and two future periods (the 2050s, the 2080s) under RCP2.6 and RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `in_shp_climcode_path` the path to the shapefile of the climate zone of the 
#' European network.
#' 
#' @export
#'
ggline_seasonality_figS3 <- function(in_dir, out_dir, in_shp_climcode_path){
  
  if(!file.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  in_hist_path <- file.path(in_dir,'historical')
  in_ssp126_nearfuture_path <- file.path(in_dir,'future/ssp126/near_future')
  in_ssp126_farfuture_path <- file.path(in_dir,'future/ssp126/far_future')
  in_ssp585_nearfuture_path <- file.path(in_dir,'future/ssp585/near_future')
  in_ssp585_farfuture_path <- file.path(in_dir,'future/ssp585/far_future')
  
  
  compute_statistics_inter_colwise <- function(in_path, in_shp_withclimcode_path, col_names){
    
    path_list <- list.files(path = in_path, pattern = '.fst', full.names = TRUE)
    date_names <- fst::read_fst(path_list[1]) %>% 
      names() %>% .[-1]
    
    reach_withclimcode_dt <- sf::read_sf(in_shp_withclimcode_path) %>%
      st_drop_geometry() %>%
      as.data.table()
    # compute streamflow intermittence colwise 
    intermittent_colwise <- lapply(seq_along(path_list), function(i, reach_withclimcode_dt){
      in_dt <- fst::read_fst(path_list[[i]]) %>%
        as.data.table()
      setnames(in_dt, 'DRYVER_RIVID', 'DRYVER_RIV')
      in_dt1 <- left_join(in_dt, reach_withclimcode_dt[,.(DRYVER_RIV, clim_code)], by='DRYVER_RIV')
      climate_zones <- reach_withclimcode_dt[,clim_code] %>% unique() %>% sort()
      intermittent_clim_zones <- lapply(seq_along(climate_zones), function(j){
        in_dt2 <- in_dt1[clim_code == climate_zones[j]]
        row_counts <- colSums(in_dt2[, lapply(.SD, function(x) x > 0)])/dim(in_dt2)[1] * 100
      })%>% 
        do.call('rbind',.) %>% .[,-c(1, 362)] %>% as.data.table()
      
    }, reach_withclimcode_dt=reach_withclimcode_dt) %>% 
      do.call('rbind',.) %>%
      `colnames<-`(date_names)  %>% 
      as.data.table()
    
    climate_zones <- reach_withclimcode_dt[,clim_code] %>% unique() %>% sort()
    intermittent_colwise[,climte_zones:=rep(climate_zones, 5)]
    intermittent_colwise[,models:=c(rep('gfdl', 6), rep('ipsl',6), rep('mpi', 6),
                                    rep('mri', 6), rep('ukesm', 6))]
    #mediterranean climate zone 
    mediter_clim_zone <- intermittent_colwise[climte_zones==5,]
    mediter_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(mediter_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    mediter_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    # subtropical climate zone
    subtropical_clim_zone <- intermittent_colwise[climte_zones==14,]
    subtropical_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(subtropical_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    subtropical_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    # temprate oceanic climate zone
    toceanic_clim_zone <- intermittent_colwise[climte_zones==15,]
    toceanic_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(toceanic_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    toceanic_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    # humid contenetial climate zone
    hcontinetal_clim_zone <- intermittent_colwise[climte_zones==25,]
    hcontinetal_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(hcontinetal_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    hcontinetal_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    # sub-arctic climate zone
    subarctic_clim_zone <- intermittent_colwise[climte_zones==27,]
    subarctic_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(subarctic_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    subarctic_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    # polar climate zone
    polar_clim_zone <- intermittent_colwise[climte_zones==30,]
    polar_clim_zone[, climte_zones:=NULL]
    melted_dt <- melt(polar_clim_zone, id.vars = 'models')
    melted_dt[, variable:=as.Date(variable)]
    melted_dt[, month:=month(variable)]
    melted_dt[,c('median', 'max', 'min'):=list(median(value), max(value), min(value)), by='month']
    polar_dt <- melted_dt[,.(median=unique(median), max=unique(max),min=unique(min)), by='month']
    
    out_dt <- cbind(mediter_dt, subtropical_dt, toceanic_dt,
                    hcontinetal_dt, subarctic_dt, polar_dt) %>% 
      `colnames<-`(rep(c('month', col_names), 6)) %>% 
      as.data.table()
    
    return(out_dt)
  }
  
  hist_dt <- compute_statistics_inter_colwise(in_path = in_hist_path,
                                              in_shp_withclimcode_path=in_shp_climcode_path,
                                              col_names = c('y1', 'y1hi', 'y1lo'))
  ssp126_nearfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp126_nearfuture_path,
                                                           in_shp_withclimcode_path=in_shp_climcode_path,
                                                           col_names = c('y2', 'y2hi', 'y2lo'))
  ssp126_farfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp126_farfuture_path,
                                                          in_shp_withclimcode_path=in_shp_climcode_path,
                                                          col_names = c('y3', 'y3hi', 'y3lo'))
  ssp585_nearfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp585_nearfuture_path,
                                                           in_shp_withclimcode_path=in_shp_climcode_path,
                                                           col_names = c('y4', 'y4hi', 'y4lo'))
  ssp585_farfuture_dt <- compute_statistics_inter_colwise(in_path = in_ssp585_farfuture_path,
                                                          in_shp_withclimcode_path=in_shp_climcode_path,
                                                          col_names = c('y5', 'y5hi', 'y5lo'))
  
  combined_dt <- cbind(hist_dt, ssp126_nearfuture_dt[,-1],
                       ssp126_farfuture_dt[,-1], ssp585_nearfuture_dt[,-1],
                       ssp585_farfuture_dt[,-1])
  
  # inner function for plotting the seasonality of non-perennial reaches 
  plot_monthly_intermittent <- function(in_hist_dt, in_ssp126_nearfuture_dt,
                                        in_ssp126_farfuture_dt, in_ssp585_nearfuture_dt,
                                        in_ssp585_farfuture_dt, xlab='Calendar month', xloc=0.2, yloc=0.7,
                                        region_name, col_region){
    combined_dt <- cbind(in_hist_dt, in_ssp126_nearfuture_dt,
                         in_ssp126_farfuture_dt, in_ssp585_nearfuture_dt,
                         in_ssp585_farfuture_dt)
    
    # Long format data frame:
    longDF = tidyr::pivot_longer(combined_dt , cols=!month , names_to="line" , values_to="y" )
    longDF$fill = NA
    longDF$fill[grep( "1" , longDF$line  )] = "y1"
    longDF$fill[grep( "2" , longDF$line  )] = "y2"
    longDF$fill[grep( "3" , longDF$line  )] = "y3"
    longDF$fill[grep( "4" , longDF$line  )] = "y4"
    longDF$fill[grep( "5" , longDF$line  )] = "y5"
    
    longDF1 <- longDF %>% 
      mutate(line = gsub("\\d", "", line)) %>% 
      tidyr::pivot_wider(id_cols = c(month, fill), names_from = line, values_from = y)
    longDF1 <- longDF1 %>% 
      mutate(fill=case_when(
        fill == 'y1'~'Reference',
        fill == 'y2'~'RCP2.6 2050s',
        fill == 'y3'~'RCP2.6 2080s',
        fill == 'y4'~'RCP8.5 2050s',
        fill == 'y5'~'RCP8.5 2080s'
      ))
    abb<- c("J","F","M","A","M","J","J","A","S","O","N","D") 
    # Define the order of legend labels
    desired_order <- c("Reference", "RCP2.6 2050s", "RCP2.6 2080s",
                       'RCP8.5 2050s', 'RCP8.5 2080s')
    
    
    ## plot with errorbar
    # Create data.frame with shading info
    shading <- data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(longDF1$month))), by = 1),
                          max = seq(from = 1.5, to = max(as.numeric(as.factor(longDF1$month))) + 0.5, by = 1),
                          col = c(0,1))
    
    ggplot(longDF1) +
      geom_rect(data = shading,
                aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf,
                    fill = factor(col), alpha = 0.1)) +
      scale_fill_manual(values = c("white", "gray85"))+
      geom_errorbar(aes(x=month,ymin=ylo, ymax=yhi, color = factor(fill, desired_order)), width=1,
                    position = position_dodge(0.8))+
      geom_line(aes(x=month, y=y, color = factor(fill, desired_order)), linewidth=1) +
      scale_color_manual(values = c('Reference' = "black", 'RCP2.6 2050s' = "#92c5de",
                                    'RCP2.6 2080s'='#0571b0', 'RCP8.5 2050s'='#f4a582',
                                    'RCP8.5 2080s'='#ca0020')) + 
      theme_bw(18)+
      scale_x_continuous(breaks=1:12, labels=abb) +
      annotation_custom(grid::textGrob(region_name, xloc, yloc, gp = gpar(fontsize=14, fontface="bold", col=col_region)))+
      theme(legend.position = 'bottom')+
      guides(color = guide_legend(ncol = 2, title=''),
             fill ='none', alpha='none')+
      labs(y='' , x=xlab)
  }
  
  ## Plot the seasonality of the non-perennial for the reaches in the climate zones 
  p_med <- plot_monthly_intermittent(hist_dt[,1:4],ssp126_nearfuture_dt[,2:4],
                                     ssp126_farfuture_dt[,2:4], ssp585_nearfuture_dt[,2:4],
                                     ssp585_farfuture_dt[,2:4], xlab='', xloc=0.2, yloc=0.7,
                                     region_name= 'mediterranean/\nsemi-arid', col_region='#A80000')
  
  p_subtropical <- plot_monthly_intermittent(hist_dt[,5:8],ssp126_nearfuture_dt[,6:8],
                                             ssp126_farfuture_dt[,6:8], ssp585_nearfuture_dt[,6:8],
                                             ssp585_farfuture_dt[,6:8], xlab='', xloc=0.3, yloc=0.7,
                                             region_name= 'humid subtropical', col_region='#FFA77F')
  
  p_toceanic <- plot_monthly_intermittent(hist_dt[,9:12],ssp126_nearfuture_dt[,10:12],
                                          ssp126_farfuture_dt[,10:12], ssp585_nearfuture_dt[,10:12],
                                          ssp585_farfuture_dt[,10:12],xlab='', xloc=0.3, yloc=0.7,
                                          region_name= 'temperate oceanic', col_region='#D1FF73')
  p_hcontinetal <- plot_monthly_intermittent(hist_dt[,13:16],ssp126_nearfuture_dt[,14:16],
                                             ssp126_farfuture_dt[,14:16], ssp585_nearfuture_dt[,14:16],
                                             ssp585_farfuture_dt[,14:16],xlab='', xloc=0.3, yloc=0.7,
                                             region_name= 'humid continental', col_region='#70A800')
  
  p_subarctic <- plot_monthly_intermittent(hist_dt[,17:20],ssp126_nearfuture_dt[,18:20],
                                           ssp126_farfuture_dt[,18:20], ssp585_nearfuture_dt[,18:20],
                                           ssp585_farfuture_dt[,18:20], xloc=0.4, yloc=0.7,
                                           region_name= 'sub-arctic', xlab = "", col_region='#73DFFF')
  
  p_polar <- plot_monthly_intermittent(hist_dt[,21:24],ssp126_nearfuture_dt[,22:24],
                                       ssp126_farfuture_dt[,22:24], ssp585_nearfuture_dt[,22:24],
                                       ssp585_farfuture_dt[,22:24], xloc=0.3, yloc=0.7,
                                       region_name= 'polar/alpine', col_region='#005CE6')
  
  pcols <- cowplot::plot_grid(p_med+theme(legend.position = 'none'),
                              p_subtropical+theme(legend.position = 'none'),
                              p_toceanic+theme(legend.position = 'none'),
                              p_hcontinetal+theme(legend.position = 'none'),
                              p_subarctic+theme(legend.position = 'none'),
                              p_polar+theme(legend.position = 'none'),
                              # align = 'vh',
                              labels=c('a', 'b', 'c', 'd', 'e', 'f'),
                              hjust = -0.5, label_size= 18,
                              ncol = 2)
  p_polar_legend <- p_polar + theme(legend.background = element_rect(fill = "white"))
  legend_b <- cowplot::get_legend(p_polar) 
  
  final_plot <- cowplot::plot_grid(pcols,legend_b, ncol = 1,rel_heights = c(1, .2))
  
  #create common x and y labels
  y_grob <- textGrob("Fraction of non-perennial reaches [%]", 
                     gp=gpar(col="black", fontsize=18), rot=90)
  
  #add to plot
  p_out <- gridExtra::grid.arrange(arrangeGrob(final_plot, left = y_grob))
  
  ggsave(p_out, filename = file.path(out_dir, 'all_climatezones_inone_errorbar_figS3.png'),
         dpi = 600, units = 'in', height = 12, width = 10)
  
}
## ------------- compute the median percentage of non-perennial reach-months figure S4----------
#' @title Compute the median of NPRMs
#' 
#' @description This function computes the median percentage of non-perennial reach-months (NPRMs)
#' in the reference period and of the five intermittence status classes in the 2080s under RCP8.5.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and the 2080s under RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' 
#' @export
#' 
create_shp_figS4 <- function(in_dir, eu_net_shp_dir, out_dir){
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  
  in_list <- list.files(path = in_dir, pattern = '.fst', full.names = TRUE)
  
  # define a costumized function to compute the number of intermittent months
  compute_inter_classes_pp <- function(in_list){
    lapply(seq_along(in_list), function(i){
      in_dt <- fst::read_fst(in_list[[i]]) %>%
        as.data.table()
      col_counts_zero <- rowSums(in_dt[, .SD[, -1, with = FALSE] == 0])
      col_counts_class1 <- rowSums(in_dt[, .SD[, -1, with = FALSE] == 1])
      col_counts_class2 <- rowSums(in_dt[, .SD[, -1, with = FALSE] == 2])
      col_counts_class3 <- rowSums(in_dt[, .SD[, -1, with = FALSE] == 3])
      col_counts_class4 <- rowSums(in_dt[, .SD[, -1, with = FALSE] == 4])
      
      (cbind(col_counts_zero, col_counts_class1, col_counts_class2,
             col_counts_class3, col_counts_class4)/360 * 100) %>% 
        `colnames<-`(c('perennial', 'class1', 'class2',
                       'class3', 'class4')) %>% as.data.table()
    }) %>% 
      do.call('cbind',.) %>% 
      as.data.table()
  }
  
  in_reach_ids <- fst::read_fst(in_list[1]) %>% 
    as.data.table() %>% .[,.(DRYVER_RIVID)]
  
  all_classes_dt <- compute_inter_classes_pp(in_list = in_list)
  all_classes_dt <- round(all_classes_dt, digits = 2)
  
  # median percentage point(P.P) of perennial months in GCMs
  perennial_cols <- all_classes_dt[, .SD,
                                   .SDcols = which(names(all_classes_dt) == "perennial")]
  perennial_cols[,median_row := apply(.SD, 1, median), .SDcols=names(perennial_cols)]
  
  perennial_cols[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  perennial_joined_st <- left_join(reach_shp,
                                   perennial_cols[,.(median_row, DRYVER_RIV)],
                                   by='DRYVER_RIV')
  sf::write_sf(perennial_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_perennial_figS4b.shp'))
  
  # median percentage point(P.P) of 1-5 no-flow days within month in GCMs
  class1_cols <- all_classes_dt[, .SD, .SDcols = which(names(all_classes_dt) == "class1")]
  class1_cols[,median_row := apply(.SD, 1, median), .SDcols=names(class1_cols)]
  
  class1_cols[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  class1_joined_st <- left_join(reach_shp,
                                class1_cols[,.(median_row, DRYVER_RIV)],
                                by='DRYVER_RIV')
  sf::write_sf(class1_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_class1_figS4c.shp'))
  
  # median percentage point(P.P) of 6-15 no-flow days within month in GCMs
  class2_cols <- all_classes_dt[, .SD, .SDcols = which(names(all_classes_dt) == "class2")]
  class2_cols[,median_row := apply(.SD, 1, median), .SDcols=names(class2_cols)]
  
  class2_cols[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  class2_joined_st <- left_join(reach_shp,
                                class2_cols[,.(median_row, DRYVER_RIV)],
                                by='DRYVER_RIV')
  sf::write_sf(class2_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_class2_figS4d.shp'))
  # median percentage point(P.P) of 16-29 no-flow days within month in GCMs
  class3_cols <- all_classes_dt[, .SD, .SDcols = which(names(all_classes_dt) == "class3")]
  class3_cols[,median_row := apply(.SD, 1, median), .SDcols=names(class3_cols)]
  class3_cols[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  class3_joined_st <- left_join(reach_shp,
                                class3_cols[,.(median_row, DRYVER_RIV)],
                                by='DRYVER_RIV')
  sf::write_sf(class3_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_class3_figS4e.shp'))
  # median percentage point(P.P) of 30-31 no-flow days within month in GCMs
  class4_cols <- all_classes_dt[, .SD, .SDcols = which(names(all_classes_dt) == "class4")]
  class4_cols[,median_row := apply(.SD, 1, median), .SDcols=names(class4_cols)]
  class4_cols[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  class4_joined_st <- left_join(reach_shp,
                                class4_cols[,.(median_row, DRYVER_RIV)],
                                by='DRYVER_RIV')
  sf::write_sf(class4_joined_st,
               dsn = file.path(out_dir, 'median_inter_mon_class4_figS4f.shp'))
  
  return(NULL)
}

##------------- compute the changes in reach-month for table 3 in paper -------------------
#' @title Compute he changes in reach-month for table 3
#' 
#' @description This function generates the values that are presented in table 3 in the paper
#' {}.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for and the 2050s and 
#' 2080s under RCP2.6 and RCP8.5. 
#' @param `out_dir` the path to the location that the final shapefiles are stored.
#' @param `eu_net_shp_dir` the path to the shapefile of the European river network
#' 
#' @export
#' 
compute_changes_table3 <- function(in_dir, out_dir, eu_net_shp_dir){
  
  # create the inputs paths
  in_hist_path <- file.path(in_dir, 'historical') 
  in_nearfuture_path_rcp2.6 <- file.path(in_dir, 'future/ssp126/near_future') 
  in_farfuture_path_rcp2.6 <- file.path(in_dir, 'future/ssp126/far_future')
  in_nearfuture_path_rcp8.5 <- file.path(in_dir, 'future/ssp585/near_future') 
  in_farfuture_path_rcp8.5 <- file.path(in_dir, 'future/ssp585/far_future')
  
  
  reach_shp <- sf::read_sf(eu_net_shp_dir)
  # define a costumized function to compute the number of intermittent months
  compute_number_inter_mon <- function(in_list){
    lapply(seq_along(in_list), function(i){
      in_dt <- fst::read_fst(in_list[[i]]) %>%
        as.data.table()
      col_counts_ratio <- rowSums(in_dt[, .SD[, -1, with = FALSE] > 0])
    }) %>% 
      do.call('cbind',.) %>% 
      `colnames<-`(c('gfdl', 'ipsl', 'mpi',
                     'mri', 'ukesm')) %>% 
      as.data.table()
  }
  
  ## total number of reach-months for 1.5 million reaches -> 552148560
  # Define breaks for classification
  breaks <- c(0, 10, 50, 500, 10000, Inf)
  labels <- c("1", "2", "3", '4', '5')
  in_reach_ids <- fst::read_fst(hist_path_list[1]) %>% 
    as.data.table() %>% .[,.(DRYVER_RIVID)]
  # reference period ------
  hist_path_list <- list.files(path = in_hist_path, pattern = '.fst', full.names = TRUE)
  ref_preiod <- compute_number_inter_mon(hist_path_list)
  ref_preiod[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  
  ref_preiod1 <- left_join(ref_preiod, reach_dt, by='DRYVER_RIV')
  ref_preiod1[, classification := cut(upa, breaks, labels = labels, include.lowest = TRUE)]
  d <- ref_preiod1[,.N, by='classification'] %>% .[order(classification)]
  reach_months_upstream_classes <- d[,N] * 360
  data_ref <- ref_preiod1[, lapply(.SD, sum), by='classification', .SDcols=-'classification']
  ref_med_allgcms <- data_ref[,lapply(.SD, function(x) round(x/reach_months_upstream_classes*100, 3)),
                              .SDcols=-'classification']
  
  ref_median <- ref_med_allgcms[,apply(.SD, 1, median), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  ref_min <- ref_med_allgcms[,apply(.SD, 1, min), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  ref_max <- ref_med_allgcms[,apply(.SD, 1, max), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  
  ref_med_eu <- round(colSums(data_ref[,-1])[1:5]/552148560 * 100,3) %>% median()
  ref_min_eu <- round(colSums(data_ref[,-1])[1:5]/552148560 * 100,3) %>% min()
  ref_max_eu <- round(colSums(data_ref[,-1])[1:5]/552148560 * 100,3) %>% max()
  
  ref_pre_combined <- rbind(c(ref_median, ref_med_eu), c(ref_min, ref_min_eu), c(ref_max, ref_max_eu))
  
  ref_eu_gcms <- round(colSums(data_ref[,-1])[1:5]/552148560 * 100,3)
  # RCP2.6 in the 2050s ------
  rcp2.6_near_list <- list.files(path = in_nearfuture_path_rcp2.6, pattern = '.fst', full.names = TRUE)
  m_near <- compute_number_inter_mon(rcp2.6_near_list)
  m_near[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  m_near_rcp2.6 <- left_join(m_near, reach_dt, by='DRYVER_RIV')
  
  m_near_rcp2.6[, classification := cut(upa, breaks, labels = labels, include.lowest = TRUE)]
  data_near <- m_near_rcp2.6[, lapply(.SD, sum), by='classification', .SDcols=-'classification']
  rcp2.6_med_allgcms <- data_near[,lapply(.SD, function(x) round(x/reach_months_upstream_classes*100, 3)),
                                  .SDcols=-'classification']
  recp2.6_dif_allgcms <- rcp2.6_med_allgcms - ref_med_allgcms
  rcp2.6_median <- recp2.6_dif_allgcms[,apply(.SD, 1, median), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp2.6_min <- recp2.6_dif_allgcms[,apply(.SD, 1, min), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp2.6_max <- recp2.6_dif_allgcms[,apply(.SD, 1, max), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  
  rcp2.6_eu_gcms <- round(colSums(data_near[,-1])[1:5]/552148560 * 100,3)
  rcp2.6_eu_dif <- rcp2.6_eu_gcms - ref_eu_gcms
  rcp2.6_med_eu <- rcp2.6_eu_dif %>% median()
  rcp2.6_min_eu <- rcp2.6_eu_dif %>% min()
  rcp2.6_max_eu <- rcp2.6_eu_dif %>% max()
  
  rcp2.6_2050s_pre_combined <- rbind(c(rcp2.6_median, rcp2.6_med_eu), c(rcp2.6_min, rcp2.6_min_eu),
                                     c(rcp2.6_max, rcp2.6_max_eu))
  ## RCP2.6 in the 2080s--------
  rcp2.6_far_list <- list.files(path = in_farfuture_path_rcp2.6, pattern = '.fst', full.names = TRUE)
  m_end <- compute_number_inter_mon(rcp2.6_far_list)
  m_end[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  m_end_rcp2.6 <- left_join(m_end, reach_dt, by='DRYVER_RIV')
  m_end_rcp2.6[, classification := cut(upa, breaks, labels = labels, include.lowest = TRUE)]
  data_far <- m_end_rcp2.6[, lapply(.SD, sum), by='classification', .SDcols=-'classification']
  rcp2.6_med_allgcms <- data_far[,lapply(.SD, function(x) round(x/reach_months_upstream_classes*100, 3)),
                                 .SDcols=-'classification']
  recp2.6_dif_allgcms <- rcp2.6_med_allgcms - ref_med_allgcms
  rcp2.6_median <- recp2.6_dif_allgcms[,apply(.SD, 1, median), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp2.6_min <- recp2.6_dif_allgcms[,apply(.SD, 1, min), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp2.6_max <- recp2.6_dif_allgcms[,apply(.SD, 1, max), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  
  rcp2.6_eu_gcms <- round(colSums(data_far[,-1])[1:5]/552148560 * 100,3)
  rcp2.6_eu_dif <- rcp2.6_eu_gcms - ref_eu_gcms
  
  rcp2.6_med_eu <- rcp2.6_eu_dif %>% median()
  rcp2.6_min_eu <- rcp2.6_eu_dif %>% min()
  rcp2.6_max_eu <- rcp2.6_eu_dif %>% max()
  
  rcp2.6_2080s_pre_combined <- rbind(c(rcp2.6_median, rcp2.6_med_eu), c(rcp2.6_min, rcp2.6_min_eu),
                                     c(rcp2.6_max, rcp2.6_max_eu))
  
  ## RCP8.5 in the 2050s -------------
  rcp8.5_near_list <- list.files(path = in_nearfuture_path_rcp8.5, pattern = '.fst', full.names = TRUE)
  m_near <- compute_number_inter_mon(rcp8.5_near_list)
  m_near[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  m_near_rcp8.5 <- left_join(m_near, reach_dt, by='DRYVER_RIV')
  m_near_rcp8.5[, classification := cut(upa, breaks, labels = labels, include.lowest = TRUE)]
  data_near <- m_near_rcp8.5[, lapply(.SD, sum), by='classification', .SDcols=-'classification']
  rcp8.5_med_allgcms <- data_near[,lapply(.SD, function(x) round(x/reach_months_upstream_classes*100, 3)),
                                  .SDcols=-'classification']
  recp8.5_dif_allgcms <- rcp8.5_med_allgcms - ref_med_allgcms
  rcp8.5_median <- recp8.5_dif_allgcms[,apply(.SD, 1, median), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp8.5_min <- recp8.5_dif_allgcms[,apply(.SD, 1, min), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp8.5_max <- recp8.5_dif_allgcms[,apply(.SD, 1, max), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  
  rcp8.5_eu_gcms <- round(colSums(data_near[,-1])[1:5]/552148560 * 100,3)
  rcp8.5_eu_dif <- rcp8.5_eu_gcms - ref_eu_gcms
  
  rcp8.5_med_eu <- rcp8.5_eu_dif %>% median()
  rcp8.5_min_eu <- rcp8.5_eu_dif %>% min()
  rcp8.5_max_eu <- rcp8.5_eu_dif %>% max()
  
  rcp8.5_2050s_pre_combined <- rbind(c(rcp8.5_median, rcp8.5_med_eu), c(rcp8.5_min, rcp8.5_min_eu),
                                     c(rcp8.5_max, rcp8.5_max_eu))
  ## RCP8.5 in the 2080s --------
  rcp8.5_far_list <- list.files(path = in_farfuture_path_rcp8.5, pattern = '.fst', full.names = TRUE)
  m_end <- compute_number_inter_mon(rcp8.5_far_list)
  m_end[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
  m_end_rcp8.5 <- left_join(m_end, reach_dt, by='DRYVER_RIV')
  m_end_rcp8.5[, classification := cut(upa, breaks, labels = labels, include.lowest = TRUE)]
  data_far <- m_end_rcp8.5[, lapply(.SD, sum), by='classification', .SDcols=-'classification']
  rcp8.5_med_allgcms <- data_far[,lapply(.SD, function(x) round(x/reach_months_upstream_classes*100, 3)),
                                 .SDcols=-'classification']
  recp8.5_dif_allgcms <- rcp8.5_med_allgcms - ref_med_allgcms
  rcp8.5_median <- recp8.5_dif_allgcms[,apply(.SD, 1, median), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp8.5_min <- recp8.5_dif_allgcms[,apply(.SD, 1, min), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp8.5_max <- recp8.5_dif_allgcms[,apply(.SD, 1, max), .SDcols=c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')]
  rcp8.5_eu_gcms <- round(colSums(data_far[,-1])[1:5]/552148560 * 100,3)
  rcp8.5_eu_dif <- rcp8.5_eu_gcms - ref_eu_gcms
  rcp8.5_med_eu <- rcp8.5_eu_dif %>% median()
  rcp8.5_min_eu <- rcp8.5_eu_dif %>% min()
  rcp8.5_max_eu <- rcp8.5_eu_dif %>% max()
  
  rcp8.5_2080s_pre_combined <- rbind(c(rcp8.5_median, rcp8.5_med_eu), c(rcp8.5_min, rcp8.5_min_eu),
                                     c(rcp8.5_max, rcp8.5_max_eu))
  
  # combine all the scenarios and the statistics into the final table ---------
  scenarios <- c(rep('reference', 3), rep('RCP2.6_2050s', 3), rep('RCP2.6_2080s', 3),
                 rep('RCP8.5_2050s', 3), rep('RCP8.5_2080s', 3))
  statistics_name <- c(rep(c('median', 'min', 'max'), 5)) 
  prefinal_table <- rbind(ref_pre_combined, rcp2.6_2050s_pre_combined,
                          rcp2.6_2080s_pre_combined, rcp8.5_2050s_pre_combined,
                          rcp8.5_2080s_pre_combined) %>% 
    round(., 1)
  final_table <- cbind(scenarios, statistics_name, prefinal_table) %>% as.data.table()
  
  colnames(final_table) <- c('scenarios', 'statistics','[0-10)', '[10-50)', '[50-500)',
                             '[500-10,000)', '>10,000', 'Europe')
  
  return(final_table)
  
}


## ------------ Compute the shifts in the reach-months table 4; table S3--------
#' @title Percentage of the reach-months shifts
#' 
#' @description This function produces a table that show the percentage of the reach-months that
#' a regim shifts (from non-perennial to perennial or vice versa) takes place due to climate change.
#' 
#' @param `in_dir` the path to the intermittence status of the five GCMs for the reference period,
#' and in the 2080s under RCP2.6 and RCP8.5. 
#' @param `in_shp_climcode_path` the path to the shapefile of the climate zone of the 
#' European network.
#' @param threshold this is a threshold that considers a reach as non-perenial if it stops to flow
#' at least one day per month above this threshold. `default=354`, which means that for at least 6 months 
#' over a 30-year period non-perennial occurs. 
#' 
#' @export
compute_status_shifts_table4 <- function(in_dir, in_shp_climcode_path,threshold=354){
  
  # Inner function to compute the shifts for the both directions ----
  analyze_shifts_status <- function(in_hist_path, in_ssp126_far_path, in_ssp585_far_path,
                                    in_shp_climcode_path, threshold=354, rcp_name = 'RCP8.5') {
    
    # function to calculate the percenatege of flow intermittence for the reaches and overall
    # for different scenarios
    compute_perennial_months <- function(in_dt){
      
      col_counts <- rowSums(in_dt[, .SD[ , with = FALSE] == 0, .SDcols = -1])
      
      
      return(col_counts)
    }
    
    # compute the changes for climate change scenarios
    
    in_list <- list(in_hist_path, in_ssp126_far_path, in_ssp585_far_path)
    
    reach_shp_withclimcode <- sf::read_sf(in_shp_climcode_path)
    climcode_dt <- reach_shp_withclimcode %>% st_drop_geometry() %>% as.data.table()
    all_num <- climcode_dt[,.N, by='clim_code'][order(clim_code),N]
    # compute the number of perennial months for all the reaches
    perennial_per_individial <- lapply(seq_along(in_list), function(i){
      in_dt <- fst::read_fst(in_list[[i]]) %>%
        as.data.table()
      compute_perennial_months(in_dt)
    }) %>% 
      do.call('cbind',.) %>% 
      `colnames<-`(c('hist', 'ssp126_far', 'ssp585_far')) %>% 
      as.data.table()
    
    in_reach_ids <- fst::read_fst(in_hist_path) %>% 
      as.data.table() %>% .[,.(DRYVER_RIVID)]
    perennial_per_individial[, 'DRYVER_RIV' := in_reach_ids$DRYVER_RIVID]
    perennial_joined_st <- left_join(perennial_per_individial,
                                     climcode_dt[,.(clim_code, DRYVER_RIV)],
                                     by='DRYVER_RIV')
    # Generate a complete list of clim_code values
    all_clim_codes <- unique(perennial_joined_st[, clim_code])
    
    if (rcp_name == 'RCP2.6') {
      # non-perennial to perennial for the reaches
      np2p_europe_ratio <- perennial_per_individial[(hist < threshold & ssp126_far >= threshold)][,.N]/1533415*100
      climatezones_np2p_per <- perennial_joined_st[, .(hist, ssp126_far, clim_code)] %>%
        .[(hist < threshold & ssp126_far >= threshold), .N,by='clim_code']
      
      climatezones_np2p_per <- merge(data.table(clim_code = all_clim_codes), 
                                     climatezones_np2p_per, by = "clim_code", all.x = TRUE)
      # Replace NA values with 0
      climatezones_np2p_per[is.na(N), N := 0]
      
      climatezones_np2p_ratio <- climatezones_np2p_per[order(clim_code),N]/all_num *100
      # perennial to non-perennial for the reaches
      p2np_europe_ratio <- perennial_per_individial[(hist > threshold & ssp126_far <=threshold)][,.N]/1533415*100
      climatezones_p2np_per <- perennial_joined_st[, .(hist,ssp126_far, clim_code)] %>%
        .[(hist > threshold & ssp126_far <=threshold), .N,by='clim_code']
      
      climatezones_p2np_per <- merge(data.table(clim_code = all_clim_codes), 
                                     climatezones_p2np_per, by = "clim_code", all.x = TRUE)
      # Replace NA values with 0
      climatezones_p2np_per[is.na(N), N := 0]
      
      climatezones_p2np_ratio <- climatezones_p2np_per[order(clim_code),N]/all_num *100
      
    } else {
      
      # non-perennial to perennial for the reaches
      np2p_europe_ratio <- perennial_per_individial[(hist < threshold & ssp585_far >= threshold)][,.N]/1533415*100
      climatezones_np2p_per <- perennial_joined_st[, .(hist,ssp585_far, clim_code)] %>%
        .[(hist < threshold & ssp585_far >= threshold), .N,by='clim_code']
      
      climatezones_np2p_per <- merge(data.table(clim_code = all_clim_codes), 
                                     climatezones_np2p_per, by = "clim_code", all.x = TRUE)
      # Replace NA values with 0
      climatezones_np2p_per[is.na(N), N := 0]
      
      climatezones_np2p_ratio <- climatezones_np2p_per[order(clim_code),N]/all_num *100
      # perennial to non-perennial for the reaches
      p2np_europe_ratio <- perennial_per_individial[(hist > threshold & ssp585_far <=threshold)][,.N]/1533415*100
      climatezones_p2np_per <- perennial_joined_st[, .(hist,ssp585_far, clim_code)] %>%
        .[(hist > threshold & ssp585_far <=threshold), .N,by='clim_code']
      
      climatezones_p2np_per <- merge(data.table(clim_code = all_clim_codes), 
                                     climatezones_p2np_per, by = "clim_code", all.x = TRUE)
      # Replace NA values with 0
      climatezones_p2np_per[is.na(N), N := 0]
      
      climatezones_p2np_ratio <- climatezones_p2np_per[order(clim_code),N]/all_num *100
      
    }
    
    out_list <- rbind(np2p=c(np2p_europe_ratio, climatezones_np2p_ratio),
                      p2np=c(p2np_europe_ratio, climatezones_p2np_ratio))
    
    return(out_list)
  }
  
  # model_names <- c('GFDL-ESM4', 'IPSL-CM6A-LR', 'MPI-ESM1-2-HR', 'MRI-ESM2-0', 'UKESM1-0-LL')
  in_list_hist <- list.files(file.path(in_dir, 'historical'), pattern = '.fst', full.names = TRUE)
  in_list_rcp2.6 <- list.files(file.path(in_dir, 'future/ssp126/far_future'), pattern = '.fst', full.names = TRUE)
  in_list_rcp8.5 <- list.files(file.path(in_dir, 'future/ssp585/far_future'), pattern = '.fst', full.names = TRUE)
  
  ## Compute the shifts for rcp2.6 in the 2080s -------------
  m_rcp2.6 <- lapply(1:5, function(i){
    in_hist_path <- in_list_hist[i] 
    in_ssp126_far_path <- in_list_rcp2.6[[i]]
    in_ssp585_far_path <- in_list_rcp8.5[[i]]
    
    analyze_shifts_status(in_hist_path, in_ssp126_far_path, in_ssp585_far_path,
                          in_shp_climcode_path=in_shp_climcode_path,
                          rcp_name = 'RCP2.6', threshold=threshold)
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  ## shift from non-perennial to perennial 
  np2p_allmodels <- m_rcp2.6[c(1,3, 5, 7, 9),]
  np2p_allmodels <- round(np2p_allmodels, 1)
  rcp2.6_np2p_med <- apply(np2p_allmodels, 2, median)
  rcp2.6_np2p_min <- apply(np2p_allmodels, 2, min)
  rcp2.6_np2p_max <- apply(np2p_allmodels, 2, max)
  ## shift from perennial to non-perennial
  p2np_allmodels <- m_rcp2.6[c(2,4, 6, 8, 10),]
  p2np_allmodels <- round(p2np_allmodels, 1)
  rcp2.6_p2np_med <- apply(p2np_allmodels, 2, median)
  rcp2.6_p2np_min <- apply(p2np_allmodels, 2, min)
  rcp2.6_p2np_max <- apply(p2np_allmodels, 2, max)
  
  ## Under RCP8.5 ------------
  m_rcp8.5 <- lapply(1:5, function(i){
    in_hist_path <- in_list_hist[i] 
    in_ssp126_far_path <- in_list_rcp2.6[[i]]
    in_ssp585_far_path <- in_list_rcp8.5[[i]]
    
    analyze_shifts_status(in_hist_path, in_ssp126_far_path, in_ssp585_far_path,
                          in_shp_climcode_path, threshold=threshold)
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  ## shift from non-perennial to perennial 
  np2p_allmodels <- m_rcp8.5[c(1,3, 5, 7, 9),]
  np2p_allmodels <- round(np2p_allmodels, 1)
  rcp8.5_np2p_med <- apply(np2p_allmodels, 2, median)
  rcp8.5_np2p_min <- apply(np2p_allmodels, 2, min)
  rcp8.5_np2p_max <- apply(np2p_allmodels, 2, max)
  ## shift from perennial to non-perennial
  p2np_allmodels <- m_rcp8.5[c(2,4, 6, 8, 10),]
  p2np_allmodels <- round(p2np_allmodels, 1)
  rcp8.5_p2np_med <- apply(p2np_allmodels, 2, median)
  rcp8.5_p2np_min <- apply(p2np_allmodels, 2, min)
  rcp8.5_p2np_max <- apply(p2np_allmodels, 2, max)
  
  scenarios <- c(rep('RCP2.6_np2p', 3), rep('RCP2.6_p2np', 3),
                 rep('RCP8.5_np2p', 3), rep('RCP8.5_p2np', 3))
  
  final_table <- rbind(rcp2.6_np2p_med, rcp2.6_np2p_min,
                       rcp2.6_np2p_max,rcp2.6_p2np_med, rcp2.6_p2np_min, rcp2.6_p2np_max,
                       rcp8.5_np2p_med, rcp8.5_np2p_min, rcp8.5_np2p_max,rcp8.5_p2np_med,
                       rcp8.5_p2np_min, rcp8.5_p2np_max) %>% data.table::as.data.table()
  final_table <- final_table %>% 
    dplyr::mutate(scenarios = scenarios, .before='V1')
  
  colnames(final_table) <- c( 'scenarios', 'Europe', 'mediterranean/semi-arid', 'humid_subtropical',
                              'temperate_oceanic', 'humid_continental', 'sub-arctic', 'polar/alpine')
  
  return(final_table)
  
}


## ------------- Compute the AOA for three sets ------------
#' @title The analysis of area of applicability
#' 
#' @description This function produces the inside AOA fraction of reaches within a 0.1 deg grid cell.
#' 
#' @param `in_dir` the path to the predictors in three sests, reanalysis-based, GCMs-based reference,
#' GCMs-based 2080s under RCP8.5. For dynamic predictors, the minimum values of the period was selected.
#' @param `in_model` the path to the step 1 RF model to obtain the weight of the predictors.
#' @param `in_data_obs` the path to the observations (station-months) of 3706 gauging stations to
#' compute the dissimilarity index for training points
#' @param `out_dir` the path to the aoa geotiffs
#' 
#' @export
#' 
analyze_aoa <- function(in_dir, in_model, in_data_obs, out_dir){
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
    
  }
  
  model_step1 <- qs::qread(in_model)
  model_data <- qs::qread(in_data_obs)
  
  model_data_edit <- model_data[,c(1, 3, 4:26)]
  
  model_data_edit[,target:=ifelse(target>0,1,0)]
  min_data <- model_data_edit[, lapply(.SD, min, na.rm = TRUE),
                              by = gaugeid, .SDcols = -c('gaugeid')]
  median_data <- model_data_edit[, lapply(.SD, median, na.rm = TRUE),
                                 by = gaugeid, .SDcols = -c('gaugeid')]
  max_data <- model_data_edit[, lapply(.SD, max, na.rm = TRUE),
                              by = gaugeid, .SDcols = -c('gaugeid')]
  
  data <- rbind(min_data, median_data, max_data)
  
  data[,gaugeid:=NULL]
  data[,target:=as.factor(target)]
  backend <- as_data_backend(data)
  
  task <- as_task_classif(backend, target = "target")
  rsmp_cv <- rsmp("cv", folds = 3L)$instantiate(task)
  
  traindi <- trainDI_editted(
    train = as.data.frame(task$data()),
    variables = task$feature_names,
    weight = data.frame(t(model_step1$model$learner$model$classif.ranger$model$variable.importance)),
    CVtest = rsmp_cv$instance[order(row_id)]$fold)
  
  
  aoa_fun <- lapply(1:12, function(i){
    
    splitted_pred <- rast(file.path(in_dir,
                                    paste0('split_',i,'_predictors_15arcsec.tif')))
    
    AOA <- aoa(splitted_pred, trainDI = traindi,
               train = as.data.frame(task$data()),
               variables = task$feature_names,
               weight = data.frame(t(model_step1$model$learner$model$classif.ranger$model$variable.importance)),
               CVtest = rsmp_cv$instance[order(row_id)]$fold)
    
    prediction_raster1 <- terra::aggregate(AOA$AOA, fact=24, fun='mean', na.rm=TRUE)
    
    writeRaster(prediction_raster1, filename = file.path(out_dir,
                                                         paste0('aoa_split_', i, '_0.1deg.tif')))
    
    print(i)
  })
  
  return(NULL)
}

#' @title Calculate the difference of GCMs-based sets to reanalysis-based set
#' 
#' @description This function produces the difference of inside AOA fraction between
#' GCMs-based sets and reanalysis-based set a 0.1 deg grid cell.
#' 
#' @param `in_dir_reanalysis` the path to the inside AOA fraction at 0.1 deg grid cell for
#' reanalysis-based set.
#' @param `in_dir_phase` the path to the inside AOA fraction at 0.1 deg grid cell for either
#' 2080s or reference sets.
#' @param `out_dir` the path to the aoa difference geotiffs between sets
#' 
#' @export
#' 
compute_aoa_dif <- function(in_dir_reanalysis, in_dir_phase, out_dir){
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  # merge reanalysis-based 
  file_list <- list.files(path = in_dir_reanalysis,
                          pattern = "0.1deg*\\.tif$", full.names = TRUE)
  
  file_list <- file_list[grepl("0.1deg", file_list)]
  
  # Read the GeoTIFF files into a list of SpatRaster objects
  rasters <- lapply(file_list, rast)
  # Merge the rasters into one SpatRaster object
  reanalysis_raster <- do.call(terra::mosaic, rasters)
  
  # writeRaster(reanalysis_raster, 
  #             filename = "aoa_prob_europe_0.1deg_1985to2014_min.tif",
  #             overwrite = TRUE)
  
  gcms_models <- c('gfdl', 'ipsl', 'mpi', 'mri', 'ukesm')
  lapply(seq_along(gcms_models), function(i){
    # List of file paths to the GeoTIFF files
    file_list <- list.files(path = file.path(in_dir_phase, gcms_models[i]),
                            pattern = "0.1deg*\\.tif$", full.names = TRUE)
    
    file_list <- file_list[grepl("0.1deg", file_list)]
    
    # Read the GeoTIFF files into a list of SpatRaster objects
    rasters <- lapply(file_list, rast)
    # Merge the rasters into one SpatRaster object
    merged_raster <- do.call(terra::mosaic, rasters)
    # Save the merged raster to a new GeoTIFF file
    dif_rasters <- merged_raster - reanalysis 
    
    writeRaster(dif_rasters,
                filename = file.path(out_dir, paste0("dif_aoa_prob_europe_0.1deg_", gcms_models[i],".tif")),
                overwrite = TRUE)
    
  })
  
  return(NULL)
  
}

## ------------- Create the figure 9 and figureS7 ------------
#' @title Produce the Geotiff for figures
#' 
#' @description This function produces the Geotiffs for figure 9 and figureS7
#' 
#' @param `in_dir` the path to the predictors in three sests, reanalysis-based, GCMs-based reference,
#' GCMs-based 2080s under RCP8.5. For dynamic predictors, the minimum values of the period was selected.
#' @param `figure9` `default`=TRUE, this present which figure is generating.
#' @param `phase` a character that either `reference` or `2080s`.
#' @param `out_dir` the path to the final aoa geotiffs of the figure
#' 
#' @export
#' 
create_fig9_and_S7 <- function(in_dir, out_dir, figure9=TRUE, phase ='reference') {
  
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  # merge reanalysis-based 
  file_list <- list.files(path = in_dir_reanalysis,
                          pattern = "0.1deg*\\.tif$", full.names = TRUE)
  
  file_list <- file_list[grepl("0.1deg", file_list)]
  
  # Read the GeoTIFF files into a list of SpatRaster objects
  rasters <- lapply(file_list, rast)
  # Merge the rasters into one SpatRaster object
  reanalysis_raster <- do.call(terra::mosaic, rasters)
  
  
  file_list <- list.files(path = in_dir,
                          pattern = ".tif$", full.names = TRUE)
  
  m <- lapply(seq_along(file_list), function(i){
    
    terra::rast(file_list[i])
  }) %>% 
    do.call(c, .)
  
  median_rast_ref <- terra::app(m, median)
  max_rast_ref <- terra::app(diff_ref, max)
  min_rast_ref <- terra::app(diff_ref, min)
  
  if (figure9) {
    writeRaster(reanalysis_raster,
                filename = file.path(out_dir, "aoa_prob_europe_0.1deg_1985to2014_min_fig9a.tif"),
                overwrite = TRUE)
    if (phase == 'referece') {
      
      writeRaster(median_rast_ref,
                  filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_ref_1985to2014_median_fiveGCM_fig9b.tif"),
                  overwrite = TRUE)
    } else {
      writeRaster(median_rast_ref,
                  filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_2080s_1985to2014_median_fiveGCM_fig9c.tif"),
                  overwrite = TRUE)
    }
  }
  
  if (!figure9 & phase == 'reference') {
    writeRaster(max_rast_ref,
                filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_ref_1985to2014_max_fiveGCM_figS7b.tif"),
                overwrite = TRUE)
    writeRaster(min_rast_ref,
                filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_ref_1985to2014_min_fiveGCM_figS7a.tif"),
                overwrite = TRUE)
  }
  
  if (!figure9 & phase == '2080s') {
    writeRaster(max_rast_ref,
                filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_2080s_1985to2014_max_fiveGCM_figS7d.tif"),
                overwrite = TRUE)
    writeRaster(min_rast_ref,
                filename = file.path(out_dir, "dif_aoa_prob_europe_0.1deg_2080s_1985to2014_min_fiveGCM_figS7c.tif"),
                overwrite = TRUE)
  }
  
  return(NULL)
}

