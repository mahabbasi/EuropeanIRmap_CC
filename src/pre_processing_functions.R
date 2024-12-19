
# --------------- Functions for high-resolution predictors -----------------
#' Iteratively extract values for a raster for a set of locations
#'
#' Takes paths a inputs, iteratively run terra::extract on a given subet of 
#' locations, and write the extracted values to fst format in output directory.
#'
#' @param in_rast_path (character) path to raster file to be read by terra::rast.
#' Note that this function is most efficient if run on multi-band rasters
#' (netcdf or tiff) rather than running separately on the bands.
#' @param lyrs (integer vector) indices of raster layers (bands) to extract.
#' @param in_coordtab_path (character) path to table of locations to use for 
#' extraction.Should have at least three columns: an ID col and two 
#' coordinate columns.
#' @param dropcol (integer vector) vector of column indices to drop when 
#' reading table of locations (for faster loading).
#' @xcol (character) name of the numeric column holding longitude coordinates
#' @ycol (character) name of the numeric column holding latitude coordinates
#' @outdir (character) path to directory where function outputs will be written
#' @out_pathroot (character) base name of the files that will be written
#' @iterstep (integer) number of locations to process at a time. Increasing this
#' number will speed up the analysis but require more memory.
#' @overwrite (boolean) whether to run extraction even if output file already exists
#' 
#' @details the output is written as {fst}. It can be read back as a data.frame
#' with read_fst(path).
#'
#' @return \link[data.table]{data.table} of paths to output fst files
#'
#' @export
iter_ras_extract <- function(in_rast_path, lyrs=1,
                             in_coordtab_path, dropcol=NULL, idcol="DRYVER_RIVID",
                             xcol='X', ycol='Y',
                             out_dir, out_pathroot, iterstep,
                             overwrite=F) {
  
  #Read raster to extract
  # ras <- terra::rast(in_rast_path, lyrs=lyrs)
  nc_file <- ncdf4::nc_open(in_rast_path)
  dis_array <- ncdf4::ncvar_get(nc_file, "dis")
  fillvalue <- ncdf4::ncatt_get(nc_file, "dis", "_FillValue")
  dis_array[dis_array == fillvalue$value] <- NA
  ras <- raster::raster(t(dis_array), xmn=-25, xmx=70,
                        ymn=12, ymx=84,
                        crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  #Read points, setting column names for coordinates to work with extract
  pts <- data.table::fread(in_coordtab_path, drop=1) %>%
    data.table::setnames(c(xcol, ycol), c('X', 'Y'))
  
  #Get number of points
  pts_n <- nrow(pts)
  
  #Create intervals of rows to select at a time
  if (pts_n > iterstep) {
    binlist <- seq(1, pts_n, iterstep) %>%
      data.table(
        l = .,
        u = c(.[2:length(.)], pts_n+1)-1
      )
  } else {
    binlist <- data.table(l = 1, u = pts_n)
  }
  
  #Iterate through ID intervals
  outputf_list <- lapply(
    seq(1, nrow(binlist)), function(i) {
      #Create output file path
      outf <- file.path(out_dir,
                        paste0(out_pathroot,
                               '_', format(binlist[i, l], scientific = FALSE),
                               "_", format(binlist[i, u], scientific=FALSE),
                               '.fst')
      )
      #Check whether it exists
      if ((!file.exists(outf)) | (overwrite == T)) {
        print(paste('Processing', outf))
        #Extract raster values by point, merge with original ID column
        raster::extract(
          x = ras,
          y = pts[binlist[i, l]:binlist[i, u], .(X, Y)],
          method="simple"
        ) %>%
          cbind(pts[binlist[i, l]:binlist[i, u], .(get(idcol))], value=.) %>%
          write_fst(path=outf) #Write it all to fast read/write format
      }
      return(outf)
    })
}


#' extract high resolution downscaled waterGAP streamflow for pour points
#'
#' Takes paths of downscaled streamflow NetCDF files and extract it using
#' `iter_ras_extract` function.
#'
#' @param ncs_dir (character) path to NetCDF file to be read by ncdf4::nc_open.
#' @param out_dir (character) path to the output file of the function
#' @param prpts_path (character) path to the csv file including the location of 
#' the most downstream point of a reach segment
#' @param model_name (character) the name of Global Circulation Model (GCM)
#' @param phase (character) the period phase of the data, such as `hist`, `ssp126`,
#' and `ssp585`.
#' @return the output is written as {fst}. It can be read back as a data.frame
#' with read_fst(path).
#' 
#' @export
#' 
extract_dsflow_hr <- function(ncs_dir, out_dir, prpts_path, model_name, phase){
  
  #list the ncfiles of downscaled streamflow
  files_ls <- list.files(ncs_dir, 
                         pattern = '.nc4', 
                         full.names = TRUE) %>% 
    gtools::mixedsort(.)
  
  # create the out dir, if it isn't existed
  if (!file.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  # extract the debi for the pourpoint at the reaches
  lapply(seq_along(files_ls), function(i){
    
    file_date <- files_ls[i] %>% basename(.) %>% gsub('.nc4', '',.) %>%
      strsplit(., '_') %>% unlist()
    
    output_name <- paste0('watergap_dis_net_', paste(file_date[3:4], collapse = "_"))
    
    iter_ras_extract(in_rast_path = files_ls[i],
                     lyrs = 1,
                     in_coordtab_path = prpts_path,
                     dropcol = 'OBJECTID',
                     idcol = 'DRYVER_RIVID',
                     xcol = 'POINT_X',
                     ycol = 'POINT_Y',
                     out_dir = out_dir,
                     out_pathroot = output_name,
                     iterstep = 1600000)
  })
  
  #merge all the fst into one and rename the cols Q_mon_year
  files_ls_fst <- list.files(out_dir,
                             pattern='.fst',
                             full.names=TRUE)
  #iterate over the fst files to have only one data.table
  q_watergap <- lapply(seq_along(files_ls_fst), function(i){
    
    file_date <- files_ls_fst[i] %>% basename(.) %>% gsub('.fst', '',.) %>%
      strsplit(., '_') %>% unlist()
    
    output_name <- paste0('Q_', paste(file_date[4:5], collapse = "_"))
    
    fst::read_fst(files_ls_fst[i]) %>% 
      as.data.table() %>% 
      .[,.(value)] %>% 
      data.table::setnames(output_name)
    
  }) %>% 
    do.call('cbind',.)
  
  #get the river IDs
  dryver_ids <- fst::read_fst(files_ls_fst[1]) %>% 
    as.data.table() %>% 
    .[,.(V1)]
  
  #remove the created fst files from directory
  lapply(seq_along(files_ls_fst), function(i){
    
    file.remove(files_ls_fst[i])
    
  })
  
  #combine the data.table and IDs of reaches and save in one file
  q_watergap_dt <- cbind(DRYVER_RIVID=dryver_ids$V1, q_watergap)
  
  q_watergap_dt %>% 
    write_fst(., path = file.path(out_dir,
                                  paste0('Q_dswatergap_',model_name, '_', phase, '.fst')))
  
  return(NULL)
}

#' Calculate the high-resolution predictors 
#' 
#' @param in_path (character) path to the time series of downscaled streamflow in {fst} file
#' produced by `extract_dsflow_hr` function.
#' @param out_path (character) path to the time series of the high-resultion predictors
#' @param n_window (numeric) a number (here: 3 or 12) indicates the window that the 
#' predictor is calculating.
#' @param Fun (character) a character between `mean` and `min` shows the function that
#' is using to compute the predictor
#' @param pred_path (character) the predictor name
#' @param out_name (character) the name of the predictor in output file
#' @param start_year (numeric) a number that shows the start year of the interesed period
#' @param end_year (numeric) a number that shows the last year of the interesed period
#' @param start_date (character) a character shows the start date of the interesed period
#' @param end_date (character) a character shows the last date of the interesed period
#' 
#' @return the output is written as {fst}. It can be read back as a data.frame
#' with read_fst(path).
#' 
#' @export
#' 
compute_hr_predictors <- function(in_path, out_path,n_window=3, Fun='mean',
                                  pred_name = 'mean_p3', out_name = 'q_mean_p3_eu',
                                  start_year=1985, end_year=2014,
                                  start_date='1985-01-01', end_date='2014-12-01') {
  file_ls <- list.files(path = in_path, pattern = '.fst', full.names = TRUE)
  
  if (!dir.exists(out_path)){
    dir.create(out_path, recursive = TRUE)
  }
  
  #read the file
  data <- read_fst(file_ls)
  
  fst::write_fst(data, path = file.path(out_path, 'watergap_dis_net_eu.fst'))
  #create a vector of date between the period
  create_dates <- function(start_year, end_year){
    all_dates <- c()
    
    # Loop through each month (from January to December)
    for (month in c(1,10,11,12, 2:9)) {
      # Loop through each year
      for (year in start_year:end_year) {
        # Create the date for the first day of the month
        date <- as.Date(paste(year, month, "01", sep = "-"))
        # Append the date to the vector
        all_dates <- c(all_dates, date)
      }
    }
    
    all_dates <- as.Date(all_dates, "1970-01-01") %>% as.character()
    return(all_dates)
  }
  #create a vector of dates within a given period
  dates <- create_dates(start_year = start_year, end_year = end_year)
  
  
  old_names <- names(data)
  new_names <- c('DRYVER_RIVID', dates)
  data.table::setnames(data,old_names, new_names)
  data.table::setDT(data)
  
  new_order <- seq.Date(as.Date("1981-01-01"),
                        as.Date("2019-12-31"),
                        "month") %>% as.character()
  new_order <- c('DRYVER_RIVID', new_order)
  data <- data[,..new_order]
  # ----------- the function for P3month ------------
  meanpmon_dt <- lapply(seq_along(1:dim(data)[1]), function(i){
    if (i %% 100000 == 0){
      cat('extracting of HR predictor was completed for ', round(i/dim(data)[1] * 100, digits = 2) , '% of reaches...', '\n')
    }
    
    dd <- frollapply(unlist(data[i, -1]), n = n_window, FUN = Fun,
                     align = "right", fill = NA)
    #to fill the NA out with the nearest values.
    # here we filled out the first three month with the corresponding values in the next year.
    near_values <- dd[13:(12+n_window)]
    dd[1:n_window] <- near_values
    dd
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  create_names <- function(start_date='1985-01-01', end_date='2014-12-01'){
    date_mon <- seq.Date(as.Date(start_date), as.Date(end_date), 'month')
    final_names <- paste0(pred_name,'_', lubridate::month(date_mon),
                          "_", lubridate::year(date_mon))
    return(final_names)
  }
  
  old_names <- names(meanpmon_dt)
  final_names <- create_names(start_date=start_date, end_date=end_date)
  setnames(meanpmon_dt, old_names, final_names)
  cbind('DRYVER_RIVID'=order_vector,meanpmon_dt) %>% 
    fst::write_fst(., path = file.path(out_path, paste0(out_name,'.fst')))
  
  return(NULL)
}


#' Calculate the high-resolution inter-annual predictors 
#' 
#' The standard deviation (sd) and coeficient of variation (cv) of
#' the downscaled streamflow are computing over three periods for
#' the five GCMs.
#' 
#' @param in_path (character)  path to the time series of downscaled streamflow in {fst} file
#' produced by `extract_dsflow_hr` function.
#' @param out_path (character) path to the time series of the high-resultion predictors
#' @param start_year (numeric) a number that shows the start year of the interesed period
#' @param end_year (numeric) a number that shows the last year of the interesed period
#' 
#' @return the output is written as {fst}. It can be read back as a data.frame
#' with read_fst(path).
#' 
#' @export
#' 
compute_hr_interannual_predictors <- function(in_path, out_path,
                                              start_year=1985, end_year=2014) {
  file_ls <- list.files(path = in_path, pattern = '.fst', full.names = TRUE)
  #read the file
  data <- read_fst(file_ls)
  #create a vector of date between the period
  create_dates <- function(start_year, end_year){
    all_dates <- c()
    
    # Loop through each month (from January to December)
    for (month in c(1,10,11,12, 2:9)) {
      # Loop through each year
      for (year in start_year:end_year) {
        # Create the date for the first day of the month
        date <- as.Date(paste(year, month, "01", sep = "-"))
        # Append the date to the vector
        all_dates <- c(all_dates, date)
      }
    }
    
    all_dates <- as.Date(all_dates, "1970-01-01") %>% as.character()
    return(all_dates)
  }
  #create a vector of dates within a given period
  dates <- create_dates(start_year = start_year, end_year = end_year)
  
  old_names <- names(data)
  new_names <- c('DRYVER_RIVID', dates)
  data.table::setnames(data,old_names, new_names)
  
  data.table::setDT(data)
  dm <- as.Date(dates, format = "%Y-%m-%d")
  new_order <- dm[order(dm)] %>% as.character()
  
  date_vector <- as.Date(new_order, format = "%Y-%m-%d") %>% lubridate::month()
  new_order <- c('DRYVER_RIVID', new_order)
  data <- data[,..new_order]
  
  # ----------- the function for CV ------------
  cvmon_dt <- lapply(seq_along(order_vector), function(i){
    if (i %% 100000 == 0){
      cat('extracting of CV (near future) of HR predictor was completed for ',
          round(i/dim(data)[1] * 100, digits = 2) , '% of reaches...', '\n')
    }
    
    dd <- data[i,-1] %>% t() %>%
      as.data.table()
    dd[, date:=date_vector] %>% .[date <= as.Date('2070-12-01')] %>% 
      dd[,cv_mon:=sd(V1, na.rm=TRUE)/mean(V1, na.rm=TRUE), by='date']
    
    dd[1:12,sd_mon]
    
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  old_names <- names(cvmon_dt)
  final_names <- paste0('cv_', 1:12)
  setnames(cvmon_dt, old_names, final_names)
  cbind('DRYVER_RIVID'=order_vector,cvmon_dt) %>% 
    fst::write_fst(., path = file.path(out_path, 'q_iav_cv_eu_nearfuture.fst'))
  
  cvmon_dt <- lapply(seq_along(order_vector), function(i){
    if (i %% 100000 == 0){
      cat('extracting of CV (far future) of HR predictor was completed for ',
          round(i/dim(data)[1] * 100, digits = 2) , '% of reaches...', '\n')
    }
    
    dd <- data[i,-1] %>% t() %>%
      as.data.table()
    dd[, date:=date_vector] %>% .[date > as.Date('2070-12-01')] %>% 
      dd[,cv_mon:=sd(V1, na.rm=TRUE)/mean(V1, na.rm=TRUE), by='date']
    
    dd[1:12,cv_mon]
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  old_names <- names(cvmon_dt)
  final_names <- paste0('cv_', 1:12)
  setnames(cvmon_dt, old_names, final_names)
  cbind('DRYVER_RIVID'=order_vector,cvmon_dt) %>% 
    fst::write_fst(., path = file.path(out_path, 'q_iav_cv_eu_farfuture.fst'))
  # ----------- the function for sdt ------------
  cvmon_dt <- lapply(seq_along(order_vector), function(i){
    if (i %% 100000 == 0){
      cat('extracting of sd (near future) of HR predictor was completed for ',
          round(i/dim(data)[1] * 100, digits = 2) , '% of reaches...', '\n')
    }
    
    dd <- data[i,-1] %>% t() %>%
      as.data.table()
    dd[, date:=date_vector] %>% .[date <= as.Date('2070-12-01')] %>% 
    .[,sd_mon:=sd(V1, na.rm=TRUE), by='date']
    
    dd[1:12,sd_mon]

  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  old_names <- names(cvmon_dt)
  final_names <- paste0('sd_', 1:12)
  setnames(cvmon_dt, old_names, final_names)
  cbind('DRYVER_RIVID'=order_vector,cvmon_dt) %>% 
    fst::write_fst(., path = file.path(out_path, 'q_iav_sd_eu_nearfuture.fst'))
  
  cvmon_dt <- lapply(seq_along(order_vector), function(i){
    if (i %% 100000 == 0){
      cat('extracting of sd (far future) of HR predictor was completed for ',
          round(i/dim(data)[1] * 100, digits = 2) , '% of reaches...', '\n')
    }
    
    
    dd <- data[i,-1] %>% t() %>%
      as.data.table()
    dd[, date:=date_vector] %>% .[date > as.Date('2070-12-01')] %>% 
      .[,sd_mon:=sd(V1, na.rm=TRUE), by='date']
    
    dd[1:12,sd_mon]
  }) %>% 
    do.call('rbind',.) %>% 
    as.data.table()
  
  old_names <- names(cvmon_dt)
  final_names <- paste0('sd_', 1:12)
  setnames(cvmon_dt, old_names, final_names)
  cbind('DRYVER_RIVID'=order_vector,cvmon_dt) %>% 
    fst::write_fst(., path = file.path(out_path, 'q_iav_sd_eu_farfuture.fst'))
  
  return(NULL)
}


## --------------- Functions for low-resolution predictors -----------------
#' Count the number of wet days within months
#' 
#' @description this function calculate the time series of number of wet days within a month
#' 
#' @param in_path (character) path to the the NetCDF files of wet days
#' @param out_path (character) path to the {fst} files indicating the 
#' number of wet days within month for European reaches
#' @param start_date (character) a character shows the start date of the interesed period
#' @param end_date (character) a character shows the last date of the interesed period
#' @param threshold (numeric) a value that represents the threshold for considering the grid-cell
#' is wet or not. `default =` 2.5 mm/day.
#' 
#' #' @return the output is written as stack {raster} similar to NetCDF.
#' 
#' @export
#' 
count_wetdays <- function(in_path, out_path, start_date = '1981-01-01', 
                          end_date='2014-12-31', threshold = 2.5){
  
  if (inherits(start_date, 'character') | inherits(end_date, 'character')){
    start_date <-as.Date(start_date)
    end_date <-as.Date(end_date)
  }
  #sequence over the start date and end date
  dates <- seq.Date(start_date, end_date, 'month') 
  
  #create a list of nc files in a directory
  files_ls <- list.files(path = in_path,
                         pattern = '.nc',
                         full.names = TRUE)
  #read the ncfiles and stack them together
  pre_data <- lapply(seq_along(files_ls), function(i){
    r <- raster::stack(files_ls[i])
  }) %>% 
    do.call('stack',.)
  #reclassify the cells based on the defined threshold
  classified_raster <- lapply(1:nlayers(pre_data), function(i){
    
    if (i %% 1000 == 0) {
      cat('Reclassifying process of layer was completed for ', round(i/nlayers(pre_data) * 100, digits = 2) , '% of layers...', '\n')
    }
    r <- pre_data[[i]] * 86400 #convert the precipitation values to mm/day
    raster::reclassify(r, c(2.5, Inf, 1, 0, 2.5, 0))
  }) %>% 
    do.call('stack',.)
  
  #reterive the number of days within a month
  num_days <- lubridate::days_in_month(dates) %>%
    as.vector() %>%
    cumsum() #sum the number of days over the period
  #count the number of wet days per each cell
  output <- lapply(seq_along(num_days), function(i){
    if (i ==1) {
      r <- classified_raster[[1:31]]
    } else {
      st <- num_days[i-1] +1
      end <- num_days[i]
      r <- classified_raster[[st:end]]
    }
    sum(r)
  }) %>% 
    do.call('stack',.)
  
  #remove the unnecessary raster and clear the memory
  rm(classified_raster)
  rm(pre_data)
  gc()
  #split the filename by underscore
  split_string <- unlist(strsplit(files_ls[1] %>% basename(), "_"))
  # Select the portion until the fifth underscore
  selected_portion <- paste(split_string[1:4], collapse = "_")
  output_filename <- paste0(selected_portion,'_wetdays.nc')
  #save the output stack into nc 
  writeRaster(output, filename = file.path(out_path,output_filename), format = "CDF",
              xname="Longitude",  yname="Latitude", zname="Time (Month)",
              varname="Number Wet Days", overwrite=TRUE)
  
  return(output)
  
  
}


#' qrdif_ql_ratio_mon
#' 
#' @description 
#' a raster of the ratio of diffuse groundwater recharge to runoff from land in 0.5*0.5 degree spatial 
#' resolution is calculated. this predictor is considered as a low resolution predictor.
#' 
#' @param path_in (character) the data source name of the nc files
#' @param path_out (character) the directory to the output file
#' @param phase_time (character) to create a sequence of time span within the period
#' 
qrdif_ql_ratio_mon <- function(path_in, path_out, phase_time='hist'){
  
  if (phase_time == 'hist') {
    dates <- seq.Date(as.Date('1985-01-01'),
                      as.Date('2014-12-01'), 
                      "month")
  }
  else if (phase_time == 'future') {
    dates <- seq.Date(as.Date('2041-01-01'),
                      as.Date('2100-12-01'), 
                      "month")
  }
  
  #import the files
  files_ls <- list.files(path = path_in, pattern = ".nc", full.names = TRUE)
  all_ncs <- lapply(files_ls, function(x) stack(x))
  #take a subset of the larger file to make the both at the same time span
  select_names <- all_ncs[[1]] %>% names()
  all_ncs[[2]] <- raster::subset(all_ncs[[2]], select_names)
  
  #overlay the diffuse groundwater recharge to runoff from land 
  qrdif_ql_ratio_mon <- raster::overlay(all_ncs[[2]], all_ncs[[1]], 
                                        fun = function(x, y) {
                                          return((x / y))
                                        })
  
  #replace NA values with zero in the entire stack, more likelihood to have intermittency
  qrdif_ql_ratio_mon[is.na(qrdif_ql_ratio_mon)] <- 0
  
  #set the names for the raster file
  names(qrdif_ql_ratio_mon) <- dates
  
  #split the filename by underscore
  split_string <- unlist(strsplit(files_ls[1] %>% basename(), "_"))
  # Select the portion until the fifth underscore
  selected_portion <- paste(split_string[1:5], collapse = "_")
  output_filename <- paste0(selected_portion,'_qrdifoverql.nc')
  #export stack file into a NetCDF file
  writeRaster(qrdif_ql_ratio_mon, filename = file.path(path_out,output_filename), format = "CDF",
              xname="Longitude",  yname="Latitude", zname="Time (Month)",
              varname="Ratio qrd over ql", overwrite=TRUE)
  
  return(qrdif_ql_ratio_mon)
}

#' Calculate accumulated low-resolution predictors
#' 
#' @description This function computes the number of wet days and the ratio of ground water to runoff
#' for pour point at the downstream of reach segments.
#' 
#' @param in_dir (character) path to the accumulated LR predictors in {csv} format
#' @param reachids_dir (character) path to the reach DRYvER IDs in order to select 
#' the European reaches from the csv file.
#' @param out_dir (character) path to the time series of the low-resultion predictors
#' @param period_phase (character) the period phase of the data, options: `hist` and `future`.
#' @param var_name (character) the name of the low-resolution predictor. options: `wetdays` and `qrdifoverql`
#' 
#' @return the output is written as {fst}. It can be read back as a data.frame
#' with read_fst(path).
#' 
#' @export
#' 
select_modify_lr_predictors <- function(in_dir, reachids_dir, out_dir,
                                        period_phase = 'hist', var_name='wetdays'){
  reach_ids <- data.table::fread(reachids_dir) %>% 
    .[,DRYVER_RIVID]
  
  lr_tb <- data.table::fread(input = in_dir)
  
  lr_tb_selected <- lr_tb[DRYVER_RIVID %in% reach_ids][,-1]
  lr_tb_selected <- unique(lr_tb_selected, by='DRYVER_RIVID')
  lr_tb_selected <- lr_tb_selected[order(DRYVER_RIVID)]
  
  set_colnames_lr <- function(lr_tb_selected, period_phase = 'hist', var_name='qrdifoverql'){
    if (period_phase == 'hist') {
      years <- 1985:2014
    }else
      years <- 2041:2100
    
    if (var_name == 'qrdifoverql') {
      col_names <- vector("character")
      for(i in seq_along(years)){
        for(j in 1:12){
          col_names[(i-1)*12 + j] <- paste0("gwr_",j,"_",years[i])
        }
      }
    }else {
      col_names <- vector("character")
      for(i in seq_along(years)){
        for(j in 1:12){
          col_names[(i-1)*12 + j] <- paste0("wetdays_",j,"_",years[i])
        }
      }
    }
    # output of the function
    lr_tb_selected <- lr_tb_selected %>% 
      `colnames<-`(c("DRYVER_RIVID", col_names))
  }
  
  if (var_name == 'wetdays' & period_phase == 'hist'){
    lr_tb_selected <- set_colnames_lr(lr_tb_selected,
                                      period_phase = 'hist', var_name='wetdays')
    lr_tb_selected[, (names(lr_tb_selected)[-1]) := lapply(.SD, `/`, 100), .SDcols = -"DRYVER_RIVID"]
    
  } else if (var_name == 'wetdays' & period_phase == 'future') {
    lr_tb_selected <- set_colnames_lr(lr_tb_selected,
                                      period_phase = 'future', var_name='wetdays')
    lr_tb_selected[, (names(lr_tb_selected)[-1]) := lapply(.SD, `/`, 100), .SDcols = -"DRYVER_RIVID"]
    
  } else if (var_name == 'qrdifoverql' & period_phase == 'hist') {
    lr_tb_selected <- set_colnames_lr(lr_tb_selected,
                                      period_phase = 'hist', var_name='qrdifoverql')
    lr_tb_selected[, (names(lr_tb_selected)[-1]) := lapply(.SD, `/`, 1000000), .SDcols = -"DRYVER_RIVID"]
  } else if (var_name == 'qrdifoverql' & period_phase == 'future'){
    lr_tb_selected <- set_colnames_lr(lr_tb_selected,
                                      period_phase = 'future', var_name='qrdifoverql')
    lr_tb_selected[, (names(lr_tb_selected)[-1]) := lapply(.SD, `/`, 1000000), .SDcols = -"DRYVER_RIVID"]
  }
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (var_name == 'qrdifoverql') {
    fst::write_fst(x = lr_tb_selected, path = file.path(out_dir, "gwr_to_runoff_ratio_eu.fst"))
  } else {
    fst::write_fst(x = lr_tb_selected, path = file.path(out_dir, "wetdays_net_eu.fst"))
  }
  
  return(NULL)
}

