############################## Trend methods time seores #################### 
##
## Functions generated through various research projects (STA) and SESYNC research support.
## Performing trend analyses of time series data to with theil sen, OLS and Mann Kendall.
##
## DATE CREATED: 08/11/2017
## DATE MODIFIED: 05/15/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis 
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

####### This script contains the following functions:
#
#1) calculate_theil_sen_time_series

############################
###Loading R library and packages                                                      
#library(gstat) #spatial interpolation and kriging methods
library(sp) # spatial/geographic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
library(readxl) #functionalities to read in excel type data
library(sf) # spatial ojbects simple feature model implementation OGC
#library(gstat)
#library(spacetime)
library(mblm)

###### Functions used in this script and sourced from other files

calculate_trend <- function(y,mod_obj=F,method="theil_sen"){
  #y,n,harmonic_val=NULL,mod_obj=F,figure=F
  #This function generates Theil Sen slope estimate using mblm package.
  #
  #INPUTS
  #1) y:  y variable (trend variable)
  #2) mod_obj: save mod obj?
  #4) method: theil sen or ols option for trend
  #
  #OUTPUTS
  #1)
  #
  #TO DO: 
  
  ############# BEGIN FUNCTION #############
  
  #setting up the data input: data.frame with 2 columns
  time_index <- 1:length(y) #

  #stores dates for later processing
  #if(inherits(data_df,"zoo")){
  #  dates_val <- date(df_mblm)
  #}else{
  #  dates_val <- NULL
  #}
  
  dates_val <- NULL # set this later
  
  df_val <- data.frame(y=y,time_index=time_index) #convert to data.frame since it was a zoo df
  #df_mblm$time_index <- time_index
  #names(df_mblm)
  
  if(method=="theil_sen"){
    df_val <- na.omit(df_val) # removes NA, might need to think about that later on
    ### automate formula input
    formula_str <- paste0(names(df_val)[1]," ~ ","time_index")
    formula_val <- as.formula(formula_str) #transform object into formula
    
    mod_obj <- mblm(formula_val,df_val)
  }
  
  if(method=="ols"){
    #df_val <- na.omit(df_val) # removes NA, might need to think about that later on
    ### automate formula input
    formula_str <- paste0(names(df_val)[1]," ~ ","time_index")
    formula_val <- as.formula(formula_str) #transform object into formula
    
    mod_obj <- lm(formula_val,df_val)
    
  }
  
  ##### Extract information from model object mblm
  
  #slope_theil_sen <- coef(mod_mblm)[2]
  #intercept_theil_sen <- coef(mod_mblm)[1]
  
  slope <- coef(mod_obj)[2]
  intercept <- coef(mod_obj)[1]
  
  #ID_ts <- subset_name
  #method_str <- "theil_sen"
  n_obs <- sum(!is.na(y))
  
  if(!is.null(dates_val)){
    nt <- length(dates_val)
    start_date <- dates_val[1]
    end_date <- dates_val[nt]
  }else{
    nt <- NA
    start_date <- NA
    end_date <- NA
  }
  
  n <- length(y)
  slope_sign <- sign(slope)
  #slope_sign <- sign(slope_theil_sen)
  
  df_trend <- data.frame(intercept=intercept,
                             slope=slope,
                             slope_sign=slope_sign,
                             method=method,
                             start_date= start_date,
                             end_date = end_date,
                             n=n,
                             n_obs=n_obs)#number of valid observation (not NA) 
  
  #### Prepare object to return
  trend_obj <- list(mod_obj,df_trend)
  names(trend_obj) <- c("mod_obj","df_trend")
  
  ##### save to disk
  
  #if(mod_obj==TRUE){
    #obj_theil_sen_filename <- file.path(out_dir,paste("theil_sen_obj_",subset_name,"_",out_suffix,".RData",sep=""))
    #save(obj_theil_sen,file= obj_theil_sen_filename)
  #} 
  
  return(trend_obj)
}

trend_reg_raster <- function(y,var_name,method="theil_sen"){
  #
  ## Function to generate outputs from harmonic regression for every pixel in a raster image.
  ## Note that the output needed for the calc function used is a vector of values.
  ## The output values can be amplitudes (A0, A1, A2), phases, p significance etc.
  #
  ### INPUTS:
  #1) y: input data, in this function the name of raster layer is expected
  #2) var_name: number of elements in the time series/sequence used for the harmonic reg.
  #4) method: theil_sen or ols
  ### OUTPUTS
  #1) value_var_name: numeric vector of values from harmonic modeling
  
  ####### Start script #######
  
  #n <- layers(y)
  
  trend_results <- calculate_trend(y,
                                   mod_obj=F,
                                   method="theil_sen")
    
  #df_in <- subset(harmonic_results$harmonic_df,harmonic==harmonic)
  df_in <- trend_results$df_trend
  
  ## Select variables to predict for every pixels
  value_var_name <- df_in[[var_name]]
  
  ## Return values for raster:
  return(value_var_name)
}

calcTrendRaster <- function(r,method="theil_sen",var_name="slope",file_format=".tif",multiband=F,num_cores=1,raster_name=NULL,out_dir=NULL){
  
  #var_name <- "A"
  #raster_name <- NULL
  #raster_name <- "NDVI_amplitude.tif"
  #file_format <- ".tif"
  #multiband=F
  
  ############ Start script ##################
  
  r_out <- try(calc(r, 
                    fun=function(y){trend_reg_raster(y,
                                                     var_name=var_name,
                                                     method=method)}))
  
  names(r_out) <- paste0(var_name,"_",method)
      
  ########## Write out
      
  if(is.null(raster_name)){
    out_raster_name <- paste(var_name,"_",i,file_format,sep="") 
  }else{
    out_raster_name <- sub(file_format,"",raster_name)
    out_raster_name <- paste(out_raster_name,"_",i,file_format,sep="") 
  }
      
      #if(multiband==TRUE){
      #  #raster_name_tmp <- basename(rast_name_var)
      #  #raster_name <- basename(sub(file_format,"",raster_name))
      #  if(out_suffix!=""){
      #    raster_name_tmp <- paste(raster_name,"_",out_suffix,file_format,sep="")
      #  }else{
      #    raster_name_tmp <- paste(raster_name,file_format,sep="")
      # }
      #  bylayer_val <- FALSE #don't write out separate layer files for each "band"
      #  rast_list <- file.path(out_dir,raster_name_tmp) #as return from function
      #}
      suffix_str <- names(r_out)
      
  if(multiband==FALSE){
    out_raster_name <- sub(file_format,"",out_raster_name)
    raster_name_tmp <- paste(out_raster_name,file_format,sep="") #don't add output suffix because in suffix_str
    bylayer_val <- TRUE #write out separate layer files for each "band"
    rast_list <- file.path(out_dir,(paste(out_raster_name,"_",suffix_str,file_format,sep=""))) 
  }
      
  #### return object
  
  return(rast_list)
}

########################## End of script #######################################

