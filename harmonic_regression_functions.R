############################## Harmonic regression #################### 
##
## Functions generetaed through various research projects and SESYNC research support.
## Performing harmonic regression time series data to evaluate amplitudes and phases for Managing Hurriance Group.
##
## DATE CREATED: 10/01/2018
## DATE MODIFIED: 05/28/2019
## AUTHORS: Benoit Parmentier
## Version: 1
## PROJECT: Time series analysis Managing Hurricanes
## ISSUE: 
## TO DO:
##
## COMMIT: exploration of estimation
##

####### This script contains the following functions:
#
#1) fit_harmonic
#2) harmonic_regression
#3) split_sequence
#4) extract_harmonic_coef
#5) harmonic_reg_raster
#

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

###### Functions used in this script and sourced from other files

#https://stats.stackexchange.com/questions/60500/how-to-find-a-good-fit-for-semi-sinusoidal-model-in-r


### This needs to be modified.
#SSTlm2 <- lm(Degrees ~ sin(2*pi*ToY)+cos(2*pi*ToY)
#             +sin(4*pi*ToY)+cos(4*pi*ToY),data=SST)
#summary(SSTlm2)

fit_harmonic <- function(p,n,y,mod_obj=F,figure=F){
  
  t <- 1:n
  omega_val=2*pi*p/n #may be more than 1
  
  cos_val <- lapply(omega_val,function(omega){cos(omega*t)})
  sin_val <- lapply(omega_val,function(omega){sin(omega*t)})
  
  cos_df <- as.data.frame(do.call(cbind,cos_val))
  names(cos_df) <- paste("cos",p,sep="")
  sin_df <- as.data.frame(do.call(cbind,sin_val))
  names(sin_df) <- paste("sin",p,sep="")
  
  in_df <- data.frame(y_var = y)
  in_df <- cbind(in_df,cos_df,sin_df)
  #View(in_df)
  #cos_val =cos(omega*t)
  #sin_val =sin(omega*t)
  
  #omega_val = lapply(p,function(p){2*pi*p/n})
  
  #plot(cos_val)
  #plot(sin_val)
  
  #y_var <- "invasion.status"
  #in_dir[[y_var]] <- as.factor(data[[y_var]]) #this is needed for randomForest to get a classification
  
  explanatory_variables <- names(in_df)[-1] #drop the first column
  
  right_side_formula <- paste(explanatory_variables,collapse = " + ")
  model_formula_str <- paste0("y_var"," ~ ",right_side_formula)
  
  #in_df <- data.frame(y=y,cos_val=cos_val,sin_val=sin_val)
  mod <- try(lm(model_formula_str ,data=in_df),silent = T)
  if(!inherits(mod,"try-error")){
    summary(mod)
    #mod2 <- lm(model_formula_str ,data=in_df,na.action = na.exclude)
    #mod$coefficients
    #mod2$coefficients
    #mod2$model
    
    ### Extract coefficients for each harmonic
    
    p_val<-2
    #p
    #debug(extract_harmonic_coef)
    test_df<- extract_harmonic_coef(p_val,n,mod)
    
    harmonic_df <- lapply(p,FUN=extract_harmonic_coef,n=n,mod=mod)
    
    ### Figure
    if(figure==TRUE){
      y_range <- range(mod$fitted.values,y,na.rm = T)
      plot(mod$fitted.values,ylim=y_range)
      points(y,col="blue",pch="+")
    }
    ###
    harmonic_obj <- harmonic_df
    
  }else{
    harmonic_df <- lapply(1:length(p),
           FUN=function(i) {harmonic_df <- data.frame(A0=NA,A=NA,a=NA,b=NA,
                              pr_A0=NA,pr_a=NA,pr_b=NA,
                              phase=NA,harmonic=NA,omega=NA)}
    )
    
    
  }
    
  
  if(mod_obj==T){
    harmonic_obj <- list(harmonic_df=harmonic_df,mod=mod)
  }else{
    harmonic_obj <- list(harmonic_df=harmonic_df)
  }
  ### 
  return(harmonic_obj)
}

harmonic_regression <- function(y,n,harmonic_val=NULL,mod_obj=F,figure=F){
  ##
  # if mod_obj is True then return the model object 
  #
  ## Default to two first harmonic:
  if(is.null(harmonic_val)){
    p <- 1:2
  }else{
    p <- 1:harmonic_val
  }
  
  #harmonic_val = 1 # pr from 1 to n/2
  
  #n<-24
  
  #l_harmonic_obj <- lapply(p,
  #       FUN=fit_harmonic,
  #       n=n,
  #       y=y,
  #       mod_obj=mod_obj,
  #       figure=figure)
  
  #debug(fit_harmonic)
  l_harmonic_obj <- try(fit_harmonic(p,n,y,mod_obj=mod_obj,figure=F),silent=TRUE)
  #handle error:
  l_df <- lapply(p,function(i){l_harmonic_obj$harmonic_df[[i]]})
  harmonic_df <- do.call(rbind,l_df)
  rownames(harmonic_df) <- NULL
    
  #View(harmonic_df)
  harmonic_results_obj  <- list(harmonic_df,l_harmonic_obj)
  names(harmonic_results_obj) <- c("harmonic_df","harmonic_obj")
  
  return(harmonic_results_obj)
}

split_sequence <- function(x,n,overlap=0){
  
  if(overlap==0){
    n_splits <- floor(length(x)/n)
    #n_modified <- n - overlap
    intervals_val <- seq(1,to=length(x),by=n)
    intervals_val <- c(intervals_val,length(x))
    n_splits
    list_intervals <- lapply(2:length(intervals_val),function(i){data.frame(start=intervals_val[[i-1]],end=intervals_val[[i]])})
    intervals_df <- do.call(rbind,list_intervals)
    intervals_df
    #lapply(2:length(intervals_val),function(i){intervals_val[[i]]-overlap})
    #length(intervals_val)
    
  }
  ##implement the other option later
  
  ## now split:
  list_y <- lapply(1:nrow(intervals_df),function(i){x[intervals_df[i,]$start:intervals_df[i,]$end-1]})
  
  split_obj <- list(list_y,intervals_df)
  names(split_obj) <- c("list_y","intervals")
  return(split_obj)
}


extract_harmonic_coef <- function(p_val,n,mod){
  
  summary(mod)
  coef_df <- (as.data.frame(t(mod$coefficients)))
  
  cos_term <- paste("cos",p_val,sep="")
  sin_term <- paste("sin",p_val,sep="")
  
  a <- coef_df[[sin_term]] #sine term
  b <- coef_df[[cos_term]] #cosine term
  A0 <- coef_df[["(Intercept)"]] #mean
  p_significance <- coef_df[["(Intercept)"]] #mean
  
  A = sqrt(a^2 + b^2)
  phase = atan(-b/a)
  ## Add p values later?
  #n <- nrow(mod$model)
  omega_val=2*pi*p_val/n #may be more than 1
  
  #class((summary(mod))$coefficients)
  results_df <- (as.data.frame(summary(mod)$coefficients))
  results_df <- results_df[rownames(results_df)%in%c(sin_term,cos_term,"(Intercept)"),]
  
  pr <- results_df$`Pr(>|t|)`
  pr_a <- pr[3] 
  pr_b <- pr[2]
  pr_A0 <- pr[1]
  
  harmonic_df <- data.frame(A0=A0,A=A,a=a,b=b,
                            pr_A0=pr_A0,pr_a=pr_a,pr_b=pr_b,
                            phase=phase,harmonic=p_val,omega=omega_val)
  
  return(harmonic_df)
}

#https://matinbrandt.wordpress.com/2013/11/15/pixel-wise-time-series-trend-anaylsis-with-ndvi-gimms-and-r/

### need to check the computation of Amplitudes!!!!

harmonic_reg_f1 <- function(y,n=24,harmonic=1){
  harmonic_results <-harmonic_regression(y,n=n,
                                         harmonic_val=NULL,
                                         mod_obj=T,figure=F)
  df_in <- subset(harmonic_results$harmonic_df,harmonic==1)
  A <- df_in$A
  
  return(A)
}

#Does work with calc
#Note that A has two outputs right now and it creates multiple outputs in raster
#harmonic_reg_raster <- function(y,var_name,n=24,harmonic_val=NULL){
harmonic_reg_raster <- function(y,var_name,n=24,harmonic_val=NULL){
  #
  ## Function to generate outputs from harmonic regression for every pixel in a raster image.
  ## Note that the output needed for the calc function used is a vector of values.
  ## The output values can be amplitudes (A0, A1, A2), phases, p significance etc.
  #
  ### INPUTS:
  #1) y: input data, in this function the name of raster layer is expected
  #2) var_name: number of elements in the time series/sequence used for the harmonic reg.
  #3) n: number of elements in the time series/sequence used for the harmonic reg.
  #4) harmonic_val: number of harmonic to consider, if NULL then use default 2
  ### OUTPUTS
  #1) value_var_name: numeric vector of values from harmonic modeling
  
  ####### Start script #######
  
  #n <- layers(y)
  
  #debug(harmonic_regression)
  harmonic_results <-harmonic_regression(y,n=n,
                                         harmonic_val=harmonic_val,
                                         mod_obj=T,figure=F)
  
  #df_in <- subset(harmonic_results$harmonic_df,harmonic==harmonic)
  df_in <- harmonic_results$harmonic_df
  
  ## Select variables to predict for every pixels
  value_var_name <- df_in[[var_name]]
  
  ## Return values for raster:
  return(value_var_name)
}

#fun=function(x) { if (is.na(x[1])){ NA } else { m = lm(x ~ time); summary(m)$coefficients[8] }}


calcHarmonicRaster <- function(r,harmonic_val=NULL,var_name="A",window_val=23,file_format=".tif",multiband=F,num_cores=1,raster_name=NULL,out_dir=NULL){
  
  ### for A1, function seems to work for A1 and A2
  #l_r_A1 <- vector("list",length=n_split)
  
  n_val <- window_val
  #harmonic_val <- NULL
  #var_name <- "A"
  #raster_name <- NULL
  #raster_name <- "NDVI_amplitude.tif"
  #file_format <- ".tif"
  #multiband=F
  
  ############ Start script ##################
  
  if(!is.null(window_val)){
      
    seq_val <- 1:nlayers(r)
    #debug(split_sequence)
    split_obj <- split_sequence(seq_val,n=n_val)
    intervals_df <- split_obj$intervals
    #lengh(split_obj$list_y[[1]])
    #split_obj$list_y[[1]]
    n_split <- nrow(intervals_df)
    
    list_rast_out <- vector("list", length=n_split)
    
    for(i in 1:n_split){
      start_val <- intervals_df$start[i]
      end_val <- intervals_df$end[i]
      n_val <- end_val - start_val + 1
      
      r_out <- try(calc(subset(r,end_val:start_val), 
                        fun=function(y){harmonic_reg_raster(y,
                                                            var_name=var_name,
                                                            n=n_val,
                                                            harmonic=harmonic_val)}))
      names(r_out) <- paste0(var_name,"_",1:nlayers(r_out))
      
      if(is.null(harmonic_val)){
        h_val <- 2
      }else{
        h_val <- harmonic_val
      }
      layer_names <- paste(var_name,1:h_val,sep="_")
      names(r_out) <- layer_names
      
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
      
      #Use compression option for tif
      writeRaster(r_out,
                  filename=file.path(out_dir,raster_name_tmp),
                  bylayer=bylayer_val,
                  suffix=suffix_str,
                  overwrite=TRUE,
                  #NAflag=NA_flag_val,
                  #datatype=data_type_str,
                  options=c("COMPRESS=LZW"))
      
      #l_r_A1[[i]] <- r_out
      list_rast_out[[i]] <- rast_list 
      
    }
    rast_list <- unlist(list_rast_out) 
  }
  
  ##### If no window defined, just use the brick/stack
  
  if(is.null(window_val)){
    r_out <- try(calc(r, 
                      fun=function(y){harmonic_reg_raster(y,
                                                          var_name=var_name,
                                                          n=nlayers(r),
                                                          harmonic=harmonic_val)}))
    names(r_out) <- paste0(var_name,"_",1:nlayers(r_out))
    suffix_str <- names(r_out)
    
    if(is.null(raster_name)){
      out_raster_name <- paste(var_name,"_",file_format,sep="") 
    }else{
      out_raster_name <- raster_name
      #out_raster_name <- sub(file_format,"",raster_name)
      #out_raster_name <- paste(out_raster_name,"_",i,file_format,se  p="") 
    }
    
    if(multiband==FALSE){
      out_raster_name <- sub(file_format,"",out_raster_name)
      raster_name_tmp <- paste(out_raster_name,file_format,sep="") #don't add output suffix because in suffix_str
      bylayer_val <- TRUE #write out separate layer files for each "band"
      rast_list <- file.path(out_dir,(paste(out_raster_name,"_",suffix_str,file_format,sep=""))) 
    }

    #Use compression option for tif
    writeRaster(r_out,
                filename=file.path(out_dir,raster_name_tmp),
                bylayer=bylayer_val,
                suffix=suffix_str,
                overwrite=TRUE,
                #NAflag=NA_flag_val,
                #datatype=data_type_str,
                options=c("COMPRESS=LZW"))
    
  }
  
  #### return object
  
  return(rast_list)
}


################################### End of script #######################################

