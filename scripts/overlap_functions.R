##R Functions to support social-ecological phenology analysis.
##Ian Breckheimer
##Last updated January 2017

####Workhorse Functions ####
inv.logit <- function(x){exp(x)/(1+exp(x))}

##Functional forms of the fit relationships for flowers and visitors
f1 <- function(x1,x2,height_int,height_slope,
               opt_int,opt_slope,
               width_int,width_slope) {
  
  exp(exp(width_int + width_slope * x2) * -1 * 
        (x1 -(opt_int + opt_slope * x2))^2 + (height_int + height_slope * x2))
  
}

f2 <- function(x1,x2,height_int2,height_slope2,
               opt_int2,opt_slope2,
               width_int2,width_slope2) {
  
  inv.logit(exp(width_int2 + width_slope2 * x2) * -1 * 
              (x1 -(opt_int2 + opt_slope2 * x2))^2 + (height_int2 + height_slope2 * x2))
  
}

##Density functions
f1_dens <- function(x1,x2,height_int,height_slope,
                    opt_int,opt_slope,
                    width_int,width_slope) { 
  
  y <- f1(x1,x2,height_int,height_slope,
          opt_int,opt_slope,
          width_int,width_slope)
  yi <- integrate(f1, -Inf, +Inf, x2,height_int=height_int,height_slope=height_slope,
                  opt_int=opt_int, opt_slope=opt_slope,
                  width_int=width_int,width_slope=width_slope)
  return(y/yi[[1]])
}


f2_dens <- function(x1,x2,height_int2,height_slope2,
                    opt_int2,opt_slope2,
                    width_int2,width_slope2) { 
  
  y <- f2(x1,x2,height_int2,height_slope2,
          opt_int2,opt_slope2,
          width_int2,width_slope2)
  yi <- integrate(f2, -Inf, +Inf, x2,height_int2=height_int2,height_slope2=height_slope2,
                  opt_int2=opt_int2, opt_slope2=opt_slope2,
                  width_int2=width_int2,width_slope2=width_slope2)
  return(y/yi[[1]])
}

##Function to find the joint minimum of these two functions.
min_f1f2_dens <- function(x1,x2, height_int,height_slope,
                          opt_int,opt_slope,
                          width_int, width_slope,
                          height_int2,height_slope2,
                          opt_int2,opt_slope2,
                          width_int2,width_slope2) {
  
  f1 <- f1_dens(x1,x2,height_int,height_slope,
                opt_int,opt_slope,
                width_int,width_slope) 
  f2 <- f2_dens(x1,x2,height_int2,height_slope2,
                opt_int2,opt_slope2,
                width_int2,width_slope2)
  return(pmin(f1, f2))
}

##Function to calculate overlap coefficients.
measure_overlap <- function(params,x2val){
  integrate(min_f1f2_dens, -Inf, Inf, x2=x2val,height_int=params[1],height_slope=params[2],
            opt_int=params[3], opt_slope=params[4],
            width_int=params[5],width_slope=params[6],
            height_int2=params[7],height_slope2=params[8],
            opt_int2=params[9], opt_slope2=params[10],
            width_int2=params[11],width_slope2=params[12])$value
}

####Function to extract the appropriate variables from the MCMC output
prep_mcmc_vars <- function(model_output1,model_output2,group_index=1,
                           mod1_param_names = c("height_int_site","height_slope",
                                                "opt_int","opt_slope",
                                                "width_int","width_slope"),
                           mod2_param_names = c("height_int","height_slope",
                                                "opt_int","opt_slope",
                                                "width_int","width_slope")){
  require(coda)
  require(dplyr)
  require(ggmcmc)
  
  ##Appends the appropriate group index.
  index <- paste("[",group_index,"]",sep="")
  output1_index <- paste(mod1_param_names,index,sep="")
  output1_param_strings <- c(mod1_param_names,output1_index)
  output2_index <- paste(mod2_param_names,index,sep="")
  output2_param_strings <- c(mod2_param_names,output2_index)
  
  ##Converts mcmc output to data frames
  output1_samples <- filter(ggs(model_output1),Parameter %in% output1_param_strings)
  stopifnot(dim(output1_samples)[1] == (dim(model_output1[[1]])[1] * 
                                          length(mod1_param_names) * length(model_output1)))
  output1_samples$Parameter <- factor(output1_samples$Parameter)
  output1_params <- tidyr::spread(output1_samples,Parameter,value)
  output1_means <- colMeans(output1_params)[3:8]
  print("Model 1 mean parameter values:")
  print(output1_means)
  
  output2_samples <- filter(ggs(model_output2),Parameter %in% output2_param_strings)
  stopifnot(dim(output2_samples)[1] == (dim(model_output2[[1]])[1] *
                                          length(mod2_param_names) * length(model_output1)))
  output2_samples$Parameter <- factor(output2_samples$Parameter)
  output2_params <- tidyr::spread(output2_samples,Parameter,value)
  output2_means <- colMeans(output2_params)[3:8]
  print("Model 2 mean parameter values:")
  print(output2_means)
  params <- cbind(output1_params[3:8],output2_params[,3:8])
  return(params)
}

###Function to measure phenological mismatch and associated uncertainty
measure_mismatch <- function(params,x2sdd=seq(110,250,by=4)){
  stopifnot(dim(params)[2] == 12)
  if(dim(params)[1] < 1000){warning("Less than 1000 samples,credible intervals will be crude.")}
  x2seq <- x2sdd / 100 - 2
  overmat <- matrix(NA,nrow=nrow(params),ncol=length(x2seq))
  for (i in 1:length(x2seq)){
    flush.console()
    print(paste("Computing match for day ",x2sdd[i],"(",i," of ",length(x2sdd),")"))
    overmat[,i] <- apply(params,FUN=measure_overlap,MARGIN=1,x2val=x2seq[i])
  }
  
  quantile_fun50 <- function(x){quantile(x,probs=0.5)}
  quantile_fun25 <- function(x){quantile(x,probs=0.25)}
  quantile_fun75 <- function(x){quantile(x,probs=0.75)}
  quantile_fun10 <- function(x){quantile(x,probs=0.10)}
  quantile_fun90 <- function(x){quantile(x,probs=0.90)}
  quantile_funlwr <- function(x){quantile(x,probs=0.025)}
  quantile_funupr <- function(x){quantile(x,probs=0.975)}
  
  overlap <- data.frame(SDD=x2sdd)
  overlap$overlap_q50 <- apply(overmat,FUN=quantile_fun50,MARGIN=2)
  overlap$overlap_q25 <- apply(overmat,FUN=quantile_fun25,MARGIN=2)
  overlap$overlap_q75 <- apply(overmat,FUN=quantile_fun75,MARGIN=2)
  overlap$overlap_q10 <- apply(overmat,FUN=quantile_fun10,MARGIN=2)
  overlap$overlap_q90 <- apply(overmat,FUN=quantile_fun90,MARGIN=2)
  overlap$overlap_lwr <- apply(overmat,FUN=quantile_funlwr,MARGIN=2)
  overlap$overlap_upr <- apply(overmat,FUN=quantile_funupr,MARGIN=2)
  return(overlap)
}

measure_mismatch_change <- function(params,x2sdd=c(130,211)){
  stopifnot(dim(params)[2] == 12)
  if(dim(params)[1] < 1000){warning("Less than 1000 samples,credible intervals will be crude.")}
  x2seq <- x2sdd / 100 - 2
  overmat <- matrix(NA,nrow=nrow(params),ncol=length(x2seq))
  for (i in 1:length(x2seq)){
    flush.console()
    print(paste("Computing match for day ",x2sdd[i],"(",i," of ",length(x2sdd),")"))
    overmat[,i] <- apply(params,FUN=measure_overlap,MARGIN=1,x2val=x2seq[i])
  }
  change <- (overmat[,2] - overmat[,1]) / overmat[,2]
  quantile_fun <- function(x){quantile(x,probs=c(0.5,0.25,0.75,0.1,0.9,0.025,0.975))}
  
  overlap_change_quants <-quantile_fun(change)

  return(overlap_change_quants)
}

##Function that does the heavy-lifting to estimate mismatch-snowmelt relationships between sites.
jags_mismatch_site <- function(visit_model,flower_model,group_index,measure_sdd,
                               site_names,access_indices,side_indices){
  access_index <- access_indices[group_index]
  side_index <- side_indices[group_index]
  
  ##Parameter names
  access_width_int_name <- paste("width_int_access",access_index,sep="")
  access_width_slope_name <- paste("width_slope_access",access_index,sep="")
  access_opt_slope_name <- paste("opt_slope_access",access_index,sep="")
  side_opt_slope_name <- paste("opt_slope_side",side_index,sep="")
  
  params <- prep_mcmc_vars(visit_model,flower_model,group_index=group_index,
                           mod1_param_names = c("height_int_site","height_slope",
                                                "opt_int_site",access_opt_slope_name,
                                                access_width_int_name,access_width_slope_name),
                           mod2_param_names = c("height_int","height_slope",
                                                "opt_int_site",side_opt_slope_name,
                                                "width_int_site","width_slope"))
  mismatch <- measure_mismatch(params=params,x2sdd=measure_sdd)
  mismatch$Site <- site_names[i]
  return(mismatch)
}

##Function that does the heavy-lifting to estimate mismatch-snowmelt relationships between sites.
jags_mismatch_site_change <- function(visit_model,flower_model,group_index,measure_sdd,
                               site_names,access_indices,side_indices){
  access_index <- access_indices[group_index]
  side_index <- side_indices[group_index]
  
  ##Parameter names
  access_width_int_name <- paste("width_int_access",access_index,sep="")
  access_width_slope_name <- paste("width_slope_access",access_index,sep="")
  access_opt_slope_name <- paste("opt_slope_access",access_index,sep="")
  side_opt_slope_name <- paste("opt_slope_side",side_index,sep="")
  
  params <- prep_mcmc_vars(visit_model,flower_model,group_index=group_index,
                           mod1_param_names = c("height_int_site","height_slope",
                                                "opt_int_site",access_opt_slope_name,
                                                access_width_int_name,access_width_slope_name),
                           mod2_param_names = c("height_int","height_slope",
                                                "opt_int_site",side_opt_slope_name,
                                                "width_int_site","width_slope"))
  mismatch <- measure_mismatch_change(params=params,x2sdd=measure_sdd)
  mismatch$Site <- site_names[i]
  return(mismatch)
}


####Function to estimate onset and end of season.
# Function to find the values that define 68.4% of the area under the curve.####
est_width <- function(time,pred,threshold=0.1586553){
  
  # Converts to numeric vectors
  time <- as.numeric(time)
  pred <- as.numeric(pred)
  
  # Gets the bin width from the first pred interval.
  bin_width <- time[2] - time[1]
  
  # Make sure that all predictions are positive and none are missing.
  stopifnot(anyNA(pred)==FALSE,
            anyNA(time)==FALSE,
            any(pred<0)==FALSE)
  
  # Total area under the curve.
  total_area <- sum(pred*bin_width,na.rm=TRUE)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(pred)*bin_width/total_area
  cumprop_down <- rev(cumsum(rev(pred))*bin_width)/total_area
  
  # Finds the indices of the first and last values greater than 0.158
  lwr_index <- min(which(cumprop_up >= threshold))
  upr_index <- max(which(cumprop_down >= threshold))
  
  # Finds the corresponding values of dss.
  lwr_bound <- time[lwr_index]
  upr_bound <- time[upr_index]
  bounds <- c(lwr_bound,upr_bound)
  names(bounds) <- c("lwr_bound","upr_bound")
  
  # Output
  return(bounds)
}

####Function to estimate upper and lower bounds from a large matrix of predictions.

