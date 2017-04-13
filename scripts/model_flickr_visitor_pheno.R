##Script to model the relationship between visitors and snow on Mt. Rainier.
##Ian Breckheimer
##Initially Created 27 October 2015
##Last updated January 2017

####Sets up workspace####
library(ggplot2)
library(mgcv)
library(dplyr)

####Loads and munges data####
d_f <- read.csv("./data/MORA_flickr_metadata_all_2009_2015_cleaned_final.csv")

##Calculates the number of unique users for each day every year.##
d_f$datetaken <- as.POSIXct(d_f$datetaken)
d_f$week <- as.numeric(format(d_f$datetaken,format="%U"))
d_f$DOW <- as.numeric(format(d_f$datetaken,format="%u"))
d_grp <- group_by(d_f,datetaken_DOY,week,year,nearest_center)
d_sum <- summarize(d_grp,
                   nphotos=length(id),
                   nusers=length(unique(owner)),
                   mean_dss=mean(days_since_snow),
                   mean_sdd=mean(sdd_pred),
                   mean_doy=mean(datetaken_DOY),
                   mean_travel_t=mean(acc_times))

##Calculates the number of unique visitors for every day since snow melt.
d_grp2 <- group_by(d_f,datetaken_DOY,DOW,week,year,nearest_center)
d_sum_snow <- summarize(d_grp2,
                        nphotos=length(id),
                        nusers=length(unique(owner)),
                        mean_sdd=mean(sdd_pred),
                        mean_travel_t=mean(acc_times))

##Does the same to get yearly totals.
d_grp3 <- group_by(d_f,nearest_center,year)
d_sum_year<- summarize(d_grp3,
                        nphotos=length(id),
                        nusers=length(unique(owner)),
                        mean_doy=mean(datetaken_DOY),
                        mean_travel_t=mean(acc_times))

##Adds back groups with zero records.
days<- expand.grid(DOY= 1:365, year = 2009:2015,nearest_center=levels(d_grp$nearest_center))
d_sum_complete <- left_join(days,d_sum, by=c("DOY" = "datetaken_DOY","year" = "year","nearest_center" = "nearest_center"))
d_sum_complete <- d_sum_complete %>% mutate(week = as.numeric(format(as.POSIXct(paste(year,DOY,sep="-"),format="%Y-%j"),format="%U")),
                                            DOW = as.numeric(format(as.POSIXct(paste(year,DOY,sep="-"),format="%Y-%j"),format="%w")),
                                            nphotos = ifelse(is.na(nphotos),0,nphotos),
                                            nusers = ifelse(is.na(nusers),0,nusers),
                                            study = "Flickr")
##Removes zeroes from last part of 2015
#d_sum_complete$nusers[d_sum_complete$year==2015 & d_sum_complete$week > 40] <- NA
#d_sum_complete$nphotos[d_sum_complete$year==2015 & d_sum_complete$week > 40] <- NA


##Melt dates for paradise SNOTEL 2009 - 2015
#para_dates <- c(202,205,240,209,199,204,151)
#para_elevs <- rep(1563,7)

#Adds predicted melt dates at parking lots from snow regression model.
para_preds <- c(177,191,210,191,180,181,139)
sunr_preds <- c(172,199,211,192,181,183,131)
tipsoo_preds <- c(175,173,207,187,171,171,107)
mow_preds <- c(177,187,203,185,178,176,149)
snow_preds <- c(para_preds,sunr_preds,tipsoo_preds,mow_preds)
snow_years <- rep(2009:2015,4)
snow_wks <- (snow_preds / 365) * 52
snow_locs <- factor(rep(c("Paradise","Sunrise","Tipsoo","Mowich"),each=7))
snow_elevs <- rep(c(1563,1950,1626,1506),4)

##Appends those values to a data frame.
d_snow <- data.frame(d_sum_complete[1:28,])
d_snow[,] <- rep(NA,nrow(d_snow)*ncol(d_snow))
d_snow$DOY <- snow_preds
d_snow$week <- snow_wks
d_snow$year <- snow_years
d_snow$nearest_center <- snow_locs
d_snow$study <- rep("Snow",28)

##Removes data from 2009 at Tipsoo lake, because no photos were taken there.
d_sum_complete <- filter(d_sum_complete, year != 2009 | nearest_center != "Tipsoo")
d_sum_complete <- filter(d_sum_complete,!is.na(nusers))

d_data_week <- group_by(filter(d_sum_complete,study=="Flickr"),week,year,nearest_center)
d_data_week <- summarize(d_data_week,
                         nphotos=sum(nphotos),
                         nusers=sum(nusers),
                         mean_dss=mean(mean_dss,na.rm=TRUE),
                         mean_sdd=mean(mean_sdd,na.rm=TRUE),
                         mean_travel_t=mean(mean_travel_t,na.rm=TRUE))

d_data_day <- group_by(filter(d_sum_complete,study=="Flickr"),year,week,DOY,DOW,nearest_center)
d_data_day <- summarize(d_data_day,
                         nphotos=sum(nphotos),
                         nusers=sum(nusers),
                         mean_dss=mean(mean_dss,na.rm=TRUE),
                         mean_sdd=mean(mean_sdd,na.rm=TRUE),
                         mean_travel_t=mean(mean_travel_t,na.rm=TRUE))


##Adds the predicted melt day of the nearest center predicted melt day
weeks_snowmelt <- left_join(d_data_week,d_snow[,1:4],by=c("year","nearest_center"))
colnames(weeks_snowmelt) <- c("week","year","nearest_center","nphotos","nusers","mean_dss",
                              "mean_sdd","mean_travel_t","melt_DOY","melt_week")

days_snowmelt <- left_join(d_data_day,d_snow[,1:4],by=c("year","nearest_center"))
colnames(days_snowmelt) <- c("year","week","DOY","DOW","nearest_center","nphotos","nusers",
                             "mean_dss","mean_sdd","mean_travel_t","melt_DOY","melt_week")

##Estimates day of year for each week and melt week.
weeks_snowmelt$DOY <- NA
weeks_snowmelt$DOY[weeks_snowmelt$year==2012] <- (weeks_snowmelt$week[weeks_snowmelt$year==2012] / (366/7)) * 366 - 3.5
weeks_snowmelt$DOY[weeks_snowmelt$year!=2012] <- (weeks_snowmelt$week[weeks_snowmelt$year!=2012] / (365/7)) * 365 - 3.5

##Removes ski season records because they have a lot of leverage.
weeks_snowmelt <- filter(weeks_snowmelt,DOY > 70 & DOY < 335)
days_snowmelt <- filter(days_snowmelt,DOY > 70 & DOY < 335)

####Uses a parametric approach to model weekly visitation with JAGS ####
library(rjags)
load.module("glm")

##Preps data.

##Sample of data for model development.
sample_ind <- runif(length(days_snowmelt$nusers),0,1)
days_samp <- days_snowmelt[sample_ind <= 1,]

obsy <- days_samp$nusers
x1 <- days_samp$DOY/100 - 2
obsx2 <- days_samp$melt_DOY/100 - 2
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
DOW <- as.numeric(as.factor(days_samp$DOW))
year <- as.numeric(as.factor(days_samp$year))
site <- as.numeric(days_samp$nearest_center)
group <- as.numeric(factor(paste(days_samp$year,days_samp$nearest_center)))
access <- as.numeric(days_samp$nearest_center %in% c("Paradise","Tipsoo")) + 1
n <- length(obsy)
nyears <- length(unique(year))
nsites <- length(unique(site))
ngroups <- length(unique(group))
nDOW <- length(unique(DOW))

# write model where day of week, site and year have random height effects, but common slopes and intercepts on opt and width.
cat("
    model{
    ## Priors
    height_slope ~ dnorm(0,1)T(-5,5)
    opt_int ~ dnorm(0,0.01)T(-20,20)
    opt_slope ~ dnorm(0,0.01)T(-5,5)
    width_int ~ dnorm(0,0.01)T(-20,20)
    width_slope ~ dnorm(0,0.01)T(-20,20)
    year_height_int_sd ~ dgamma(1,0.1)
    year_height_int_tau <- pow(year_height_int_sd,-2)
    site_height_int_sd ~ dgamma(1,0.1)
    site_height_int_tau <- pow(site_height_int_sd,-2)
    dow_height_int_sd ~ dgamma(1,0.1)
    dow_height_int_tau <- pow(dow_height_int_sd,-2)
    #group_opt_int_sd ~ dgamma(1,0.1)
    #group_opt_int_tau <- pow(group_opt_int_sd,-2)
    err_sd ~ dunif((sd_obs - 0.1),(sd_obs + 0.1))
    tau_obs <- 1 / (err_sd * err_sd)
    r ~ dgamma(1,0.1)
    
    ##Random effects priors.
    #for (k in 1:ngroup){
    #  opt_int_group[k] ~ dnorm(0,group_opt_int_tau)
    #}
    for (j in 1:nyr){
      height_int_yr[j] ~ dnorm(0,year_height_int_tau)
    }
    for (l in 1:nsite){
      height_int_site[l] ~ dnorm(0,site_height_int_tau)
    }
    for (m in 1:ndow){
      height_int_dow[m] ~ dnorm(0,dow_height_int_tau)
    }
    
    ## Likelihood
    for (i in 1:n){
      x2true[i] ~ dnorm(0,pop_taux)
      x2[i] ~ dnorm(x2true[i],tau_obs)
      height[i] <- height_slope * x2true[i] + height_int_yr[year[i]] + height_int_site[site[i]] + height_int_dow[dow[i]]
      opt[i] <- opt_int + opt_slope * x2true[i]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      log(mu[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      p[i] <- r/(r+mu[i])
      y[i] ~ dnegbin(p[i],r)
    }
    }
    ",
    fill=TRUE, file="./scratch/xyerror_day_combined_nb.txt")

# write model where day of week, site and year have random height effects, 
# but common fixed slopes and intercepts on opt and width.
cat("
    model{
    ## Priors
    height_slope ~ dnorm(0,1)T(-5,5)
    opt_int ~ dnorm(0,0.01)T(-20,20)
    opt_slope ~ dnorm(0,0.01)T(-5,5)
    width_int ~ dnorm(0,0.01)T(-20,20)
    width_slope ~ dnorm(0,0.01)T(-20,20)
    year_height_int_sd ~ dgamma(1,0.1)
    year_height_int_tau <- pow(year_height_int_sd,-2)
    site_height_int_sd ~ dgamma(1,0.1)
    site_height_int_tau <- pow(site_height_int_sd,-2)
    dow_height_int_sd ~ dgamma(1,0.1)
    dow_height_int_tau <- pow(dow_height_int_sd,-2)
    #group_opt_int_sd ~ dgamma(1,0.1)
    #group_opt_int_tau <- pow(group_opt_int_sd,-2)
    err_sd ~ dunif((sd_obs - 0.1),(sd_obs + 0.1))
    tau_obs <- 1 / (err_sd * err_sd)
    r ~ dgamma(1,0.1)
    
    ##Random effects priors.
    #for (k in 1:ngroup){
    #  opt_int_group[k] ~ dnorm(0,group_opt_int_tau)
    #}
    for (j in 1:nyr){
    height_int_yr[j] ~ dnorm(0,year_height_int_tau)
    }
    for (l in 1:nsite){
    height_int_site[l] ~ dnorm(0,site_height_int_tau)
    opt_int_site[l] ~ dnorm(0,0.001)
    opt_slope_site[l] ~ dnorm(0,0.001)
    width_int_site[l] ~ dnorm(0,0.001)
    }
    for (m in 1:ndow){
    height_int_dow[m] ~ dnorm(0,dow_height_int_tau)
    }
    
    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_slope * x2true[i] + height_int_yr[year[i]] + height_int_site[site[i]] + height_int_dow[dow[i]]
    opt[i] <- opt_int_site[site[i]] + opt_slope_site[site[i]] * x2true[i]
    width[i] <- exp(width_int_site[site[i]] + width_slope * x2true[i]) * -1
    log(mu[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    p[i] <- r/(r+mu[i])
    y[i] ~ dnegbin(p[i],r)
    }
    }
    ",
    fill=TRUE, file="./scratch/xyerror_day_site_nb.txt")


# write model where day of week, site and year have random height effects, 
# but fixed slopes and intercepts on opt and width that vary by access.
cat("
    model{
    ## Priors
    height_slope ~ dnorm(0,1)T(-5,5)
    opt_int ~ dnorm(0,0.01)T(-20,20)
    opt_slope ~ dnorm(0,0.01)T(-5,5)
    width_int ~ dnorm(0,0.01)T(-20,20)
    width_slope ~ dnorm(0,0.01)T(-20,20)
    year_height_int_sd ~ dgamma(1,0.1)
    year_height_int_tau <- pow(year_height_int_sd,-2)
    site_height_int_sd ~ dgamma(1,0.1)
    site_height_int_tau <- pow(site_height_int_sd,-2)
    site_opt_int_sd ~ dgamma(1,0.1)
    site_opt_int_tau <- pow(site_opt_int_sd,-2)
    dow_height_int_sd ~ dgamma(1,0.1)
    dow_height_int_tau <- pow(dow_height_int_sd,-2)
    #group_opt_int_sd ~ dgamma(1,0.1)
    #group_opt_int_tau <- pow(group_opt_int_sd,-2)
    err_sd ~ dunif((sd_obs - 0.1),(sd_obs + 0.1))
    tau_obs <- 1 / (err_sd * err_sd)
    r ~ dgamma(1,0.1)
    
    ##Random effects priors.
    for (j in 1:nyr){
    height_int_yr[j] ~ dnorm(0,year_height_int_tau)
    }
    for (l in 1:nsite){
    height_int_site[l] ~ dnorm(0,site_height_int_tau)
    opt_int_site[l] ~ dnorm(0,site_opt_int_tau)
    }
    for (m in 1:ndow){
    height_int_dow[m] ~ dnorm(0,dow_height_int_tau)
    }
    for (p in 1:naccess){
    opt_slope_access[p] ~ dnorm(0,0.001)
    width_int_access[p] ~ dnorm(0,0.001)
    width_slope_access[p] ~ dnorm(0,0.001)
    }
    
    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_slope * x2true[i] + height_int_yr[year[i]] + height_int_site[site[i]] + height_int_dow[dow[i]]
    opt[i] <- opt_int_site[site[i]] + opt_slope_access[access[i]] * x2true[i]
    width[i] <- exp(width_int_access[access[i]] + width_slope_access[access[i]] * x2true[i]) * -1
    log(mu[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    p[i] <- r/(r+mu[i])
    y[i] ~ dnegbin(p[i],r)
    }
    }
    ",
    fill=TRUE, file="./scratch/xyerror_day_access_nb.txt")


# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,year = year,site=site, access=access,dow=DOW,
               nyr=nyears, nsite=nsites, naccess=2,ndow=nDOW, n = n)

# initiate common model
mod2 <- jags.model("./scratch/xyerror_day_combined_nb.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)

# Save to disk
save(mod2,file="./scratch/visitors_jags_common_2009_2015_model.Rdata",compress=TRUE)


# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_slope","height_int_dow",
                                      "height_int_yr","height_int_site",
                                      "opt_int","opt_slope", 
                                      "width_int","width_slope",
                                      "year_height_int_sd","site_height_int_sd",
                                      "dow_height_int_sd","r"))
gelman.diag(out2)
save(out2,file="./scratch/visitors_jags_common_2009_2015.Rdata",compress=TRUE)

# Collects DIC samples.
mod2_dic <- dic.samples(mod2,thin=10,n.iter=10000,type="pD")
save(mod2_dic,file="./scratch/visitors_jags_common_2009_2015_dic.Rdata",compress=TRUE)

# initiate site-specific model
mod3 <- jags.model("./scratch/xyerror_day_site_nb.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod3, n.iter=10000)

# Save to disk
save(mod3,file="./scratch/visitors_jags_site_commonwidthslope_2009_2015_model.Rdata",compress=TRUE)

# simulate posterior
out3 <- coda.samples(mod3, n.iter=10000, thin=10,
                     variable.names=c("height_slope","height_int_dow",
                                      "height_int_yr","height_int_site",
                                      "opt_int_site","opt_slope_site", 
                                      "width_int_site","width_slope",
                                      "year_height_int_sd","site_height_int_sd",
                                      "dow_height_int_sd","r"))
gelman.diag(out3)
save(out3,file="./scratch/visitors_jags_site_commonwidthslope_2009_2015.Rdata",compress=TRUE)

# Collects DIC samples.
mod3_dic <- dic.samples(mod3,thin=10,n.iter=10000,type="pD")
save(mod3_dic,file="./scratch/visitors_jags_site_commonwidthslope_2009_2015_dic.Rdata",compress=TRUE)

# initiate common model
mod4 <- jags.model("./scratch/xyerror_day_access_nb.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod4, n.iter=10000)

# Save to disk
save(mod4,file="./scratch/visitors_jags_access_2009_2015_model.Rdata",compress=TRUE)

# simulate posterior
out4 <- coda.samples(mod4, n.iter=10000, thin=10,
                     variable.names=c("height_slope","height_int_dow",
                                      "height_int_yr","height_int_site",
                                      "opt_int_site","opt_slope_access", 
                                      "width_int_access","width_slope_access",
                                      "year_height_int_sd","site_height_int_sd",
                                      "dow_height_int_sd","r"))
gelman.diag(out4)
save(out4,file="./scratch/visitors_jags_access_2009_2015.Rdata",compress=TRUE)

# Collects DIC samples.
mod4_dic <- dic.samples(mod4,thin=10,n.iter=10000,type="pD")
save(mod4_dic,file="./scratch/visitors_jags_access_2009_2015_dic.Rdata",compress=TRUE)

###Compares DIC
mod2_dic
mod3_dic
mod4_dic

####Visualizes the JAGS fit ####
jags_pred_data <- expand.grid(DOY=seq(1,365,by=1),
                              melt_DOY=seq(100,220,by=1),
                              DOW=seq(1,7))

jags_fit_fun <- function(height_int,height_int_DOW,
                         height_slope,
                         opt_int,opt_slope,
                         width_int,width_slope,
                         x1vec,x2vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  
  #Calculates the value of the response.
  height <- height_int + height_int_DOW + height_slope * x2vec
  opt <- opt_int + opt_slope * x2vec
  width<- exp(width_int + width_slope * x2vec) * -1
  mu <- exp(width * (x1vec - opt)^2 + height)
  return(mu)
}

jags_meds <- summary(out2)$quantiles[,3]
jags_pred_data$predmu <- NA

for (i in 1:7){
  jags_pred_data$predmu[jags_pred_data$DOW==i] <- jags_fit_fun(height_int=jags_meds["height_int_site[2]"],
                                                               height_int_DOW=jags_meds[paste("height_int_dow[",i,"]",sep="")],
                                                               height_slope=jags_meds["height_slope"],
                                                               width_int=jags_meds["width_int"],
                                                               width_slope=jags_meds["width_slope"],
                                                               opt_int=jags_meds["opt_int"],
                                                               opt_slope=jags_meds["opt_slope"],
                                                               x1vec=jags_pred_data$DOY[jags_pred_data$DOW==i]/100 - 2,
                                                               x2vec=jags_pred_data$melt_DOY[jags_pred_data$DOW==i]/100 - 2)
}

ggplot(filter(jags_pred_data,DOW==7))+
  geom_raster(aes(y=DOY,x=melt_DOY,fill=predmu),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  scale_color_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  geom_point(aes(y=DOY,x=melt_DOY),size=1,color="white",
             position=position_jitter(width=5),
             data=filter(days_snowmelt,nusers==0),alpha=0.1)+
  geom_point(aes(y=DOY,x=melt_DOY,size=nusers),color="black",
             position=position_jitter(width=5),
             data=filter(days_snowmelt,nusers>0),alpha=0.15)+
  scale_x_continuous(limits=c(110,220),expand=c(0,0))+
  scale_y_continuous(limits=c(110,335),expand=c(0,0))+
  xlab("Day of Snow Melt")+
  ylab("Day of Year")+
  guides(fill=guide_colorbar("Predicted Mean"),
         size=guide_legend("Number of Users"))+
  theme_bw()+
  theme(text=element_text(family="sans"))

##Visualizes the fit against data for individual sites and years.
week_labels <- c("Sunday","Monday","Tuesday","Wednesday", "Thursday","Friday","Saturday")
ggplot(filter(jags_pred_data,melt_DOY==180))+
  geom_line(aes(x=DOY,y=predmu,color=factor(DOW)))+
  xlim(c(70,335))+
  ylab("Number of Unique Users")+
  xlab("Day of Year")+
  scale_color_discrete("Day of Week",labels=week_labels)+
  geom_point(aes(x=DOY,y=nusers,color=factor(DOW+1)),position=position_jitter(height=0.2),
            data=subset(days_snowmelt,nearest_center=="Paradise" & year==2013))+
  theme_bw()
