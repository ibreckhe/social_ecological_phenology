##Script to model wildflower phenology in Flickr data##
##Ian Breckheimer
##Initially Created 10 November 2015
##Last updated January 2017

####Sets up workspace####
library(Hmisc)
library(dplyr)
library(ggplot2)
library(mgcv)
library(rjags)
setwd("~/code/MORA_pheno")

####Brings in data####
f_all <- read.csv("./data/MORA_flickr_classified_2009_2015_final.csv")

###Removes photos with more than 20 photos per unique coordinate####
f_all$locfact <- as.factor(paste(f_all$long,f_all$lat,sep="by"))
f_all_grp <- group_by(f_all,locfact)
nphotos <- summarize(f_all_grp,nphotos=n())
small_locs <- nphotos$locfact[which(nphotos$nphotos < 20)]

###Creates a new dataset indicating whether there are any flowers at all.
f_all <- filter(f_all,Phenophase!="np"&
                      Species != "other" &
                      days_since_snow > -100 &
                       datetaken_DOY > 90 &
                       datetaken_DOY < 330 &
#                      year != 2012 &
#                      year != 2014 &
                      locfact %in% small_locs)

###Abundance of the 10 focal species in the dataset.
f_flowering <- filter(f_all,Phenophase=="flowering")

###Relative abundance of Erythronium and Anemone in each year.

f_grp <- group_by(f_all,id,owner,year,datetaken_DOY,
                  nearest_center,days_since_snow, relev_30m,
                  canopy_pct,srad_noc,sdd_pred,X_UTM,Y_UTM)
f_flwr <- summarise(f_grp,FlwrYN = as.numeric(any(Phenophase=="flowering")),
                          EarlyYN = as.numeric(any(Species %in% c("erythronium montanum","anemone occidentalis") & 
                                          Phenophase == "flowering")),
                          LateYN = as.numeric(any(Species %in% c("ligusticum grayi","polygonum bistortoides") & 
                                    Phenophase == "flowering")),
                          NoEarly = as.numeric(any(Species %nin% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering")),
                          NoLate = as.numeric(any(Species %nin% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering")),
                          nSpp = sum(as.numeric(Phenophase=="flowering")))
f_flwr$Flwr <- as.numeric(f_flwr$FlwrYN)
f_flwr$Early <- as.numeric(f_flwr$EarlyYN)
f_flwr$Late <- as.numeric(f_flwr$LateYN)
f_flwr$notEarly <- as.numeric(f_flwr$NoEarly)
f_flwr$notLate <- as.numeric(f_flwr$NoLate)
f_flwr$yearfact <- as.factor(f_flwr$year)

##Computes the proportion of photos with early and late flowers in every year.
f_yr <- group_by(f_flwr,year)
f_prop <- summarise(f_yr,pFlwr=mean(Flwr),
                         nFlwr=sum(Flwr),
                         nEarly=sum(Early),
                         nLate=sum(Late),
                         nNotEarly=sum(notEarly),
                         nNotLate=sum(notLate),
                         pEarly=nEarly/nFlwr,
                         pLate=nLate/nFlwr)
f_prop_2009_2014 <- f_prop[f_prop$year %in% c(2009:2014),c("nEarly","nLate","nNotEarly","nNotLate")]
sums_2009_2014 <- colSums(f_prop_2009_2014)

f_prop_2015 <- f_prop[f_prop$year == 2015,c("nEarly","nLate","nNotEarly","nNotLate")]
sums_2015 <- colSums(f_prop_2015)

##Chi-squared test of the proportion of photos with early and late flowering species.
early_ctable <- matrix(NA,ncol=2,nrow=2)
early_ctable[1,1] <- sums_2009_2014["nEarly"]
early_ctable[1,2] <- sums_2015["nEarly"]
early_ctable[2,1] <- sums_2009_2014["nNotEarly"]
early_ctable[2,2] <- sums_2015["nNotEarly"]
colnames(early_ctable) <- c("Typical Melt","Early Melt")
rownames(early_ctable) <- c("Early Flowering","Not Early Flowering")

early_ctable
chisq.test(early_ctable,simulate.p.value=TRUE,B=10000)

late_ctable <- matrix(NA,ncol=2,nrow=2)
late_ctable[1,1] <- sums_2009_2014["nLate"]
late_ctable[1,2] <- sums_2015["nLate"]
late_ctable[2,1] <- sums_2009_2014["nNotLate"]
late_ctable[2,2] <- sums_2015["nNotLate"]
colnames(late_ctable) <- c("Typical Melt","Early Melt")
rownames(late_ctable) <- c("Late Flowering","Not Late Flowering")

late_ctable
chisq.test(late_ctable,simulate.p.value=TRUE,B=10000)

####Plots data.####
ggplot(f_flwr)+
  geom_point(aes(x=sdd_pred,y=datetaken_DOY,color=as.factor(FlwrYN)),size=0.1)+
  geom_abline(slope=1,intercept=0,linetype=1)+
  scale_color_grey(start=0.8,end=0.2,guide=FALSE)+
  geom_smooth(aes(x=sdd_pred,y=datetaken_DOY,color=as.factor(FlwrYN)),
              method="lm",formula=y~x,se=TRUE)+
  xlab("Predicted Snow Disappearance Day")+
  ylab("Day of Year")+
#  facet_grid(facets=year~.)+
  theme_bw()+
  theme(text=element_text(family="sans"))

ggplot(f_flwr,groups=yearfact)+
  geom_point(aes(x=days_since_snow,y=Flwr,color=nearest_center),
             position=position_jitter(height=0.1),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr,color=nearest_center),method="glm",formula=y~poly(x,2),se=TRUE,
              method.args=list(family="binomial"))+
#  facet_grid(facets=year~.)+
  theme_bw()+
  theme(text=element_text(family="sans"))

ggplot(f_flwr,groups=yearfact)+
  geom_point(aes(x=days_since_snow,y=nSpp),
             position=position_jitter(height=0.1),size=0.1)+
  geom_smooth(aes(x=days_since_snow,y=Flwr),method="gam",formula=y~s(x),se=TRUE,
              method.args=list(family="nb",optimizer="perf"))+
#  facet_grid(facets=nearest_center~.)+
  ylim(c(0,7))+
  theme_bw()+
  theme(text=element_text(family="sans"))

####Fits a parametric model using JAGS####

## Subsamples data for model development
rand <- runif(dim(f_flwr)[1],0,1)
f_samp <- f_flwr[rand <= 1,]

## Prep data.
obsy <- f_samp$Flwr
x1 <- f_samp$datetaken_DOY/100 - 2
obsx2 <- f_samp$sdd_pred/100 - 2
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
groups <- as.numeric(factor(f_samp$owner))
ngroups <- length(unique(groups))
year <- as.numeric(factor(f_samp$year))
nyears <- length(unique(year))
sites <- as.numeric(factor(f_samp$nearest_center))
nsites <- length(unique(sites))
sides <- as.numeric(factor(f_samp$nearest_center %in% c("Sunrise","Tipsoo")))
nsides <- length(unique(sides))

## specify model with common parameters across all sites except height.
cat("
    model {
    ## Priors
    height_int ~ dnorm(0,0.01)
    height_slope ~ dnorm(0,0.01)
    opt_int ~ dnorm(0,0.01)
    opt_slope ~ dnorm(1,0.1)
    width_int ~ dnorm(0,0.1)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    height_int_yr_sd ~ dunif(0,20)
    height_int_yr_tau <- pow(height_int_yr_sd,-2)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)
   
    ## Random effects priors
    for (j in 1:ngroups){
      height_grp[j] ~ dnorm(0,group_tau)
    }
    for (k in 1:nyears){
      height_int_yr[k] ~ dnorm(0,height_int_yr_tau)
    }
    
    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_int + height_slope * x2true[i] + 
                 height_grp[group[i]] + height_int_yr[year[i]]
    opt[i] <- opt_int + opt_slope * x2true[i]
    width[i] <- exp(width_int + width_slope * x2true[i]) * -1
    logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_common.txt")


## specify model with a separate opt intercept for each site.
cat("
    model {
    ## Priors
    height_int ~ dnorm(0,0.01)
    height_slope ~ dnorm(0,0.01)
    opt_int ~ dnorm(0,0.01)
    opt_slope ~ dnorm(1,0.1)
    width_int ~ dnorm(0,0.1)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    height_int_yr_sd ~ dunif(0,20)
    height_int_yr_tau <- pow(height_int_yr_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)

    ## Random effects priors
    for (j in 1:ngroups){
      height_grp[j] ~ dnorm(0,group_tau)
    }
    for (k in 1:nsites){
      opt_int_site[k] ~ dnorm(0,0.001)
    }
    for (l in 1:nyears){
      height_int_yr[l] ~ dnorm(0,height_int_yr_tau)
    }

    ## Likelihood
    for (i in 1:n){
      x2true[i] ~ dnorm(0,pop_taux)
      x2[i] ~ dnorm(x2true[i],tau_obs)
      height[i] <- height_int + height_slope * x2true[i] +
                   height_grp[group[i]] + height_int_yr[year[i]]
      opt[i] <- opt_slope * x2true[i] + opt_int_site[site[i]]
      width[i] <- exp(width_int + width_slope * x2true[i]) * -1
      logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
      y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_siteopt.txt")

## specify model with a separate opt intercept and width intercept for each site.
cat("
    model {
    ## Priors
    height_int ~ dnorm(0,0.01)
    height_slope ~ dnorm(0,0.01)
    opt_int ~ dnorm(0,0.01)
    opt_slope ~ dnorm(1,0.1)
    width_int ~ dnorm(0,0.1)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    height_int_yr_sd ~ dunif(0,20)
    height_int_yr_tau <- pow(height_int_yr_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)
    
    ## Random effects priors
    for (j in 1:ngroups){
      height_grp[j] ~ dnorm(0,group_tau)
    }
    for (k in 1:nsites){
      opt_int_site[k] ~ dnorm(0,0.01)
      width_int_site[k] ~ dnorm(3,0.1)
    }
    for (l in 1:nyears){
      height_int_yr[l] ~ dnorm(0,height_int_yr_tau)
    }

    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_int + height_slope * x2true[i] +
                  height_grp[group[i]] + height_int_yr[year[i]]
    opt[i] <- opt_slope * x2true[i] + opt_int_site[site[i]]
    width[i] <- exp(width_int_site[site[i]] + width_slope * x2true[i]) * -1
    logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_siteopt_sitewidth.txt")

## specify model with a separate opt intercept and width intercept for each site.
cat("
    model {
    ## Priors
    height_int ~ dnorm(0,0.01)
    height_slope ~ dnorm(0,0.01)
    opt_int ~ dnorm(0,0.01)
    opt_slope ~ dnorm(1,0.1)
    width_int ~ dnorm(0,0.1)
    width_slope ~ dnorm(0,0.1)T(-10,10)
    group_sd ~ dunif(0,20)
    group_tau <- pow(group_sd,-2)
    height_int_yr_sd ~ dunif(0,20)
    height_int_yr_tau <- pow(height_int_yr_sd,-2)
    err_sd ~ dunif((sd_obs - 0.01),(sd_obs + 0.01))
    tau_obs <- 1 / (err_sd * err_sd)
    
    ## Random effects priors
    for (j in 1:ngroups){
    height_grp[j] ~ dnorm(0,group_tau)
    }
    for (k in 1:nsites){
    opt_int_site[k] ~ dnorm(0,0.01)
    width_int_site[k] ~ dnorm(0,0.1)
    }
    for (l in 1:nyears){
    height_int_yr[l] ~ dnorm(0,height_int_yr_tau)
    }
    for (m in 1:nsides){
    opt_slope_side[m] ~ dnorm(1,0.1)
    }
    
    ## Likelihood
    for (i in 1:n){
    x2true[i] ~ dnorm(0,pop_taux)
    x2[i] ~ dnorm(x2true[i],tau_obs)
    height[i] <- height_int + height_slope * x2true[i] +
    height_grp[group[i]] + height_int_yr[year[i]]
    opt[i] <- opt_slope_side[side[i]] * x2true[i] + opt_int_site[site[i]]
    width[i] <- exp(width_int_site[site[i]] + width_slope * x2true[i]) * -1
    logit(alpha[i]) <- width[i] * (x1[i] - opt[i])^2 + height[i] 
    y[i] ~ dbern(alpha[i])
    }
    }
    ", fill=T, file="./scratch/xerror_siteopt_sitewidth_sideslope.txt")

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = groups, site = sites, year=year,side=sides,
               ngroups = ngroups,nsites=nsites,nyears=nyears,nsides=nsides,n = length(x1))

# initiate model
mod2 <- jags.model("./scratch/xerror_common.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod2, n.iter=10000)
save(mod2,file="./scratch/jags_flower_output_common_2009_2015_model.Rdata",compress=TRUE)
   
# simulate posterior
out2 <- coda.samples(mod2, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int",
                                      "width_int","width_slope",
                                      "group_sd"))
gelman.diag(out2)
save(out2,file="./scratch/jags_flower_output_common_2009_2015.Rdata",compress=TRUE)
load("./scratch/jags_flower_output_common_2011_2015.Rdata")

# Compute DIC
out2_dic <- dic.samples(mod2,n.iter=10000,thin=10,type="pD")
save(out2_dic,file="./scratch/jags_flower_output_common_2009_2015_dic.Rdata",compress=TRUE)

# initiate model
mod3 <- jags.model("./scratch/xerror_siteopt.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod3, n.iter=10000)
save(mod3,file="./scratch/jags_flower_output_siteopt_2009_2015_model.Rdata",compress=TRUE)
#load("./scratch/jags_flower_output_siteopt_2009_2015_model.Rdata")

# simulate posterior
out3 <- coda.samples(mod3, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int_site",
                                      "width_int","width_slope",
                                      "site_opt_int_sd","group_sd"))
gelman.diag(out3)
save(out3,file="./scratch/jags_flower_output_siteopt_2009_2015.Rdata",compress=TRUE)

# compute DIC
out3_dic <- dic.samples(mod3,n.iter=10000,thin=10,type="pD")
save(out3_dic,file="./scratch/jags_flower_output_siteopt_2009_2015_dic.Rdata",compress=TRUE)

# initiate model
mod4 <- jags.model("./scratch/xerror_siteopt_sitewidth.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod4, n.iter=20000)
save(mod4,file="./scratch/jags_flower_output_siteoptwidth_2009_2015_model.Rdata",compress=TRUE)

# simulate posterior
out4 <- coda.samples(mod4, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int_site",
                                      "width_int_site","width_slope",
                                      "site_opt_int_sd","group_sd"))
gelman.diag(out4)
save(out4,file="./scratch/jags_flower_output_siteoptwidth_2009_2015.Rdata",compress=TRUE)

# Compute DIC
out4_dic <- dic.samples(mod4,n.iter=10000,thin=10,type="pD")
save(out4_dic,file="./scratch/jags_flower_output_siteoptwidth_2009_2015_dic.Rdata",compress=TRUE)

# initiate model
mod5 <- jags.model("./scratch/xerror_siteopt_sitewidth_sideslope.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod5, n.iter=20000)
save(mod5,file="./scratch/jags_flower_output_siteoptwidth_sideslope_2009_2015_model.Rdata",compress=TRUE)

# simulate posterior
out5 <- coda.samples(mod5, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope_side", "opt_int_site",
                                      "width_int_site","width_slope",
                                      "site_opt_int_sd","group_sd"))
gelman.diag(out5)
save(out5,file="./scratch/jags_flower_output_siteoptwidth_sideslope_2009_2015.Rdata",compress=TRUE)

# Compute DIC
out5_dic <- dic.samples(mod5,n.iter=10000,thin=10,type="pD")
save(out5_dic,file="./scratch/jags_flower_output_siteoptwidth_sideslope_2009_2015_dic.Rdata",compress=TRUE)


##Compares DIC
out2_dic
out3_dic
out4_dic
out5_dic

####Checks if snow error affects the results####
sdobs <- 0.08

# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = groups, site = sites, year=year,side=sides,
               ngroups = ngroups,nsites=nsites,nyears=nyears,nsides=nsides,n = length(x1))

# initiate model
mod6 <- jags.model("./scratch/xerror_common.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod6, n.iter=10000)
save(mod6,file="./scratch/jags_flower_output_common_2009_2015_model_snow2x.Rdata",compress=TRUE)

# simulate posterior
out6 <- coda.samples(mod6, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int",
                                      "width_int","width_slope",
                                      "group_sd","r"))
gelman.diag(out6)
save(out6,file="./scratch/jags_flower_output_common_2009_2015_snow2x.Rdata",compress=TRUE)

####Checks if species selection affects the results####
f_nolup <- filter(f_all,Species!="Lupinus arcticus")
f_grp_nolup <- group_by(f_nolup,id,owner,year,datetaken_DOY,
                  nearest_center,days_since_snow, relev_30m,
                  canopy_pct,srad_noc,sdd_pred,X_UTM,Y_UTM)
f_flwr_nolup <- summarise(f_grp_nolup,FlwrYN = any(Phenophase=="flowering"),
                    EarlyYN = any(Species %in% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering"),
                    LateYN = any(Species %in% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering"),
                    NoEarly = any(Species %nin% c("erythronium montanum","anemone occidentalis") & 
                                    Phenophase == "flowering"),
                    NoLate = any(Species %nin% c("ligusticum grayi","polygonum bistortoides") & 
                                   Phenophase == "flowering"),
                    nSpp = sum(as.numeric(Phenophase=="flowering")))
f_flwr_nolup$Flwr <- as.numeric(f_flwr_nolup$FlwrYN)
f_flwr_nolup$Early <- as.numeric(f_flwr_nolup$EarlyYN)
f_flwr_nolup$Late <- as.numeric(f_flwr_nolup$LateYN)
f_flwr_nolup$notEarly <- as.numeric(f_flwr_nolup$NoEarly)
f_flwr_nolup$notLate <- as.numeric(f_flwr_nolup$NoLate)
f_flwr_nolup$yearfact <- as.factor(f_flwr_nolup$year)

## Subsamples data for model development
rand <- runif(dim(f_flwr_nolup)[1],0,1)
f_samp <- f_flwr_nolup[rand <= 1,]

## Prep data.
obsy <- f_samp$Flwr
x1 <- f_samp$datetaken_DOY/100 - 2
obsx2 <- f_samp$sdd_pred/100 - 2
sdobs <- 0.04
pop_taux2 <- 1/(sd(obsx2)^2)
groups <- as.numeric(factor(f_samp$owner))
ngroups <- length(unique(groups))
year <- as.numeric(factor(f_samp$year))
sites <- as.numeric(factor(f_samp$nearest_center))
nsites <- length(unique(sites))

ggdat <- data.frame(FlwrYN=obsy,DOY=x1,SDD=obsx2,year=factor(f_samp$year),sites=factor(f_samp$nearest_center))
ggplot(ggdat)+
  geom_point(aes(x=SDD,y=DOY,color=as.factor(FlwrYN)),size=0.1)+
  geom_abline(slope=1,intercept=0,linetype=1)+
  scale_color_grey(start=0.8,end=0.2)+
  geom_smooth(aes(x=SDD,y=DOY,color=as.factor(FlwrYN)),
              method="lm",formula=y~x,se=TRUE)+
  facet_grid(facets=.~sites)+
  theme_bw()


# bundle data
jags_d <- list(x1 = x1, x2 = obsx2, y = obsy, sd_obs = sdobs, 
               pop_taux = pop_taux2,group = groups, site = sites, year=year,side=sides,
               ngroups = ngroups,nsites=nsites,nyears=nyears,nsides=nsides,n = length(x1))

# initiate model
mod7 <- jags.model("./scratch/xerror_common.txt", data=jags_d,
                   n.chains=3, n.adapt=1000)
update(mod7, n.iter=10000)
save(mod7,file="./scratch/jags_flower_output_common_2009_2015_model_nolup.Rdata",compress=TRUE)

# simulate posterior
out7 <- coda.samples(mod7, n.iter=10000, thin=10,
                     variable.names=c("height_int","height_slope", 
                                      "opt_slope", "opt_int",
                                      "width_int","width_slope",
                                      "group_sd","r"))
gelman.diag(out7)
save(out7,file="./scratch/jags_flower_output_common_2009_2015_nolup.Rdata",compress=TRUE)

####Plots JAGS Predictions for the Phenology model####

##Creates a data frame with all covariate combinations where we want predictions
jags_pred_data <- expand.grid(sdd_pred=seq(110,230,by=1),datetaken_DOY=seq(100,350,by=1),
                             nearest_center=c("Mowich","Paradise","Sunrise","Tipsoo"),yearfact=2014)
##Inverse logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}

##Function for the fit
jags_fit_fun <- function(height_int,height_slope,
                         opt_int,opt_slope,
                         width_int,width_slope,
                         x1vec,x2vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  
  #Calculates the value of the response.
  height <- height_int + height_slope * x2vec
  opt <- opt_int + opt_slope * x2vec
  width<- exp(width_int + width_slope * x2vec) * -1
  mu <- inv.logit(width * (x1vec - opt)^2 + height)
  return(mu)
}

##Extracts median values of parameters from JAGS output
jags_meds <- summary(out5)$quantiles[,3]
jags_pred_data$predmu <- jags_fit_fun(height_int=jags_meds["height_int"],
                                      height_slope=jags_meds["height_slope"],
                                      width_int=jags_meds["width_int"],
                                      width_slope=jags_meds["width_slope"],
                                      opt_int=jags_meds["opt_int"],
                                      opt_slope=jags_meds["opt_slope"],
                                      x1vec=jags_pred_data$datetaken_DOY/100 - 2,
                                      x2vec=jags_pred_data$sdd_pred/100 - 2)

pdf("./figs/jags_fit_flower_common_snow2x.pdf",width=6,height=4)
ggplot(jags_pred_data)+
  geom_raster(aes(y=datetaken_DOY,x=sdd_pred,fill=predmu),alpha=1)+
  scale_fill_gradientn(colours = rev(rainbow(8,start=0,end=0.6)))+
  geom_abline(slope=1,intercept=0)+
  geom_point(aes(y=datetaken_DOY,x=sdd_pred),
                 shape=16,color="white",alpha=0.2,
             data=filter(f_samp,Flwr==0))+
  geom_point(aes(y=datetaken_DOY,x=sdd_pred),
             shape=17,color="black",alpha=0.5,
             data=filter(f_samp,Flwr==1))+
  scale_color_grey(start=1,end=0)+
  scale_x_continuous(limits=c(110,220),expand=c(0,0))+
  scale_y_continuous(limits=c(110,350),expand=c(0,0))+
  xlab("Day of Snow Melt")+
  ylab("Day of Year")+
  guides(fill=guide_colorbar("Predicted Prob."),
         size=guide_legend("Flowering"))+
  theme_bw()+
  theme(text=element_text(family="sans"))
dev.off()
