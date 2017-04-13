##Script to estimate overlap coefficients for visitor and flower phenology
##Ian Breckheimer
##Initally Created November 13th, 2015
##Last updated January 2017

####Sets up workspace####
library(dplyr)
library(rjags)
library(ggmcmc)
library(doParallel)
library(foreach)
setwd("~/code/MORA_pheno/")
source("./scripts/overlap_functions.R")

####Analysis and figures using the common-parameter model####
##Brings in fit model objects and converts them to data frames.
load("./scratch/visitors_jags_common_2009_2015.Rdata")
visit_model <- out2
load("./scratch/jags_flower_output_common_2009_2015.Rdata")
flower_model <- out2

##Processes model output to measure phenological mismatch for each site.

##Names and days for measurement
site_names <- c("Mowich","Paradise","Sunrise","Chinook Pass")
measure_sdd <- seq(110,240,by=10)

##Function that does the heavy-lifting.
jags_mismatch <- function(visit_model,flower_model,group_index,measure_sdd,
                          site_names){
  params <- prep_mcmc_vars(visit_model,flower_model,group_index=group_index,
                           mod1_param_names = c("height_int_site[2]","height_slope",
                                                "opt_int","opt_slope",
                                                "width_int","width_slope"),
                           mod2_param_names = c("height_int","height_slope",
                                                "opt_int","opt_slope",
                                                "width_int","width_slope"))
  mismatch <- measure_mismatch(params=params,x2sdd=measure_sdd)
  mismatch$Site <- site_names[i]
  return(mismatch)
}

##Registers a parallel backend.
cl <- makePSOCKcluster(2)
registerDoParallel(cl)

##Runs the computation in parallel.
mismatch <- foreach(i=1:length(site_names),.combine='rbind') %dopar% {
                  jags_mismatch(visit_model,flower_model,group_index=i,measure_sdd=measure_sdd,
                                site_names=site_names)
}
stopCluster(cl)

##Graphs the results.
month_breaks <- c(121,152,182,213,244)
month_labels <- c("May","June","July","Aug","Sept")

ggplot(data=mismatch)+
  geom_ribbon(aes(x=SDD,ymax=overlap_upr,ymin=overlap_lwr),fill="grey20",alpha=0.2)+
#  geom_ribbon(aes(x=SDD,ymax=overlap_q10,ymin=overlap_q90),fill="grey20",alpha=0.4)+
  geom_ribbon(aes(x=SDD,ymax=overlap_q25,ymin=overlap_q75),fill="grey20",alpha=0.6)+
  geom_line(aes(x=SDD,y=overlap_q50))+
  scale_x_continuous(breaks = month_breaks,labels = month_labels )+
  xlab("Last Snow Melt")+
  ylab("Phenological Match")+
  facet_wrap(facets=~Site)+
  ylim(c(0,1))+
  theme_bw()+
  theme(text=element_text(family="sans"))

##Graphs functions and overlap at 2011 and 2015 Snow Disappearance Dates
paramelt2015 <- 139
paramelt2011 <- 211
xseq <- seq(1,365,by=1) / 100 - 2
x2_2015 <- paramelt2015 / 100 - 2
x2_2011 <- paramelt2011 / 100 - 2

params <- prep_mcmc_vars(model_output1=visit_model,model_output2=flower_model,group_index=2,
                         mod1_param_names = c("height_int_site[2]","height_slope",
                                              "opt_int","opt_slope",
                                              "width_int","width_slope"),
                         mod2_param_names = c("height_int","height_slope",
                                              "opt_int","opt_slope",
                                              "width_int","width_slope"))
visit_params <- params[,1:6]
flower_params <- params[,7:12]

visit_fun_2011 <- function(x) { f1_dens(x1=xseq, x2=x2_2011, 
                              height_int=x[1],height_slope=x[2],
                              opt_int=x[3],opt_slope=x[4],
                              width_int=x[5],width_slope=x[6])}
visits_2011 <- apply(visit_params, FUN=visit_fun_2011,MARGIN = 1)

visit_fun_2015 <- function(x) { f1_dens(x1=xseq, x2=x2_2015, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
visits_2015 <- apply(visit_params, FUN=visit_fun_2015,MARGIN = 1)

flower_fun_2011 <- function(x) { f2_dens(x1=xseq, x2=x2_2011, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
flowers_2011 <- apply(flower_params, FUN=flower_fun_2011,MARGIN = 1)

flower_fun_2015 <- function(x) { f2_dens(x1=xseq, x2=x2_2015, 
                                   height_int=x[1],height_slope=x[2],
                                   opt_int=x[3],opt_slope=x[4],
                                   width_int=x[5],width_slope=x[6])}
flowers_2015 <- apply(flower_params, FUN=flower_fun_2015,MARGIN = 1)

visits_quant_2011 <- data.frame(t(apply(visits_2011,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
visits_quant_2011$group <- "Visitors"
visits_quant_2011$Year <- 2011


visits_quant_2015 <- data.frame(t(apply(visits_2015,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
visits_quant_2015$group <- "Visitors"
visits_quant_2015$Year <- 2015

flowers_quant_2011 <- data.frame(t(apply(flowers_2011,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
flowers_quant_2011$group <- "Flowers"
flowers_quant_2011$Year <- 2011

flowers_quant_2015 <- data.frame(t(apply(flowers_2015,MARGIN = 1,
                           FUN=function(x){quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))})))
flowers_quant_2015$group <- "Flowers"
flowers_quant_2015$Year <- 2015

visit_means <- colMeans(visit_params)
flower_means <- colMeans(flower_params)

match_2011 <- min_f1f2_dens(x1=xseq,x2=x2_2011,
                            height_int=visit_means[1],height_slope=visit_means[2],
                            opt_int=visit_means[3],opt_slope=visit_means[4],
                            width_int=visit_means[5], width_slope=visit_means[6],
                            height_int2=flower_means[1],height_slope2=flower_means[2],
                            opt_int2=flower_means[3],opt_slope2=flower_means[4],
                            width_int2=flower_means[5],width_slope2=flower_means[6])
match_2015 <- min_f1f2_dens(x1=xseq,x2=x2_2015,
                            height_int=visit_means[1],height_slope=visit_means[2],
                            opt_int=visit_means[3],opt_slope=visit_means[4],
                            width_int=visit_means[5], width_slope=visit_means[6],
                            height_int2=flower_means[1],height_slope2=flower_means[2],
                            opt_int2=flower_means[3],opt_slope2=flower_means[4],
                            width_int2=flower_means[5],width_slope2=flower_means[6])

match_2011 <- data.frame(doy=(xseq + 2) * 100,match=match_2011)
match_2011$year <- 2011
match_2015 <- data.frame(doy=(xseq + 2) * 100,match=match_2015)
match_2015$year <- 2015
match_2011_2015 <- rbind(match_2011,match_2015)

flowers_visitors <- rbind(visits_quant_2015,visits_quant_2011,
                          flowers_quant_2015,flowers_quant_2011)
flowers_visitors$DOY <- rep((xseq + 2) * 100,4)
colnames(flowers_visitors) <- c("lwr","lwr25","median","upr75","upr","group","year","doy")

##Graphs the results.
month_breaks <- c(1,32,60,91,121,152,182,213,244,274,305,335)
month_labels <- c("Jan","","Mar","","May","",
                  "July","","Sept","","Nov","")

ggplot(data=flowers_visitors)+
  geom_ribbon(aes(x=doy,ymin=lwr,ymax=upr,fill=group,alpha=factor(year)))+
  geom_ribbon(aes(x=doy,ymax=match,ymin=0,alpha=factor(year)),fill="grey60",
              data=match_2011_2015)+
  geom_line(aes(x=doy,y=median,color=group,alpha=factor(year)))+
  scale_x_continuous(limits=c(0,366),breaks=month_breaks,labels=month_labels)+
  scale_alpha_discrete(range=c(0.3,0.9))+
#  facet_wrap(facets=~year)+
  guides(alpha="none",fill="none",color="none")+
  xlab("Day of Year")+
  ylab("Density")+
  theme_bw()+
  theme(panel.grid=element_blank())

###Graphs the estimates of the phenological peak for people and flowers.

params <- prep_mcmc_vars(visit_model,flower_model,group_index=2,
                         mod1_param_names = c("height_int_site[2]","height_slope",
                                              "opt_int","opt_slope",
                                              "width_int","width_slope"),
                         mod2_param_names = c("height_int","height_slope",
                                              "opt_int","opt_slope",
                                              "width_int","width_slope"))
visit_params <- params[,1:6]
flower_params <- params[,7:12]


##Creates dataset with quantiles for peak in each model.
x1seq_unscaled <- seq(1,365,by=2)
x2seq_unscaled <- seq(80,300,by=2)

x2seq <- x2seq_unscaled / 100 - 2
opt_fun <- function(x,x2seq=x2seq){x[3] + x[4] * x2seq}
visit_opts <- apply(visit_params,FUN=opt_fun,MARGIN = 1,x2seq=x2seq)
flower_opts <- apply(flower_params,FUN=opt_fun,MARGIN = 1,x2seq=x2seq)
quantile_fun <- function(x) {quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))}
visit_quants <- data.frame(t((apply(visit_opts, FUN=quantile_fun, MARGIN=1) + 2) * 100))
colnames(visit_quants) <- c("lwr","lwr25","median","upr75","upr")
visit_quants$SDD <- (x2seq + 2) * 100
visit_quants$group <- "visitors"
flower_quants <- data.frame(t((apply(flower_opts, FUN=quantile_fun, MARGIN=1) + 2) * 100))
colnames(flower_quants) <- c("lwr","lwr25","median","upr75","upr")
flower_quants$SDD <- (x2seq + 2) * 100
flower_quants$group <- "flowers"
opt_quants <- rbind(flower_quants,visit_quants)
opt_quants$stat <- "optimum"

##Creates dataset with quantiles for width in each model.
xseq_unscaled <- expand.grid(DOY=x1seq_unscaled,melt_DOY=x2seq_unscaled)
xseq <- xseq_unscaled / 100 - 2
pred_fun_flwr <- function(x) { f2(x1=xseq[,1], x2=xseq[,2], 
                                        height_int=x[1],height_slope=x[2],
                                        opt_int=x[3],opt_slope=x[4],
                                        width_int=x[5],width_slope=x[6])}
pred_fun_visit <- function(x) { f1(x1=xseq[,1], x2=xseq[,2], 
                                  height_int=x[1],height_slope=x[2],
                                  opt_int=x[3],opt_slope=x[4],
                                  width_int=x[5],width_slope=x[6])}
flower_preds <- apply(flower_params,FUN=pred_fun_flwr,MARGIN=1)
visit_preds <- apply(visit_params,FUN=pred_fun_visit,MARGIN=1)

flower_start_mat <- matrix(NA,nrow=dim(flower_preds)[2],ncol=length(x2seq_unscaled))
flower_end_mat <- matrix(NA,nrow=dim(flower_preds)[2],ncol=length(x2seq_unscaled))
for (i in 1:dim(flower_preds)[2]){
  print(paste("Rep. ",i))
  flower_frame <- data.frame(DOY=xseq$DOY,melt_DOY=xseq$melt_DOY,pred=flower_preds[,i])
  flower_pred_grp <- group_by(flower_frame,melt_DOY)
  width <- summarise(flower_pred_grp,daystart=est_width(DOY,pred,threshold=0.25)[1],
                     dayend=est_width(DOY,pred,threshold=0.25)[2])
  flower_start_mat[i,] <- as.numeric(width$daystart)
  flower_end_mat[i,] <- as.numeric(width$dayend)
}

visit_start_mat <- matrix(NA,nrow=dim(visit_preds)[2],ncol=length(x2seq_unscaled))
visit_end_mat <- matrix(NA,nrow=dim(visit_preds)[2],ncol=length(x2seq_unscaled))
for (i in 1:dim(visit_preds)[2]){
  print(paste("Rep. ",i))
  visit_frame <- data.frame(DOY=xseq$DOY,melt_DOY=xseq$melt_DOY,pred=visit_preds[,i])
  visit_pred_grp <- group_by(visit_frame,melt_DOY)
  width <- summarise(visit_pred_grp,daystart=est_width(DOY,pred,threshold=0.25)[1],
                     dayend=est_width(DOY,pred,threshold=0.25)[2])
  visit_start_mat[i,] <- as.numeric(width$daystart)
  visit_end_mat[i,] <- as.numeric(width$dayend)
}

quantile_fun <- function(x) {quantile(x,probs=c(0.025,0.25,0.5,0.75,0.975))}

visit_start_quants <- (data.frame(t(apply(visit_start_mat, FUN=quantile_fun, MARGIN=2))) + 2) * 100
colnames(visit_start_quants) <- c("lwr","lwr25","median","upr75","upr")
visit_start_quants$SDD <- x2seq_unscaled
visit_start_quants$group <- "visitors"
visit_end_quants <- (data.frame(t(apply(visit_end_mat, FUN=quantile_fun, MARGIN=2))) + 2) * 100
colnames(visit_end_quants) <- c("lwr","lwr25","median","upr75","upr")
visit_end_quants$SDD <- x2seq_unscaled
visit_end_quants$group <- "visitors"


flower_start_quants <- (data.frame(t(apply(flower_start_mat, FUN=quantile_fun, MARGIN=2))) + 2) * 100
colnames(flower_start_quants) <- c("lwr","lwr25","median","upr75","upr")
flower_start_quants$SDD <- x2seq_unscaled
flower_start_quants$group <- "flowers"
flower_end_quants <- (data.frame(t(apply(flower_end_mat, FUN=quantile_fun, MARGIN=2))) + 2) * 100
colnames(flower_end_quants) <- c("lwr","lwr25","median","upr75","upr")
flower_end_quants$SDD <- x2seq_unscaled
flower_end_quants$group <- "flowers"

start_quants <- rbind(flower_start_quants,visit_start_quants)
end_quants <- rbind(flower_end_quants,visit_end_quants)

x_month_breaks <- c(91,121,152,182,213,244)
x_month_labels <- c("April","May","June","July","Aug","Sept")
y_month_breaks <- c(1,32,60,91,121,152,182,213,244,274,305,335)
y_month_labels <- c("Jan","Feb","Mar","April","May","June",
                  "July","Aug","Sept","Oct","Nov","Dec")

ggplot()+
#  geom_ribbon(aes(x=SDD,ymin=lwr,ymax=upr,fill=group),alpha=0.2,data=opt_quants)+
#  geom_ribbon(aes(x=SDD,ymin=lwr25,ymax=upr75,fill=group),alpha=0.5,data=opt_quants)+
#  geom_line(aes(x=SDD,y=median,color=group),data=opt_quants)+
  geom_ribbon(aes(x=start_quants$SDD,ymin=start_quants$median,ymax=end_quants$median,
                  fill=start_quants$group),alpha=0.2)+
  geom_line(aes(x=start_quants$SDD,y=start_quants$median,
                  color=start_quants$group),linetype="dotted")+
  geom_line(aes(x=end_quants$SDD,y=end_quants$median,
                color=end_quants$group),linetype="dotted")+
  geom_line(aes(x=SDD,y=median,color=group),linetype="solid",data=opt_quants)+
  scale_fill_discrete(name="Peak Season")+
  scale_color_discrete(name="Peak Season")+
  scale_x_continuous(limits=c(90,240),breaks=x_month_breaks,labels=x_month_labels)+
  scale_y_continuous(limits=c(70,280),breaks=y_month_breaks,labels=y_month_labels)+
  geom_abline(aes(intercept=0,slope=1))+
  xlab("Snow Disappearance Date")+
  ylab("Day of Year")+
  theme_bw()

####Analysis of site-specific overlap.

##Brings in fit model objects and converts them to data frames.
load("./scratch/visitors_jags_access_2009_2015.Rdata")
visit_model <- out4
load("./scratch/jags_flower_output_siteoptwidth_sideslope_2009_2015.Rdata")
flower_model <- out5

##Processes model output to measure phenological mismatch for each site.

##Names and days for measurement
site_names <- c("Mowich","Paradise","Sunrise","Chinook Pass")
access_indices <- c("[2]","[1]","[2]","[1]")
side_indices <- c("[1]","[1]","[2]","[2]")
measure_sdd <- seq(110,240,by=10)

##Registers a parallel backend.
cl <- makePSOCKcluster(2)
registerDoParallel(cl)

##Runs the computation in parallel.
mismatch_site <- foreach(i=1:length(site_names),.combine='rbind') %dopar% {
  jags_mismatch_site(visit_model,flower_model,group_index=i,measure_sdd=measure_sdd,
                site_names=site_names,access_indices=access_indices,
                side_indices=side_indices)
}
stopCluster(cl)

cl <- makePSOCKcluster(2)
registerDoParallel(cl)
mismatch_change_site <- foreach(i=1:length(site_names),.combine='rbind') %dopar% {
  jags_mismatch_site_change(visit_model,flower_model,group_index=i,measure_sdd=c(130,211),
                     site_names=site_names,access_indices=access_indices,
                     side_indices=side_indices)
}
stopCluster(cl)

##Graphs the results.
month_breaks <- c(121,152,182,213,244)
month_labels <- c("May","June","July","Aug","Sept")

ggplot(data=mismatch_site)+
  geom_ribbon(aes(x=SDD,ymax=overlap_upr,ymin=overlap_lwr),fill="grey20",alpha=0.2)+
  #  geom_ribbon(aes(x=SDD,ymax=overlap_q10,ymin=overlap_q90),fill="grey20",alpha=0.4)+
  geom_ribbon(aes(x=SDD,ymax=overlap_q25,ymin=overlap_q75),fill="grey20",alpha=0.6)+
  geom_line(aes(x=SDD,y=overlap_q50))+
  scale_x_continuous(breaks = month_breaks,labels = month_labels )+
  xlab("Last Snow Melt")+
  ylab("Phenological Match")+
  facet_wrap(facets=~Site)+
  ylim(c(0,1))+
  theme_bw()+
  theme(text=element_text(family="sans"))