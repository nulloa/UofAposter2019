library(rstan)
library(shinystan)
library(ggplot2)
library(latex2exp)
library(tidyr)
library(xtable)
library(dplyr)

# Load in the forecasts mcmc
file_names <- as.list(dir(path="~/Downloads/bfdaRes",pattern="FC",full.names=TRUE))
file_names <- file_names[grep("5|10|15|20", file_names)]

for(files in file_names){
  load(as.character(files))
  assign(sub(".*bfdaRes\\/ *(.*?) *\\.RData.*", "\\1", as.character(files)), mcmc.res)
  rm(mcmc.res)
}



# Load in Overall Data
W <- readRDS("../../bfdaFluScripts/Forecasts/Data/FCfull.rds")
colmeans <- readRDS("../../bfdaFluScripts/Forecasts/Data/colmeansFCfull.rds")
W2 <- W + colmeans
W = reshape2::melt(W)
W2 = reshape2::melt(W2)
names(W) <- c("subj", "week", "centered_ili")
names(W2) <- c("subj", "week", "orig_ili")
season <- readRDS("../../bfdaFluScripts/Forecasts/Data/seasonGrpFCfull.rds")
season$subj <- as.integer(rownames(season))
names(season)[1:2] <- c("Season", "Num_Seas")
region <- readRDS("../../bfdaFluScripts/Forecasts/Data/regionGrpFCfull.rds")
region$subj <- as.integer(rownames(region))
names(region)[1:2] <- c("Region", "Num_Reg")
W <- dplyr::left_join(W, season)
W <- dplyr::left_join(W, region)
W <- dplyr::left_join(W, W2)


# Load in the forecasting data
W5 <- readRDS("../../bfdaFluScripts/Forecasts/Data/FC5.rds")
E <- readRDS("../../bfdaFluScripts/Data/eigenfnReg.rds")

colMeans5 <- readRDS("../../bfdaFluScripts/Forecasts/Data/colmeansFC5.rds")
season5 <- readRDS("../../bfdaFluScripts/Forecasts/Data/seasonGrpFC5.rds")
region5 <- readRDS("../../bfdaFluScripts/Forecasts/Data/regionGrpFC5.rds")

colMeans10 <- readRDS("../../bfdaFluScripts/Forecasts/Data/colmeansFC10.rds")
season10 <- readRDS("../../bfdaFluScripts/Forecasts/Data/seasonGrpFC10.rds")
region10 <- readRDS("../../bfdaFluScripts/Forecasts/Data/regionGrpFC10.rds")

colMeans15 <- readRDS("../../bfdaFluScripts/Forecasts/Data/colmeansFC15.rds")
season15 <- readRDS("../../bfdaFluScripts/Forecasts/Data/seasonGrpFC15.rds")
region15 <- readRDS("../../bfdaFluScripts/Forecasts/Data/regionGrpFC15.rds")

colMeans20 <- readRDS("../../bfdaFluScripts/Forecasts/Data/colmeansFC20.rds")
season20 <- readRDS("../../bfdaFluScripts/Forecasts/Data/seasonGrpFC20.rds")
region20 <- readRDS("../../bfdaFluScripts/Forecasts/Data/regionGrpFC20.rds")



N_subj <- dim(W5)[1]
N_obs <- dim(W5)[2]
N_region <- length(unique(region5$orig_nm))
N_season <- length(unique(season5$orig_nm))

ext.indep.5  <- rstan::extract(FC5_Indep)
ext.region.season.5 <- rstan::extract(FC5_Region_Season)
ext.region.5 <- rstan::extract(FC5_RegionHier)
ext.season.5 <- rstan::extract(FC5_SeasonHier)

ext.indep.10  <- rstan::extract(FC10_Indep)
ext.region.season.10 <- rstan::extract(FC10_Region_Season)
ext.region.10 <- rstan::extract(FC10_RegionHier)
ext.season.10 <- rstan::extract(FC10_SeasonHier)

ext.indep.15  <- rstan::extract(FC15_Indep)
ext.region.season.15 <- rstan::extract(FC15_Region_Season)
ext.region.15 <- rstan::extract(FC15_RegionHier)
ext.season.15 <- rstan::extract(FC15_SeasonHier)

ext.indep.20  <- rstan::extract(FC20_Indep)
ext.region.season.20 <- rstan::extract(FC20_Region_Season)
ext.region.20 <- rstan::extract(FC20_RegionHier)
ext.season.20 <- rstan::extract(FC20_SeasonHier)




# Remove raw stan results
for(files in file_names){
  rm(list=as.character(sub(".*bfdaRes\\/ *(.*?) *\\.RData.*", "\\1", as.character(files))))
}





# Forecasting Results Function
forecastres <- function(mcmc, season, region, E, colmeans, iter, model, forecastlvl){
  if(exists(as.character(substitute(mcmc)))!=TRUE){
    message(paste(as.character(substitute(mcmc))," is missing… :-("))
  }
  
  numseas  <- length(unique(season))
  numreg   <- length(unique(region))
  numgroup <- numseas*numreg
  weeks    <- 1:32
  outdf    <- NULL
  preds    <- centered_preds <- matrix(data=NA, ncol=length(weeks), nrow=iter)
  mpred    <- upred <- lpred <- rep(NA, length(weeks))
  
  for(g in 1:numgroup){
    for(i in 1:iter){
      for(week in weeks){
        centered_preds[i,week] <- rnorm(1, sum(mcmc$beta[i, g, 1:30]*E[week,]), sd=mcmc$eps[i])
        preds[i,week] <- rnorm(1, sum(mcmc$beta[i, g, 1:30]*E[week,]), sd=mcmc$eps[i]) + colmeans[g, week]
      }
    }
    mpred <- colMeans(preds)
    upred <- apply(preds, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    lpred <- apply(preds, 2, quantile, probs = c(0.025),  na.rm = TRUE)
    
    centered_mpred <- colMeans(centered_preds)
    centered_upred <- apply(centered_preds, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    centered_lpred <- apply(centered_preds, 2, quantile, probs = c(0.025),  na.rm = TRUE) 
    
    tmp <- data.frame(week  = weeks,
                      region = region[g],
                      season = season[g],
                      meanForecast=mpred,
                      upperForecast=upred,
                      lowerForecast=lpred,
                      centered_meanForecast=centered_mpred,
                      centered_upperForecast=centered_upred,
                      centered_lowerForecast=centered_lpred,
                      model=model,
                      forecastlvl = forecastlvl)
    outdf <- rbind(outdf, tmp)
  }


  return(outdf)
}

# Function for Beta Results

beta_res <- function(mcmc, N_subj, region, season, model, forecastlvl){
  
  out_betares <- data.frame(est = rep(NA, 30*N_subj),
                            upper = rep(NA, 30*N_subj),
                            lower = rep(NA, 30*N_subj),
                            eignfn = rep(NA, 30*N_subj),
                            season = rep(NA, 30*N_subj),
                            region = rep(NA, 30*N_subj),
                            model = model,
                            forecastlvl = forecastlvl
                            
  )
  
  j = 1
  for(s in 1:N_subj){
    for(e in 1:30){
      out_betares$est[j] <- mean(mcmc$beta[,s,e])
      out_betares$upper[j] <- quantile(mcmc$beta[,s,e], prob=0.975)
      out_betares$lower[j] <- quantile(mcmc$beta[,s,e], prob=0.025)
      out_betares$eignfn[j] <- e
      out_betares$region[j] <- region$orig_nm[s]
      out_betares$season[j] <- season$orig_nm[s]
      j = j+1
    }
  }
  
  return(out_betares)
}






###### Five Week Forecasts #####
indep.fc5 <- forecastres(ext.indep.5, season5$orig_nm, region5$orig_nm, E, colMeans5, 2000, "Indep", "FC5")
indep.fc5$region  <- factor(indep.fc5$region, levels=paste("Region", 1:10))

region.season.fc5 <- forecastres(ext.region.season.5, season5$orig_nm, region5$orig_nm, E, colMeans5, 2000, "Region-Season", "FC5")
region.season.fc5$region  <- factor(region.season.fc5$region, levels=paste("Region", 1:10))

region.fc5 <- forecastres(ext.region.5, season5$orig_nm, region5$orig_nm, E, colMeans5, 2000, "Region", "FC5")
region.fc5$region  <- factor(region.fc5$region, levels=paste("Region", 1:10))

season.fc5 <- forecastres(ext.season.5, season5$orig_nm, region5$orig_nm, E, colMeans5, 2000, "Season", "FC5")
season.fc5$region  <- factor(season.fc5$region, levels=paste("Region", 1:10))

comb_fc5 <- rbind(indep.fc5, region.season.fc5, region.fc5, season.fc5)

comb_fc5 = merge(comb_fc5, W, by.x=c("week","season", "region"), by.y=c("week","Season", "Region"))

fres <- comb_fc5 %>% filter(season %in% paste("16-17") & !(model %in% c("Indep"))) %>% data.frame
ggplot(fres, aes(x=week)) + geom_point(aes(y=centered_ili)) + 
  geom_line(aes(y=centered_meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=centered_lowerForecast, ymax=centered_upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="Centered ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("centeredfc5.png", width=9, height=6, units="in")

ggplot(fres, aes(x=week)) + geom_point(aes(y=orig_ili)) + 
  geom_line(aes(y=meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=lowerForecast, ymax=upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("fc5.png", width=9, height=6, units="in")


# Get Beta Results
indep_betares <- beta_res(ext.indep.5, N_subj, region5, season5, "Indep", "FC5")
indep_betares$region  <- factor(indep_betares$region, levels=paste("Region", 1:10))
region_season_betares <- beta_res(ext.region.season.5, N_subj, region5, season5, "Region-Season", "FC5")
region_season_betares$region  <- factor(region_season_betares$region, levels=paste("Region", 1:10))
region_betares <- beta_res(ext.region.5, N_subj, region5, season5, "Region", "FC5")
region_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))
season_betares <- beta_res(ext.season.5, N_subj, region5, season5, "Season", "FC5")
season_betares$region <- factor(season_betares$region, levels=paste("Region", 1:10))

comb_betares <- rbind(indep_betares, region_season_betares, region_betares, season_betares)

# All regions and seasons
fbeta_overall <- comb_betares %>% filter(eignfn %in% paste(1:10) & !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_overall, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_overall_fc5.png", width=9, height=9, units="in")

# Same Region different Seasons
fbeta_reg <- comb_betares %>% filter(region %in% paste("Region", c(7,9)) &
                                     eignfn %in% paste(1:10) &
                                     !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_reg, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_region_fc5.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_seas <- comb_betares %>% filter(season %in% paste0(16,"-",17) &
                                      eignfn %in% paste(1:10) &
                                      !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_seas, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) + 
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(season~region) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_seas_fc5.png", width=9, height=6, units="in")


















###### Ten Week Forecasts #####
indep.fc10 <- forecastres(ext.indep.10, season10$orig_nm, region10$orig_nm, E, colMeans10, 2000, "Indep", "FC10")
indep.fc10$region  <- factor(indep.fc10$region, levels=paste("Region", 1:10))

region.season.fc10 <- forecastres(ext.region.season.10, season10$orig_nm, region10$orig_nm, E, colMeans10, 2000, "Region-Season", "FC10")
region.season.fc10$region  <- factor(region.season.fc10$region, levels=paste("Region", 1:10))

region.fc10 <- forecastres(ext.region.10, season10$orig_nm, region10$orig_nm, E, colMeans10, 2000, "Region", "FC10")
region.fc10$region  <- factor(region.fc10$region, levels=paste("Region", 1:10))

season.fc10 <- forecastres(ext.season.10, season10$orig_nm, region10$orig_nm, E, colMeans10, 2000, "Season", "FC10")
season.fc10$region  <- factor(season.fc10$region, levels=paste("Region", 1:10))

comb_fc10 <- rbind(indep.fc10, region.season.fc10, region.fc10, season.fc10)

comb_fc10 = merge(comb_fc10, W, by.x=c("week","season", "region"), by.y=c("week","Season", "Region"))

fres <- comb_fc10 %>% filter(season %in% paste("16-17") & !(model %in% c("Indep"))) %>% data.frame
ggplot(fres, aes(x=week)) + geom_point(aes(y=centered_ili)) + 
  geom_line(aes(y=centered_meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=centered_lowerForecast, ymax=centered_upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="Centered ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("centeredfc10.png", width=9, height=6, units="in")

ggplot(fres, aes(x=week)) + geom_point(aes(y=orig_ili)) + 
  geom_line(aes(y=meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=lowerForecast, ymax=upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("fc10.png", width=9, height=6, units="in")


# Get Beta Results
indep_betares <- beta_res(ext.indep.10, N_subj, region10, season10, "Indep", "FC10")
indep_betares$region  <- factor(indep_betares$region, levels=paste("Region", 1:10))
region_season_betares <- beta_res(ext.region.season.10, N_subj, region10, season10, "Region-Season", "FC10")
region_season_betares$region <- factor(region_season_betares$region, levels=paste("Region", 1:10))
region_betares <- beta_res(ext.region.10, N_subj, region10, season10, "Region", "FC10")
region_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))
season_betares <- beta_res(ext.season.10, N_subj, region10, season10, "Season", "FC10")
season_betares$region <- factor(season_betares$region, levels=paste("Region", 1:10))

comb_betares <- rbind(indep_betares, region_season_betares, region_betares, season_betares)

# All regions and seasons
fbeta_overall <- comb_betares %>% filter(eignfn %in% paste(1:10) & !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_overall, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_overall_fc10.png", width=9, height=9, units="in")

# Same Region different Seasons
fbeta_reg <- comb_betares %>% filter(region %in% paste("Region", c(7,9)) &
                                     eignfn %in% paste(1:10) &
                                     !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_reg, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_region_fc10.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_seas <- comb_betares %>% filter(season %in% paste0(16,"-",17) &
                                      eignfn %in% paste(1:10) &
                                      !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_seas, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) + 
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(season~region) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_seas_fc10.png", width=9, height=6, units="in")




















###### Fifteen Week Forecasts #####
indep.fc15 <- forecastres(ext.indep.15, season15$orig_nm, region15$orig_nm, E, colMeans15, 2000, "Indep", "FC15")
indep.fc15$region  <- factor(indep.fc15$region, levels=paste("Region", 1:10))

region.season.fc15 <- forecastres(ext.region.season.15, season15$orig_nm, region15$orig_nm, E, colMeans15, 2000, "Region-Season", "FC15")
region.season.fc15$region  <- factor(region.season.fc15$region, levels=paste("Region", 1:10))

region.fc15 <- forecastres(ext.region.15, season15$orig_nm, region15$orig_nm, E, colMeans15, 2000, "Region", "FC15")
region.fc15$region  <- factor(region.fc15$region, levels=paste("Region", 1:10))

season.fc15 <- forecastres(ext.season.15, season15$orig_nm, region15$orig_nm, E, colMeans15, 2000, "Season", "FC15")
season.fc15$region  <- factor(season.fc15$region, levels=paste("Region", 1:10))

comb_fc15 <- rbind(indep.fc15, region.season.fc15, region.fc15, season.fc15)

comb_fc15 = merge(comb_fc15, W, by.x=c("week","season", "region"), by.y=c("week","Season", "Region"))

fres <- comb_fc15 %>% filter(season %in% paste("16-17") & !(model %in% c("Indep"))) %>% data.frame
ggplot(fres, aes(x=week)) + geom_point(aes(y=centered_ili)) + 
  geom_line(aes(y=centered_meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=centered_lowerForecast, ymax=centered_upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="Centered ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("centeredfc15.png", width=9, height=6, units="in")

ggplot(fres, aes(x=week)) + geom_point(aes(y=orig_ili)) + 
  geom_line(aes(y=meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=lowerForecast, ymax=upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("fc15.png", width=9, height=6, units="in")


# Get Beta Results
indep_betares <- beta_res(ext.indep.15, N_subj, region15, season15, "Indep", "FC15")
indep_betares$region  <- factor(indep_betares$region, levels=paste("Region", 1:10))
region_season_betares <- beta_res(ext.region.season.15, N_subj, region15, season15, "Region-Season", "FC15")
region_season_betares$region <- factor(region_season_betares$region, levels=paste("Region", 1:10))
region_betares <- beta_res(ext.region.15, N_subj, region15, season15, "Region", "FC15")
region_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))
season_betares <- beta_res(ext.season.15, N_subj, region15, season15, "Season", "FC15")
season_betares$region <- factor(season_betares$region, levels=paste("Region", 1:10))

comb_betares <- rbind(indep_betares, region_season_betares, region_betares, season_betares)

# All regions and seasons
fbeta_overall <- comb_betares %>% filter(eignfn %in% paste(1:10) & !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_overall, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_overall_fc15.png", width=9, height=9, units="in")

# Same Region different Seasons
fbeta_reg <- comb_betares %>% filter(region %in% paste("Region", c(7,9)) &
                                     eignfn %in% paste(1:10) &
                                     !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_reg, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_region_fc15.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_seas <- comb_betares %>% filter(season %in% paste0(16,"-",17) &
                                      eignfn %in% paste(1:10) &
                                      !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_seas, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) + 
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(season~region) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_seas_fc15.png", width=9, height=6, units="in")











###### Twenty Week Forecasts #####
indep.fc20 <- forecastres(ext.indep.20, season20$orig_nm, region20$orig_nm, E, colMeans20, 2000, "Indep", "FC20")
indep.fc20$region  <- factor(indep.fc20$region, levels=paste("Region", 1:10))

region.season.fc20 <- forecastres(ext.region.season.20, season20$orig_nm, region20$orig_nm, E, colMeans20, 2000, "Region-Season", "FC20")
region.season.fc20$region  <- factor(region.season.fc20$region, levels=paste("Region", 1:10))

region.fc20 <- forecastres(ext.region.20, season20$orig_nm, region20$orig_nm, E, colMeans20, 2000, "Region", "FC20")
region.fc20$region  <- factor(region.fc20$region, levels=paste("Region", 1:10))

season.fc20 <- forecastres(ext.season.20, season20$orig_nm, region20$orig_nm, E, colMeans20, 2000, "Season", "FC20")
season.fc20$region  <- factor(season.fc20$region, levels=paste("Region", 1:10))

comb_fc20 <- rbind(indep.fc20, region.season.fc20, region.fc20, season.fc20)

comb_fc20 = merge(comb_fc20, W, by.x=c("week","season", "region"), by.y=c("week","Season", "Region"))

fres <- comb_fc20 %>% filter(season %in% paste("16-17") & !(model %in% c("Indep"))) %>% data.frame
ggplot(fres, aes(x=week)) + geom_point(aes(y=centered_ili)) + 
  geom_line(aes(y=centered_meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=centered_lowerForecast, ymax=centered_upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="Centered ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("centeredfc20.png", width=9, height=6, units="in")

ggplot(fres, aes(x=week)) + geom_point(aes(y=orig_ili)) + 
  geom_line(aes(y=meanForecast, color=model)) + 
  geom_ribbon(aes(ymin=lowerForecast, ymax=upperForecast, fill=model), alpha=0.35) + 
  facet_grid(season~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("fc20.png", width=9, height=6, units="in")


# Get Beta Results
indep_betares <- beta_res(ext.indep.20, N_subj, region20, season20, "Indep", "FC20")
indep_betares$region  <- factor(indep_betares$region, levels=paste("Region", 1:10))
region_season_betares <- beta_res(ext.region.season.20, N_subj, region20, season20, "Region-Season", "FC20")
region_season_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))
region_betares <- beta_res(ext.region.20, N_subj, region20, season20, "Region", "FC20")
region_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))
season_betares <- beta_res(ext.season.20, N_subj, region20, season20, "Season", "FC20")
season_betares$region <- factor(season_betares$region, levels=paste("Region", 1:10))

comb_betares <- rbind(indep_betares, region_season_betares, region_betares, season_betares)

# All regions and seasons
fbeta_overall <- comb_betares %>% filter(eignfn %in% paste(1:10) & !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_overall, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_overall_fc20.png", width=9, height=9, units="in")

# Same Region different Seasons
fbeta_reg <- comb_betares %>% filter(region %in% paste("Region", c(7,9)) &
                                     eignfn %in% paste(1:10) &
                                     !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_reg, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_region_fc20.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_seas <- comb_betares %>% filter(season %in% paste0(16,"-",17) &
                                      eignfn %in% paste(1:10) &
                                      !(model %in% c("Indep"))) %>% data.frame
ggplot(fbeta_seas, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) + 
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(season~region) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_seas_fc20.png", width=9, height=6, units="in")



















##### Forecasting MSE #####
comb_fc <- rbind(comb_fc5, comb_fc10, comb_fc15, comb_fc20)

mse = comb_fc %>% 
  group_by(forecastlvl, model) %>% 
  summarise(
    MSE = mean((orig_ili - meanForecast)^2)
  ) %>% 
  mutate(
    weeks = substr(forecastlvl, 3, 5)
  ) %>%
  select(
    Model=model,
    `Weeks Used` = weeks,
    MSE
  ) %>% data.frame


mse$forecastlvl <- NULL
names(mse)[2] <- "Weeks Used"
print(
  xtable(mse, caption="This table shows the Mean Square Error of all four models at all four times. 
         The the hierarchical models outpreform the independent model in all cases and the region 
         preforms best with the exception of the forecast using 20 weeks. At 20 weeks, the region-season
         model preforms best though it only slightly outpreforms the region model.",
         digits=c(0,3,3,3), label = "tab:fc_mse"
  ), 
  include.rownames=FALSE, file = "fc_mse_tab.tex")

print(
  xtable(t(mse), caption="This table shows the Mean Square Error of all four models at all four times. 
         The the hierarchical models outpreform the independent model in all cases and the region 
         preforms best with the exception of the forecast using 20 weeks. At 20 weeks, the region-season
         model preforms best though it only slightly outpreforms the region model.",
         digits=c(0,rep(3, nrow(mse))), label = "tab:fc_mse"
  ), 
  include.colnames=FALSE, file = "fc_mse_tab_t.tex")

