library(rstan)
library(shinystan)
library(ggplot2)
library(latex2exp)
library(tidyr)
library(xtable)
library(dplyr)
library(coda)



# Geweke functions to be used
stan_to_coda <- function(fit){
  t <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  mcmc <- mcmc.list(lapply(1:ncol(t), function(x) mcmc(t[,x,])))
  return(mcmc)
}

geweke_stan <- function(x){
  x <- stan_to_coda(x)
  geweke <- geweke.diag(x)
  unlistgd <- unlist(geweke)
  return(list(beta = unlistgd[grep("z.beta",names(unlistgd))],
              tau = unlistgd[grep("z.tau",names(unlistgd))],
              lambda = unlistgd[grep("z.lamb",names(unlistgd))],
              eps = unlistgd[grep("z.eps",names(unlistgd))])
  )
}

# Fitted Results Function
fitres <- function(mcmc, season, region, W, E, colmeans, iter, model){
  if(exists(as.character(substitute(mcmc)))!=TRUE){
    message(paste(as.character(substitute(mcmc))," is missing… :-("))
  }
  
  numseas  <- length(unique(season))
  numreg   <- length(unique(region))
  numgroup <- numseas*numreg
  weeks    <- 1:32
  outdf    <- NULL
  fit      <- centered_mfit <- matrix(data=NA, ncol=length(weeks), nrow=iter)
  pred     <- centered_mpred <- matrix(data=NA, ncol=length(weeks), nrow=iter)
  mfit     <- ufit <- lfit <- rep(NA, length(weeks))
  mpred    <- upred <- lpred <- rep(NA, length(weeks))
  
  for(g in 1:numgroup){
    for(i in 1:iter){
      for(week in weeks){
        centered_mfit[i,week] <- mcmc$m[i, g, week]
        fit[i,week] <- mcmc$m[i, g, week] + colmeans[g, week]
        
        tmppred <- rnorm(1, mcmc$m[i, g, week], sd=mcmc$eps[i])
        centered_mpred[i,week] <- tmppred
        pred[i,week] <- tmppred + colmeans[g, week]
      }
    }
    # Mean Fit
    mfit <- colMeans(fit)
    ufit <- apply(fit, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    lfit <- apply(fit, 2, quantile, probs = c(0.025),  na.rm = TRUE)
    
    centered_fit <- colMeans(centered_mfit)
    centered_ufit <- apply(centered_mfit, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    centered_lfit <- apply(centered_mfit, 2, quantile, probs = c(0.025),  na.rm = TRUE)
    
    # Predicted Fit
    mpred <- colMeans(pred)
    upred <- apply(pred, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    lpred <- apply(pred, 2, quantile, probs = c(0.025),  na.rm = TRUE)
    
    centered_pred <- colMeans(centered_mpred)
    centered_upred <- apply(centered_mpred, 2, quantile, probs = c(0.975),  na.rm = TRUE) 
    centered_lpred <- apply(centered_mpred, 2, quantile, probs = c(0.025),  na.rm = TRUE) 
    
    tmp <- data.frame(week  = weeks,
                      ili = W[g,] + colmeans[g,],
                      centered_ili = W[g,],
                      region = region[g],
                      season = season[g],
                      meanFit=mfit,
                      upperFit=ufit,
                      lowerFit=lfit,
                      centered_meanFit=centered_fit,
                      centered_upperFit=centered_ufit,
                      centered_lowerFit=centered_lfit,
                      meanPred=mpred,
                      upperPred=upred,
                      lowerPred=lpred,
                      centered_meanPred=centered_pred,
                      centered_upperPred=centered_upred,
                      centered_lowerPred=centered_lpred,
                      model=model)
    outdf <- rbind(outdf, tmp)
  }
  
  
  return(outdf)
}

# Function for Beta Results
beta_res <- function(mcmc, N_subj, region, season, model){
  
  out_betares <- data.frame(est = rep(NA, 30*N_subj),
                            upper = rep(NA, 30*N_subj),
                            lower = rep(NA, 30*N_subj),
                            eignfn = rep(NA, 30*N_subj),
                            season = rep(NA, 30*N_subj),
                            region = rep(NA, 30*N_subj),
                            model = model)
  
  j = 1
  for(s in 1:N_subj){
    for(e in 1:30){
      out_betares$est[j] <- mean(mcmc$beta[,s,e])
      out_betares$upper[j] <- quantile(mcmc$beta[,s,e], prob=0.975)
      out_betares$lower[j] <- quantile(mcmc$beta[,s,e], prob=0.025)
      out_betares$eignfn[j] <- e
      out_betares$region[j] <- region[s]
      out_betares$season[j] <- season[s]
      j = j+1
    }
  }
  
  return(out_betares)
}

# Function to look at postior prob of beta
beta_postprob <- function(mcmc, N_subj, region, season, model){
  
  out_betares <- data.frame(prob   = rep(NA, 30*N_subj),
                            eignfn = rep(NA, 30*N_subj),
                            season = rep(NA, 30*N_subj),
                            region = rep(NA, 30*N_subj),
                            model  = model)
  
  j = 1
  for(s in 1:N_subj){
    for(e in 1:30){
      out_betares$prob[j] <- mean(mcmc$beta[,s,e] < -1 | mcmc$beta[,s,e] > 1)
      out_betares$eignfn[j] <- e
      out_betares$region[j] <- region[s]
      out_betares$season[j] <- season[s]
      j = j+1
    }
  }
  
  return(out_betares)
}



# Load in the data
files = list.files(path = "~/Downloads/bfdaRes/OGfit", pattern = ".RData")
for(f in files){
  cat("Loading",f,"\n")
  load(paste0("~/Downloads/bfdaRes/OGfit/", f))
  assign(paste0(sapply(strsplit(f, "\\."), `[`, 1), ".mcmc"), mcmc.res)
  
  # if(f %in% c("LASSO_Region.RData", "LASSO_Season.RData")){}
  # else{
    cat("Est Geweke Diag for",f,"\n")
    assign(paste0("geweke.", sapply(strsplit(f, "\\."), `[`, 1)), geweke_stan(mcmc.res))
  # }
  
  cat("Extracting",f,"\n")
  ext.mcmc <- rstan::extract(mcmc.res)
  assign(paste0("ext.", sapply(strsplit(f, "\\."), `[`, 1)), ext.mcmc)
  cat("Removing mcmc.res\n")
  rm(mcmc.res, ext.mcmc)
}

W <- readRDS("~/Github/bfdaFluScripts/Data/matdatWReg.rds")
E <- readRDS("~/Github/bfdaFluScripts/Data/eigenfnReg.rds")

colMeans <- readRDS("~/Github/bfdaFluScripts/Data/colmeans.rds")
season <- readRDS("~/Github/bfdaFluScripts/Data/seasonGrp.rds")
region <- readRDS("~/Github/bfdaFluScripts/Data/regionGrp.rds")

N_subj <- dim(W)[1]
N_obs <- dim(W)[2]
N_region <- length(unique(region$orig_nm))
N_season <- length(unique(season$orig_nm))





##### Create data frame of actual data and estimates to plot #####

##### Create data frame of actual data and estimates to plot #####
indep <- fitres(ext.Indep, season$orig_nm, region$orig_nm, W, E, colMeans, 2000, "Indep")
indep$region  <- factor(indep$region, levels=paste("Region", 1:10))

RegionSeason <- rbind(fitres(ext.LASSO_RegionSeason, season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "LASSO - RegSeas"),
                      fitres(ext.HSt_RegionSeason,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HSt - RegSeas"),
                      fitres(ext.HS_RegionSeason,    season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HS - RegSeas"),
                      fitres(ext.FHS_RegionSeason,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "FHS - RegSeas")
                      )
RegionSeason$region  <- factor(RegionSeason$region, levels=paste("Region", 1:10))

Region <- rbind(fitres(ext.LASSO_Region, season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "LASSO - Region"),
                fitres(ext.HSt_Region,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HSt - Region"),
                fitres(ext.HS_Region,    season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HS - Region"),
                fitres(ext.FHS_Region,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "FHS - Region")
                )
Region$region  <- factor(Region$region, levels=paste("Region", 1:10))

Season <- rbind(fitres(ext.LASSO_Season, season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "LASSO - Season"),
                fitres(ext.HSt_Season,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HSt - Season"),
                fitres(ext.HS_Season,    season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "HS - Season"),
                fitres(ext.FHS_Season,   season$orig_nm, region$orig_nm, W, E, colMeans, 3000, "FHS - Season")
                )
Season$region  <- factor(Season$region, levels=paste("Region", 1:10))

comb_dfres <- rbind(indep, RegionSeason, Region, Season)

fres <- comb_dfres %>% filter(season %in% paste("15-16")) %>% data.frame

ggplot(fres, aes(x=week)) + geom_point(aes(y=ili)) + 
  facet_grid(~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("ILI.png", width=9, height=6, units="in")

ggplot(fres, aes(x=week)) + geom_point(aes(y=ili)) + 
  geom_line(aes(y=meanFit, color=model)) + geom_ribbon(aes(ymin=lowerFit, ymax=upperFit, fill=model), alpha=0.15) + 
  facet_grid(~region) + labs(x="Week", y="ILI Percentage", color="Model", fill="Model") + 
  theme_bw()
ggsave("ILIfit.png", width=9, height=6, units="in")



###### Original Data vs Demeaned Data ######
mfres <- indep %>% filter(season %in% paste("10-11") & region %in% "Region 1") %>% data.frame
mfres <- reshape2::melt(mfres[,c("ili", "week", "centered_ili")], id="week")
levels(mfres$variable) <- c("Original", "Centered")
mfres$variable <- factor(mfres$variable, levels = c("Original", "Centered"))
ggplot(mfres, aes(x=week)) + geom_point(aes(y=value)) + 
  facet_grid(~variable) + labs(x="Week", y="ILI Percentage") + 
  theme_bw()
ggsave("origvscent.png", width=7, height=4, units="in")




###### Create dataframe of betas for plotting ######

indep_betares <- beta_res(ext.Indep, N_subj, region$orig_nm, season$orig_nm, "Indep")
indep_betares$region  <- factor(indep_betares$region, levels=paste("Region", 1:10))

region_season_betares <- rbind(beta_res(ext.LASSO_RegionSeason, N_subj, region$orig_nm, season$orig_nm, "LASSO - RegSeas")
                               )
region_season_betares$region  <- factor(region_season_betares$region, levels=paste("Region", 1:10))

region_betares <- rbind(beta_res(ext.LASSO_Region, N_subj, region$orig_nm, season$orig_nm, "LASSO - Region"),
                        beta_res(ext.HSt_Region,   N_subj, region$orig_nm, season$orig_nm, "HSt - Region"),
                        beta_res(ext.HS_Region,    N_subj, region$orig_nm, season$orig_nm, "HS - Region"),
                        beta_res(ext.FHS_Region,   N_subj, region$orig_nm, season$orig_nm, "FHS - Region")
                        )
region_betares$region <- factor(region_betares$region, levels=paste("Region", 1:10))

season_betares <- rbind(beta_res(ext.LASSO_Season, N_subj, region$orig_nm, season$orig_nm, "LASSO - Season"),
                        beta_res(ext.HSt_Season,   N_subj, region$orig_nm, season$orig_nm, "HSt - Season"),
                        beta_res(ext.HS_Season,    N_subj, region$orig_nm, season$orig_nm, "HS - Season"),
                        beta_res(ext.FHS_Season,   N_subj, region$orig_nm, season$orig_nm, "FHS - Season")
                        )
season_betares$region <- factor(season_betares$region, levels=paste("Region", 1:10))


comb_betares <- rbind(indep_betares, region_season_betares, region_betares, season_betares)

# All regions and seasons
fbeta_overall <- comb_betares %>% filter(eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_overall, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_overall.png", width=9, height=9, units="in")

# Same Region different Seasons
fbeta_reg <- comb_betares %>% filter(region %in% paste("Region", c(7,9)),
                                     eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_reg, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) +
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(region~season) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_region.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_seas <- comb_betares %>% filter(season %in% paste0(15,"-",16),
                                      eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_seas, aes(x=est, color=model)) + geom_point(aes(y=eignfn)) + 
  geom_segment(aes(x=lower, xend=upper, y=eignfn, yend=eignfn)) + facet_grid(~region) +
  labs(x=TeX("$\\hat{\\beta}$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_sing_seas.png", width=9, height=6, units="in")




##### Look at posterior shrinkage with m_eff #####
n = N_obs
s2 = apply(E, 2, var)

meff <- function(mcmc, n, s2, region, season, model, hier){
  
  tmp <- data.frame(kappa=rep(NA, 3000),
                    region=rep(NA, 3000),
                    season=rep(NA, 3000))
  kappa <- matrix(NA, nrow=3000, ncol=30)
  out <- NULL
  sigma = mcmc$eps^2
  
  for(s in 1:80){
    for(e in 1:30){
      if(hier=="Region"){
        if(model=="FHS"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,region$numeric[s]]^2) * s2[e] * (mcmc$lamb_tilde[,region$numeric[s],e]^2))
        }else if(model=="LASSO"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,region$numeric[s]]^2) * s2[e] * (1))
        }else{
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,region$numeric[s]]^2) * s2[e] * (mcmc$lamb[,region$numeric[s],e]^2))
        }
      }else if(hier=="Season"){
        if(model=="FHS"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,season$numeric[s]]^2) * s2[e] * (mcmc$lamb_tilde[,season$numeric[s],e]^2))
        }else if(model=="LASSO"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,season$numeric[s]]^2) * s2[e] * (1))
        }else{
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,season$numeric[s]]^2) * s2[e] * (mcmc$lamb[,season$numeric[s],e]^2))
        }
      }else if(hier=="Region-Season"){
        if(model=="FHS"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,s]^2) * s2[e] * (mcmc$lamb_tilde[,s,e]^2))
        }else if(model=="LASSO"){
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,s]^2) * s2[e] * (1))
        }else{
          kappa[,e] = 1/(1 + n * sigma * (mcmc$tau[,s]^2) * s2[e] * (mcmc$lamb[,s,e]^2))
        }
      }
    }
    
    tmp$kappa = rowSums(1-kappa)
    tmp$region = region$orig_nm[s]
    tmp$season = season$orig_nm[s]
    out <- rbind(out, tmp)
    
  }
  
  out$model = model
  out$hier = hier
  return(out)
}

m_eff_fhs_reg      <- meff(ext.FHS_Region, n, s2, region, season, "FHS", "Region")
m_eff_fhs_seas     <- meff(ext.FHS_Season, n, s2, region, season, "FHS", "Season")
m_eff_fhs_reg_seas <- meff(ext.FHS_RegionSeason, n, s2, region, season, "FHS", "Region-Season")

m_eff_hs_reg      <- meff(ext.HS_Region, n, s2, region, season, "HS", "Region")
m_eff_hs_seas     <- meff(ext.HS_Season, n, s2, region, season, "HS", "Season")
m_eff_hs_reg_seas <- meff(ext.HS_RegionSeason, n, s2, region, season, "HS", "Region-Season")

m_eff_hst_reg      <- meff(ext.HSt_Region, n, s2, region, season, "HSt", "Region")
m_eff_hst_seas     <- meff(ext.HSt_Season, n, s2, region, season, "HSt", "Season")
m_eff_hst_reg_seas <- meff(ext.HSt_RegionSeason, n, s2, region, season, "HSt", "Region-Season")

m_eff_lasso_reg      <- meff(ext.LASSO_Region, n, s2, region, season, "LASSO", "Region")
m_eff_lasso_seas     <- meff(ext.LASSO_Season, n, s2, region, season, "LASSO", "Season")
m_eff_lasso_reg_seas <- meff(ext.LASSO_RegionSeason, n, s2, region, season, "LASSO", "Region-Season")

m_eff_comb <- rbind(m_eff_fhs_reg, m_eff_fhs_seas, m_eff_fhs_reg_seas,
                    m_eff_hs_reg, m_eff_hs_seas, m_eff_hs_reg_seas,
                    m_eff_hst_reg, m_eff_hst_seas, m_eff_hst_reg_seas,
                    m_eff_lasso_reg, m_eff_lasso_seas, m_eff_lasso_reg_seas)

ggplot(m_eff_comb) + geom_density(aes(x=kappa, color=model)) + 
  facet_grid(region~season) + theme_bw()

ggplot(m_eff_comb) + geom_boxplot(aes( x=season, y=kappa, color=model)) + 
  facet_grid(~region) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Season", y=TeX("$M_{eff}$"), color="Model")
ggsave("meff.png", width=9, height=6, units="in")

ggplot(m_eff_comb) + geom_boxplot(aes( x=season, y=kappa, color=model)) + 
  facet_grid(hier~region, scales="free_y") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x="Season", y=TeX("$M_{eff}$"), color="Model")
ggsave("meff.png", width=9, height=6, units="in")




###### Check number of beta parameters via posterior prob ######

indep_postprob <- beta_postprob(ext.Indep, N_subj, region$orig_nm, season$orig_nm, "Indep")
indep_postprob$region  <- factor(indep_postprob$region, levels=paste("Region", 1:10))

region_season_postprob <- rbind(beta_postprob(ext.LASSO_RegionSeason, N_subj, region$orig_nm, season$orig_nm, "LASSO - RegSeas")
)
region_season_postprob$region  <- factor(region_season_postprob$region, levels=paste("Region", 1:10))

region_postprob <- rbind(beta_postprob(ext.LASSO_Region, N_subj, region$orig_nm, season$orig_nm, "LASSO - Region"),
                         beta_postprob(ext.HSt_Region,   N_subj, region$orig_nm, season$orig_nm, "HSt - Region"),
                         beta_postprob(ext.HS_Region,    N_subj, region$orig_nm, season$orig_nm, "HS - Region"),
                         beta_postprob(ext.FHS_Region,   N_subj, region$orig_nm, season$orig_nm, "FHS - Region")
)
region_postprob$region <- factor(region_postprob$region, levels=paste("Region", 1:10))

season_postprob <- rbind(beta_postprob(ext.LASSO_Season, N_subj, region$orig_nm, season$orig_nm, "LASSO - Season"),
                        beta_postprob(ext.HSt_Season,    N_subj, region$orig_nm, season$orig_nm, "HSt - Season"),
                        beta_postprob(ext.HS_Season,     N_subj, region$orig_nm, season$orig_nm, "HS - Season"),
                        beta_postprob(ext.FHS_Season,    N_subj, region$orig_nm, season$orig_nm, "FHS - Season")
)
season_postprob$region <- factor(season_postprob$region, levels=paste("Region", 1:10))


comb_postprob <- rbind(indep_postprob, region_season_postprob, region_postprob, season_postprob)

# All regions and seasons
fbeta_postprob <- comb_postprob %>% filter(eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_postprob, aes(x=prob, color=model)) + geom_point(aes(y=eignfn)) + 
  facet_grid(region~season) +
  labs(x=TeX("$P(\\beta_{i} \\notin (-1,1))$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_postprob.png", width=9, height=9, units="in")


# Same Region different Seasons
fbeta_postprob <- comb_postprob %>% filter(region %in% paste("Region", c(7,9)),
                                           eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_postprob, aes(x=prob, color=model)) + geom_point(aes(y=eignfn)) +
  facet_grid(region~season) +
  labs(x=TeX("$P(\\beta_{i} \\notin (-1,1))$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_postprob_sing_region.png", width=9, height=6, units="in")

# Same Season, different Regions
fbeta_postprob <- comb_postprob %>% filter(season %in% paste0(10,"-",11),
                                           eignfn %in% paste(1:10)) %>% data.frame
ggplot(fbeta_postprob, aes(x=prob, color=model)) + geom_point(aes(y=eignfn)) + 
  facet_grid(season~region) +
  labs(x=TeX("$P(\\beta_{i} \\notin (-1,1))$"), y=TeX("$\\phi (x)$"), color="Model") + 
  theme_bw() + coord_flip() + scale_y_continuous(breaks=seq(1,10,1))
ggsave("beta_postprob_sing_seas.png", width=9, height=6, units="in")



specific_counts = comb_postprob %>% 
  group_by(season, region, model) %>%
  summarize(nshrunk = sum(prob < .1),
            prop_shrunk = nshrunk/n())

model_counts = comb_postprob %>% 
  group_by(model) %>%
  summarize(nshrunk = sum(prob < .1),
            prop_shrunk = nshrunk/n())

names(model_counts) <- c("Model", "Num. Shrunk", "\\% Shrunk")
print(
  xtable(t(model_counts), caption="This table shows the number of beta parameters whose posterior lie mainly within -1 and 1.",
         digits=c(0,rep(3, length(unique(model_counts$Model)))), label = "tab:mse"
  ), 
  include.colnames=FALSE,sanitize.text.function=function(x){x}, file = "prob_shrunk_tab.tex")








###### Correlation plot ######
fbeta_corr <- comb_betares %>% filter(eignfn %in% paste(1:6)) %>% select(-c(lower, upper)) %>% data.frame
melted_corr <- spread(fbeta_corr, eignfn, est)
pairs(melted_corr[,4:9], col=melted_corr$model)


makePairs <- function(data){
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
               x = data[, xcol], y = data[, ycol], data)
  }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
  }))
  list(all=all, densities=densities)
}

# expand beta data frame for pairs plot
gg1 = makePairs(melted_corr[,4:9])

# new data frame mega iris
mega_beta = data.frame(gg1$all, Model=rep(melted_corr$model, length=nrow(gg1$all)))

# pairs plot
ggplot(mega_beta, aes_string(x = "x", y = "y")) + 
  facet_grid(xvar ~ yvar, scales = "free") + 
  geom_point(aes(colour=Model), na.rm = TRUE, alpha=0.8) + 
  stat_density(aes(x = x, y = ..scaled.. * diff(range(x)) + min(x)), 
               data = gg1$densities, position = "identity", 
               colour = "grey20", geom = "line") +
  labs(x="", y="") + theme_bw()

ggsave("beta_correlation.png", width=9, height=6, units="in")


