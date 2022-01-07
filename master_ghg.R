# Master file for NHC ghg analyses
# A Carter
# 2020-02-14
setwd('C:/Users/Alice Carter/git/ghg_patterns_nhc/')

# Compile ghg dataframe and all drivers
  # this script also fills in missing driver points using other data
  # and by discharge interpolation, rating curves for depth, K600 and by
  # allowing instantaneous data and daily averages to fill when appropriate.
  source('src/compile_dataframe_gasfluxdata.R')

#plot spatial and temporal ghg patterns
  source('src/ghg_spatial_temporal_plots.R')

# Build linear mixed effects models for gas conc and flux
  # attempt at PCA's
  # Plots of all NHC UNHC water chem data, pca for discharge
  source('src/ghg_analysis.R')

# Mass balance in each reach along the continuum
  source('src/ghg_mass_balance.R')


g(fit_nreg_sd7_ss05)Visualizeplot_metab_preds(predict_metab(fit_nreg_sd7_ss05))plot_metab_preds(predict_metab(fit_nreg_fixedK))K600 vs ERKvER <- get_fit(fit_raymond_ss24)KvER <- get_fit(fit_raymond_ss05)KvER <- get_fit(fit_nreg_ss24)KvER <- get_fit(fit_nreg_sd7_ss05)KvER <- get_fit(fit_nreg_fixedK)plot(KvER$daily$K600_daily_mean, KvER$daily$ER_daily_mean)cor(KvER$daily$K600_daily_mean,KvER$daily$ER_daily_mean,use = "na.or.complete")Write Fileswritefiles <- function(mod){data <- get_fit(mod)for (i in seq_along(data)) {filename = paste(names(data)[i], ".csv")write.csv(data[[i]], filename)  }write.csv(unlist(get_specs(mod)),"specs.csv")write.csv(get_data_daily(mod), "datadaily.csv")write.csv(get_data(mod),"mod_and_obs_DO.csv")}Create new folder for site and write csv infoReset working directory!!getwd()writefiles(fit_207_lim)Fixed K model ####the only place that Q enters into the model is to predict the K-Q relationshipfor this reason, we can use a fake Q as nodes and with our data to simplifyusing a fixed K/Q relationship as a prior. In this case, we are binning thedata and then replacing everything within each bin with a node, assigned tothe corresponding K derived from the Hall 1972 paperprep_fake_Q <- function(dat, bayes_specs_fixedK){dat <- dat %>%mutate(date = as.Date(solar.time))daily <- dat %>%group_by(date) %>%summarize(discharge = mean(discharge, na.rm = T)) %>%mutate(logQ = log(discharge))nbins <-  length(bayes_specs_fixedK$K600_lnQ_nodes_centers)Qbreaks <- c(bayes_specs_fixedK$K600_lnQ_nodes_centers[1] -diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2,bayes_specs_fixedK$K600_lnQ_nodes_centers +diff(bayes_specs_fixedK$K600_lnQ_nodes_centers)[1]/2)rQ <- range(daily$logQ, na.rm = T)if(Qbreaks[1] > rQ[1]){Qbreaks[1] <- rQ[1]}if(Qbreaks[nbins] < rQ[2]){Qbreaks[nbins] <- rQ[2]}daily <- daily %>%mutate(Q = cut(daily$logQ, Qbreaks, labels = 1:nbins)) %>%select(date, Q)daily$Q <- as.numeric(daily$Q)dat <- left_join(dat, daily, by = "date") %>%select(-discharge) %>%mutate(discharge = exp(Q)) %>%select(-Q,-date)return(dat)}
