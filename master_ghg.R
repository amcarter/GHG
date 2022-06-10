# Master file for NHC ghg analyses
# A Carter
# 2020-02-14
setwd('C:/Users/alice.carter/git/ghg_patterns_nhc/')

# Compile ghg dataframe and all drivers
  # this script also fills in missing driver points using other data
  # and by discharge interpolation, rating curves for depth, K600 and by
  # allowing instantaneous data and daily averages to fill when appropriate.
  source('src/compile_dataframe_gasfluxdata.R')

#plot spatial and temporal ghg patterns
  source('src/plots/ghg_spatial_temporal_plots.R')

# Build linear mixed effects models for gas conc and flux
  # attempt at PCA's
  # Plots of all NHC UNHC water chem data, pca for discharge
  source('src/GHG_linear_models.R')

# Mass balance in each reach along the continuum
  source('src/ghg_mass_balance.R')


