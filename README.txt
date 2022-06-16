READ ME for GHG repository

Created by: AM Carter
Feb 2022

Summary:
Code for processing and analyzing New Hope Creek greenhouse gas data.


DATA:
1. NHCsite_metadata.csv: Contains full site names and sitecodes, coordinates, geomorphic data, and calculated from the natrional hydrography dataset: watershed area, percent impervious surface cover (isc.percent), NHD comid, sltream bed slope 

2. NHC_2019-2020_processed_GHGdata.csv: Contains sitecodes, dates and times in EST, water chemistry, and dissolved greenhouse gas concentrations for each site x sample time in the dataset. error_pct is given in cases where sample volume had to be estimated.

3. all_nhc_ysi_data.csv: Contains sitecodes, dates and times in UTC, and water chemistry data collected with a ysi probe, and water stages.

4. NOAA_airpres.csv: 15 minute air pressure data pulled from the National Oceanic and Atmospheric Administration database from the location nearest to New Hope Creek.

5. metabolism.csv: Raw 15 minute data from stream sensors that was used to estimate metabolism as well as calculated discharge values.

6. met_preds_stream_metabolizer_O2.rds: Daily metabolism and gas exchange estimates from all study sites for the duration of the study estimated using streamMetabolizer. See Carter 2021 for more information. Metabolism estimates are in gC/m2/d

7. water_chemistry_2019-2020_compiled.csv: grab sample ion data for sample time points in dataset.

8. StreampulseWQDec2020.csv: grab sample ion data for the NHC and UNHC sites

SCRIPTS:
1. compile_dataframe_gasfluxdata.R
	- Combines greenhouse gas concentrations with water chemistry, metabolism estimates, discharge measurements, air pressure data, and site metadata from files 1-8. 
	- For water chemistry parameters where there are both instantaneous ysi samples as well as in stream 15 minute sensor data, the ysi probe measurements were kept unless they were missing, in which case they were filled in with the matching 15 minute sensor data point.
	- Samples were placed in sample groups based on the sample dates. In most cases, all samples in a group were taken on the same day. In two cases, there were some samples taken one day after the first samples - these were grouped together for longitudinal analysis.
	- Gas flux rates are calculated based on the estimated gas exchange rates and the measured dissolved gas concentrations using the marecalc package.

2. ghg_spatial_temporal_plots.R
	- Takes the compiled dataframes generated in script 1 and makes plots

3. Spatial_Alice.R
	- uses a DEM file to extract stream slopes for each study site using the whitebox package.

4. GHG_linear_models.R
	- scales and centers data for model building
	- contains functions for model selection for both linear models and linear mixed effects models.
	- contains functions for calculating model fits and information criterions.

5. plot_excess_gas_ratios.R
	- calculates and plots ratios of gas departures from equillibrium

6. calc_numbers_for_results.R
	- various summary functions to generate the numbers reported in the results section

7. plot_CO2_NEP_fraction.R
	- calculates the fraction of CO2 flux attributable to NEP
	- plots it for various RQs

8. Plot_correlation.R
	- plots correlations and significance for model predictor variables

