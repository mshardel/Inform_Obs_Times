# Inform_Obs_Times
R code for simulation study for "Time-Varying Latent Effect Models for Repeated Measurements to Address Informative Observation Times in the U.S. Medicare Minimum Data Set."

Simulations for:
1. Semi-parametric longitudinal model with time-varying latent effects and time-varying covariates.
2. Inverse conditional intensity rate ratio (ICIRR)-weighted semi-parametric longitudinal model with time-varying latent effects and time-varying covariates.

sims_generate_proposed_Lin-Ying.R: simulates data, analyzes simulated data using proposed and Lin-Ying method, and writes CSV files of estimates to sims_results.

sims_tables_figures_proposed_Lin-Ying.R: creates Figure 1, Figure 2, Figure S.1, Figure S.2, Table 1, and Table S.1, and saves them to sims_results/sims_postprocess.

sims_generate_ICIRR-weighted-proposed_Buzkova-Lumley.R: simulates data, analyzes simulated data using proposed ICIRR-weighted and Buzkova-Lumley method, and writes CSV files of estimates to wt_sims_results.

sims_tables_ICIRR-weighted-proposed_Buzkova-Lumley.R: creates Table S.2 and Table S.3, and saves them to wt_sims_results/wt_sims_postprocess.



