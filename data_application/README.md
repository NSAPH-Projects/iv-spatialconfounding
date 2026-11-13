# Data application: PM2.5 and all-cause mortality

The exposure is the 2001-2010 average of fine particulate matter (PM2.5) at the
ZIP-code level, and the outcome is the all-cause mortality rate among Medicare
enrollees aged 65+ over 2011-2016 (n = 33,254 ZIP codes; adapted from Wu et
al., 2020). The Medicare data are not publicly
available; the analysis was run in Harvard University's Regulated Data
Environment (ReD).

Run from `data_application/`:

1. `0_preprocessing.R`: merges exposure, confounders and outcomes, builds the
   thin plate spline and graph Laplacian bases (`../R/basis_construction.R`),
   and saves them to `preprocessed/`. Also writes `images/PM_ivs_20pct_2010.png`
   (exposure and estimated instrumental / non-instrumental variation).
2. `1_run_analysis.R`: baseline, spatial coordinates, IV-TPS and IV-GL estimates
   of the truncated exposure effect for cutoffs 6-12 µg/m³, with the sequential
   stability rule (`../R/basis_selection.R`). Writes CSVs to `results/` and
   `images/selected_estimates_across_cutoffs_*.png`.
3. `2_maps_layout.R`: rearranges `images/PM_ivs_20pct_2010.png` into the
   two-row layout used in the paper (`images/PM_ivs_20pct_2010_2x3.png`).

`tables/` holds the exported Hausdorff distances between each method's confidence
interval and the oracle's (Supplementary Section 5.2), copied from
`results/full_sample_basis_selection/`.
