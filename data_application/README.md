# Data Application

This folder contains the code for the data application section of the paper. 

- `PM_mortality.R` contains the code to load the data, preprocess the data, estimate the truncated exposure effects using seven methods, and create plots of the estimates.
- The plots are saved in the `images` folder.

### Sources of Data
The exposure, confounders, and mortality data used in the data application were adapted from:
[X. Wu, D. Braun, J. Schwartz, M. A. Kioumourtzoglou, F. Dominici, Evaluating the impact of long-term exposure to fine particulate matter on mortality among the elderly. Sci. Adv. 6, eaba5692 (2020).](https://www.science.org/doi/full/10.1126/sciadv.aba5692). 

- Medicare claims data were obtained from the Centers for Medicare and Medicaid services. 
- Exposure data was obtained from a [well-validated ensemble-based prediction model](https://www.sciencedirect.com/science/article/pii/S0160412019300650).
- Confounder data was obtained from the U.S. Census, American Community Survey, Behavioral Risk Factor Surveillance System, and Gridmet via Google Earth Engine.

This data is not publicly available.