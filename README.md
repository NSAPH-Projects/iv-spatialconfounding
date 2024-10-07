# An Instrumental Variables Framework to Unite Spatial Confounding Methods

This repository contains the code for the simulation and data application sections of the paper. The workflow is:

- `funcs.R` contains all utility functions for the simulation and data applicaton. The main function is 'ctseff' adapted from the [npcausal](https://github.com/ehkennedy/npcausal/tree/master) package with minor modifications.
- `simulation` contains the code for the simulation section of the paper.
- `data_application` contains the code for the data application section of the paper.

### Acknowledgements
The computations in this paper were run on the FASRC cluster supported by the FAS Division of Science Research Computing Group at Harvard University. Medicare mortality data are stored at a Level-3 secured data platform on Research Computing Environment, supported by the Institute for Quantitative Social Science in the Faculty of Arts and Sciences at Harvard University. 

The exposure, confounders, and mortality data used in the data application were adapted from:
[X. Wu, D. Braun, J. Schwartz, M. A. Kioumourtzoglou, F. Dominici, Evaluating the impact of long-term exposure to fine particulate matter on mortality among the elderly. Sci. Adv. 6, eaba5692 (2020).](https://www.science.org/doi/full/10.1126/sciadv.aba5692). We would like to thank the authors for sharing the data.

### Contact us
If you have any questions, please contact us at [swoodward@g.harvard.edu](mailto:swoodward@g.harvard.edu) or [fdominic@hsph.harvard.edu](mailto:fdominic@hsph.harvard.edu).