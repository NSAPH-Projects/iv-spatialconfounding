# An Instrumental Variables-Like Framework to Unite Spatial Confounding Methods

Code for the simulation study and data application of the paper.

```
R/                  shared functions
  basis_construction.R   thin plate spline and graph Laplacian spatial bases
  basis_selection.R      sequential stability rule for the number of basis functions
  truncated_effect.R     doubly robust estimator of the truncated exposure effect
                         (adapted from npcausal::ctseff, Kennedy et al. 2017)
simulation/         simulation study (Section 5 and Supplementary Section 4)
data_application/   PM2.5 and mortality application (Section 6)
```

See `simulation/README.md` for how to reproduce the simulation results, tables
and figures, and `data_application/README.md` for the data application. The
Medicare data used in the data application are not publicly available; the
analysis was run in Harvard University's Regulated Data Environment (ReD).

## R packages

sf, dplyr, ggplot2, mgcv, RSpectra, Matrix, SuperLearner (with gam and earth),
KernSmooth, RANN; the data application also uses readr, stringr, xtable, here
and gridExtra.

## Citation

```bibtex
@article{woodward2024instrumental,
  title   = {An Instrumental Variables Framework to Unite Spatial Confounding Methods},
  author  = {Woodward, Sophie M. and Tec, Mauricio and Dominici, Francesca},
  journal = {arXiv preprint arXiv:2411.10381},
  year    = {2024},
  url     = {https://arxiv.org/abs/2411.10381}
}
```

### Acknowledgements
The computations in this paper were run on the FASRC Cannon cluster, supported by the FAS Division of Science Research Computing Group at Harvard University, as well as the Regulated Data (ReD) Environment, operated by Harvard University Information Technology and supported by the Office of the Vice Provost for Research at Harvard University.

The exposure, confounders, and mortality data used in the data application were adapted from:
[X. Wu, D. Braun, J. Schwartz, M. A. Kioumourtzoglou, F. Dominici, Evaluating the impact of long-term exposure to fine particulate matter on mortality among the elderly. Sci. Adv. 6, eaba5692 (2020).](https://www.science.org/doi/full/10.1126/sciadv.aba5692). We would like to thank the authors for sharing the data.

### Contact us
If you have any questions, please contact us at [swoodward@g.harvard.edu](mailto:swoodward@g.harvard.edu) or [fdominic@hsph.harvard.edu](mailto:fdominic@hsph.harvard.edu).
