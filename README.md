
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Zooplankton Model of Size Spectrum (ZooMSS)

<!-- # planktonr <a href='https://github.com/MathMarEcol/zoomss'><img src='man/figures/planktonr.png' align="right" width="139px" /></a> -->

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Windows](https://github.com/MathMarEcol/zoomss/actions/workflows/Windows.yaml/badge.svg)](https://github.com/MathMarEcol/zoomss/actions/workflows/Windows.yaml)
[![Linux](https://github.com/MathMarEcol/zoomss/actions/workflows/Linux.yaml/badge.svg)](https://github.com/MathMarEcol/zoomss/actions/workflows/Linux.yaml)
[![MacOS](https://github.com/MathMarEcol/zoomss/actions/workflows/MacOS.yaml/badge.svg)](https://github.com/MathMarEcol/zoomss/actions/workflows/MacOS.yaml)
[![issues -
zoomss](https://img.shields.io/github/issues/MathMarEcol/zoomss)](https://github.com/MathMarEcol/zoomss/issues)
[![Codecov test
coverage](https://codecov.io/gh/MathMarEcol/zoomss/graph/badge.svg)](https://app.codecov.io/gh/MathMarEcol/zoomss)
<!-- badges: end -->

## Overview of ZooMSS

The Zooplankton Model of Size Spectra (ZooMSS) is a functional
size-spectrum model of the marine ecosystem (following Heneghan et
al. 2016) to resolve phytoplankton, nine zooplankton functional groups
(heterotrophic flagellates and ciliates, omnivorous and carnivorous
copepods, larvaceans, euphausiids, salps, chaetognaths and jellyfish)
and three size-based fish groups. Zooplankton functional groups are
resolved using their body-size ranges, size-based feeding
characteristics and carbon content, and the zooplankton community
emerges from the model across global environmental gradients, depending
on the functional traits of the different groups.

We developed the Zooplankton Model of Size Spectra (ZooMSSv2) based on
the prototype of Heneghan et al. (2016). ZooMSS uses the functional
size-spectrum framework (Blanchard et al., 2017) to resolve the body
size ranges, size-based feeding characteristics and carbon content of
nine zooplankton groups and three fish groups. The model supports
time-varying environmental conditions enabling studies of seasonal
cycles, climate change scenarios, and ecosystem responses to
environmental variability.

ZooMSS represents the marine ecosystem as three communities:
phytoplankton, zooplankton and fish. The zooplankton community consists
of nine of the most abundant zooplankton groups, and the fish community
was made up of a small, medium and large group. Dynamics of the
phytoplankton are not explicitly resolved in the model, rather the mean
size structure of the phytoplankton community is estimated directly from
satellite chlorophyll a observations (Brewin et al., 2010; Barnes et
al., 2011; Hirata et al., 2011). Abundances of the zooplankton and fish
communities are driven by size-dependent processes of growth and
mortality, with the temporal dynamics of each functional group governed
by separate second-order McKendrick-von Foerster equations.

## Installation

You can install the development version of zoomss from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("MathMarEcol/zoomss")
```

## Publications

1.  Heneghan, R.F., Everett, J.D., Blanchard, J.L., Richardson,
    A.J., 2016. Zooplankton Are Not Fish: Improving Zooplankton Realism
    in Size-Spectrum Models Mediates Energy Transfer in Food Webs.
    Front. Mar. Sci. 3, 1–15. <https://doi.org/10.3389/fmars.2016.00201>

2.  Heneghan, R.F., Everett, J.D., Sykes, P., Batten, S.D., Edwards, M.,
    Takahashi, K., Suthers, I.M., Blanchard, J.L., Richardson, A.J., in
    review, A global size-spectrum model of the marine ecosystem that
    resolves zooplankton composition. Ecological Modelling

## Getting Help

If you encounter problems running the model, raise an issue on GitHub:
<https://github.com/MathMarEcol/ZoopSizeSpectraModel/issues>

If you find errors or want to improve the model, we’d love you to make
the changes and submit a pull request for us to review and approve.
