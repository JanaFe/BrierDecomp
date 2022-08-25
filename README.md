## 3-component Decomposition of Brier Score
### Numpy-implementation for 3-component decomposition of Brier Score into Reliability, Resolution and Uncertainty, and their estimated standard deviations

This code was derived from the R-package 'Specs Verification - Forecast Verification Routines for Ensemble Forecasts of Weather and Climate' (Function: BrierDecomp) 
by Stefan Siegert, Jonas Bhend, Igor Kroener and Matteo De Felice.
Dept. Mathematics and Statistics, University of Exeter, UK
- CRAN Repository: https://rdrr.io/cran/SpecsVerification/src/R/BrierDecomp.R (Date/Publication: 2020-02-26 15:40:06 UTC)
- Github: https://github.com/sieste/SpecsVerification

The function BrierDecomp() has the following input and output.

INPUT:
- p: 1d array of predicted probabilities
- y: 1d array of binary observations (0 or 1)
- bins: to estimate the calibration function, default: 10
- bias_corrected: logical, default=False, if false, the standard (biased) decomposition of Murphy (1973) is used. If true, the bias-corrected decomposition of Ferro (2012) (See References in Readme)

OUTPUT:
- dictionary with estimates of Reliability, Resolution and Uncertainty
- dictionary with corresponding std. deviations

References:
- Murphy (1973): A New Vector Partition of the Probability Score. J. Appl. Met. https://doi.org/10.1175/1520-0450(1973)012<0595:ANVPOT>2.0.CO;2 
- Ferro and Fricker (2012): A bias-corrected decomposition of the Brier score. QJRMS. https://rmets.onlinelibrary.wiley.com/doi/10.1002/qj.1924 
- Siegert (2013): Variance estimation for Brier Score decomposition. QJRMS. https://arxiv.org/abs/1303.6182 
