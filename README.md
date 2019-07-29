# trends_harmonic_patterns
Extraction of periodic patterns (e.g. seasonality) is often useful in time series analysi. This repository contains time series functionalities for trends analysis and pattern extraction using harmonic regression.

The goal is to provide easy functions to operate on univariate time series as well as time series raster stack (e.g. satellite image of vegetation).

Several functions are provided to extract generate amplitudes and phases for any frequency given a time series. The default is two harmonics (annual and bi-annual if 12 months periods are considered) but the package also provide ability to extract any harmonic.

In addition, short-term Fourier transform (also known as Window Fourier Transform) is provided using any size window desired. The current implementation allows for missing values without the requirement for filling of NA.
