# Optimizing-Coupon-Strike-for-an-Autocallable-Note-with-Local-Volatility-Model
This repository contains the methodology and implementation for optimizing coupon strike parameters in autocallable structured notes using a local volatility framework. The model prices a 2-year autocallable note linked to the Nikkei 225, S&P 500, and HSI indices, with features including knock-in (50%) and knock-out (110%) barriers.

# Key Features
- **Risk-free Rates Interpolation**: Generate continuous OIS rate using cubic spline interpolation
- **Local Volatility Modelling**: Implements Dupire's equation with linear interpolation to fit volatility surfaces
- **Multi-Asset Monte Carlo**: Simulates correlated paths for three indices using Cholesky decomposition
- **Bisection Optimization**: Solves coupon strike that prices the note at 98% of the issue price
