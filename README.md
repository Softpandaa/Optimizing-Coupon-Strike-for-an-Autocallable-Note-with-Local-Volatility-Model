# Optimizing-Coupon-Strike-for-an-Autocallable-Note-with-Local-Volatility-Model
This repository contains the methodology and implementation for optimizing coupon strike parameters in autocallable structured notes using a local volatility framework. The model prices a 2-year autocallable note linked to the Nikkei 225, S&P 500, and HSI indices, with features including knock-in (50%) and knock-out (110%) barriers.

# Key Features
- **Risk-free Rates Interpolation**: Generate continuous OIS rate using cubic spline interpolation
- **Local Volatility Modelling**: Implements Dupire's equation with linear interpolation to fit volatility surfaces
- **Multi-Asset Monte Carlo**: Simulates correlated paths for three indices using Cholesky decomposition
- **Bisection Optimization**: Solves coupon strike that prices the note at 98% of the issue price

# Note Structure
- **Trade Date**: 11 December 2023
- **Maturity**: 2 years
- **Underlying Indices**: Nikkei 225 (NKY), S&P 500 (SPX), Hang Seng Index (HSI)
- **Issue Price**: 100% of $10,000 denomination

1. **Variable Interest** (Paid semi-annually)
   - 2% p.a. if reference price ≥ coupon strike
   - 0.01% p.a. otherwise
   - Reference price: Closing price of laggard index (lowest S<sub>t</sub>/S<sub>0</sub>)

2. **Knock-Out Event**
   - Triggered if laggard index ≥ 110% of initial spot
   - Pays 100% principal + accrued interest
   - Note immediately terminates

3. **Knock-In Event**
   - Triggered if laggard index ≤ 50% of initial spot at any time
   - Affects final redemption

4. **Final Redemption** (At maturity if not knocked out)
   - No knock-in: 100% principal
   - Knock-in occurred: Principal × min⁡(100%, S<sub>T</sub>/S<sub>0</sub> of laggard)

# Code Structure
1. main.ipynb: The main script to run pricing examples and configure input settings for the autocallable note
2. VolSurfaceFunctions.py:
   - Process the data
   - Define the local volatility model
   - Fit the volatility surface
