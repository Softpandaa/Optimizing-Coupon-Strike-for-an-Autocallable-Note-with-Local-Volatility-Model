import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline, make_interp_spline
from functools import partial
from scipy.misc import derivative

import warnings
warnings.filterwarnings("ignore")

def get_continuous_r(r, t):
    """
    This function turns daily compounded r into continuously compounded r
    r: daily compounded r
    t: number of days compounded
    
    Return: continuously compounded r
    """
    return t*np.log(1+r/t)

def get_r_func(file, t):
    """
    This function returns a function that get r for each T
    file: file name containing all the daily compounded r
    t: number of days compounded
    
    Return: a function that get r for each T
    """
    temp = pd.read_excel(file, index_col=0)
    r = np.array([get_continuous_r(_, t)/100 for _ in temp.iloc[0]])
    T = temp.columns
    return CubicSpline(T, r) 

def get_ATMF(chain):
    """
    This function estimates the values of C_ATMF and K_ATMF based on the given chain, the methodology can be found in the word doc
    chain: raw option chain
    
    Return: C_ATMF, K_ATMF
    """
    C = chain["Call"]
    P = chain["Put"]
    K = chain["Strike"]
    n = len(C)
    for i in range(1, n):
        if C[i] >= P[i] and C[i+1] <= P[i+1]:
            C_K1, P_K1, K1 = C[i], P[i], K[i]
            C_K2, P_K2, K2 = C[i+1], P[i+1], K[i+1]
            break
    C_ATMF = (C_K1+P_K1)/2
    K_ATMF = (C_ATMF-C_K1)/(C_K2-C_K1) * (K2-K1) + K1
    return C_ATMF, K_ATMF

def get_ATM_Call(chain, S):
    """
    This function estimates the values of C_ATM based on the given chain, the methodology can be found in the word doc
    chain: raw option chain
    
    Return: C_ATM
    """
    C = chain["Call"]
    K = chain["Strike"]
    n = len(C)
    for i in range(1, n):
        if S >= K[i] and S <= K[i+1]:
            C_K1, K1 = C[i], K[i]
            C_K2, K2 = C[i+1], K[i+1]
            break
    C_ATM = (C_K2-C_K1)/(K2-K1) * (S-K1) + C_K1
    return C_ATM

def get_implied_div(S, K, r, T):
    """
    This functions estimated the continuous dividend yield based on a forward price
    S: Spot price of underlying
    K: Forward price of underlying
    r: risk-free rate
    T: Year to maturity
    
    Return: implied dividend yield
    """
    q = r - np.log(K/S) / T
    return max(q, 0) # Prevent implied dividend yield < 0, which doesn't make sense

def BS_Euro_Pricer_F(F, K, r, T, sig):
    """
    This function computes the European option price based on a forward price
    F: forward price at T
    K: Strike
    r: risk-free rate (in decimal)
    T: Year to maturity
    sig: Volatility of the underlying (in decimal)
    
    Return: European option price
    """
    d1 = (np.log(F/K) + ((sig**2)/2)*T) / (sig*np.sqrt(T))
    d2 = d1 - sig * np.sqrt(T)
    return np.exp(-r*T) * (F * norm.cdf(d1) - K * norm.cdf(d2))

def BS_Euro_Pricer(isCall, S, K, r, q, T, sig):
    """
    This function computes the European option price
    isCall: True for Call; False for Put
    S: Spot price of underlying
    K: Strike
    r: risk-free rate (in decimal)
    q: continuous dividend yield (in decimal)
    T: Year to maturity
    sig: Volatility of the underlying (in decimal)
    
    Return: European option price
    """
    d1 = (np.log(S/K) + (r - q + (sig**2) / 2) * T) / (sig * np.sqrt(T))
    d2 = d1 - sig * np.sqrt(T)
    if isCall:
        return S * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    else:
        return K * np.exp(-r * T) * norm.cdf(-d2) - S * np.exp(-q * T) * norm.cdf(-d1)

def IV_Bisect(P, tol, isCall, S, K, r, q, T):
    """
    This function estimates the implied volatility of the given option
    P: Observed price of the option
    tol: Maximum error of approximation
    isCall: True for Call; False for Put
    S: Spot price of underlying
    K: Strike
    r: risk-free rate (in decimal)
    q: continuous dividend yield (in decimal)
    T: Year to maturity
    
    Return: implied volatility (in sigma)
    """
    sig_left = 0.0001
    sig_right = 0.0001
    while BS_Euro_Pricer(isCall, S, K, r, q, T, sig_right) - P < 0:
        sig_right += 0.2
    sig_temp = (sig_left + sig_right) / 2
    while abs(BS_Euro_Pricer(isCall, S, K, r, q, T, sig_temp) - P) > tol:
        if BS_Euro_Pricer(isCall, S, K, r, q, T, sig_temp) - P > 0:
            sig_right = sig_temp
        else:
            sig_left = sig_temp
        sig_temp = (sig_left + sig_right) / 2
    return sig_temp

def get_dC_dT(F, K, r, T, sig):
    """
    Compute dC/dT
    F: forward price at T
    K: Strike
    r: risk-free rate
    T: Year to maturity
    sig: implied volatility
    
    Return: dC/dT
    """
    f = partial(BS_Euro_Pricer_F, F, K, r, sig=sig)
    return derivative(f, T, dx=1/252) # dT is 1 day

def get_d2C_dK2(F, K, r, T, sig, step):
    """
    Compute d2C/dK2
    F: forward price at T
    K: Strike
    r: risk-free rate
    T: Year to maturity
    sig: implied volatility
    step: difference between 2 strike
    
    Return: d2C/dK2
    """
    f = partial(BS_Euro_Pricer_F, F, r=r, T=T, sig=sig)
    return derivative(f, K, dx=step*.01, n=2) # dK is 1% of strike price step size

def get_local_vol_Dupire(F, K, r, T, sig, step):
    """
    Compute local vol according to Dupire Equation
    F: forward price at T
    K: Strike
    r: risk-free rate
    T: Year to maturity
    sig: implied volatility
    step: difference between 2 strike
    
    Return: local vol (in sigma^2)
    """
    dC_dT = get_dC_dT(F, K, r, T, sig)
    d2C_dK2 = get_d2C_dK2(F, K, r, T, sig, step)
    return dC_dT / ((K**2)*d2C_dK2/2)

def chain_processor(chain, S, r, q, T, K_ATMF, step):
    """
    This function process the raw option chain and return a processed chain for further analysis
    chain: unprocessed option chain (From BBG OMON)
    S: Spot price of underlying
    r: risk-free rate (in decimal)
    q: continuous dividend yield (in decimal)
    T: Year to maturity
    K_ATMF: Strike of ATM forward (forward price)
    step: difference between 2 strike
    
    Return: processed chain containing Moneyness, IV, Call price etc.
    """
    Processed_Chain = pd.DataFrame()
    C = chain["Call"]
    P = chain["Put"]
    K = chain["Strike"]
    n = len(C)
    for i in range(1, n+1):
        try:
            x = np.log(K[i] / K_ATMF) # Forward Moneyness
            X = np.log(K[i] / S) # Spot Moneyness
            if x > 0: # Use OTM Call/ Put to calculate implied volatility
                if C[i] == 0:
                    continue
                else:
                    IV = IV_Bisect(C[i], .0001, True, S, K[i], r, q, T)
            else:
                if P[i] == 0:
                    continue
                else:
                    IV = IV_Bisect(P[i], .0001, False, S, K[i], r, q, T)
            LV = np.sqrt(get_local_vol_Dupire(K_ATMF, K[i], r, T, IV, step))
            if np.isnan(LV) or np.isinf(LV):
                continue 
            else:
                processed_chain = pd.DataFrame({"Spot Moneyness": X,
                                                "IV": IV,
                                               "Call": C[i],
                                               "Put": P[i],
                                               "Strike": K[i],
                                               "LV": LV},
                                              index=(x, ))
                Processed_Chain = pd.concat([Processed_Chain, processed_chain])
        except:
            continue
    Processed_Chain.index.name = "Forward Moneyness"
    return Processed_Chain


def get_processed_data(file, S, r_func, step):
    """
    This function returns a list containing the processed option chains of different maturity, and a table containing rf, q etc.
    file: Excel file name containing all the raw option chains
    S: Spot price of the underlying
    r_func: a function that gets r for each t
    step: difference between 2 strike
    
    Return: list of processed chains, info table
    """
    processed_chains = []
    processed_tables = pd.DataFrame()
    xl = pd.ExcelFile(file)
    for sheet_name in xl.sheet_names:
        raw_data = pd.read_excel(file, sheet_name=sheet_name)
        raw_chain = raw_data.loc[1:].copy() # Get rid of the info row
        raw_chain.rename(columns={"Last":"Call", "Last.1":"Put"}, inplace=True)
        chain_des = raw_data.loc[0][0].split() # First row of each sheet contains the info 
        T = int(chain_des[1][1:-3])/365
        r = r_func(T)
        C_ATMF, K_ATMF = get_ATMF(raw_chain)
        implied_div = get_implied_div(S, K_ATMF, r, T)
        sig_ATMF = IV_Bisect(C_ATMF, .00001, True, S, K_ATMF, r, implied_div, T)
        LV_ATMF = np.sqrt(get_local_vol_Dupire(K_ATMF, K_ATMF, r, T, sig_ATMF, step))
        processed_chains.append(chain_processor(raw_chain, S, r, implied_div, T, K_ATMF, step))
        processed_table = pd.DataFrame({"rf": r,
                                       "q": implied_div,
                                       "C_ATMF": C_ATMF,
                                       "K_ATMF": K_ATMF,
                                       "sig_ATMF": sig_ATMF,
                                       "M_ATMF": np.log(K_ATMF/S),
                                        "LV_ATMF": LV_ATMF},
                                      index=(T,))
        processed_tables = pd.concat([processed_tables, processed_table])
        print(f"T={T} processed")
    processed_tables.index.name = "Time to maturity"
    return processed_chains, processed_tables

def vol_sur_M(M, sig_ATMF, params):
    """
    This function estimates the implied volatility given the moneyness of the option at time t
    M: Moneyness, ln(K/S)
    sig_ATMF: implied volatility (in sigma^2) of options at M = ln(K/K_ATMF) = 0
    params: [delta, gamma, kappa]
    
    Return: implied volatility (in sigma^2)
    """
    return sig_ATMF + params[0]*np.tanh((M)*params[2])/params[2] + params[1]/2*(np.tanh((M)*params[2])/params[2])**2

def get_RSS_func(M, Vol, sig_ATMF):
    """
    Get the squared distance between expected implied vol and actual implied vol (RSS)
    M: Moneyness of the observed vol
    Vol: Observed implied vol (in sigma^2)
    sig_ATMF: implied volatility (in sigma^2) of options at M = ln(K/K_ATMF) = 0
    
    Return: a distance function based on the params
    """
    def RSS_func(params):
        error = vol_sur_M(M, sig_ATMF, params) - Vol
        return error.dot(error)
    
    return RSS_func


def fit_vol_sur_M(M, Vol, sig_ATMF, initial_guess=[-1, .25, .25]):
    """
    Return the vol surface function which the error is minimized
    M: Moneyness of the observed vol 
    Vol: Observed implied vol (in sigma^2)
    sig_ATMF: implied volatility (in sigma^2) of options at M = ln(K/K_ATMF) = 0
    initial_guess: initial guess of the parameters
    
    Return: fitted vol surface function at time t
    """
    rss_func = get_RSS_func(M, Vol, sig_ATMF)
    res = minimize(rss_func, initial_guess, method="SLSQP")
    return partial(vol_sur_M, sig_ATMF=sig_ATMF, params=res.x)

def get_vol_sur_M_all(chain, table, initial_guess=[-1, .25, .25]):
    """
    Fit vol surface of all t in the info table
    chain: processed option chain
    table: processed info table
    initial_guess: initial guess of the parameters
    
    Return: list containing the fitted vol surface and corresponding t
    """
    surface = []
    for i, T in enumerate(table.index):
        M = np.array(chain[i].index)
        Vol = np.square(np.array(chain[i]["IV"]))
        sig_ATMF = table.loc[T, "sig_ATMF"]**2
        surface.append([fit_vol_sur_M(M, Vol, sig_ATMF, initial_guess), T])
    return surface

def get_loc_vol_sur_M_all(chain, table, initial_guess=[-1, .25, .25]):
    """
    Fit local vol surface of all t in the info table
    chain: processed option chain
    table: processed info table
    initial_guess: initial guess of the parameters
    
    Return: list containing the fitted vol surface (in forward moneyness) and corresponding t
    """
    surface = []
    for i, T in enumerate(table.index):
        M = np.array(chain[i].index)
        Vol = np.square(np.array(chain[i]["LV"]))
        LV_ATMF = table.loc[T, "LV_ATMF"]**2
        surface.append([fit_vol_sur_M(M, Vol, LV_ATMF, initial_guess), T])
    return surface

def fit_vol_sur_t_Cubic(vol_surfaces, M):
    """
    M: log forward Moneyness
    Vol_sur_M: array of fitted vol sur function at t
    
    Return: CubicSpline function that return the implied volatility of a specific moneyness at time t
    """
    T = []
    Sig = []
    for vol_sur, t in vol_surfaces:
        T.append(t)
        Sig.append(vol_sur(M))
    return CubicSpline(T, Sig)

def get_IV_M_T_Cubic(vol_surfaces, M, T):
    """
    M: log forward Moneyness
    Vol_sur_M: array of fitted vol sur function at t
    T: time at which you want to know the implied volatility
    
    Return: implied/ local volatility of a specific moneyness at time t fitted by Cubic Spline
    """
    cs = fit_vol_sur_t_Cubic(vol_surfaces, M)
    return cs(T)

def get_IV_M_T_Linear(vol_surfaces, M, T_):
    """
    M: log forward Moneyness
    Vol_sur_M: array of fitted vol sur function at t
    T_: time at which you want to know the implied volatility
    
    Return: implied/ local volatility of a specific moneyness at time t fitted by linear interpolation
    """
    T = []
    Sig = []
    for vol_sur, t in vol_surfaces:
        T.append(t)
        Sig.append(vol_sur(M))
    bs = make_interp_spline(T, Sig, k=1, check_finite=False)
    return bs(T_)

def get_w(vol_surfaces, y, T, Cubic=True):
    """
    Compute total vol
    vol_surfaces: fitted implied vol surface functions
    y: Log forward moneyness
    T: Year to maturity
    Cubic: Use CubicSpline if True, use linear if False
    
    Return: total vol
    """
    if Cubic:
        return get_IV_M_T_Cubic(vol_surfaces, y, T)*T
    else:
        return get_IV_M_T_Linear(vol_surfaces, y, T)*T

def get_dw_dy(vol_surfaces, y, T, Cubic=True):
    """
    Compute dw/dy at T and y
    vol_surfaces: fitted vol surface functions
    y: Log forward moneyness
    T: Year to maturity
    Cubic: Use CubicSpline if True, use linear if False
    
    Return: dw/dy
    """
    f = partial(get_w, vol_surfaces, T=T, Cubic=Cubic)
    return derivative(f, y, dx=1e-3)

def get_d2w_dy2(vol_surfaces, y, T, Cubic=True):
    """
    Compute d2w/dy2 at T and y
    vol_surfaces: fitted vol surface functions
    y: Log forward moneyness
    T: Year to maturity
    Cubic: Use CubicSpline if True, use linear if False
    
    Return: d2w/dy2
    """
    f = partial(get_w, vol_surfaces, T=T, Cubic=Cubic)
    return derivative(f, y, dx=1e-3, n=2)

def get_dw_dT(vol_surfaces, y, T, Cubic=True):
    """
    Compute dw/dT at T and y
    vol_surfaces: fitted vol surface functions
    y: Log forward moneyness
    T: Year to maturity
    Cubic: Use CubicSpline if True, use linear if False
    
    Return: dw/dT
    """
    f = partial(get_w, vol_surfaces, y, Cubic=Cubic)
    return derivative(f, T, dx=1/252)

def get_local_vol(vol_surfaces, y, T, Cubic=True):
    """
    Compute local vol at T and y
    vol_surfaces: fitted vol surface functions
    y: Log forward moneyness
    T: Year to maturity
    Cubic: Use CubicSpline if True, use linear if False
    
    Return: local vol (in sigma^2)
    """
    dw_dT = get_dw_dT(vol_surfaces, y, T, Cubic=Cubic)
    dw_dy = get_dw_dy(vol_surfaces, y, T, Cubic=Cubic)
    d2w_dy2 = get_d2w_dy2(vol_surfaces, y, T, Cubic=Cubic)
    w = get_w(vol_surfaces, y, T, Cubic=Cubic)
    local_vol = dw_dT / (1-y/w*dw_dy+1/4*(-1/4-1/w+(y**2)/(w**2))*(dw_dy)**2+d2w_dy2/2)
    return local_vol