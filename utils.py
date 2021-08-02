#!/usr/bin/env python

import numpy as np
from scipy import special
import pylab as plt
from scipy.optimize import curve_fit


def sample_from_powerlaw(logminval, logmaxval, power):
    """Draw a value randomly from the specified power law"""

    logmin = logminval #np.log10(minval)
    logmax = logmaxval #np.log10(maxval)

    c = -1.0 * logmax * power
    nmax = 10.0**(power*logmin + c)

    # slightly worried about inf loops here...
    while True:
        log = np.random.uniform(logmin, logmax)
        n = 10.0**(power*log + c)

        if nmax*np.random.random() <= n:
            break

    return 10.0**log


def gauss_integral(fstart, fend, mu_f, sig_f):
    # https://www.wolframalpha.com/input/
    # ?i=integral+of+%281%2F%28sqrt%282*pi%29*sigma%29%29*e%5E%28-
    # %28%28x-mu%29%2F%28sigma%29%29%5E2%2F2%29+dx
    
    # Verify if this works!!
    inp_start = (fstart - mu_f)/(np.sqrt(2*np.pi) * sig_f)
    inp_end = (fend - mu_f)/(np.sqrt(2*np.pi) * sig_f)
    integral = 0.5 * (special.erf(inp_end) - special.erf(inp_start))
    return integral


def energy_to_fluence(E, sig_f, distance=972):
    # E in ergs
    # sig_f in MHz
    D_l = distance * 10**6 * 3 * 10**18 # cm 
    fac = 10**(-23)
    return E / (2.355*sig_f*10**6) / (4*np.pi*D_l**2) / fac # Jy s 


def fluence_to_energy(S, freq_width, distance=972):
    # S in Jy s
    # freq_width in MHz
    D_l = distance * 10**6 * 3 * 10**18 # cm 
    fac = 10**(-23)
    return S * (freq_width*10**6) * (4*np.pi*D_l**2) * fac # ergs


def single_power_law(x, a, b):
    return a*x**(b)


def fit_single_power_law(Ns, Es):
    popt, perr = curve_fit(single_power_law, Es, Ns)
    slope = popt[1]
    scale = popt[0]
    return scale, slope


def double_power_law(x, amplitude, x_break, alpha_1, alpha_2):
    """One dimensional broken power law model function"""

    alpha = np.where(x < x_break, alpha_1, alpha_2)
    xx = x / x_break
    return amplitude * xx ** (alpha)


def fit_double_power_law(Ns, Es):
    popt, perr = curve_fit(double_power_law, Es, Ns, 
                           bounds=([-np.inf, 10**(0), -5, -5],
                                   [np.inf, 10**(4), 5, 5]))
    slope1 = popt[-2]
    slope2 = popt[-1]
    scale = popt[0]
    break_e = popt[1] 
    return scale, slope1, slope2, break_e
