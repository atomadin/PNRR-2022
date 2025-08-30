#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import pickle
import multiprocessing as mp

from .graphene_constants import GrapheneConstants


_warnings = []

################################################################################
# General parameters and formulas.

# Physical parameters.
p = GrapheneConstants().c.copy()

# Energy mesh [eV].
# enMesh = np.linspace(0.001, 3.0, 1001)
enMesh = np.logspace(-5.0, np.log10(6.0), 2001)


# Fermi-Dirac distribution.
def fermi_dirac_func(e, mu, tempK, range=15.0):
    temp = p["kB"] * tempK
    x = (e - mu) / temp
    if (type(x) == np.ndarray):
        rA = np.where(x < -range, 1.0, 0.0)
        x2 = np.where(np.abs(x) < range, x, 0.0)
        rB = np.where(np.abs(x) < range, 1.0 / (np.exp(x2) + 1.0), 0.0)
        r = rA + rB
    else:
        if (x > range):
            r = 0.0
        elif (x < -range):
            r = 1.0
        else:
            r = 1.0 / (np.exp(x) + 1.0)
    return r


# Derivative of the Fermi-Dirac distribution.
def fermi_dirac_deriv_func(e, mu, tempK, range=15.0):
    temp = p["kB"] * tempK
    x = (e - mu) / temp
    if (type(x) == np.ndarray):
        x2 = np.where(np.abs(x) < range, x, 0.0)
        u = np.where(np.abs(x) < range, np.exp(x2), 0.0)
        c = 1 + u
        r = - u / (c * c * temp)
    else:
        if (x > range):
            r = 0.0
        elif (x < -range):
            r = 0.0
        else:
            u = np.exp(x)
            c = 1 + u
            r = - u / (c * c * temp)
    return r


# Bose-Einstein distribution.
def bose_einstein_func(e, mu, tempK):
    temp = p["kB"] * tempK
    r = 1.0 / (np.exp((e - mu) / temp) - 1.0)
    return r


# Density of states.
def dos_func(e):
    r = p["Nf"] * np.abs(e) / (2.0 * np.pi * p["hbarvf"] * p["hbarvf"])
    return r


# Evaluate the density of states on the mesh.
nuMesh = dos_func(enMesh)


################################################################################
# One band.
# The carrier density in the band is conserved.

# Carrier density [nm^-2] from Fermi energy [eV].
def oneband_dens_ef_func(ef):
    l = np.sign(ef)
    if (l < 0.0):
        r = 0.0
    else:
        r = ef * ef / (np.pi * p["hbarvf"] * p["hbarvf"])
    return r


# Fermi energy [eV] from carrier density [nm^-2].
def oneband_ef_dens_func(dens):
    if (np.any(dens < 0.0)):
        print("Error: Density in a band must be non negative.")
        ef = np.nan * dens
    else:
        ef = p["hbarvf"] * np.sqrt(np.pi * dens)
    return ef


# Carrier density [nm^-2].
def oneband_dens_func(mu, tempK):
    # fMesh = np.array([fermi_dirac_func(en, mu, tempK) for en in enMesh])
    fMesh = fermi_dirac_func(enMesh, mu, tempK)
    dens = scipy.integrate.simps(nuMesh * fMesh, enMesh)
    return dens


# Minimum one-band density that can reliably be represented on the mesh
# as temperature varies and chemical potential decreases.
# It depends on the mesh and the cutoff that we impose on the exponentials
# of the Fermi-Dirac distribution.
# Use the function below to evaluate it "by eye" once the mesh is fixed.
_densmin = 1.0e-11

def exploreDensityMinimum():
    # Calculate the minimum non-zero density that we obtain as the
    # chemical potential varies, as a function of the temperature.
    # The density eventually goes to exactly zero because of the cutoff
    # that we impose on the tails of the Fermi-Dirac distribution
    # (see below).
    tt = np.logspace(np.log10(10.0), np.log10(3000.0), 101)
    mm = np.linspace(-5.0, 2.0, 10001)
    densmin = []
    for tempK in tt:
        dd = np.array([gt.oneband_dens_func(mu, tempK) for mu in mm])
        densmin.append(np.min(dd[np.nonzero(dd)]))
    densmin = np.array(densmin)
    return densmin


# Energy density [eV nm^-2].
def oneband_energy_func(mu, tempK):
    # fMesh = np.array([fermi_dirac_func(en, mu, tempK) for en in enMesh])
    fMesh = fermi_dirac_func(enMesh, mu, tempK)
    en = scipy.integrate.simps(nuMesh * fMesh * enMesh, enMesh)
    return en


# Chemical potential [eV] at fixed temperature [K] from carrier density
# [nm^-2].
def oneband_mu_func_bisect(dens, tempK, muMin=-6.0, muMax=2.0, quiet=True):
    try:
        mu = scipy.optimize.bisect(
            lambda mu1: (dens - oneband_dens_func(mu1, tempK)),
            muMin, muMax)
    except ValueError:
        if not quiet:
            print("Error in oneband_mu_func_bisect")
            print("    dens = %f\n    tempK = %f" % (dens, tempK))
        _warnings.append("Error in oneband_mu_func_bisect: dens = %f, tempK = %f" % (dens, tempK))
        raise
    return mu

# By default, calculate chemical potential by bisection.
oneband_mu_func = oneband_mu_func_bisect

    
# Heat capacity [eV K^-1 nm^-2].
def oneband_heat_capacity_func(dens, tempK, dtK=10.0):
    # Make sure that the step is large enough to avoid numerical noise
    # in the value of the energy.
    en0 = oneband_energy_func(oneband_mu_func(dens, tempK), tempK)
    large_enough = False
    while (not large_enough):
        tempK1 = tempK + dtK
        en1 = oneband_energy_func(oneband_mu_func(dens, tempK1), tempK1)
        if np.abs((en1 - en0)/en0) < 0.05:
            dtK = dtK * 2.0
        else:
            large_enough = True
    # Make sure that the step for the finite difference does not
    # produce negative temperatures.
    dtK = np.minimum(dtK, tempK / 4.0)
    # Calculate the derivative dE/dT with a fifth-order formula.
    tK_l = [tempK - 2.0 * dtK, tempK - dtK, tempK + dtK, tempK + 2.0 * dtK]
    en_l = [oneband_energy_func(oneband_mu_func(dens, tK), tK) for tK in tK_l]
    cv = (en_l[0] - 8.0 * en_l[1] + 8 * en_l[2] - en_l[3]) / 12.0 / dtK
    return cv


def oneband_heat_capacity_undoped(tempK):
    temp = p["kB"] * tempK
    r = 21.6 * p["kB"] * temp * temp / (np.pi * p["hbarvf"] * p["hbarvf"])
    return r
    

################################################################################
# One electron and one hole band.
# The difference of the carrier densities is conserved.

# In the following, dens represents the e/h (i.e. electron when positive or
# hole when negative) carrier density at zero temperature, in the absence of
# photoexcitation, i.e. when the electron distribution is a Fermi step.

# In other words, dens is the density of electrons in the system, minus the
# the density of electrons at the charge neutrality point.

# The electron density in conduction band and the hole density in valence band
# at finite temperature are denoted nE and nH, respectively.

# Because of particle conservation, at any temperature, and even after
# photoexcitation, it holds that:
# dens = nE - nH

# Fermi energy from the e/h density.
def twobands_ef_dens_func(dens):
    ef = np.sign(dens) * p["hbarvf"] * np.sqrt(np.pi * np.abs(dens))
    return ef


# e/h density from the Fermi energy.
def twobands_dens_ef_func(ef):
    dens = np.sign(ef) * oneband_dens_ef_func(np.abs(ef))
    return dens

# e/h density [nm^-2].
def twobands_dens_func(mu, tempK):
    nE = oneband_dens_func(mu, tempK)
    nH = oneband_dens_func(-mu, tempK)
    dens = nE - nH
    return dens

# Energy density in both bands.
def twobands_energy_func(mu, tempK):
    en = oneband_energy_func(mu, tempK) + oneband_energy_func(-mu, tempK)
    return en


# Chemical potential [eV] at fixed temperature [K] from e/h density [nm^-2] at
# zero temperature.
def twobands_mu_func_bisect(dens, tempK, muMin=-1.0, muMax=1.0, quiet=True):
    try:
        mu = scipy.optimize.bisect(
            lambda mu1: (
                dens - (oneband_dens_func(mu1, tempK)
                        - oneband_dens_func(-mu1, tempK))),
            muMin, muMax)
    except ValueError:
        if not quiet:
            print("Error in twobands_mu_func_bisect")
            print("    dens = %f\n    tempK = %f" % (dens, tempK))
        _warnings.append("Error in twobands_mu_func_bisect: dens = %f, tempK = %f" % (dens, tempK))
        raise
    return mu
            
# By default, calculate chemical potentials by bisection.
twobands_mu_func = twobands_mu_func_bisect


# Electron and hole densities in the two bands at finite temperature.
def twobands_nE_nH_func(dens, tempK):
    mu = twobands_mu_func(dens, tempK)
    nE = oneband_dens_func(mu, tempK)
    nH = oneband_dens_func(-mu, tempK)
    return (nE, nH)


# e/h density [nm^-2] at zero temperature given the carrier density in either
# band at finite temperature [K].
# In this function we defined the carrier density as nE if positive
# or -nH if negative.
# WARNING: This function was called twobands_dens_func in a previous
# version of the code.
def twobands_dens_carrier_func(nC, tempK):
    if nC >= 0.0:
        nE = nC
        # Chemical potential at finite temperature.
        mu = oneband_mu_func(nE, tempK)
        # Hole density at finite temperature.
        nH = oneband_dens_func(-mu, tempK)
    else:
        nH = -nC
        mu = -oneband_mu_func(nH, tempK)
        nE = oneband_dens_func(mu, tempK)
    # e/h density at zero temperature.
    dens = nE - nH
    return dens


# Heat capacity [eV K^-1 nm^-2].
def twobands_heat_capacity_func(dens, tempK, dtK=10.0): 
    # Make sure that the step is large enough to avoid numerical noise
    # in the value of the energy.
    en0 = twobands_energy_func(twobands_mu_func(dens, tempK), tempK)
    large_enough = False
    while (not large_enough):
        tempK1 = tempK + dtK
        en1 = twobands_energy_func(twobands_mu_func(dens, tempK1), tempK1)
        if np.abs((en1 - en0)/en0) < 0.05:
            dtK = dtK * 2.0
        else:
            large_enough = True
    # Make sure that the step for the finite difference does not
    # produce negative temperatures.
    dtK = np.minimum(dtK, tempK / 4.0)
    # Calculate the derivative dE/dT with a fifth-order formula.
    tK_l = [tempK - 2.0 * dtK, tempK - dtK, tempK + dtK, tempK + 2.0 * dtK]
    en_l = [twobands_energy_func(twobands_mu_func(dens, tK), tK) for tK in tK_l]
    cv = (en_l[0] - 8.0 * en_l[1] + 8 * en_l[2] - en_l[3]) / 12.0 / dtK
    return cv


# Analytical approximations to the heat capacity at equilibrium.
def heat_capacity_undoped(tempK):
    temp = p["kB"] * tempK
    r = 21.6 * p["kB"] * temp * temp / (np.pi * p["hbarvf"] * p["hbarvf"])
    return r


def heat_capacity_doped(dens, tempK):
    ef = twobands_ef_dens_func(dens)
    temp = p["kB"] * tempK
    r = np.pi * np.pi * dos_func(ef) * p["kB"] * temp / 3.0
    return r


################################################################################
# Two bands out of equilibrium, with finite photoexcited density.

# Energy density in both bands.
def photoexc_energy_func(muC, muV, tempK):
    en = oneband_energy_func(muC, tempK) + oneband_energy_func(-muV, tempK)
    return en


# Electron and hole density in the two bands after photoexcitation.
def photoexc_nE_nH_func(dens, dnE, tempK):
    # Densities in the two bands at finite temperature before photoexcitation.
    nE0, nH0 = twobands_nE_nH_func(dens, tempK)
    # Increase the density of both electrons and holes by the photoexcited
    # density.  It holds:
    # dens = nE0 - nH0 = nE - nH
    # Photoexcitation and temperature increase "commute," it does not matter
    # how one reaches a certain state given tempK and dnE, because an increase
    # in the temperature alone cannot change dnE.
    nE = nE0 + dnE
    nH = nH0 + dnE
    return (nE, nH)


# Chemical potentials of a photoexcited system.
def photoexc_mu_func(dens, dnE, tempK):
    # Calculate the electron and hole densities.
    nE, nH = photoexc_nE_nH_func(dens, dnE, tempK)
    # Calculate the chemical potentials in each band, corresponding to
    # the temperature tempK and the densities just calculated.
    
    # Here we may encounter the problem that if one of the two densities is
    # too low the chemical potential cannot be found.  This happens e.g. if the
    # temperature is very low and the Fermi energy is away from the Dirac point.
    # In this case we assume that the two chemical potentials are equal.
    # We assume that at least one of the two densities is sizable, otherwise
    # there is a problem of physics in the calling code.
    if (np.abs(nE) < _densmin):
        muV = -oneband_mu_func(nH, tempK)
        muC = muV
        _warnings.append("photoexc_mu_func: nE too small.")
    elif (np.abs(nH) < _densmin):
        muC = oneband_mu_func(nE, tempK)
        muV = muC
        _warnings.append("photoexc_mu_func: nH too small.")
    else:
        muC = oneband_mu_func(nE, tempK)
        muV = -oneband_mu_func(nH, tempK)
    return (muC, muV)


# Heat capacity of a photoexcited system.
# The photoexcited density is kept fixed.
def photoexc_heat_capacity_func(dens, dnE, tempK, dtK=10.0):
    # Make sure that the step is large enough to avoid numerical noise
    # in the value of the energy.
    en0 = photoexc_energy_func(*photoexc_mu_func(dens, dnE, tempK), tempK)
    large_enough = False
    while (not large_enough):
        tempK1 = tempK + dtK
        en1 = photoexc_energy_func(*photoexc_mu_func(dens, dnE, tempK1), tempK1)
        if np.abs((en1 - en0)/en0) < 0.05:
            dtK = dtK * 2.0
        else:
            large_enough = True
    # Make sure that the step for the finite difference does not
    # produce negative temperatures.
    dtK = np.minimum(dtK, tempK / 4.0)
    # Calculate the derivative dE/dT with a fifth-order formula.
    tK_l = [tempK - 2.0 * dtK, tempK - dtK, tempK + dtK, tempK + 2.0 * dtK]
    en_l = [photoexc_energy_func(*photoexc_mu_func(dens, dnE, tK), tK) for tK in tK_l]
    cv = (en_l[0] - 8.0 * en_l[1] + 8 * en_l[2] - en_l[3]) / 12.0 / dtK
    return cv


# Sum of the heat capacity of electrons and holes.
# The electron and hole densities are kept fixed.
def photoexc_heat_capacity_func_2(dens, dnE, tempK, dtK=10.0):    
    # Calculate the electron and hole densities.
    nE, nH = photoexc_nE_nH_func(dens, dnE, tempK)
    # Calculate the heat capacity separately in the two bands.
    # If the density in a band is too low, just use vanishing heat capacity.
    # Its value will be negligible in the final summation.
    if (nE > _densmin):
        cvE = oneband_heat_capacity_func(nE, tempK, dtK=10.0)
    else:
        cvE = 0.0
    if (nH > _densmin):
        cvH = oneband_heat_capacity_func(nH, tempK, dtK=10.0)
    else:
        cvH = 0.0
    cv = cvE + cvH
    return cv


# Derivatives of the chemical potentials in the two bands.

# Used in the calculation of linear Boltzmann transport coefficients.
# The chemical potentials are seen as functions of:
# - the Fermi energy, which is a 1-to-1 to the total electron density;
# - the electron temperature, which is common to the two bands;
# - the photoexcited density in the conduction band.

def saturation(x_val, x_sat):
    y_val = np.abs(x_val)
    s = x_val / y_val
    y_sat = np.abs(x_sat) * 2.0 / np.pi
    r = s * y_sat * np.arctan(y_val / y_sat)
    return r


def photoexc_mu_deriv(ef, tempK, dnE, d_ef=0.010, d_tempK=10.0, d_dnE=1.0e-5):
    # Make sure that increments are not too large.
    # The photoexcited density in conduction band has to be positive,
    # and it cannot be too small or it is negligible.  Of course this
    # could be problematic if at the same time the Fermi energy and the
    # temperature are very small.  But we are not interested in this
    # degenerate scenario.  Remember that 0.01nm^-2 is 10^12cm^-2.
    # So a change in photoexcited density 10^-5nm^-2, i.e. 10^9cm^-2 is
    # already at the order of disorder fluctuations.
    # The temperatures that we deal with are always at least room temperature,
    # so there is no need to modify d_tempK.
    # The Fermi energy can go through zero.  In this case, let us maximize
    # the derivative step using the temperature.  Moreover, in the following
    # derivatives wrt the Fermi energy will be taken symmetric.
    temp = p["kB"] * tempK
    d_ef = saturation(d_ef, 0.2 * np.sqrt(ef**2 + temp**2))

    # Incremented variables.
    dens = twobands_dens_ef_func(ef)
    densL = twobands_dens_ef_func(ef - 0.5 * d_ef)
    densR = twobands_dens_ef_func(ef + 0.5 * d_ef)
    tempK1 = tempK + d_tempK
    dnE1 = dnE + d_dnE

    # Calculate discrete derivatives.
    dmu = {}

    # With respect to the Fermi energy.
    # Notice that the function is called on the corresponding densities.
    muCL, muVL = photoexc_mu_func(dens=densL, dnE=dnE, tempK=tempK)
    muCR, muVR = photoexc_mu_func(dens=densR, dnE=dnE, tempK=tempK)
    dmu["muC_ef"] = (muCR - muCL) / d_ef
    dmu["muV_ef"] = (muVR - muVL) / d_ef

    # With respect to the temperature.
    muC0, muV0 = photoexc_mu_func(dens=dens, dnE=dnE, tempK=tempK)
    muC1, muV1 = photoexc_mu_func(dens=dens, dnE=dnE, tempK=tempK1)
    dmu["muC_tempK"] = (muC1 - muC0) / d_tempK
    dmu["muV_tempK"] = (muV1 - muV0) / d_tempK

    # With respect to the photoexcited density.
    muC1, muV1 = photoexc_mu_func(dens=dens, dnE=dnE1, tempK=tempK)
    dmu["muC_dnE"] = (muC1 - muC0) / d_dnE
    dmu["muV_dnE"] = (muV1 - muV0) / d_dnE

    return dmu



################################################################################
# Prepare the calculation of the chemical potentials by interpolation.
# When this function is called, the functions that calculate the chemical
# potential by e.g. bisection are used to populate a lookup table.
# An interpolating function is then built and substituted as the default
# calculation method for the chemical potentials.

# Function to calculate the loookup table of chemical potentials using
# parallelization.

def _one_mu_val(x):
    err = False
    try:
        mu_one = oneband_mu_func_bisect(
            x["dens"], x["tempK"], muMin=x["muMin"], muMax=x["muMax"], quiet=True)
    except ValueError:
        mu_one = np.NaN
        err = True
    try:
        mu_two = twobands_mu_func_bisect(
            x["dens"], x["tempK"], muMin=x["muMin"], muMax=x["muMax"])
    except ValueError:
        mu_two = np.NaN
        err = True
    if err:
        _warnings.append("Value errors in the calculation of the "
                         "lookup table for the chemical potentials.")
    return ({"i_dens": x["i_dens"], "i_tempK": x["i_tempK"],
             "mu_one": mu_one, "mu_two": mu_two})

# Save values when interpolation does not work.
_out_of_bounds_one = []
_out_of_bounds_two = []


def mu_func_use_2d_interpolation(dens_min=1.0e-5, dens_max=1.0, dens_num=30,
                                 tempK_min=100.0, tempK_max=4000.0, tempK_num=30,
                                 muMin=-2.0, muMax=2.0,
                                 save_to=None,
                                 load_from=None,
                                 proc_num=1):
    # Create or modify a few global functions.
    global oneband_mu_func_interp2d
    global twobands_mu_func_interp2d
    global oneband_mu_func
    global twobands_mu_func
    
    # First step: calculate or read the lookup table.
    if load_from is not None:
        # The data for the interpolation has already been calculated.
        with open(load_from, "rb") as f:
            mu_func_interpolation_data = pickle.load(f)
        tempK_log_list = mu_func_interpolation_data["tempK_log_list"]
        dens_log_list = mu_func_interpolation_data["dens_log_list"]
        oneband_mu_mesh = mu_func_interpolation_data["oneband_mu_mesh"]
        twobands_mu_mesh_n = mu_func_interpolation_data["twobands_mu_mesh_n"]
    else:
        # Define two logarithmic meshes for density and temperature.
        tempK_log_list = np.linspace(np.log10(tempK_min), np.log10(tempK_max), tempK_num)
        dens_log_list = np.linspace(np.log10(dens_min), np.log10(dens_max), dens_num)
        # Generate a list of parameters for the parallelization.
        xx = []
        for i_dens,dens_log in enumerate(dens_log_list):
            for i_tempK,tempK_log in enumerate(tempK_log_list):
                dens = 10**dens_log
                tempK = 10.0**tempK_log
                xx.append({"i_dens": i_dens, "i_tempK": i_tempK,
                           "dens": dens, "tempK": tempK,
                           "muMin": muMin, "muMax": muMax})
        # Calculate the chemical potentials for all the parameters.
        with mp.Pool(proc_num) as p:
            yy = p.map(_one_mu_val, xx)
        # Copy the results to matrices.
        oneband_mu_mesh = np.zeros((dens_num, tempK_num))
        twobands_mu_mesh_n = np.zeros((dens_num, tempK_num))
        for y in yy:
            oneband_mu_mesh[y["i_dens"],y["i_tempK"]] = y["mu_one"]
            twobands_mu_mesh_n[y["i_dens"],y["i_tempK"]] = y["mu_two"]
        # Save interpolation data.
        mu_func_interpolation_data = {
            "dens_log_list": dens_log_list,
            "tempK_log_list": tempK_log_list,
            "oneband_mu_mesh": oneband_mu_mesh,
            "twobands_mu_mesh_n": twobands_mu_mesh_n
            # "twobands_mu_mesh_p": twobands_mu_mesh_p,
        }
        if save_to is not None:
            with open(save_to, "wb") as f:
                pickle.dump(mu_func_interpolation_data, f)

    # Second step: define the interpolating functions.

    # Prepare the interpolator for the values of mu on the mesh.
    oneband_mu_interpolator = scipy.interpolate.RegularGridInterpolator(
        (dens_log_list, tempK_log_list), oneband_mu_mesh)
    twobands_mu_interpolator_n = scipy.interpolate.RegularGridInterpolator(
        (dens_log_list, tempK_log_list), twobands_mu_mesh_n)

    # The function that calculates the chemical potential now must first
    # calculate the logarithms of the arguments and then call the interpolator.

    def oneband_mu_func_interp2d(dens, tempK, muMin=muMin, muMax=muMax):
        try:
            mu = float(oneband_mu_interpolator((np.log10(dens), np.log10(tempK))))
        except ValueError:
            mu = oneband_mu_func_bisect(dens, tempK, muMin=muMin, muMax=muMax, quiet=True)
            _out_of_bounds_one.append([dens, tempK, mu])
        return mu
    
    def twobands_mu_func_interp2d(dens, tempK, muMin=muMin, muMax=muMax):
        if (dens >= 0.0):
            try:
                mu = float(twobands_mu_interpolator_n((np.log10(dens), np.log10(tempK))))
            except ValueError:
                mu = twobands_mu_func_bisect(dens, tempK, muMin=muMin, muMax=muMax)
                _out_of_bounds_two.append([dens, tempK, mu])
        else:
            try:
                mu = -float(twobands_mu_interpolator_n((np.log10(-dens), np.log10(tempK))))
            except ValueError:
                mu = twobands_mu_func_bisect(dens, tempK, muMin=muMin, muMax=muMax)
                _out_of_bounds_two.append([dens, tempK, mu])
        return mu
            
    # Once we prepare the interpolation we assume that the default function
    # to calculate the chemical potential is the interpolated one.
    oneband_mu_func = oneband_mu_func_interp2d
    twobands_mu_func = twobands_mu_func_interp2d
    
    # Return the interpolation data for inspection.
    return mu_func_interpolation_data

def mu_func_use_bisection():
    global oneband_mu_func_bisect
    global twoband_mu_func_bisect
    global oneband_mu_func
    global twobands_mu_func
    oneband_mu_func = oneband_mu_func_bisect
    twobands_mu_func = twobands_mu_func_bisect

    