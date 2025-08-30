#/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.integrate import simps

from . import graphene_thermodynamics_v2 as gt
from .graphene_constants import GrapheneConstants

# Refractive indices.

n_r = {
    "silica": 1.44, # fused silica, SiO2
    "ionic": 1.42, # ionic gate
    "sapphire": 1.73
}


class GrapheneOptics:

    def __init__(self):
        # Physical parameters.
        self.p = GrapheneConstants().c.copy()

    # Optical conductivity [nm/fs], at zero wave vector, given the photon
    # energy [eV].
    def optical_conductivity(self, eph):
        p = self.p
        # The electron and hole energy are half the photon energy.
        en = eph / 2.0
        sigma = self.p["sigma_uni"] * (
            1 - gt.fermi_dirac_func(en, p["muC"], p["tempK"])
            - gt.fermi_dirac_func(en, -p["muV"], p["tempK"]))
        return sigma


    # Absorbance given the photon energy [eV].
    def absorbance(self, eph):
        sigma = self.optical_conductivity(eph)
        # Expression valid for |sigma / c| << 1.
        alpha = self.p["alpha_coeff"] * sigma
        return alpha


    # Differential transmission given the photon energy [eV].
    def transmittance(self, eph):
        p = self.p
        sigma = self.optical_conductivity(eph)
        # Expression valid for |sigma / c| << 1.
        transm = p["transm_fresnel"] * (1.0 - p["transm_coeff"] * sigma)
        return transm


    # Set the chemical potentials [eV] and the temperature [K] of the electron
    # distribution.
    def set_electron_thermo(self, muC, muV, tempK):
        p = self.p
        # Chemical potential in conduction band [eV].
        p["muC"] = muC
        # Chemical potential in valence band [eV].
        p["muV"] = muV
        # Temperature [K].
        p["tempK"] = tempK


    # Refractive indices.  The incident and reflected wave are in the top
    # region and the transmitted wave is in the substrate.
    def set_refractive_indices(self, nsub=n_r["silica"], ntop=n_r["ionic"]):
        p = self.p
        # Refractive indices.
        p["nsub"] = nsub
        p["ntop"] = ntop
        # Calculate derived coefficients.
        # Coefficient in the absorbance.
        p["alpha_coeff"] = (
            (4.0 * np.pi / self.p["c"])
            * 4.0 * p["ntop"]
            / ((p["ntop"] + p["nsub"]) * (p["ntop"] + p["nsub"])))
        # Coefficient in the transmittance.
        p["transm_coeff"] = (
            (4.0 * np.pi / self.p["c"]) * 2.0 / (p["ntop"] + p["nsub"]))
        # Transmittance for vanishing conductivity.
        p["transm_fresnel"] = (
            4.0 * (p["ntop"] * p["nsub"])
            / ((p["ntop"] + p["nsub"]) * (p["ntop"] + p["nsub"])))


        
class GrapheneTHG():
    
    def __init__(self, nsub=n_r["silica"], ntop=n_r["ionic"], eph=0.32,
                 gammaConst=None, gammaProp=None, gammaInv=None,
                 maldagueMin=0.001, maldagueMax=2.0, maldagueNum=301):
        # Physical parameters.
        self.p = GrapheneConstants().c.copy()        
        p = self.p
        # Refractive indices.
        p["nsub"] = nsub
        p["ntop"] = ntop
        # Energy of probe photon [eV]. THG will be 3x.
        p["eph"] = eph
        # Coefficients for the calculation of the scattering rate.
        self.coeffGamma(gammaConst=gammaConst, gammaProp=gammaProp, gammaInv=gammaInv)
        # Coefficients for the calculation of eta.
        # Coefficient K [eV^6 nm^4 fs^2].
        p["K"] = (np.pi * p["hbar"] * p["e_sq"]**2 / 24.0)**2 * (p["v_F"] / p["c"])**4
        # Ratio of refractive indices.
        p["nfact"] = nsub / (ntop + nsub)**2 / ntop**3
        # Mesh for the Maldague integral of S.
        self.eMesh = np.linspace(maldagueMin, maldagueMax, maldagueNum)
        # self.eMesh = np.logspace(np.log10(maldagueMin), np.log10(maldagueMax), maldagueNum)

        
    def coeffGamma(self, gammaConst=None, gammaProp=None, gammaInv=None):
        # Coefficients for the calculation of the scattering rate Gamma.
        # Define as a separate function in case that it is convenient to reuse
        # the same GrapheneTHG object with different coefficients, e.g. if they
        # depend on time.  In this case, just give the correct value of
        # gamma as gammaConst, and set the other parameters to zero.
        ## The default values: gammaProp=0.0, gammaInv=0.0013 have been
        ## obtained by Rostami and Soavi using a fit to the conductivity.
        if gammaConst is None:
            gammaConst = 0.0
        if gammaProp is None:
            gammaProp = 0.0
        if gammaInv is None:
            gammaInv = 0.0013
        # Constant [eV].
        self.p["gammaConst"] = gammaConst
        # Proportional to the Fermi energy [dimensionless].
        self.p["gammaProp"] = gammaProp
        # Inversely proportional to the Fermi energy [eV^2].
        self.p["gammaInv"] = gammaInv
        
        
    def funcGamma(self, muC=0.2, muV=0.2):
        muAv = 0.5 * (np.abs(muC) + np.abs(muV))
        g = self.p["gammaConst"] + self.p["gammaProp"] * muAv + self.p["gammaInv"] / muAv
        return g

        
    def calculate_eta(self, powdens=1.0e-2, muC=0.2, muV=0.2, tempK=300.0):
        p = self.p
        # Probe power density [eV / fs nm^2].
        p["powdens"] = powdens
        # Chemical potentials in the two bands [eV].
        p["muC"] = muC
        p["muV"] = muV
        # Electron temperature [K].
        p["tempK"] = tempK
        # Convert temperature to [eV].
        p["temp"] = p["kB"] * p["tempK"]
        # The relaxation rate [eV].
        p["Gamma"] = self.funcGamma(p["muC"], p["muV"])
        # Calculate the function S [eV^-4] by performing the Maldague integral.
        self.calculate_S()
        # The THG efficiency [dimensionless].
        p["eta"] = p["K"] * p["nfact"] * p["powdens"]**2 * np.abs(p["S"])**2
        
    
    def funcS0(self, eF=0.2):
        p = self.p
        # eF is a Fermi energy [eV].
        # Dimensionless function.
        G = lambda z: np.log((z + 1.0) / (z - 1.0))
        # Dimensionless complex variable.
        # Here the values of the chemical potentials used to calculate Gamma
        # are *not* related to eF, because this function is used as the
        # integrand in the Maldague integral, where eF is an integration
        # variable, while the chemical potentials (and hence Gamma) refer to
        # the actual system at finite temperature.
        w = p["eph"] + 1.0j * p["Gamma"]
        # When eF=0, assume that the limit is taken from the direction of
        # positive Gamma, i.e. negative imaginary part of eF/w.
        eps = 1.0e-12
        zF = np.abs(eF) / w - 1.0j * eps
        s = ( 17.0 * G(2.0 * zF) - 64.0 * G(zF) + 45.0 * G(2.0 * zF / 3.0) )
        s = s / w**4
        return s
    
    
    def calculate_S(self):
        p = self.p
        # Evaluate the contributions to the Maldague integrand on the mesh.
        zm = (self.eMesh - p["muC"]) / p["temp"]
        zp = (self.eMesh + p["muV"]) / p["temp"]
        emax = 15.0
        zmk = np.where(np.abs(zm) < emax, zm, 0.0)
        kC = np.where(np.abs(zm) < emax, 1.0 / np.cosh(0.5 * zmk)**2, 4.0 * np.exp(-2.0 * np.abs(zm)))
        zpk = np.where(np.abs(zp) < emax, zp, 0.0)
        kV = np.where(np.abs(zp) < emax, 1.0 / np.cosh(0.5 * zpk)**2, 4.0 * np.exp(-2.0 * np.abs(zp)))
        ss = self.funcS0(self.eMesh)
        yy = (kC + kV) * ss
        self.ss = ss  # save integrand for inspection
        self.yy = yy  # save integrand for inspection
        r1 = simps(yy, self.eMesh) / 4.0 / p["temp"]
        # Evaluate the correction depending on the value at eF = 0.
        rC = gt.fermi_dirac_func(0.0, p["muC"], p["tempK"])
        rV = gt.fermi_dirac_func(0.0, p["muV"], p["tempK"])
        s0 = self.funcS0(0.0)
        r2 = (rC - rV) * s0
        # The result is the difference of the two values.
        p["S"] = r1 - r2
