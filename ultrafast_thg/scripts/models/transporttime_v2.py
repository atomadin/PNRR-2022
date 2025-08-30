#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy

from .graphene_constants import GrapheneSample
from . import lindhard_v2 as lh
from . import graphene_thermodynamics_v2 as gt


#%% Several transport models.

class RelaxationTime():

    def tauFunc(self, e):
        pass



class ConstantTime(RelaxationTime):

    def __init__(self, tau=1000.0):
        p = {}
        # Constant relaxation time [fs].
        p["tau"] = tau
        self.p = p

    def tauFunc(self, e):
        r = self.p["tau"] + 0.0 * e
        return r



class PowerTime(RelaxationTime):

    # Calculates the relaxation time as
    # tau(energy) = a (energy)^b

    def __init__(self, a=100.0, b=1.0):
        p = {}
        p["a"] = a
        p["b"] = b
        self.p = p


    def tauFunc(self, e):
        r = self.p["a"] * np.power(e, self.p["b"])
        return r



class CoulombImpuritiesTime(RelaxationTime):

    # Temperature-dependent RPA screening.
    # Follow:
    # E.H. Hwang and S. Das Sarma, Phys. Rev. B 79, 165404 (2009).

    def __init__(self, nimpCM=0.5, dimp=0.0, epst=1.0, epsb=1.0,
                 theory="rpa"):

        # Physical parameters.
        p = GrapheneSample(epst=epst, epsb=epsb, qcutoff=0.0).p

        # Impurity density [10^12 cm^-2].
        p["nimpCM"] = nimpCM
        # Impurity density [nm^-2].
        p["nimp"] = nimpCM * 0.01
        # Distance between the impurity and the graphene sheet [nm].
        p["dimp"] = dimp
        # Theory for the relaxation time.
        p["theory"] = theory

        # Maximum absolute energy of the conical band dispersion [eV].
        p["enCutoff"] = 3.0

        # Interpolation mesh in the energy space.
        p["ee"] = np.logspace(-3, np.log10(p["enCutoff"]), 51)

        # Integral for the scattering time.  Perform the integration on a
        # dimensionless mesh from 0 to 2, as in the notes by Tallarico.
        # The variable y corresponds to the ratio (hbar v_F q) / energy.
        p["yy"] = np.linspace(0.001, 2.0, 21)
        p["yySq"] = p["yy"] * p["yy"]
        p["coeff"] = p["yySq"] * np.sqrt(1.0 - 0.25 * p["yySq"])

        self.p = p

        # Select the appropriate screening function.
        if (p["theory"] == "rpa"):
            self.lindhard = lh.RPAStatic(epst=p["epst"], epsb=p["epsb"])
        elif (p["theory"] == "thomasfermi"):
            self.lindhard = lh.ThomasFermiStatic(epst=p["epst"], epsb=p["epsb"])

        # Initialize the function to constant.
        self._tauFunc = lambda e: 0.0


    def tauFunc(self, e):
        return self._tauFunc(e)


    # Redefine the interpolation mesh.
    def interpMeshInit(self, eMin, eMax, eNum):
        self.p["ee"] = np.logspace(np.log10(eMin), np.log10(eMax), eNum)


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
        p["temp"] = p["tempK"] * p["kB"]

        # Carrier density in the two bands [nm^-2].
        nE = gt.oneband_dens_func(mu=muC, tempK=tempK)
        nH = gt.oneband_dens_func(mu=-muV, tempK=tempK)
        # Electron density [nm^-2].
        dens = nE - nH
        # Fermi energy [eV].
        p["eF"] = gt.twobands_ef_dens_func(dens=dens)

        # Prepare the calculation of the Lindhard function.
        if (p["theory"] == "rpa"):
            # The maximum wave vector qmax where the dielectric function epsFunc
            # is calculated is qmax = 2 emax / (hbar v_F), where emax is the
            # maximum energy where the relaxation time is calculated.
            qMax = 2.0 * p["ee"][-1] / p["hbarvf"]
            qNum = int(np.floor(1.0 * p["ee"][-1] / p["temp"]))
            self.lindhard.interpMeshInit(qMax, qNum)
            # The Maldague mesh needs to be large enough to allow for high
            # temperatures and Fermi energies.  This expressions seems to work
            # fine in a reasonable range.
            eMax = 10.0 * (np.abs(p["eF"]) + p["temp"])
            eNum = 301
            self.lindhard.maldagueMeshInit(eMax, eNum)

        # Initialize the Lindhard function at the thermodynamic parameters.
        self.lindhard.set_electron_thermo(muC=muC, muV=muV, tempK=tempK)

        # Calculate the transport time on the energy mesh and interpolate.
        tau = np.array([self.tauIntegral(e) for e in p["ee"]])
        self._tauFunc = scipy.interpolate.interp1d(p["ee"], tau, kind='linear')


    def tauIntegral(self, e):
        p = self.p
        # Calculate the static screened interaction potential on the mesh.
        q = e / p["hbarvf"]
        ww = np.array([self.coulombScreened(y, q) for y in p["yy"]])
        wwSq = np.abs(ww * ww)
        # Perform the integration.
        ical = (
            (1.0 / (2.0 * np.pi))
            *  scipy.integrate.simps(p["coeff"] * wwSq, p["yy"]))
        # Calculate the scattering time.
        tau = e * p["hbar"] / (p["nimp"] * p["hbarvf"] * p["hbarvf"] * ical)
        return tau


    # Static screened Coulomb potential between electrons and impurities,
    # dimensionless.
    def coulombScreened(self, y, q):
        # y: dimensionless, (hbar v_F q) / energy
        # dimp: distance between impurities and electrons
        p = self.p
        w = (
            np.exp(-y * q * p["dimp"]) * (2.0 * np.pi * p["alphaee"] / y) /
            self.lindhard.epsFunc(y * q))
        return w
