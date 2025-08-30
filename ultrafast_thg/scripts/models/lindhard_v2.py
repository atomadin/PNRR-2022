#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy
import scipy.interpolate

from .graphene_constants import GrapheneSample
from . import graphene_thermodynamics_v2 as gt


#%% Several screening models.


class EpsStatic():

    def epsFunc(self, q):
        pass


    def set_electron_thermo(self, muC, muV, tempK):
        pass



class NoScreeningStatic(EpsStatic):

    def epsFunc(self, q):
        r = 1.0
        return r



class ThomasFermiStatic(EpsStatic):

    def __init__(self, epst=1.0, epsb=1.0):

        # Physical parameters.
        p = GrapheneSample(epst=epst, epsb=epsb, qcutoff=0.0).p

        self.p = p


    def set_electron_thermo(self, muC, muV, tempK):
        p = self.p
        p["muC"] = muC
        p["muV"] = muV
        p["tempK"] = tempK

        # Calculate the Thomas-Fermi wave vector.

        # Carrier density in the two bands [nm^-2].
        nE = gt.oneband_dens_func(muC, tempK)
        nH = gt.oneband_dens_func(-muV, tempK)
        # Electron density [nm^-2].
        dens = nE - nH
        # Fermi energy [eV].
        p["eF"] = gt.twobands_ef_dens_func(dens)

        # Thomas-Fermi wave vector [nm^-1].
        p["qTF"] = 4.0 * p["alphaee"] * np.abs(p["eF"]) / p["hbarvf"]


    def epsFunc(self, q):
        r = 1.0 + self.p["qTF"] / q
        return r



class RPAStatic(EpsStatic):

    def __init__(self, epst=1.0, epsb=1.0):

        # Physical parameters.
        self.gs = GrapheneSample(epst=epst, epsb=epsb, qcutoff=0.0)
        p = self.gs.p

        # Maximum absolute energy of the conical band dispersion [eV].
        p["enCutoff"] = 3.0

        # Energy mesh for the Maldague integration [eV].
        p["ee"] = np.linspace(0.0, p["enCutoff"], 101)
        # We found robust results with maximum value set to
        # 10 * eFermi and 301 points.
        # Use extent 10 * (eFermi + temp) if the temperature is much larger
        # than eFermi.

        # Wavevector mesh for the interpolation of the polarization function.
        p["qq"] = np.linspace(0.0, 3.0, 101)

        self.p = p

        # Interpolation function, initialized to a constant.
        self._chi0Func = lambda x : 0.0 * x


    # Redefine the interpolation mesh.
    def interpMeshInit(self, qMax, qNum):
        self.p["qq"] = np.linspace(0.0, qMax, qNum)


    # Redefine the Maldague mesh.
    def maldagueMeshInit(self, eMax, eNum):
        self.p["ee"] = np.linspace(0.0, eMax, eNum)


    # Static dielectric function.
    def epsFunc(self, q):
        r = 1.0 - self.gs.coulomb(q) * self.chi0Func(q)
        return r


    # Static polarization function.
    def chi0Func(self, q):
        return self._chi0Func(q)


    # Calculate the polarization function on a mesh, interpolate and then
    # initialize the interpolation function.
    def set_electron_thermo(self, muC, muV, tempK):
        p = self.p
        p["muC"] = muC
        p["muV"] = muV
        p["tempK"] = tempK

        # Calculate the polarization function on the mesh of wave vectors.
        chi = np.array([self.chi0Thermo(q, muC, muV, tempK) for q in p["qq"]])

        # Interpolate the polarization function.
        self._chi0Func = scipy.interpolate.interp1d(
            p["qq"], chi, kind='linear')


    # Calculate the static polarization function at finite chemical potentials
    # and temperature.
    # Perform the Maldague integral.
    # Define the function with arguments for the chemical potentials and
    # temperature, so it might be called directly if neeeded.
    def chi0Thermo(self, q, muC, muV, tempK):
        p = self.p
        temp = tempK * p["kB"]
        de = p["ee"][1] - p["ee"][0]
        if (temp < de):
            # If the temperature is too small compared to the energy mesh
            # step, the integration becomes the summation of two Dirac
            # delta functions.
            # To obtain continuous results as the temperature is varied,
            # the Delta functions are evaluated on the point of the energy
            # mesh which is closer to the desired value.
            if (muC > 0.0):
                ie0 = np.argmin(np.abs(p["ee"] - muC))
                y1 = self.chi0Zero(q, p["ee"][ie0])
            else:
                y1 = 0.0
            if (-muV > 0.0):
                ie0 = np.argmin(np.abs(p["ee"] + muV))
                y2 = self.chi0Zero(q, p["ee"][ie0])
            else:
                y2 = 0.0
            chi = y1 + y2
        else:
            # Value of the correlation function at zero temperature on the
            # Maldague mesh.
            chi0s = 0.0 * p["ee"]
            for ie,e in enumerate(p["ee"]):
                chi0s[ie] = self.chi0Zero(q, e)
            # Integration coefficients on the Maldague mesh.
            coeff = 0.0 * p["ee"]
            # Calculate the coefficient on the energy mesh.
            for ie,e in enumerate(p["ee"]):
                x1 = (e - muC)/(2.0 * temp)
                x2 = (e + muV)/(2.0 * temp)
                if (np.abs(x1) < 10.0):
                    y1 = 1.0 / (4.0 * temp * np.power(np.cosh(x1), 2.0))
                else:
                    y1 = 0.0
                if (np.abs(x2) < 10.0):
                    y2 = 1.0 / (4.0 * temp * np.power(np.cosh(x2), 2.0))
                else:
                    y2 = 0.0
                coeff[ie] = y1 + y2
            # Perform the integration.
            chi = scipy.integrate.simps(coeff * chi0s, p["ee"])
        # Add the residual contribution.
        y1 = gt.fermi_dirac_func(0.0, muC, tempK)
        y2 = gt.fermi_dirac_func(0.0, muV, tempK)
        chi = chi - (y1 - y2) * self.chi0Zero(q, 0.0)
        return chi


    # Non-interacting polarization, static limit,
    # at zero temperature and no photoexcited density.
    # See Eq.(2.16) in:
    # Kotov et al., Rev. Mod. Phys. 84, 1067 (2012)

    def chi0Zero(self, q, eF):
        p = self.p
        kF = np.abs(eF) / p["hbarvf"]
        kF2 = 2.0 * kF
        chi = - kF2 / (np.pi * p["hbarvf"])
        if (q > kF2):
            r = kF2 / q
            chi = chi + q / (2.0 * np.pi * p["hbarvf"]) * (
                r * np.sqrt(1.0 - r*r) + np.arcsin(r) - 0.5 * np.pi)
        return chi
