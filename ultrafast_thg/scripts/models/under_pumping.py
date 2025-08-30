#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy
import scipy.integrate
import scipy.interpolate
import scipy.optimize

from .graphene_constants import GrapheneConstants
from . import graphene_thermodynamics_v2 as gt
from .graphene_optics import GrapheneOptics
from .graphene_optics import n_r


class UnderPumping:

    def __init__(self, fluenceSI=20.0, dtp=150.0, tlatt=300.0, eph=0.8,
                 eF=0.200, tau=200.0, alpha_res=0.001,
                 nsub=n_r["silica"], ntop=n_r["ionic"]):

        # Physical constants.
        p = GrapheneConstants().c.copy()

        # Parameters of the system.

        # Fluence [uJ / cm^2].
        p["fluenceSI"] = fluenceSI
        # Fluence [eV / nm^2].
        p["fluence"] = self.unit_convert_fluence_from_SI(fluenceSI)
        # Pulse duration [fs].
        p["dtp"] = dtp
        # Lattice temperature [K].
        p["tlatt"] = tlatt
        # Photon energy [eV].
        p["eph"] = eph
        # Fermi energy [eV].
        p["eF"] = eF
        # Temperature relaxation time [fs].
        p["tau"] = tau
        # Residual absorption.
        p["alpha_res"] = alpha_res
        # Refractive indices.  The incident and reflected waves are in the top
        # region and the transmitted wave is in the substrate.
        p["nsub"] = nsub
        p["ntop"] = ntop

        # Derived quantities.

        # Electron density measured from the charge neutrality point.
        p["dens"] = gt.twobands_dens_ef_func(p["eF"])

        self.p = p

        # Optics calculations.
        self.go = GrapheneOptics()
        self.go.set_refractive_indices(nsub=p["nsub"], ntop=p["ntop"])

        # Dictionary with the thermodynamic quantities in the steady state.
        self.steady_state = {
            "tempK": None,
            "dnE": None,
            "muC": None,
            "muV": None
        }


    # Calculate thermodynamic quantities in the steady state.
    def calculate_steady_state(self, tempKMin=200.0, tempKMax=2000.0,
                               dnEMin=0.000, dnEMax=0.100):
        p = self.p
        # Temperature.
        tempK = self.steady_state_temperature(
            tempKMin=tempKMin, tempKMax=tempKMax,
            dnEMin=dnEMin, dnEMax=dnEMax)
        # Photoexcited density.
        dnE = self.steady_state_pe_density(
            tempK, dnEMin=dnEMin, dnEMax=dnEMax)
        # Chemical potentials.
        muC, muV = gt.photoexc_mu_func(dens=p["dens"], dnE=dnE, tempK=tempK)
        # Save to state dictionary.
        self.steady_state["tempK"] = tempK
        self.steady_state["dnE"] = dnE
        self.steady_state["muC"] = muC
        self.steady_state["muV"] = muV


    # Steady-state photoexcited density [nm^-2].
    def steady_state_pe_density_root(self, dnE, tempK):
        p = self.p
        # Calculate the absorption coefficient.
        muC, muV = gt.photoexc_mu_func(dens=p["dens"], dnE=dnE, tempK=tempK)
        self.go.set_electron_thermo(muC=muC, muV=muV, tempK=tempK)
        alpha = self.go.absorbance(eph=p["eph"])
        # Take the absorbance only if positive.
        alpha = alpha if (alpha > 0) else 0.0
        # Quantity to nullify.
        # Here we do not add alpha_res to alpha, because the photoexcited
        # density is due to inter-band transitions, while the residual
        # absorption is supposed to represent intra-band transitions.
        r = dnE - (p["tau"] / p["eph"]) * alpha * (p["fluence"] / p["dtp"])
        return r


    def steady_state_pe_density(self, tempK, dnEMin=0.000, dnEMax=0.100):
        dnE_ss = scipy.optimize.bisect(
            lambda dnE: self.steady_state_pe_density_root(dnE, tempK),
            dnEMin, dnEMax)
        return dnE_ss


    # Steady-state temperature [K].
    def steady_state_temperature_root(self, tempK, dnEMin=0.000, dnEMax=0.100):
        p = self.p
        # Calculate the photoexcited density.
        dnE = self.steady_state_pe_density(tempK, dnEMin=dnEMin, dnEMax=dnEMax)
        # Calculate the absorption coefficient.
        muC, muV = gt.photoexc_mu_func(dens=p["dens"], dnE=dnE, tempK=tempK)
        self.go.set_electron_thermo(muC=muC, muV=muV, tempK=tempK)
        alpha = self.go.absorbance(eph=p["eph"])
        # Calculate the heat capacity.
        cv = gt.photoexc_heat_capacity_func(
            dens=p["dens"], dnE=dnE, tempK=tempK)
        # Quantity to nullify.
        r = tempK - (
            p["tlatt"] + p["tau"] * (alpha + p["alpha_res"]) 
            / cv * (p["fluence"] / p["dtp"]))
        return r


    def steady_state_temperature(self, tempKMin=100.0, tempKMax=3000.0,
                                 dnEMin=0.000, dnEMax=0.100):
        tempK_ss = scipy.optimize.bisect(
            lambda tempK: self.steady_state_temperature_root(
                tempK, dnEMin=dnEMin, dnEMax=dnEMax),
            tempKMin, tempKMax)
        return tempK_ss


    # Convert fluence from SI [uJ/cm^-2] to [eV / nm^2].
    def unit_convert_fluence_from_SI(self, fluenceSI):
        return (fluenceSI * 0.062)


