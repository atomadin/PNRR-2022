#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import pickle
import numpy as np
import pandas as pd
from scipy import signal

from . import graphene_thermodynamics_v2 as gt
from .under_pumping import UnderPumping
from .cooling_phonons_v2 import CoolingPhonons
from .graphene_optics import GrapheneOptics


# Calculate the steady state at the end of the pump with the module
# under_pumping, and then calculate the relaxation due to phonon emission
# with the module cooling_phonons_v2.


class EvasExperiment:

    def __init__(self):
        pass
        

    def to_pickle(self, fileName):
        d = {
            "p1": self.p1,
            "p2": self.p2,
            "r2": self.r2,
            "steady_state": self.steady_state,
            "dynamics": self.dynamics
            }
        with open(fileName, "wb") as f:
            pickle.dump(d, f)        


    def read_pickle(self, fileName):
        with open(fileName, "rb") as f:
            d = pickle.load(f)
        self.p1 = d["p1"]
        self.p2 = d["p2"]
        self.r2 = d["r2"]
        self.steady_state = d["steady_state"]
        self.dynamics = d["dynamics"]


    def set_parameters(
        self,
        # Fluence [uJ / cm^2].
        fluenceSI=20.0,
        # Pulse duration [fs].
        dtp=150.0,
        # Lattice temperature [K].
        tlatt=300.0,
        # Fermi energy [eV].
        eF=0.200,
        # Photon energy [eV].
        eph=0.8,
        # Temperature relaxation time [fs].
        tau=100.0,
        # Residual absorption.
        alpha_res=0.0,
        # Substrate refraction index.
        nsub=1.44,  # fused silica
        ntop=1.42,   # ionic liquid
        # Optical phonons lifetime [fs].
        tau_ph=1200.0,
        # Parameter that determines the EPCs [eV/nm].
        dtdb=45.0,
        # Parameter that fixes the phonon density [eV].
        # If None, it is estimated in time.
        enmax=None,
        # Integration time [fs].
        tmax=15.0e+3,
        # Number of stroboscopic times.
        tnum=201,
        # Approximate time step [fs].
        dt_approx=10.0):

        # Parameters of the dynamics under pumping.
        self.p1 = {
            "fluenceSI": fluenceSI,
            "dtp": dtp,
            "tlatt": tlatt,
            "eF": eF,
            "eph": eph,
            "tau": tau,
            "alpha_res": alpha_res,
            "nsub": nsub,
            "ntop": ntop
        }

        # Parameters of the cooling to phonons.
        self.p2 = {
            "tau_ph": tau_ph,
            "dtdb": dtdb,
            "enmax": enmax
        }

        # Integration parameters.
        self.r2 = {
            "tmax": 15.0e+3,
            "tnum": 201,
            "dt_approx": 10.0,
            "print_time": False
        }
        

    def solve_dynamics(self, verbose=False):
        # Solve the dynamics under pumping.
        up = UnderPumping(**self.p1)
        if verbose:
            print("pump_probe/solve_dynamics: Calculating steady state.")
        up.calculate_steady_state(
            tempKMin=200.0, tempKMax=2000.0,
            dnEMin=0.000, dnEMax=0.100)
        if verbose:
            print("pump_probe/solve_dynamics: Done steady state.")
            print("pump_probe/solve_dynamics: tempK = %f, dnE = %f" % (up.steady_state["tempK"], up.steady_state["dnE"]))
        
        # Solve the cooling dynamics.
        # Results from the previous integration are needed as parameters to
        # define the initial state of the system.
        p2_init = {
            "eF": self.p1["eF"],
            "tempK_eq": self.p1["tlatt"],
            "tempK_0": up.steady_state["tempK"],
            "dnE_0": 0.0
        }
        p2_tot = {
            "tau_ph": self.p2["tau_ph"],
            **p2_init
        }

        co = CoolingPhonons(
            dtdb=self.p2["dtdb"],
            enmax=self.p2["enmax"])
        co.set_params(**p2_tot)
        co.run(**self.r2)

        # Save results to instance variables.
        self.steady_state = up.steady_state.copy()
        self.dynamics = co.dynamics.copy(deep=True)

        
    def calculate_thermo(self):
        # Return the thermodynamic variables in time.
        tt = self.dynamics["t"].to_numpy()
        tempK = self.dynamics["tempK"].to_numpy()
        muC = self.dynamics["muC"].to_numpy()
        muV = self.dynamics["muV"].to_numpy()
        return (tt, tempK, muC, muV)


    def calculate_signal(self, probe_en_list):
        # Object to perform optical calculations.
        go = GrapheneOptics()
        go.set_refractive_indices(
            nsub=self.p1["nsub"], ntop=self.p1["ntop"])

        # Calculate the transmittance in time.
        transm_t = []
        # Iterate over the times.
        for i_t,val_t in self.dynamics.iterrows():
            # Calculate the chemical potentials.
            dens = gt.twobands_dens_ef_func(self.p1["eF"])
            muC, muV = gt.photoexc_mu_func(dens, val_t["dnE"], val_t["tempK"])
            # Calculate the transmittance.
            go.set_electron_thermo(muC, muV, val_t["tempK"])
            transm_t.append([
                go.transmittance(probe_en) for probe_en in probe_en_list])
        transm_t = np.array(transm_t)
        n_time, n_probe = transm_t.shape

        # Calculate the transmittance at equilibrium.
        # Calculate the chemical potentials.
        dens_0 = gt.twobands_dens_ef_func(self.p1["eF"])
        muC_0, muV_0 = gt.photoexc_mu_func(dens_0, 0.0, self.p1["tlatt"])
        # Calculate the transmittance.  Copy the array for each time.
        go.set_electron_thermo(muC_0, muV_0, self.p1["tlatt"])
        transm_0 = np.tile(np.array([
            go.transmittance(probe_en) for probe_en in probe_en_list]),
            [n_time,1])

        # Calculate the differential transmission.
        # Use the definition, not the approximate relation to the conductivity.
        diff_transm_t = (transm_t - transm_0) / transm_0

        # Convolve the differential transmission with a Gaussian pulse.
        tt = self.dynamics["t"].to_numpy()
        gaussian = np.exp(-np.power((tt - 0.5 * tt[-1]) / self.p1["dtp"], 2.0))
        bleaching_t = np.zeros((n_time,n_probe))
        for i_probe in np.arange(n_probe):
            bleaching_t[:,i_probe] = (signal.convolve(
                diff_transm_t[:,i_probe], gaussian, mode='same')
                / sum(gaussian))

        return (tt, bleaching_t)

    