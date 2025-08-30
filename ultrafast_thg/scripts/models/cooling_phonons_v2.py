#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Temperature and phonon dynamics.
# Equations of motion as in:
# H. Wang et al., Appl. Phys. Lett. 96, 081917 (2010).
# For the electron-phonon emission rates, see also:
# F. Rana et al., Phys. Rev. B 79, 115447 (2009);


import numpy as np
import pandas as pd
import scipy
from scipy.interpolate import interp1d

from . import graphene_thermodynamics_v2 as gt
from .graphene_constants import GrapheneConstants


class CoolingPhonons:

    def __init__(self, dtdb=45.0, enmax=None):

        # Physical parameters.
        p = GrapheneConstants().c.copy()

        # Coefficients for the recombination rates.
        # Derivative of the hopping integral wrt the Carbon-Carbon
        # distance [eV / nm].
        p["dt/db"] = dtdb  # default: 45.0
        # Mass density of the graphene sheet [eV fs^2 / nm^4].
        # p["rho"] = (7.6 * 10e-7) * 0.63 * 1.0e+13
        p["rho"] = 2.0 * p["mC"] / p["A0"]
        # Optical phonon energy at the Gamma point [eV].
        p["en_G"] = p["hc"] * 1580.0 * 1.0e-7
        p["w_G"] = p["en_G"] / p["hbar"]
        # Optical phonon energy at the K point [eV].
        p["en_K"] = p["hc"] * 1300.0 * 1.0e-7
        p["w_K"] = p["en_K"] / p["hbar"]
        
        # Cutoff energy for phonon emission [eV].
        # If None, it is estimated based on the chemical potential
        # and electron temperature.
        p["enmax"] = enmax  # default: 0.8
        
        # EPC coefficients in the form as in Piscanec PRL 2004 [eV^2].
        p["<g_G>^2"] = (
            p["hbar"] / (2.0 * p["mC"] * p["w_G"])
            * (9.0 / 4.0) * p["dt/db"]**2)
        p["<g_K>^2"] = (
            p["hbar"] / (2.0 * p["mC"] * p["w_K"])
            * (9.0 / 2.0) * p["dt/db"]**2)

        # Upper integration energy [eV].
        p["en_inf"] = 3.0

        self.p = p
                
        # Derived coefficient for the recombination rates.
        w = {}
        w["cr_G"] = (
            9.0 / (np.pi * p["rho"] * p["w_G"] * p["hbarvf"]**4) * p["dt/db"]**2)
        w["cr_K"] = (
            9.0 / (np.pi * p["rho"] * p["w_K"] * p["hbarvf"]**4) * p["dt/db"]**2)
        
        # Integration meshes, and shifted meshes [eV].
        dE_approx = 0.010  # approximate integration step [eV]
        neG = int(np.ceil(p["en_inf"] - p["en_G"]) / dE_approx))
        w["enMeshG"] = np.linspace(p["en_G"], p["en_inf"], neG)
        w["enMeshG_shift"] = w["enMeshG"] - p["en_G"]
        neK = int(np.ceil(p["en_inf"] - p["en_K"]) / dE_approx))
        w["enMeshK"] = np.linspace(p["en_K"], p["en_inf"], neK)
        w["enMeshK_shift"] = w["enMeshK"] - p["en_K"]

        self.w = w
                

    def estimate_enmax(self, muC, muV, tempK):
        # Calculate phonon energy cutoff.
        p = self.p
        
        if (p["enmax"] is None):
            # "Typical" maximal electron energy [eV].
            temp = p["kB"] * tempK
            mu_eff = np.sqrt((muC**2 + muV**2)/2.0)
            ek_max = mu_eff + 2.0 * temp  # phenomenological estimate
            enmax = 2.0 * ek_max
        else:
            enmax = p["enmax"]

        return enmax

    
    def minimum_smooth(self, x, xmin, x0):
        # Returns xmin for x smaller than x0, for large x returns s,
        # smoothly joins the two asymptotes.
        # x0 must be less than xmin, or the interpolation is not continuous

        # Raise an error if the phonon density is negative.
        if ((xmin < x0)):
            raise ValueError("cooling_phonons_v2/minimum_smooth: x0 must be less than xmin.")
        
        if (x < x0):
            r = xmin
        else:
            r = np.sqrt((x-x0)**2 + (xmin-x0)**2) + x0
        return r
        
        
    def dens_phonon_modes(self, muC, muV, tempK):
        # Calculate the density [nm^-2] of phonon modes that take part in
        # electronic transitions.
        p = self.p
        
        # Estimate the phonon energy cutoff.
        # We use a formula which is the same for the two modes.
        enmax = self.estimate_enmax(muC, muV, tempK)
        
        # Density of phonon modes [nm^-2].
        mm_G = (
            2.0 / (4.0 * np.pi)
            * ( (enmax / p["hbarvf"])**2 - (p["en_G"] / p["hbarvf"])**2 ))
        mm_K = (
            2.0 / (4.0 * np.pi)
            * ( (enmax / p["hbarvf"])**2 - (p["en_K"] / p["hbarvf"])**2 ))
        
        # Minimum allowed value for the phonon densities.
        mm_G_min = 2.0 / (4.0 * np.pi) * (p["en_G"] / p["hbarvf"])**2
        mm_K_min = 2.0 / (4.0 * np.pi) * (p["en_K"] / p["hbarvf"])**2

        # Smooth between the value calculated above and the minimum.
        mm_G = self.minimum_smooth(mm_G, mm_G_min, 0.5*mm_G_min)
        mm_K = self.minimum_smooth(mm_K, mm_K_min, 0.5*mm_K_min)

        # Raise an error if the phonon density is negative.
        if ((mm_G < 0.0) or (mm_K < 0.0)):
            raise ValueError("cooling_phonons_v2/dens_phonon_modes: Phonon space space is negative.")
        
        return(mm_G, mm_K)

    
    def set_params(self, tau_ph=1200.0, eF=0.200, tempK_eq=300.0,
                   dnE_0=0.0, tempK_0=1000.0,
                   intra_band=True,
                   inter_band=True,
                   phonons_at_tempK_0=False):
        p = self.p
        w = self.w
        
        # Optical phonons lifetime [fs].
        p["tau_ph"] = tau_ph
        # Fermi energy [eV].
        p["eF"] = eF
        # Equilibrium temperature [K].
        p["tempK_eq"] = tempK_eq
        # Photoexcited density [nm^-2].
        p["dnE_0"] = dnE_0
        # Initial temperature [K].
        p["tempK_0"] = tempK_0
        # Which transitions to include.
        # It also affects the calculation of the heat capacity.
        p["intra_band"] = intra_band
        p["inter_band"] = inter_band
        # Whether phonons and electrons are thermalized at t=0.
        p["phonons_at_tempK_0"] = phonons_at_tempK_0

        # Derived parameters.
        # Equilibrium density [nm^-2].
        w["dens"] = gt.twobands_dens_ef_func(ef=p["eF"])

        # Equilibrium phonon mode occupation.
        w["nG_eq"] = gt.bose_einstein_func(
            e=p["en_G"], mu=0.0, tempK=p["tempK_eq"])
        w["nK_eq"] = gt.bose_einstein_func(
            e=p["en_K"], mu=0.0, tempK=p["tempK_eq"])

    
    def calculate_monitor_values(self, dens, dnE, tempK, nG, nK):
        # Chemical potentials [eV].
        muC, muV = gt.photoexc_mu_func(
            dens=dens, dnE=dnE, tempK=tempK)
        # Phonon energy cutoff [eV].
        enmax = self.estimate_enmax(muC, muV, tempK)

        return (enmax, muC, muV)
        

    # Calculate the time-evolution of the variables.
    def run(self, tmax=1.0e+4, tnum=101, dt_approx=10.0, print_time=True):
        p = self.p
        w = self.w
        
        # Time mesh where the variables are saved.
        tt = np.linspace(0.0, tmax, tnum)
        # Time steps between two save times.
        snum = np.int(np.ceil((tt[1] - tt[0]) / dt_approx))
        dt = (tt[1] - tt[0]) / snum
        if (print_time):
            print("Using dt = %.3f fs." % dt)
        # Arrays with the results.
        dnE_t = np.zeros(tnum)
        tempK_t = np.zeros(tnum)
        nG_t = np.zeros(tnum)
        nK_t = np.zeros(tnum)
        # Auxiliary quantities that are monitored in time.
        enmax_t = np.zeros(tnum)
        muC_t = np.zeros(tnum)
        muV_t = np.zeros(tnum)
        
        
        # Variables at initial time.
        dnE_t[0] = p["dnE_0"]
        tempK_t[0] = p["tempK_0"]
        # Initial phonon mode occupation.
        if p["phonons_at_tempK_0"]:
            # At time t=0 electrons and phonons are thermalized.
            nG_t[0] = gt.bose_einstein_func(
                e=p["en_G"], mu=0.0, tempK=p["tempK_0"])
            nK_t[0] = gt.bose_einstein_func(
                e=p["en_K"], mu=0.0, tempK=p["tempK_0"])
        else:
            # At time t=0 phonons are cool, at equilibrium.
            nG_t[0] = w["nG_eq"]
            nK_t[0] = w["nK_eq"]
        # Calculate monitor quantities.
        enmax_t[0], muC_t[0], muV_t[0] = self.calculate_monitor_values(
            w["dens"], dnE_t[0], tempK_t[0], nG_t[0], nK_t[0])
            
        # Calculate variables on the time mesh.
        for it,t in enumerate(tt[:-1]):
            # Copy the variables at this save time.
            dnE_1, tempK_1, nG_1, nK_1 = (
                dnE_t[it], tempK_t[it], nG_t[it], nK_t[it])
            # Evolve the variables between save times.
            for ls in np.arange(snum):
                dnE_1, tempK_1, nG_1, nK_1 = self.model_one_step(
                    w["dens"], dnE_1, tempK_1, nG_1, nK_1, dt)
            # Copy the variables to the next save time.
            dnE_t[it+1], tempK_t[it+1], nG_t[it+1], nK_t[it+1] = (
                dnE_1, tempK_1, nG_1, nK_1)
            # Calculate monitor quantities.
            enmax_t[it+1], muC_t[it+1], muV_t[it+1] = self.calculate_monitor_values(
                w["dens"], dnE_t[it+1], tempK_t[it+1], nG_t[it+1], nK_t[it+1])
            # Print progress every 10 save times.
            if (print_time):
                if (np.mod(it + 1, tnum // 10) == 0):
                    print("Done %d." % (it + 1))
        # Define a matrix with the results.
        self.dynamics_m =  np.c_[tt, tempK_t, nG_t, nK_t, dnE_t, enmax_t, muC_t, muV_t]
        # Define a DataFrame with the results.
        self.dynamics = pd.DataFrame(self.dynamics_m, columns=[
            "t", "tempK", "nG", "nK", "dnE", "enmax", "muC", "muV"])


    # Calculate the optical phonon emission rates [nm^-2 fs^-1].
    def optical_phonon_emission_rates(self, muC, muV, tempK, nG, nK):
        p = self.p
        w = self.w

        # Electron distributions in conduction band.
        fe_meshG = gt.fermi_dirac_func(
            e=w["enMeshG"], mu=muC, tempK=tempK)
        fe_meshG_shift = gt.fermi_dirac_func(
            e=w["enMeshG_shift"], mu=muC, tempK=tempK)
        fe_meshK = gt.fermi_dirac_func(
            e=w["enMeshK"], mu=muC, tempK=tempK)
        fe_meshK_shift = gt.fermi_dirac_func(
            e=w["enMeshK_shift"], mu=muC, tempK=tempK)

        # Hole distributions in valence band.
        fh_meshG = gt.fermi_dirac_func(
            e=w["enMeshG"], mu=-muV, tempK=tempK)
        fh_meshG_shift = gt.fermi_dirac_func(
            e=w["enMeshG_shift"], mu=-muV, tempK=tempK)
        fh_meshK = gt.fermi_dirac_func(
            e=w["enMeshK"], mu=-muV, tempK=tempK)
        fh_meshK_shift = gt.fermi_dirac_func(
            e=w["enMeshK_shift"], mu=-muV, tempK=tempK)

        # Phonon emission rates due to electrons in conduction band.
        R_Ge = w["cr_G"] * scipy.integrate.simps(
            w["enMeshG"] * w["enMeshG_shift"] * (
                fe_meshG * (1.0 - fe_meshG_shift) * (1.0 + nG)
                - fe_meshG_shift * (1.0 - fe_meshG) * nG),
            w["enMeshG"])
        R_Ke = w["cr_K"] * scipy.integrate.simps(
            w["enMeshK"] * w["enMeshK_shift"] * (
                fe_meshK * (1.0 - fe_meshK_shift) * (1.0 + nK)
                - fe_meshK_shift * (1.0 - fe_meshK) * nK),
            w["enMeshK"])

        # Phonon emission rates due to holes in valence band.
        R_Gh = w["cr_G"] * scipy.integrate.simps(
            w["enMeshG"] * w["enMeshG_shift"] * (
                fh_meshG * (1.0 - fh_meshG_shift) * (1.0 + nG)
                - fh_meshG_shift * (1.0 - fh_meshG) * nG),
            w["enMeshG"])
        R_Kh = w["cr_K"] * scipy.integrate.simps(
            w["enMeshK"] * w["enMeshK_shift"] * (
                fh_meshK * (1.0 - fh_meshK_shift) * (1.0 + nK)
                - fh_meshK_shift * (1.0 - fh_meshK) * nK),
            w["enMeshG"])

        # Total emission rates.
        R_G = R_Ge + R_Gh
        R_K = R_Ke + R_Kh
        return (R_G, R_K)


    # Calculate the derivative of the variables of the model.
    def model_derivative(self, dens, dnE, tempK, nG, nK):
        p = self.p
        w = self.w

        # Heat capacity.
        if p["heat_capacity_is_sum"]:
            # Calculate the heat capacity without recombination.
            cv = gt.photoexc_heat_capacity_func_2(dens=dens, dnE=dnE, tempK=tempK)
        else:
            # Calculate the heat capacity allowing for recombination.
            cv = gt.photoexc_heat_capacity_func(dens=dens, dnE=dnE, tempK=tempK)

        # Chemical potentials.
        muC, muV = gt.photoexc_mu_func(dens=dens, dnE=dnE, tempK=tempK)

        # Optical phonon emission rates.
        R_G, R_K = self.optical_phonon_emission_rates(muC, muV, tempK, nG, nK)

        # Density of optical phonon modes involved in the dynamics.
        mm_G, mm_K = self.dens_phonon_modes(muC, muV, tempK)
        
        # Phenomenological relaxation of photoexcited density.
        # d_dnE = -dnE / 200.0
        d_dnE = 0.0 # identically zero, instantaneous relaxation

        # Cooling.
        d_tempK = - (R_G * p["en_G"] + R_K * p["en_K"]) / cv

        # Phonon dissipation.
        d_nG = R_G / mm_G - (nG - w["nG_eq"]) / p["tau_ph"]
        d_nK = R_K / mm_K - (nK - w["nK_eq"]) / p["tau_ph"]
        return (d_dnE, d_tempK, d_nG, d_nK)


    # Evolve the variables by a time dt.  First-order Euler.
    def model_one_step_euler(self, dens, dnE_0, tempK_0, nG_0, nK_0, dt):
        x0 = np.array((dnE_0, tempK_0, nG_0, nK_0))
        dx = np.array(self.model_derivative(dens, *x0))
        x1 = x0 + dt * dx
        dnE, tempK, nG, nK = x1
        return (dnE, tempK, nG, nK)


    # Evolve the variables by a time dt.  Fourth-order Runge-Kutta.
    def model_one_step(self, dens, dnE_0, tempK_0, nG_0, nK_0, dt):
        dt2 = dt / 2.0
        y_n = np.array((dnE_0, tempK_0, nG_0, nK_0))
        k1 = np.array(self.model_derivative(dens, *(y_n           ) ))
        k2 = np.array(self.model_derivative(dens, *(y_n + dt2 * k1) ))
        k3 = np.array(self.model_derivative(dens, *(y_n + dt2 * k2) ))
        k4 = np.array(self.model_derivative(dens, *(y_n + dt  * k3)) )
        y_n1 = y_n + (dt/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        dnE, tempK, nG, nK = y_n1
        return (dnE, tempK, nG, nK)


    # Calculate the temperature of a bosonic mode from its occupation and energy.
    def temp_bose(self, n, en):
        r = en / (self.p["kB"] * np.log(1.0 + 1.0 / n))
        return r
