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
from scipy.integrate import simps

from . import graphene_thermodynamics_v2 as gt
from .graphene_constants import GrapheneConstants
from .graphene_optics import GrapheneOptics
from .graphene_optics import n_r
from .graphene_optics import GrapheneTHG

# Equations of motion for:
# - electron temperature
# - photoexcited density
# - phonon mode occupation
# taking into account phonon emission and an initial pump pulse.
# Valid for long pump pulses, such that the dynamics of interest is much
# slower than the time-scale for electron thermalization.

class CoolingPhonons:

    def __init__(self,
                 # Phonon parameters.
                 tempK_eq = 300.0, dtdb=45.0, enmax=None, tau_ph=1200.0,
                 # Pump pulse parameters.
                 pu_fluenceSI=20.0, pu_dt=150.0, pu_eph=0.8,
                 # Probe pulse parameters.
                 pr_delay=1000.0, pr_fluenceSI=1.0, pr_dt=150.0, pr_eph=0.3,
                 # Optical parameters.
                 nsub=n_r["silica"], ntop=n_r["ionic"], alpha_res=0.001,
                 # Calculation of the THG.
                 gammaConst=None, gammaProp=None, gammaInv=None,
                 # Sample parameters.
                 eF=0.200,
                 # Parameter tweaks.
                 ph_en_frac=1.0,  # factor of phonon energy
                 abs_frac=1.0,  # factor of phonon absorption
                 dtdb_cutoff_dens=None,  # cutoff of electron-phono coupling
                 thg_no_pe=False,  # calculate THG without photoexcitation
                 pe_relax_rate=0.0  # extra relaxation of photoexcited density
                ):
        
        # Physical parameters.
        p = GrapheneConstants().c.copy()

        # Phonons.
        # --------

        # Equilibrium temperature [K].
        p["tempK_eq"] = tempK_eq
        # Derivative of the hopping integral wrt the Carbon-Carbon
        # distance [eV / nm].
        p["dt/db"] = dtdb  # default: 45.0
        # Cutoff energy for phonon emission [eV].
        # If None, it is estimated based on the chemical potential
        # and electron temperature.
        p["enmax"] = enmax  # default: 0.8
        # Optical phonons lifetime [fs].
        p["tau_ph"] = tau_ph
        # Phenomenological factor of phonon energies.
        p["ph_en_frac"] = ph_en_frac
        # Phenomenological factor of phonon absorption.
        p["abs_frac"] = abs_frac
        # Phenomenological density cutoff of electron-phonon coupling.
        p["dtdb_cutoff_dens"] = dtdb_cutoff_dens
            
        # Mass density of the graphene sheet [eV fs^2 / nm^4].
        # p["rho"] = (7.6 * 10e-7) * 0.63 * 1.0e+13
        p["rho"] = 2.0 * p["mC"] / p["A0"]
        # Optical phonon energy at the Gamma point [eV].
        p["en_G"] = p["hc"] * 1580.0 * 1.0e-7 * p["ph_en_frac"]
        p["w_G"] = p["en_G"] / p["hbar"]
        # Optical phonon energy at the K point [eV].
        p["en_K"] = p["hc"] * 1300.0 * 1.0e-7 * p["ph_en_frac"]
        p["w_K"] = p["en_K"] / p["hbar"]
        
        # EPC coefficients in the form as in Piscanec PRL 2004 [eV^2].
        p["<g_G>^2"] = (
            p["hbar"] / (2.0 * p["mC"] * p["w_G"])
            * (9.0 / 4.0) * p["dt/db"]**2)
        p["<g_K>^2"] = (
            p["hbar"] / (2.0 * p["mC"] * p["w_K"])
            * (9.0 / 2.0) * p["dt/db"]**2)

        # Upper integration energy [eV].
        p["en_inf"] = 3.0

                
        # Derived coefficient for the recombination rates.
        w = {}
        w["cr_G"] = (
            9.0 / (np.pi * p["rho"] * p["w_G"] * p["hbarvf"]**4) * p["dt/db"]**2)
        w["cr_K"] = (
            9.0 / (np.pi * p["rho"] * p["w_K"] * p["hbarvf"]**4) * p["dt/db"]**2)
        
        # Integration meshes, and shifted meshes, for intra- and inter-band
        # processes [eV].
        dE_approx = 0.005  # approximate integration step [eV]

        ne = int(np.ceil(p["en_inf"] - p["en_G"]) / dE_approx)
        w["enIntraG"] = np.linspace(p["en_G"], p["en_inf"], ne)
        w["enIntraG_shift"] = w["enIntraG"] - p["en_G"]  # e - hw

        ne = int(np.ceil(p["en_G"]) / dE_approx)
        w["enInterG"] = np.linspace(0.0, p["en_G"], ne)
        w["enInterG_shift"] = p["en_G"] - w["enInterG"]  # hw - e
        
        ne = int(np.ceil(p["en_inf"] - p["en_K"]) / dE_approx)
        w["enIntraK"] = np.linspace(p["en_K"], p["en_inf"], ne)
        w["enIntraK_shift"] = w["enIntraK"] - p["en_K"]  # e - hw

        ne = int(np.ceil(p["en_K"]) / dE_approx)
        w["enInterK"] = np.linspace(0.0, p["en_K"], ne)
        w["enInterK_shift"] = p["en_K"] - w["enInterK"]  # hw - e
        
        w["mesh_labels"] = [
            "enIntraG", "enIntraG_shift", "enInterG", "enInterG_shift",
            "enIntraK", "enIntraK_shift", "enInterK", "enInterK_shift"]
        
        # Equilibrium phonon mode occupation.
        w["nG_eq"] = gt.bose_einstein_func(
            e=p["en_G"], mu=0.0, tempK=p["tempK_eq"])
        w["nK_eq"] = gt.bose_einstein_func(
            e=p["en_K"], mu=0.0, tempK=p["tempK_eq"])        

        
        # Pump pulse.
        # -----------

        # Fluence [uJ / cm^2].
        p["pu_fluenceSI"] = pu_fluenceSI
        # Fluence [eV / nm^2].
        p["pu_fluence"] = self.unit_convert_fluence_from_SI(pu_fluenceSI)
        # Pulse duration [fs].
        p["pu_dt"] = pu_dt
        # To smooth the pulse, assume a ramp-up and ramp-down of 10%.
        p["pu_ra"] = 0.1 * pu_dt
        # Photon energy [eV].
        p["pu_eph"] = pu_eph
        # Pump power density, calculated as fluence divided by duration.
        p["pu_powdens"] = p["pu_fluence"] / p["pu_dt"]

        
        # Probe pulse.
        # ------------
        
        # Delay [fs].
        p["pr_delay"] = pr_delay
        # Fluence [uJ / cm^2].
        p["pr_fluenceSI"] = pr_fluenceSI
        # Fluence [eV / nm^2].
        p["pr_fluence"] = self.unit_convert_fluence_from_SI(pr_fluenceSI)
        # Pulse duration [fs].
        p["pr_dt"] = pr_dt
        # To smooth the pulse, assume a ramp-up and ramp-down of 10%.
        p["pr_ra"] = 0.1 * pr_dt
        # Photon energy [eV].
        p["pr_eph"] = pr_eph
        # Probe power density, calculated as fluence divided by duration.
        p["pr_powdens"] = p["pr_fluence"] / p["pr_dt"]
        
        
        # Optical parameters of the sample.
        # ---------------------------------
        
        # Refractive indices.  The incident and reflected waves are in the top
        # region and the transmitted wave is in the substrate.
        p["nsub"] = nsub
        p["ntop"] = ntop
        # Residual absorption.
        p["alpha_res"] = alpha_res
        
        # Optics calculations.
        self.go = GrapheneOptics()
        self.go.set_refractive_indices(nsub=p["nsub"], ntop=p["ntop"])

        # Calculation of the THG from the probe pulse.
        p["thg_data"] = {}
        p["thg_data"]["init"] = {
            "nsub": p["nsub"],
            "ntop": p["ntop"],
            "eph": p["pr_eph"],
            "maldagueMin": 0.001,
            "maldagueMax": 3.0,
            "maldagueNum": 301,
            "gammaConst": gammaConst,
            "gammaProp": gammaProp,
            "gammaInv": gammaInv
        }
        self.thg = GrapheneTHG(**p["thg_data"]["init"])
        
        
        # Electronic parameters of the sample.
        # ------------------------------------
        
        # Fermi energy [eV].
        p["eF"] = eF
        # Electron density measured from the charge neutrality point.
        p["dens"] = gt.twobands_dens_ef_func(p["eF"])
        
        # Phenomenological relaxation rate of the photoexcited density.
        p["pe_relax_rate"] = pe_relax_rate
        
        # THG calculation options.
        # ------------------------
        
        # Ignore the photoexcited density in the calculation of the THG,
        # but keep it in the dynamics.
        p["thg_no_pe"] = thg_no_pe
        
        self.p = p
        self.w = w

        
    # Convert fluence from SI [uJ/cm^-2] to [eV / nm^2].
    def unit_convert_fluence_from_SI(self, fluenceSI):
        return (fluenceSI * 0.062)


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
        # Returns xmin for x smaller than x0, for large x returns x,
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


    # Calculate the optical phonon emission rates [nm^-2 fs^-1].
    def optical_phonon_emission_rates(self, muC, muV, tempK, nG, nK):
        p = self.p
        w = self.w

        # Electron distribution on the meshes.
        ff = {}
        for ml in w["mesh_labels"]:
            ff[ml] = gt.fermi_dirac_func(e=w[ml], mu=muC, tempK=tempK)
        
        # Hole distributionon the meshes.
        hh = {}
        for ml in w["mesh_labels"]:
            hh[ml] = gt.fermi_dirac_func(e=w[ml], mu=-muV, tempK=tempK)

        # Phenomenological reduction of electron-phonon coupling.
        if (p["dtdb_cutoff_dens"] is not None):
            nE = gt.oneband_dens_func(muC, tempK)  # electron density
            nH = gt.oneband_dens_func(muV, tempK)  # hole density
            rC = np.exp(- nE / np.abs(p["dtdb_cutoff_dens"]))
            rV = np.exp(- nH / np.abs(p["dtdb_cutoff_dens"]))
        else:
            rC = 1.0
            rV = 1.0
            
        # Phonon emission rates due to electrons in conduction band.
        e = "enIntraG"
        e_w = "enIntraG_shift"
        R_Gi = w["cr_G"] * scipy.integrate.simps(
            w[e] * w[e_w] * (
                ff[e] * (1.0 - ff[e_w]) * (nG + 1.0) * rC
                - ff[e_w] * (1.0 - ff[e]) * nG * p["abs_frac"] * rC
                + hh[e] * (1.0 - hh[e_w]) * (nG + 1.0) * rV
                - hh[e_w] * (1.0 - hh[e]) * nG * p["abs_frac"] * rV),
            w[e])
        e = "enIntraK"
        e_w = "enIntraK_shift"
        R_Ki = w["cr_K"] * scipy.integrate.simps(
            w[e] * w[e_w] * (
                ff[e] * (1.0 - ff[e_w]) * (nK + 1.0) * rC
                - ff[e_w] * (1.0 - ff[e]) * nK * p["abs_frac"] * rC
                + hh[e] * (1.0 - hh[e_w]) * (nK + 1.0) * rV
                - hh[e_w] * (1.0 - hh[e]) * nK * p["abs_frac"] * rV),
            w[e])
        e = "enInterG"
        w_e = "enInterG_shift"
        R_Gx = w["cr_G"] * scipy.integrate.simps(
            w[e] * w[w_e] * (
                ff[e] * hh[w_e] * (nG + 1.0)
                - (1.0 - ff[e]) * (1.0 - hh[w_e]) * nG * p["abs_frac"]),
            w[e]) * np.sqrt(rC * rV)
        e = "enInterK"
        w_e = "enInterK_shift"
        R_Kx = w["cr_K"] * scipy.integrate.simps(
            w[e] * w[w_e] * (
                ff[e] * hh[w_e] * (nK + 1.0)
                - (1.0 - ff[e]) * (1.0 - hh[w_e]) * nK * p["abs_frac"]),
            w[e]) * np.sqrt(rC * rV)
        
        return (R_Gi, R_Gx, R_Ki, R_Kx)
    
    
    def temp_bose(self, n, en):
        # Calculate the temperature of a bosonic mode from its occupation and energy.
        r = en / (self.p["kB"] * np.log(1.0 + 1.0 / n))
        return r

    
    def calculate_monitor_values(self, dens, dnE, tempK, nG, nK):
        # Chemical potentials [eV].
        muC, muV = gt.photoexc_mu_func(
            dens=dens, dnE=dnE, tempK=tempK)
        # Phonon energy cutoff [eV].
        enmax = self.estimate_enmax(muC, muV, tempK)

        return (enmax, muC, muV)


    # Calculate the time-evolution of the variables.
    def run(self, tmax=1.0e+4, tnum=101, dt_approx=10.0, print_time=True, step="rk4"):
        p = self.p
        w = self.w
        
        # Define the method for the integration step.
        if (step == "euler"):
            self.model_one_step = self.model_one_step_euler
        elif (step == "rk4"):
            self.model_one_step = self.model_one_step_rk4
            
        # Time mesh where the variables are saved.
        tt = np.linspace(0.0, tmax, tnum)
        # Time steps between two save times.
        snum = np.int(np.ceil((tt[1] - tt[0]) / dt_approx))
        dt = (tt[1] - tt[0]) / snum
        # Smaller time step to be used during the pulses.
        snum1 = 10
        dt1 = dt / snum1
        if (print_time):
            print("Using dt = %.3f fs." % dt)
        # Arrays with the results.
        dnE_t = np.zeros(tnum)  # photoexcited density
        tempK_t = np.zeros(tnum)  # temperature
        nG_t = np.zeros(tnum)  # phonon occupation of G mode
        nK_t = np.zeros(tnum)  # phonon occupation of K mode
        # Auxiliary quantities that are monitored in time.
        enmax_t = np.zeros(tnum)
        muC_t = np.zeros(tnum)  # chemical potential in conduction band
        muV_t = np.zeros(tnum)  # chemical potential in valence band
        # THG efficiency, calculated during the probe pulse on a finer mesh.
        pr_tt = []
        pr_thg = []
        pr_thg_args = []
                
        # Variables at initial time.
        # # Exactly zero photoexcited density might lead to exponentially small
        # # densities in valence band at very low temperatures, which then
        # # create problems in the calculation of the chemical potential.
        # # Assume the photoexcited density is always finite - it might just
        # # represent fluctuations due to disorder.
        # dnE_t[0] = 1.0e-6  # nm^-2
        dnE_t[0] = 0.0
        tempK_t[0] = p["tempK_eq"]
        nG_t[0] = w["nG_eq"]
        nK_t[0] = w["nK_eq"]
        # Calculate monitor quantities.
        enmax_t[0], muC_t[0], muV_t[0] = self.calculate_monitor_values(
            p["dens"], dnE_t[0], tempK_t[0], nG_t[0], nK_t[0])
        
        # Iterate over save times.
        for i_save,t_save in enumerate(tt[:-1]):  # save mesh counter
            # Copy the variables from the latest save time.
            dnE_1, tempK_1, nG_1, nK_1 = (
                dnE_t[i_save], tempK_t[i_save], nG_t[i_save], nK_t[i_save])
            # Evolve the variables between save times using snum steps.
            for ls in np.arange(snum):
                t = t_save + ls * dt
                # The step is broken into snum1 smaller steps if we are within
                # a pulse or if a pulse starts before the end of the step.
                if (self.during_pump_pulse(t) or self.during_probe_pulse(t)
                    or self.during_pump_pulse(t + dt) or self.during_probe_pulse(t + dt)):
                    for ls1 in np.arange(snum1):
                        t1 = t + ls1 * dt1
                        dnE_1, tempK_1, nG_1, nK_1 = self.model_one_step(
                            p["dens"], dnE_1, tempK_1, nG_1, nK_1, dt1, t1)
                        # Calculate the THG efficiency during the probe only.
                        if self.during_probe_pulse(t1):
                            if (p["thg_no_pe"]):
                                dnE_thg = 0.0
                            else:
                                dnE_thg = dnE_1
                            enmax_1, muC_1, muV_1 = self.calculate_monitor_values(
                                p["dens"], dnE_thg, tempK_1, nG_1, nK_1)
                            # There is a minor inconsistency in the treatment of the pulses.
                            # In the absorption part of the dynamics, a smooth profile is
                            # is used, with a ramp-up and a higher central part to preserve
                            # the total fluence.  Here instead we use a constant power.
                            # These are just numerical tricks and the final results should
                            # not be affected by these details.
                            thg_args = {
                                "powdens": p["pr_powdens"],
                                "muC": muC_1,
                                "muV": muV_1,
                                "tempK": tempK_1
                            }
                            self.thg.calculate_eta(**thg_args)
                            pr_tt.append(t1)
                            pr_thg.append(self.thg.p["eta"])
                            pr_thg_args.append(thg_args)
                else:
                    dnE_1, tempK_1, nG_1, nK_1 = self.model_one_step(
                        p["dens"], dnE_1, tempK_1, nG_1, nK_1, dt, t)

            # Copy the variables to the next save time.
            dnE_t[i_save+1], tempK_t[i_save+1], nG_t[i_save+1], nK_t[i_save+1] = (
                dnE_1, tempK_1, nG_1, nK_1)
            # Calculate monitor quantities.
            enmax_t[i_save+1], muC_t[i_save+1], muV_t[i_save+1] = self.calculate_monitor_values(
                p["dens"], dnE_t[i_save+1], tempK_t[i_save+1], nG_t[i_save+1], nK_t[i_save+1])
            # Print progress every 10 save times.
            if (print_time):
                if (np.mod(i_save + 1, tnum // 10) == 0):
                    print("Done %d." % (i_save + 1))

        # Define a matrix with the results.
        self.dynamics_m =  np.c_[tt, tempK_t, muC_t, muV_t, dnE_t,
                                 nG_t, nK_t, enmax_t]
        # Define a DataFrame with the results.
        self.dynamics = pd.DataFrame(self.dynamics_m, columns=[
            "t", "tempK", "muC", "muV", "dnE", "nG", "nK", "enmax"])

        # Save the THG efficiency.
        p["pr_tt"] = np.array(pr_tt)
        p["etaTHG"] = np.array(pr_thg)
        p["thg_data"]["time"] = [t for t in pr_tt]
        p["thg_data"]["args"] = pr_thg_args
        p["etaTHG_max"] = p["etaTHG"].max()
        p["etaTHG_avg"] = simps(p["etaTHG"], p["pr_tt"]) / (p["pr_tt"][-1] - p["pr_tt"][0])

        
    # Functions to check if the pump or probe pulse is active.        
    def during_pump_pulse(self, t):
        if (t < self.p["pu_dt"]):
            r = True
        else:
            r = False
        return r

    def during_probe_pulse(self, t):
        if (t > self.p["pr_delay"] and t < (self.p["pr_delay"] + self.p["pr_dt"])):
            r = True
        else:
            r = False
        return r

    # Model intensity profile of the pulses, normalized to have the same
    # integral of the constant unity between tmin and tmax (i.e. tmax-tmin).
    def profile_pulse(self, t, tmin, tmax, tramp):
        if (t < tmin or t > tmax):
            r = 0.0
        elif (t < (tmin + tramp) and t < (tmax + tmin)/2.0):
            r = 1.0 - 0.5 * (np.cos(np.pi * (t - tmin) / tramp) + 1.0)
        elif (t > (tmax - tramp) and t > (tmax + tmin)/2.0):
            r = 1.0 - 0.5 * (np.cos(np.pi * (tmax - t) / tramp) + 1.0)
        else:
            r = 1.0
        s = 1.0 * (tmax - tmin - 2.0 * tramp) + tramp
        r = r * (tmax - tmin) / s
        return r
    
    
    # Calculate the derivative of the variables of the model.
    def model_derivative(self, dens, dnE, tempK, nG, nK, t):
        p = self.p
        w = self.w

        # Heat capacity, calculated for electrons and holes separately.
        cv = gt.photoexc_heat_capacity_func_2(dens=dens, dnE=dnE, tempK=tempK)

        # Chemical potentials.
        muC, muV = gt.photoexc_mu_func(dens=dens, dnE=dnE, tempK=tempK)
       
        # Absorbed power during initial pulse.
        if self.during_pump_pulse(t):
            # Absorbance.
            self.go.set_electron_thermo(muC=muC, muV=muV, tempK=tempK)
            self.go.set_refractive_indices(nsub=p["nsub"], ntop=p["ntop"])
            alpha = self.go.absorbance(eph=p["pu_eph"])
            # Take the absorbance only if positive.
            alpha = alpha if (alpha > 0) else 0.0
            # Add the residual absorbance.
            alpha = alpha + p["alpha_res"]
            # Absorbed power density.
            power_abs_pu = alpha * p["pu_powdens"] * self.profile_pulse(
                t, 0.0, p["pu_dt"], p["pu_ra"])
            # Photoexcited electrons per unit of time.
            d_n_abs_pu = power_abs_pu / p["pu_eph"]
        else:
            power_abs_pu = 0.0
            d_n_abs_pu = 0.0

        # Absorbed power during the second pulse.
        # We know that we will not have more than two pulses in these
        # experiments, so we just repeat the above code twice.
        if self.during_probe_pulse(t):
            # Absorbance.
            self.go.set_electron_thermo(muC=muC, muV=muV, tempK=tempK)
            self.go.set_refractive_indices(nsub=p["nsub"], ntop=p["ntop"])
            alpha = self.go.absorbance(eph=p["pr_eph"])
            # Take the absorbance only if positive.
            alpha = alpha if (alpha > 0) else 0.0
            # Add the residual absorbance.
            alpha = alpha + p["alpha_res"]
            # Absorbed power density.
            power_abs_pr = alpha * p["pr_powdens"] * self.profile_pulse(
                t, p["pr_delay"], p["pr_delay"] + p["pr_dt"], p["pr_ra"])
            # Photoexcited electrons per unit of time.
            d_n_abs_pr = power_abs_pr / p["pr_eph"]
        else:
            power_abs_pr = 0.0
            d_n_abs_pr = 0.0

        # Sum the contributions of the two pulses.
        power_abs = power_abs_pu + power_abs_pr
        d_n_abs = d_n_abs_pu + d_n_abs_pr
    
        # Optical phonon emission rates, intraband and interband.
        R_Gi, R_Gx, R_Ki, R_Kx = self.optical_phonon_emission_rates(muC, muV, tempK, nG, nK)
        R_G = R_Gi + R_Gx
        R_K = R_Ki + R_Kx
        
        # Density of optical phonon modes involved in the dynamics.
        mm_G, mm_K = self.dens_phonon_modes(muC, muV, tempK)

        # Derivatives of the variables.
        d_dnE = d_n_abs - R_Gx - R_Kx - dnE * p["pe_relax_rate"]
        d_tempK = (power_abs - R_G * p["en_G"] - R_K * p["en_K"]) / cv
        d_nG = R_G / mm_G - (nG - w["nG_eq"]) / p["tau_ph"]
        d_nK = R_K / mm_K - (nK - w["nK_eq"]) / p["tau_ph"]

        return (d_dnE, d_tempK, d_nG, d_nK)


    # Evolve the variables by a time dt.  First-order Euler.
    def model_one_step_euler(self, dens, dnE_0, tempK_0, nG_0, nK_0, dt, t):
        x0 = np.array((dnE_0, tempK_0, nG_0, nK_0))
        dx = np.array(self.model_derivative(dens, *x0, t))
        x1 = x0 + dt * dx
        dnE, tempK, nG, nK = x1
        return (dnE, tempK, nG, nK)


    # Evolve the variables by a time dt.  Fourth-order Runge-Kutta.
    def model_one_step_rk4(self, dens, dnE_0, tempK_0, nG_0, nK_0, dt, t):
        dt2 = dt / 2.0
        y_n = np.array((dnE_0, tempK_0, nG_0, nK_0))
        k1 = np.array(self.model_derivative(dens, *(y_n           ), t       ))
        k2 = np.array(self.model_derivative(dens, *(y_n + dt2 * k1), t + dt2 ))
        k3 = np.array(self.model_derivative(dens, *(y_n + dt2 * k2), t + dt2 ))
        k4 = np.array(self.model_derivative(dens, *(y_n + dt  * k3), t + dt  ))
        y_n1 = y_n + (dt/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        dnE, tempK, nG, nK = y_n1
        return (dnE, tempK, nG, nK)
