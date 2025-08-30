#/usr/bin/env python3
# -*- coding: utf-8 -*-

# Calculate the transport coefficients using the linearized Boltzmman equation.
# Follow:
# E.H. Hwang, E. Rossi, and S. Das Sarma, PRB 80, 235415 (2009).
# See also:
# J. Ziman, "Electrons and Phonon", Eqs. (7.5.1), (9.9.6)
# For a short summary of the derivation of the coefficients from Boltzmann:
# A. Cantarero and F.X. Alvarez, doi:10.1007/978-3-319-02012-9_1

import numpy as np
import scipy
import scipy.integrate

from .graphene_constants import GrapheneSample
from . import graphene_thermodynamics_v2 as gt
from .transporttime_v2 import CoulombImpuritiesTime


#%% Linear transport coefficients with corrections due to photoexcitation.


class ThermopowerPhotoexcited():

    def __init__(self, epst=1.0, epsb=1.0, dgate=100.0,
                 nimpCM=0.5, dimp=0.0,
                 rel_length_ratio=1.0, en_max=1.5):

        s = GrapheneSample(epst=epst, epsb=epsb)
        # Set the gate layer thickness d [nm].
        s.set_bottom_gate(d=dgate)
        self.s = s
        p = {}
        # Impurity density [10^12 cm^-2].
        p["nimpCM"] = nimpCM
        # Distance between the impurities and the graphene sheet [nm].
        p["dimp"] = dimp
        # Ratio between the decay length of the temperature to
        # photoexcited density.
        p["alpha"] = rel_length_ratio
        # Upper cutoff of the energy mesh for the integration [eV].
        p["en_max"] = en_max
        # Dimensionful scales for the conductivity and the thermopower.
        p["unit_cond"] = s.p["e_sq"] / s.p["h_planck"]
        p["unit_thpo"] = s.p["kB"] / s.p["e"]
        self.p = p
        # Dictionary for the values of the transport coefficients.
        self.val = {}


    def set_electron_distro(self, ef, tempK, dnECM):
        s = self.s
        p = self.p

        p["ef"] = ef
        # Convert the temperature from K to eV.
        p["tempK"] = tempK
        temp = s.p["kB"] * tempK
        p["temp"] = temp
        # Convert photoexcited density from 10^12cm^-2 to nm^-2.
        dnE = 0.01 * dnECM
        p["dnE"] = dnE

        # Equilibrium density at zero temperature [nm^-2].
        p["dens"] = gt.twobands_dens_ef_func(ef=p["ef"])
        # Chemical potentials in the two bands [eV].
        # Notice: these are the chem.pot. for the electrons, not the holes.
        p["muC"], p["muV"] = gt.photoexc_mu_func(
            dens=p["dens"], dnE=p["dnE"], tempK=p["tempK"])
        # Chemical potential of the system at equilibrium, with the same
        # temperature and total density.
        p["mu0"] = gt.twobands_mu_func(dens=p["dens"], tempK=p["tempK"])
        # Density of states at the Fermi energy [eV^-1 nm^-2].
        # Includes valley and spin degeneracy.
        p["nu"] = s.dos(p["ef"])
        # The derivatives of the chemical potentials wrt ef, tempK, dnE.
        p["dmu"] = gt.photoexc_mu_deriv(
            ef=p["ef"], tempK=p["tempK"], dnE=p["dnE"])

        # Calculate the thermopower.
        self.calculate_photoexcited()
        self.calculate_equilibrium()


    def calculate_photoexcited(self):
        s = self.s
        p = self.p

        # On the "theoretical" and "phenomenological" definitions of the
        # the thermopower, see:
        # G.D. Mahan, J. Appl. Phys. 87, 7326 (2000)
        # J. Cai and G.D. Mahan, Phys. Rev. B 74, 075201 (2006)
        #
        # Here we calculate the "theoretical" values, because the
        # "phenomenological" ones are not defined in the presence of
        # photoexcitation.
        #
        # In the function "calculate_equilibrium" below, we calculate the
        # "phenomenological" ones corresponding to a system with the same
        # temperature and Fermi energy, but with a single chemical potential,
        # i.e. no photoexcitation.

        # Parameters for the linear transport coefficients, including
        # corrections due to the photoexcitation.
        # Parameter \tilde{\epsilon} in the conductivity.
        epsTildeCM1 = (
            1.0 + (s.p["gate_cap"] / s.p["e_sq"])
            * p["dmu"]["muC_ef"] / p["nu"])
        epsTildeVM1 = (
            1.0 + (s.p["gate_cap"] / s.p["e_sq"])
            * p["dmu"]["muV_ef"] / p["nu"])
        # Parameter \tilde{\mu} in the thermopower.
        muTildeC = (
            p["muC"] - p["alpha"] * p["dmu"]["muC_dnE"] * p["dnE"]
            - p["dmu"]["muC_tempK"] * p["tempK"])
        muTildeV = (
            p["muV"] - p["alpha"] * p["dmu"]["muV_dnE"] * p["dnE"]
            - p["dmu"]["muV_tempK"] * p["tempK"])

        # Save the coefficients of the bar-red transport coefficients.
        self.barc = {
            "epsTildeCM1": epsTildeCM1,
            "epsTildeVM1": epsTildeVM1,
            "muTildeC": muTildeC,
            "muTildeV": muTildeV
        }

        # Integration over the energy.
        en_num = int(np.floor(5.0 * p["en_max"] / p["temp"]))
        en_mesh = np.logspace(-4, np.log10(p["en_max"]), en_num)

        # Define the derivatives of the electron and hole distributions.
        dfel = np.array([
            gt.fermi_dirac_deriv_func(en, p["muC"], p["tempK"])
            for en in en_mesh])
        dfho = np.array([
            gt.fermi_dirac_deriv_func(en, -p["muV"], p["tempK"])
            for en in en_mesh])

        # Calculate the transport scattering time (due to Coulomb impurities)
        # on the energy mesh.
        # Note: the value is the same for electrons and holes.
        coulImp = CoulombImpuritiesTime(
            nimpCM=p["nimpCM"], dimp=p["dimp"],
            epst=s.p["epst"], epsb=s.p["epsb"])
        coulImp.interpMeshInit(0.0001, p["en_max"], 51)
        coulImp.set_electron_thermo(p["muC"], p["muV"], p["tempK"])
        tau_mesh = np.array([coulImp.tauFunc(en) for en in en_mesh])

        # ? ? ? ? ? ? ? ? ? ? ? ? ?
        # IS TAU_MESH THE SAME
        # FOR ELECTRONS AND HOLES
        # ? ? ? ? ? ? ? ? ? ? ? ? ?

        # Calculate the transport coefficients for the various bands.

        # Dimensionful coefficients.
        cq = self.p["unit_cond"]
        ct = self.p["unit_thpo"]
        c2 = 2.0 / s.p["hbar"]

        # Conductivity in conduction band.
        coeff = - cq * c2 * epsTildeCM1
        integrand = en_mesh * tau_mesh * dfel
        barCondC = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Conductivity in valence band.
        coeff = - cq * c2 * epsTildeVM1
        integrand = en_mesh * tau_mesh * dfho
        barCondV = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Product of conductivity and thermopower in conduction band.
        coeff = ct * cq * c2 / p["temp"]
        integrand = en_mesh * tau_mesh * (en_mesh - muTildeC) * dfel
        barProdC = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Product of conductivity and thermopower in valence band.
        coeff = - ct * cq * c2 / p["temp"]
        integrand = en_mesh * tau_mesh * (en_mesh + muTildeV) * dfho
        barProdV = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Final expressions for the transport coefficients.
        barCond = barCondC + barCondV
        barThPo = (barProdC + barProdV) / barCond

        self.val["barCond"] = barCond
        self.val["barThPo"] = barThPo


    def calculate_equilibrium(self):
        s = self.s
        p = self.p

        # Integration over the energy.
        en_num = int(np.floor(5.0 * p["en_max"] / p ["temp"]))
        en_mesh = np.logspace(-4, np.log10(p["en_max"]), en_num)

        # Define the derivatives of the electron and hole distributions.
        dfel = np.array([
            gt.fermi_dirac_deriv_func(en, p["muC"], p["tempK"])
            for en in en_mesh])
        dfho = np.array([
            gt.fermi_dirac_deriv_func(en, -p["muV"], p["tempK"])
            for en in en_mesh])

        # Calculate the transport scattering time (due to Coulomb impurities)
        # on the energy mesh.
        # Note: the value is the same for electrons and holes.
        coulImp = CoulombImpuritiesTime(
            nimpCM=p["nimpCM"], dimp=p["dimp"],
            epst=s.p["epst"], epsb=s.p["epsb"])
        coulImp.interpMeshInit(0.0001, p["en_max"], 51)
        coulImp.set_electron_thermo(p["muC"], p["muV"], p["tempK"])
        tau_mesh = np.array([coulImp.tauFunc(en) for en in en_mesh])

        # ? ? ? ? ? ? ? ? ? ? ? ? ?
        # IS TAU_MESH THE SAME
        # FOR ELECTRONS AND HOLES
        # ? ? ? ? ? ? ? ? ? ? ? ? ?

        # Calculate the transport coefficients for the various bands.

        # Dimensionful coefficients.
        cq = self.p["unit_cond"]
        ct = self.p["unit_thpo"]
        c2 = 2.0 / s.p["hbar"]

        # Conductivity in conduction band.
        coeff = - cq * c2
        integrand = en_mesh * tau_mesh * dfel
        condC = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Conductivity in valence band.
        coeff = - cq * c2
        integrand = en_mesh * tau_mesh * dfho
        condV = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Product of conductivity and thermopower in conduction band.
        coeff = ct * cq * c2 / p["temp"]
        integrand = en_mesh * tau_mesh * (en_mesh - p["muC"]) * dfel
        prodC = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Product of conductivity and thermopower in valence band.
        coeff = - ct * cq * c2 / p["temp"]
        integrand = en_mesh * tau_mesh * (en_mesh + p["muV"]) * dfho
        prodV = coeff * scipy.integrate.simps(integrand, en_mesh)

        # Final expressions for the transport coefficients.
        cond = condC + condV
        thpo = (prodC + prodV) / cond

        self.val["Cond"] = cond
        self.val["ThPo"] = thpo




