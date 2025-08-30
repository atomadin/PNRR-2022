#/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class GrapheneConstants:

    def __init__(self):
        c = {}
        # Reduced Planck's constant [eV fs].
        c["hbar"] = 0.66
        # Planck's constant [eV fs].
        c["h_planck"] = 2.0 * np.pi * c["hbar"]
        # Speed of light [nm/fs].
        c["c"] = 299.0
        # Useful parameter [eV nm].
        c["hc"] = c["h_planck"] * c["c"]
        # Graphene cone slope [eV nm].
        c["hbarvf"] = 0.66  # or 0.6
        # Fermi velocity in graphene [nm / fs].
        c["v_F"] = c["hbarvf"] / c["hbar"]
        # Boltzmann constant [eV / K].
        c["kB"] = 1.0 / 11600.0
        # Spin and valley degeneracy.
        c["Nf"] = 4
        # Electron charge squared [eV nm].
        c["e_sq"] = 1.43
        # Electron charge [sqrt(eV nm)].
        c["e"] = np.sqrt(c["e_sq"])
        # Carbon-Carbon distance [nm].
        c["aCC"] = 0.142
        # Primitive vector length [nm].
        c["a_latt"] = np.sqrt(3.0) * c["aCC"]
        # Primitive cell area [nm^2].
        c["A0"] = 1.5 * np.sqrt(3.0) * c["aCC"] * c["aCC"]
        # Carbon mass [eV fs^2 / nm^2].
        c["mC"] = 1.99e-26 * (1.0e+31/1.6)
        # Universal conductivity [nm / fs].
        c["sigma_uni"] = c["e_sq"] / (4.0 * c["hbar"])
        # Quantum of conductance [nm / fs].
        c["cond_quant"] = c["e_sq"] / (2.0 * np.pi * c["hbar"])

        c["epsSiO"] = 3.9  # if using hbarvf=0.66 this leads to alphaee=0.88

        self.c = c


# The following is an example of how to use the constants defined above.
# Instantiate the GrapheneConstants class and copy the dictionary.
# (The copy operation should allow the garbage collector to release the rest
# of the object, I guess.)
# The GrapheneConstant class is not intended to produce "hot" objects where
# constants are added or changed.  For example, if a constant is changed
# its value is not reflected on the other constants that depend on it.
# The idea is that these constants should be set in the source code, on a
# per-project basis.  If one of these constants has to be changed more than
# once (for example, calculation something for several values of v_F)
# new ad-hoc variables should be defined in the calling code.

class GrapheneSample:

    def __init__(self, epst=1.0, epsb=1.0, qcutoff=0.0):
        p = GrapheneConstants().c.copy()
        # Dielectic constants of top and bottom substrates.
        p["epst"] = epst
        p["epsb"] = epsb
        # Average dielectric constant.
        p["epsAv"] = 0.5 * (epst + epsb)
        # Dimensionless coupling constant.
        p["alphaee"] = p["e_sq"] / (p["epsAv"] * p["hbarvf"])
        # Cutoff for the Coulomb potential [nm^-1].
        p["qcutoff"] = qcutoff
        p["qcutoff_sq"] = qcutoff * qcutoff
        # Coefficient of the energy in the DOS.  Includes spin and valley.
        p["dos_coeff"] = 2.0 * 2.0 / (2.0 * np.pi * p["hbarvf"] * p["hbarvf"])
        self.p = p


    # 2D Fourier transform of the Coulomb potential on the graphene sheet.
    def coulomb(self, q):
        r = (2.0 * np.pi * (self.p["e_sq"] / self.p["epsAv"])
             / np.sqrt(q*q + self.p["qcutoff_sq"]))
        return r


    # Density of states at given energy.
    def dos(self, en):
        r = self.p["dos_coeff"] * np.abs(en)
        return r


    # Bottom gate.
    def set_bottom_gate(self, d=100.0):
        # Gate layer thickness [nm].
        self.p["gate_d"] = d
        # Capacitance per unit area [nm^-1].
        self.p["gate_cap"] = self.p["epsb"] / (4.0 * np.pi * d)

