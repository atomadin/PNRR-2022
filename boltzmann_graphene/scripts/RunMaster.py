#/usr/bin/env python3
# -*- coding: utf-8 -*-

## Set up script.

import os
import sys
import re
import datetime
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt
import json


## Functions to manage the run.

def plot_energy_density(runDir, fileName):
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")
    ee = 100.0/6.242*np.loadtxt("%s/out-evol-energy.csv" % runDir, delimiter=",")
    et = ee[:,0] + ee[:,1]
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    plt.plot(times, ee[:,0], "-r", label="electron")
    plt.plot(times, ee[:,1], "-b", label="phonon")
    plt.plot(times, et, "-k", label="total")
    plt.ylim([0.0, 1.05*np.max(et)])
    plt.xlabel(r"$t~[{\rm fs}]$")
    plt.ylabel(r"${\cal E}~[\mu{\rm J}/{\rm cm}^{2}]$")
    plt.legend(loc="center left", frameon=False, fontsize="small")
    plt.savefig("%s/%s.png" % (runDir, fileName), dpi=300)


def plot_differential_transmission(runDir, fileName):
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")
    wp = np.loadtxt("%s/out-diff-transm-freqs.csv" % runDir, delimiter=",")
    dtt = np.loadtxt("%s/out-evol-diff-transm.csv" % runDir, delimiter=",")
    cmap_name = "gnuplot_r"
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(wp))]
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    for iw,(w,c) in enumerate(zip(wp,mycolors)):
        plt.plot(times, dtt[:,iw] / np.max(dtt[:,iw]), "-", lw=1.0, c=c, label=r"%.0f~{\rm nm}" % (np.round(1240.0/w*0.01)*100.0))
    plt.xlabel(r"$t~[{\rm fs}]$")
    plt.ylabel(r"$\Delta T/T\vert_{\rm norm}$")
    plt.ylim([-0.05, 1.05])
    plt.legend(loc="lower right", frameon=False, fontsize="small", labelspacing=0.3)
    plt.savefig("%s/%s.png" % (runDir, fileName), dpi=300)


def params2namelist(params):
    s = ""
    for nm in params:
        s += "&%s\n" % nm
        for p in params[nm]:
            s += "%s = %s\n" % (p, params[nm][p])
        s += "/\n\n"
    return s


def merge_params(base_params, params_mod):
    # Create a new set of parameters.
    # Do it explicitly to avoid issues with nested copies.
    params = {}
    for nm in base_params:
        params[nm] = {}
        for p in base_params[nm]:
            params[nm][p] = base_params[nm][p]
    # Update the base parameters with the run parameters.
    # Assume that run_params does not contain new keys!
    # base_params must be a valid, complete set of parameters for the run.
    for nm in params_mod:
        for p in params_mod[nm]:
            params[nm][p] = params_mod[nm][p]
    return params


def runlist2paramlist(run_list):
    # Check that the same number of run values is given for all parameters.
    run_num = len(run_list[0][2])
    for s in run_list:
        # s[0] is the label of the namelist
        # s[1] is the label of the parameter
        # s[2] is the list of the values of the parameter, one for each run
        if (len(s[2]) != run_num):
            sys.exit("runlist2paramlist: List of run parameters not consistent.")
    # Create a list of dictionaries with the correct structure.
    param_list = []
    for ir in range(run_num):
        # Create an empty dictionary of parameters for this run.
        params = {}
        # Iterate over all the parameters.
        for s in run_list:
            # Create the dictionary corresponding to the namelist if necessary.
            if s[0] not in params:
                params[s[0]] = {}
            # Make sure that the parameter is saved as a string.
            params[s[0]][s[1]] = str(s[2][ir])
        param_list.append(params)
    return param_list


def execute_runs(runDir, base_params, run_list=None, param_list=None):
    # If a compact "run list" is given instead of the more flexible
    # "param list", convert the former to the latter.
    if run_list is not None:
        param_list = runlist2paramlist(run_list)
    if (run_list is None and param_list is None):
        sys.exit("execut_runs: Missing run parameters!")
    # Save the param_list to the run folder.
    with open("%s/param_list.txt" % runDir, "w") as f:
        f.write(json.dumps(param_list, indent=4))
    # Iterate over the runs.
    for ir, params_mod in enumerate(param_list):
        # Produce the parameters for this run.
        params = merge_params(base_params, params_mod)
        # Produce the namelist for the input file.
        input_txt = params2namelist(params)
        # Write the input file.
        with open("%s/input.txt" % runDir, "w") as f:
            f.write(input_txt)
        # Execute the run.
        print("\n\nRUN NUMBER %02d\n-------------\n\n" % ir)
        os.system("%s/boltzmann_graphene.x" % runDir)
        # Create relevant plots.
        plot_energy_density(runDir, "plot-energy-density-%02d" % ir)
        plot_differential_transmission(runDir, "plot-diff-transm-%02d" % ir)
        # Copy run data to a folder.
        h = datetime.datetime.now()
        saveDirName = "%s_%02d" % (h.strftime("%Y-%m-%d_%H%M"), ir)
        os.system("cd %s && mkdir %s" % (runDir, saveDirName))
        os.system("cd %s && mv out* input.txt %s" % (runDir, saveDirName))


## Default set of parameters.

eva_1 = {
    "band_structure_params": {
        "factor_wavevector_mesh": "13",
        "define_symmetric_bz": "T",
        "min_level_num_per_bin": "3",
        "absorption_M_point_eV": "4.65"
    },
    "phonon_modes_params": {
        "phonon_wavevector_cutoff_use": "T",
        "phonon_wavevector_cutoff_ratio": "0.5"
    },
    "time_evolution_params": {
        "time_max_fs": "501.0",
        "desired_timestep_fs": "5.0",
        "num_save_times": "51",
        "snapshot_times_fs": "10.0 50.0 300. 500.0",
        "snapshot_all_times": "F"
    },
    "electron_distribution_params": {
        "initial_state": "3",
        "fermi_energy_eV": "0.200",
        "temperature_eq_K": "300.0"
    },
    "distro_fermi_dirac_params": {
        "temperature_K": "300.0"
    },
    "distro_hot_electrons_params": {
        "photoexc_density_cm": "5.0",
        "temperature_K": "2000.0"
    },
    "distro_photoexcited_params": {
        "pump_wavelength_nm": "265.0",
        "pump_broadening_eV": "0.150",
        "photoexc_density_cm": "0.5",
        "use_fluence": "T",
        "abs_fluence_uJ_cm": "0.3"
    },
    "phonon_distribution_params": {
        "temperature_ph_K": "300.0",
        "substrate_relaxation_fs": "1200.0"
    },
    "coulomb_scattering_params": {
        "relative_eps_top": "1.0",
        "relative_eps_bottom": "2.1",
        "maximum_scatt_wavevector_nm_m_1": "5.0",
        "maximum_scatt_energy_eV": "2.0",
        "neglect_electronic_screening": "F",
        "broadening_factor": "1.1",
        "broadening_is_constant": "F"
    },
    "electron_phonon_scatt_params": {
        "broadening_eV": "0.032",
        "broadening_use_level_spacing": "F",
        "broadening_local_use": "T",
        "broadening_local_factor": "0.05",
        "electron_wavevec_cutoff_use": "T",
        "electron_wavevec_cutoff_ratio": "0.7"
    },
    "diff_transm_params": {
        "probe_wavelengths_nm": "500.0 600.0 700.0 900.0 1200.0",
        "broadening_eV": "0.050"
    }
}


## Run lists.

scan_photoexc = [
    ["distro_photoexcited_params", "photoexc_density_cm", ["0.1", "0.3", "0.5", "1.0"]]
]

scan_probe_freq = [
    ["distro_photoexcited_params", "pump_wavelength_nm", ["260.0", "280.0", "300.0", "350.0" ]]
]

scan_broadening_coulomb = [
    ["coulomb_scattering_params", "broadening_factor", ["0.1", "0.5", "1.0", "2.0" ]]
]

scan_timestep = [
    ["time_evolution_params", "desired_timestep_fs", ["5.0", "8.0", "11.0", "14.0" ]]
]

scan_broadening_phonons = [
    ["electron_phonon_scatt_params", "broadening_eV", ["0.005", "0.010", "0.020", "0.050" ]]
]

timestep_and_broadening = [
    {
        "time_evolution_params": {
            "desired_timestep_fs": "1.0"
        },
        "coulomb_scattering_params": {
            "broadening_is_constant": "T"
        }
    },
    {
        "time_evolution_params": {
            "desired_timestep_fs": "1.0"
        },
        "coulomb_scattering_params": {
            "broadening_is_constant": "F"
        }
    },
    {
        "time_evolution_params": {
            "desired_timestep_fs": "5.0"
        },
        "coulomb_scattering_params": {
            "broadening_is_constant": "T"
        }
    },
    {
        "time_evolution_params": {
            "desired_timestep_fs": "5.0"
        },
        "coulomb_scattering_params": {
            "broadening_is_constant": "F"
        }
    }
]

scan_mesh_size = [
    ["band_structure_params", "factor_wavevector_mesh", ["6", "7", "8", "9" ]]
]

mesh_size_and_phonon_broadening = [
    {
        "band_structure_params": {
            "factor_wavevector_mesh": "9"
        },
        "electron_phonon_scatt_params": {
            "broadening_eV": "0.020"
        }
    },
    {
        "band_structure_params": {
            "factor_wavevector_mesh": "9"
        },
        "electron_phonon_scatt_params": {
            "broadening_eV": "0.010"
        }
    },
    {
        "band_structure_params": {
            "factor_wavevector_mesh": "11"
        },
        "electron_phonon_scatt_params": {
            "broadening_eV": "0.020"
        }
    },
    {
        "band_structure_params": {
            "factor_wavevector_mesh": "11"
        },
        "electron_phonon_scatt_params": {
            "broadening_eV": "0.010"
        }
    }
]

## Execute task.

execute_runs(".", eva_1, run_list=scan_mesh_size)

#execute_runs(".", eva_1, param_list=mesh_size_and_phonon_broadening)

