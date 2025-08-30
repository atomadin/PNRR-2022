#/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import datetime
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt

################################################################################
## UTILITY FUNCTIONS

def plot_cells():
    '''Function to plot the primitive cell and the BZ.'''
    ba = np.loadtxt("%s/out-geo-rec-lat-basis.csv" % runDir, delimiter=",")
    a2 = ba[0,:]
    a6 = ba[1,:]
    a1 = a2 + a6
    a3 = -a6
    a4 = -a1
    a5 = -a2
    cor = np.loadtxt("%s/out-geo-wavevectors-K-K1.csv" % runDir, delimiter=",")
    k1 = cor[0,:]
    k6 = cor[1,:]
    k2 = k6 - a6
    k3 = -k6
    k4 = -k1
    k5 = -k2
    plt.plot([0,a6[0],a1[0],a2[0],0],[0,a6[1],a1[1],a2[1],0],"-c", lw=0.5)
    plt.plot([k1[0],k2[0],k3[0],k4[0],k5[0],k6[0],k1[0]],[k1[1],k2[1],k3 [1],k4[1],k5[1],k6[1],k1[1]],"-m", lw=0.5)


def trunc_f(ff,ib,fmin,fmax):
    '''Truncate values for plotting two bands.'''
    if (ib == 0):  # valence, holes
        rr = np.array([min(max(1.0-f,fmin),fmax) for f in ff])
    else:  # conductance, electrons
        rr = np.array([min(max(f,fmin),fmax) for f in ff])
    return np.log10(rr)


def trunc_g(ff,fmin,fmax):
    '''Truncate to a given interval and calculate the logarithm.'''
    rr = np.array([min(max(f,fmin),fmax) for f in ff])
    return np.log10(rr)


def trunc_log(ff,lmin,lmax):
    '''Truncate and map [10^lmin,10^lmax] to [0,lmax-lmin] for plotting.'''
    return np.array([(np.log10(min(max(f,10.0**lmin),10.0**lmax))-lmin) for f in ff])


################################################################################
## STANDARD ANALYSES (DYNAMICS)

def electron_distro_azimuthal():
    '''Electron distribution, averaged over the angle.'''

    # Save times.
    times = np.loadtxt("%s/out-snap-times.csv" % runDir, delimiter=",")

    # Extreme of the bins for the average distribution.
    bins = np.loadtxt("%s/out-geo-energy-bins.csv" % runDir, delimiter=",")
    # Center point of the bins, for easier plotting.
    bins_plot = 0.5 * (bins[:-1] + bins[1:])

    # Load distribution (averaged over the angle) for several times.
    distro_avg_t = []
    for it,t in enumerate(times):
        distro_avg_t.append(np.loadtxt("%s/out-snap-%04d-el-avg-e.csv" % (runDir, it+1), delimiter=","))

    # Plot the average distribution at three times.
    it_plot = [0,1,-1]
    cmap_name = "gnuplot_r"
    # mycolors = ["red", "blue", "black"]
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(it_plot))]
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    for it,c in zip(it_plot, mycolors):
        t = times[it]
        f = distro_avg_t[it][:,0]
        # Plot separately the curves above and below the Dirac point.
        nb2 = int(len(bins_plot)/2)+1
        plt.plot(bins_plot[:nb2], f[:nb2], "-o", color=c, lw=0.5, mew=None, ms=1.0, label=r"$%.0f~{\rm fs}$" % t)
        plt.plot(bins_plot[-nb2:], f[-nb2:], "-o", color=c, lw=0.5, mew=None, ms=1.0,)
    plt.xlim([-10.0, 10.0])
    plt.ylim([-0.05, 1.05])
    plt.xlabel(r"$\varepsilon~[{\rm eV}]$")
    plt.ylabel(r"$f(\varepsilon)$")
    plt.legend(loc="lower left", fontsize="small", frameon=False, fancybox=False)
    plt.savefig("%s/plot-distro-avg.png" % runDir, dpi=300)


def electron_distro_primitive_cell():
    '''Electron distribution in the primitive cell.'''

    # Save times.
    times = np.loadtxt("%s/out-snap-times.csv" % runDir, delimiter=",")

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Distance between k points.
    dk = np.sqrt(np.vdot(kk[1,:] - kk[0,:],kk[1,:] - kk[0,:]))

    # Load the electron distribution in time.
    distro_t = []
    for it,t in enumerate(times):
        distro_t.append(np.loadtxt("%s/out-snap-%04d-el-distro-k.csv" % (runDir, it+1), delimiter=","))

    # Logarithmic scale for the colorbar.
    yticks = [-5, -4, -3, -2, -1]
    yticklabels = [r"$<\!10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$>\!10^{-1}$"]
    fmin = 10**(yticks[0])
    fmax = 10**(yticks[-1])

    # Plot the distribution in each band in time.
    cmap_name = "gnuplot_r"
    for ib,bl in zip([0,1],["va", "co"]):
        for it,(t,ff) in enumerate(zip(times,distro_t)):
            plt.figure(figsize=(3.5,3.5), frameon=True)
            ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
            # plt.scatter(kk[:,0], kk[:,1], s=4.0, c=ff[:,1], cmap=mpl.colormaps[cmap_name])
            p = mpl.collections.PatchCollection([plt.Circle((k[0],k[1]), radius=0.4*dk, linewidth=0.0) for k in kk], cmap=mpl.colormaps[cmap_name])
            p.set_array(trunc_f(ff[:,ib],ib,fmin,fmax))
            plt.gca().add_collection(p)
            plt.xlim([-1.0, 30.0])
            plt.ylim([-30.0, 30.0])
            plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
            plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
            plt.title(r"\scriptsize $t = %.0f~{\rm fs}$" % t)
            ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
            mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
            ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
            ax2.set_yticklabels(yticklabels)
            plt.savefig("%s/plot-distro-%s-%04d.png" % (runDir, bl, it), dpi=300)


def electron_distro_path():
    '''Electron distribution on a high-symmetry path.'''

    # Save times.
    times = np.loadtxt("%s/out-snap-times.csv" % runDir, delimiter=",")

    # Load the electron distribution in time.
    distro_t = []
    for it,t in enumerate(times):
        distro_t.append(np.loadtxt("%s/out-snap-%04d-el-distro-k.csv" % (runDir, it+1), delimiter=","))

    # Load the indices of the wave vectors in the path.
    path = np.loadtxt("%s/out-geo-path.csv" % runDir, delimiter=",", dtype="int32")

    # Subtract 1 to have indices for numpy arrays instead of Fortran arrays.
    path = path - 1

    n_pri = (len(path) - 1)
    n_K = (1 + n_pri/3) - 1
    n_M = (n_K + n_pri/6) - 1

    lmax = 0  # max exp
    lmin = -5  # min exp

    plt.figure(figsize=(3.5,3.0), frameon=True)
    cmap_name = "gnuplot_r"
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(times))]
    plt.axes([0.2,0.2,0.7,0.7])
    plt.axvline(n_K, c="gray", lw=0.5)
    plt.axvline(n_M, c="gray", lw=0.5)
    for it,(t,c) in enumerate(zip(times,mycolors)):
        # conduction, electrons, plot upwards
        plt.plot(trunc_log(distro_t[it][path,1],lmin,lmax), "-o", color=c, lw=0.5, mew=None, ms=1.0, label=(r"$t=%.0f~{\rm fs}$" % t))
        # valence, holes, plot downwards
        plt.plot(-trunc_log(1.0-distro_t[it][path,0],lmin,lmax), "-o", color=c, lw=0.5, mew=None, ms=1.0)
    plt.xlabel(r"${\bm k}$")
    plt.ylabel(r"$f_{{\rm h},{\bm k}}(t) \quad\quad\quad f_{{\rm e},{\bm k}}(t)$")
    plt.xticks([0, n_K, n_M, n_pri], [r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"])
    plt.xlim([0,n_pri])
    yticks = np.linspace(-(lmax-lmin),lmax-lmin,2*(lmax-lmin)+1)
    yticklabels = [r"$1$"] + [(r"$10^{%d}$" % l) for l in range(lmax-1, lmin, -1)] + [r"$0$" ] + [(r"$10^{%d}$" % l) for l in range(lmin+1, lmax)] + [r"$1$"]
    plt.gca().set_yticks(yticks)
    plt.gca().set_yticklabels(yticklabels)
    plt.ylim([-(lmax-lmin), lmax-lmin])
    plt.legend(loc="upper right", frameon=False, fontsize="small")
    plt.savefig("%s/plot-path-distribution.png" % runDir, dpi=300)


def electron_density_in_time():
    '''Electron density conservation.'''

    # Save times.
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")

    # Load the electron (hole) density in the conduction (valence) band.
    # Convert the density from [nm^-2] to [10^12 cm^-2].
    rr = 100.0*np.loadtxt("%s/out-evol-el-dens.csv" % runDir, delimiter=",")

    # Calculate the difference between electrons in conduction band and holes
    # in valence band.  (Sum over spins.)  This has to be constant.
    re = rr[:,1] + rr[:,3]
    rh = rr[:,0] + rr[:,2]
    rp = re - rh

    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    if (rp[0] > 0.0):
        plt.plot(times, np.log10(rp), "-k", lw=1.0, label="eq. electrons")
    else:
        plt.plot(times, np.log10(-rp), "--k", lw=1.0, label="eq. holes")
    plt.plot(times, np.log10(re), "-b", lw=1.0, label="electrons")
    plt.plot(times, np.log10(rh), "--r", lw=1.0, label="holes")
    plt.xlabel(r"$t~[{\rm fs}]$")
    plt.ylabel(r"$n~[{\rm cm}^{-2}]$")
    #ymax = np.ceil(max(rp.max(),re.max(),rh.max()) / 10.0) * 10.0
    #ymin = (np.floor(min(rp.min(),re.min(),rh.min()) / 10.0)) * 10.0
    yticks = [-2,-1,0,1,2]
    yticklabels = [r"$10^{%d}$" % (12+l) for l in yticks]
    plt.gca().set_yticks(yticks)
    plt.gca().set_yticklabels(yticklabels)
    plt.ylim([-2, 2])
    plt.legend(loc="upper left", frameon=False, fontsize="small")
    plt.savefig("%s/plot-density.png" % runDir, dpi=300)


def energy_density_in_time():
    '''Energy density.'''
    # The total must be conserved in the absence of losses to the substrate.

    # Save times.
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")

    # Load the particle density in the two bands.
    # Convert the density from [nm^-2] to [uJ cm^-2].
    ee = 100.0/6.242*np.loadtxt("%s/out-evol-energy.csv" % runDir, delimiter=",")

    # Sum of electron and phonon energy density.
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
    plt.savefig("%s/plot-energy-density.png" % runDir, dpi=300)


def differential_transmission_in_time():
    '''Differential transmission.'''

    # Save times.
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")

    # Load the probe frequencies.
    wp = np.loadtxt("%s/out-opt-diff-transm-freqs.csv" % runDir, delimiter=",")

    # Load the differential transmission in time for several probe frequencies.
    dtt = np.loadtxt("%s/out-evol-diff-transm.csv" % runDir, delimiter=",")

    # Plot the differential transmission in time.
    cmap_name = "gnuplot_r"
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(wp))]
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.53,0.7])
    for iw,(w,c) in enumerate(zip(wp,mycolors)):
        if (iw < 10):
            dt_max = np.max(dtt[:,iw])
            plt.plot(times, dtt[:,iw] / dt_max, "-", lw=1.0, c=c, label=r"%.0f~{\rm nm}" % (np.round(1240.0/w*0.1)*10.0))
    plt.xlabel(r"$t~[{\rm fs}]$")
    plt.ylabel(r"$\Delta T/T\vert_{\rm norm}$")
    plt.xlim([0.0, 120.0])  # time range to compare with experiments
    # plt.ylim([-0.25, 1.05])
    plt.legend(loc="lower left", bbox_to_anchor=(1.0,0.0), frameon=False, fontsize="small", labelspacing=0.3)
    plt.savefig("%s/plot-diff-transm.png" % runDir, dpi=300)

    # Plot the differential transmission for several times as a function of
    # frequency without normalizing to the max.
    cmap_name = "gnuplot_r"
    tt_target = np.array([10.0, 25.0, 50.0])
    itt0 = [np.abs(t - times).argmin() for t in tt_target]
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(itt0))]
    mycolors = ["g", "r", "b"]
    wmin = 0.300
    wmax = 5.000
    # plt.figure(figsize=(3.5,3.0), frameon=True)
    # ax1 = plt.axes([0.2,0.2,0.7,0.65])
    plt.figure(figsize=(3.5,1.75), frameon=True)
    ax1 = plt.axes([0.15,0.23,0.77,0.65])
    for it,c in zip(itt0,mycolors):
        plt.plot(wp, dtt[it,:] * 100.0, "o-", lw=1.0, ms=2.0, mec=None, c=c, label=r"%.0f~{\rm fs}" % times[it])
    plt.xticks([0.3, 1.0, 2.0, 3.0, 4.0, 5.0])
    #plt.xlim([wmin, wmax])
    plt.xlim([0.87, 3.8])
    plt.ylim([-0.05, 0.55])
    plt.xlabel(r"$\hbar \omega_{\rm pr}~[{\rm eV}]$")
    plt.ylabel(r"$\Delta T/T~[\%]$")
    plt.legend(loc="upper right", frameon=False, fontsize="small", labelspacing=0.3)
    ax2 = ax1.twiny()
    lticks = np.array([1200.0, 700.0, 500.0, 400.0, 300.0])
    ax2.set_xticks(1240.0/lticks)
    ax2.set_xticklabels([r"$%.0f$" % l for l in lticks])
    plt.xlim([wmin, wmax])
    plt.xlabel(r"$\lambda_{\rm pr}~[{\rm nm}]$")
    plt.savefig("%s/plot-diff-transm-freq.png" % runDir, dpi=300)

    # Plot the differential transmission for several times as a function of
    # frequency.  Each set of data at a given frequency is normalized to the
    # maximum of the data at that frequency.
    cmap_name = "gnuplot_r"
    tt_target = np.array([10.0, 25.0, 50.0])
    itt0 = [np.abs(t - times).argmin() for t in tt_target]
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(itt0))]
    wmin = 0.300
    wmax = 5.000
    dtt_max = np.array([np.max(dtt[:,iw]) for iw,w in enumerate(wp)])
    plt.figure(figsize=(3.5,3.0), frameon=True)
    ax1 = plt.axes([0.2,0.2,0.7,0.65])
    for it,c in zip(itt0,mycolors):
        plt.plot(wp, dtt[it,:] / dtt_max, "o-", lw=1.0, ms=2.0, mec=None, c=c, label=r"%.0f~{\rm fs}" % times[it])
    plt.xticks([0.3, 1.0, 2.0, 3.0, 4.0, 5.0])
    plt.xlim([wmin, wmax])
    #plt.ylim([0.00, 0.05])
    plt.xlabel(r"$\hbar \omega_{\rm pr}~[{\rm eV}]$")
    plt.ylabel(r"$\Delta T/T\vert_{\rm norm}$")
    plt.legend(loc="upper right", frameon=False, fontsize="small", labelspacing=0.3)
    ax2 = ax1.twiny()
    lticks = np.array([1200.0, 700.0, 500.0, 400.0, 300.0])
    ax2.set_xticks(1240.0/lticks)
    ax2.set_xticklabels([r"$%.0f$" % l for l in lticks])
    plt.xlim([wmin, wmax])
    plt.xlabel(r"$\lambda_{\rm pr}~[{\rm nm}]$")
    plt.savefig("%s/plot-diff-transm-freq-norm.png" % runDir, dpi=300)


def phonon_occupation_in_time():
    '''Phonon average occupation in time.'''

    # Save times.
    times = np.loadtxt("%s/out-evol-times.csv" % runDir, delimiter=",")

    # Load the phonon average occupation [dimensionless].
    ua = np.loadtxt("%s/out-evol-ph-avg-occ.csv" % runDir, delimiter=",")

    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    plt.plot(times, ua[:,0], "-k", label=r"${\rm LO}\Gamma$")
    plt.plot(times, ua[:,1], "-r", dashes=[1.0,1.0], label=r"${\rm TO}\Gamma$")
    plt.plot(times, ua[:,2], "-g", dashes=[2.0,1.0], label=r"$K$")
    plt.plot(times, ua[:,3], "-b", dashes=[3.0,1.0,1.0,1.0], label=r"$K'$")
    plt.xlabel(r"$t~[{\rm fs}]$")
    plt.ylabel(r"$n_{\nu}$")
    plt.legend(loc="lower right", frameon=False, fontsize="small")
    plt.savefig("%s/plot-phonon-occupation.png" % runDir, dpi=300)


################################################################################
## AUXILIARY ANALYSES (GEOMETRY)

def electron_geometry_valleys():
    '''Electron valleys.'''

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Load the valley label.
    vv = np.loadtxt("%s/out-geo-valley-labels.csv" % runDir, dtype="int32")

    # Plot wave vectors in the two valleys.
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    for k,v in zip(kk,vv):
        if (v == 1):
            c = "blue"
        else:
            c = "red"
        plt.scatter(k[0], k[1], marker="x", lw=0.5, s=4.0, c=c)
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-valley-labels.png" % runDir, dpi=300)

    # Load the K and K1 points.
    cor = np.loadtxt("%s/out-geo-wavevectors-K-K1.csv" % runDir, delimiter=",")

    # Load the wavevectors with respect to the K point.
    kk_K = np.loadtxt("%s/out-geo-valley-K.csv" % runDir, delimiter=",")
    kk_K1 = np.loadtxt("%s/out-geo-valley-K1.csv" % runDir, delimiter=",")

    # Plot wavevectors with respect to the center of each valley.
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    for k1,k2,v in zip(kk_K,kk_K1,vv):
        k0 = cor[v-1,:]
        if (v == 1):
            dk = k1
        else:
            dk = k2
        plt.arrow(k0[0], k0[1], dk[0], dk[1], alpha=0.2)
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-valley-K-K1.png" % runDir, dpi=300)


def band_dispersion_path():
    '''Band dispersion on a high-symmetry path.'''

    # Load the band energy.
    bands = np.loadtxt("%s/out-geo-bands.csv" % runDir, delimiter=",")

    # Load the indices of the wave vectors in the path.
    path = np.loadtxt("%s/out-geo-path.csv" % runDir, delimiter=",", dtype="int32")

    # Subtract 1 to have indices for numpy arrays instead of Fortran arrays.
    path = path - 1

    n_pri = (len(path) - 1)
    n_K = (1 + n_pri/3) - 1
    n_M = (n_K + n_pri/6) - 1

    en_M = 4.65 / 2.0

    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    # Phonon energy range, for comparison.
    eph = 0.200  # eV
    plt.axhspan(-eph/2.0, eph/2.0, color="orange", lw=0, alpha=0.5)
    plt.axhspan(1.0-eph/2.0, 1.0+eph/2.0, color="orange", lw=0, alpha=0.5)
    # High-symmetry points.
    plt.axvline(n_K, c="gray", lw=0.5)
    plt.axvline(n_M, c="gray", lw=0.5)
    plt.axhline(en_M, c="cyan", lw=0.5)
    plt.axhline(-en_M, c="cyan", lw=0.5)
    plt.plot(bands[path,0], "-o", color="red", lw=0.5, mew=None, ms=1.0)
    plt.plot(bands[path,1], "-o", color="blue", lw=0.5, mew=None, ms=1.0)
    plt.xlabel(r"${\bm k}$")
    plt.ylabel(r"$\varepsilon_{{\bm k},\lambda}~[{\rm eV}]$")
    plt.xticks([0, n_K, n_M, n_pri], [r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"])
    plt.xlim([0,n_pri])
    plt.ylim([-7.5, 7.5])
    plt.savefig("%s/plot-path-bands.png" % runDir, dpi=300)

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Distance between k points.
    dk = np.sqrt(np.vdot(kk[1,:] - kk[0,:],kk[1,:] - kk[0,:]))

    cc = np.zeros(len(kk))
    cc[path] = 1.0

    # Plot the high-symmetry path in reciprocal space.
    cmap_name = "gnuplot_r"
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
    # plt.scatter(kk[:,0], kk[:,1], s=4.0, c=ff[:,1], cmap=mpl.colormaps[cmap_name])
    p = mpl.collections.PatchCollection([plt.Circle((k[0],k[1]), radius=0.4*dk, linewidth=0.0) for k in kk], cmap=mpl.colormaps[cmap_name])
    p.set_array(cc)
    plt.gca().add_collection(p)
    plt.xlim([-1.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-path.png" % runDir, dpi=300)


def band_dispersion_brillouin():
    '''Band dispersion.'''

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Distance between k points.
    dk = np.sqrt(np.vdot(kk[1,:] - kk[0,:],kk[1,:] - kk[0,:]))

    # Load the band dispersion.
    ee = np.loadtxt("%s/out-geo-bands.csv" % runDir, delimiter=",")

    # Plot the arrows corresponding to the gradient.
    cmap_name = "gnuplot"
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    p = mpl.collections.PatchCollection([plt.Circle((k[0],k[1]), radius=0.4*dk, linewidth=0.0) for k in kk], cmap=mpl.colormaps[cmap_name])
    p.set_array(ee[:,1])
    plt.gca().add_collection(p)
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-bands.png" % runDir, dpi=300)


def band_gradient_brillouin():
    '''Gradient of the band dispersion.'''

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Load the gradient.
    gg = np.loadtxt("%s/out-geo-gradient.csv" % runDir, delimiter=",")

    # Load the band dispersion.
    ee = np.loadtxt("%s/out-geo-bands.csv" % runDir, delimiter=",")

    # Gradient length.
    ll = np.sqrt(gg[:,0]**2 + gg[:,1]**2)

    # How to color the arrows.
    nn = ll / ll.max()
    nn = (ee[:,1] - ee[:,1].min()) / (ee[:,1].max() - ee[:,1].min())

    # Plot the arrows corresponding to the gradient.
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    plt.quiver(kk[:,0], kk[:,1], gg[:,2], gg[:,3], nn, cmap="gnuplot")
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-band-gradient.png" % runDir, dpi=300)


def density_of_states():
    '''Density of states.'''

    # Load the energy bins.
    ee = np.loadtxt("%s/out-geo-energy-bins.csv" % runDir, delimiter=",")

    # Load the density of states.  This is spin-resolved, need to multiply by 2.
    nu = np.loadtxt("%s/out-geo-dos.csv" % runDir, delimiter=",")

    hbarvf = 0.50  # eV nm

    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    plt.bar(ee[:-1], 2.0*nu, color="gray", align="edge", width=(ee[1]-ee[0]))
    # Linear DOS around the Dirac points.
    xx = np.array([-2.0, 0.0, 2.0])
    yy = np.abs((2.0*2.0/(2.0*np.pi)/hbarvf**2) * xx)
    plt.plot(xx, yy, "-r", lw=0.5)
    # M points.
    plt.axvline(4.65/2.0, color="red", lw=0.5)
    plt.axvline(-4.65/2.0, color="red", lw=0.5)
    plt.xlabel(r"$\varepsilon~[{\rm eV}]$")
    plt.ylabel(r"$\nu(\varepsilon)~[{\rm eV}^{-1}\,{\rm nm}^{-2}]$")
    #plt.ylim([0.0, 30.0])
    plt.xticks(np.linspace(-8.0,8.0,9))
    plt.xlim([-8.0,8.0])
    plt.savefig("%s/plot-dos.png" % runDir, dpi=300)


def coulomb_scattering_mesh():
    '''Coulomb scattering wave vector mesh.'''

    # Coulomb scattering wave vectors.
    kk = np.loadtxt("%s/out-geo-coulomb-scatt-mesh.csv" % runDir, delimiter=",")

    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    plt.scatter(kk[:,0], kk[:,1], marker="x", lw=0.5, s=4.0, c="black")
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$q_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$q_{y}~[{\rm nm}^{-1}]$")
    plt.title(r"\scriptsize Coulomb scattering wave vector mesh")
    plt.savefig("%s/plot-scattering-mesh.png" % runDir, dpi=300)


def electron_phonon_scattering_support():
    '''Support of phonon modes and electrons included in the scattering.'''

    # Wave vectors in the BZ.
    kkbz = np.loadtxt("%s/out-geo-brillouin-zone.csv" % runDir, delimiter=",")
    kkpc = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Support of the four phonon modes.
    supp = np.loadtxt("%s/out-geo-phonon-support.csv" % runDir, delimiter=",", dtype="int32")

    # Select wave vectors for each modes.
    kk_LOG = kkbz[supp[:,0]==1,:]
    kk_TOG = kkbz[supp[:,1]==1,:]
    kk_A1K = kkpc[supp[:,2]==1,:]
    kk_A1K1 = kkpc[supp[:,3]==1,:]

    # Mark the wave vectors where the several modes have support.
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    plot_cells()
    plt.scatter(kk_LOG[:,0], kk_LOG[:,1], marker="+", lw=0.2, s=1.0, c="black")
    plt.scatter(kk_TOG[:,0], kk_TOG[:,1], marker="x", lw=0.2, s=1.0, c="red")
    plt.scatter(kk_A1K[:,0], kk_A1K[:,1], marker="+", lw=0.2, s=1.0, c="green")
    plt.scatter(kk_A1K1[:,0], kk_A1K1[:,1], marker="x", lw=0.2, s=1.0, c="blue")
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-phonon-modes-support.png" % runDir, dpi=300)

    # Support of the electrons included in the scattering.
    # Since the supports are centered around the K and K' points, plot in
    # the primitive cell instead of the Brillouin zone.
    supp_el = np.loadtxt("%s/out-geo-phonon-scatt-el-supp.csv" % runDir, delimiter=",", dtype="int32")
    kk_el = kkpc[supp_el[:]==1,:]

    # Mark the wave vectors corresponding to the included electrons.
    plt.figure(figsize=(3.5,3.5), frameon=True)
    ax1 = plt.axes([0.25,0.2,0.7*(46.0/60.0),0.7])
    plt.scatter(kk_el[:,0], kk_el[:,1], marker="x", lw=0.5, s=5.0, c="white", ec="black")
    plot_cells()
    plt.xlim([-16.0, 30.0])
    plt.ylim([-30.0, 30.0])
    plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
    plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
    plt.savefig("%s/plot-phonon-scatt-el-supp.png" % runDir, dpi=300)


################################################################################
## DEVELOPMENT ANALYSES (CHECKS)

def check_electron_thermodynamics():
    '''Electron thermodynamics.'''

    # 2024-06-24
    # OK but needs large number of wave vectors (factor >10) to achieve good
    # results at Fermi energies of ~200meV.

    thrm = np.loadtxt("%s/out-thrm-el-1.csv" % runDir, delimiter=",")

    # Density vs. Fermi energy.
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    xx = np.linspace(0.0, 3.0, 101)
    # Calculate n = eF^2 / (pi hbar), and convert nm^-2 to cm^-2.
    yy = 100.0 * xx*xx / (np.pi * 0.6**2.0)
    plt.plot(xx, yy, "-r", lw=0.5)
    plt.plot(thrm[:,0], thrm[:,1]*100.0, "-ok", lw=0.5, ms=3.0)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim([0.1, 3.0])
    plt.ylim([1.0, 1000.0])
    plt.xlabel(r"$\varepsilon_{\rm F}~[{\rm eV}]$")
    plt.ylabel(r"$n~[10^{12}\,{\rm cm}^{-2}]$")
    plt.savefig("%s/plot-thrm-dens-vs-eF.png" % runDir, dpi=300)

    # Chemical potential vs. Fermi energy for several temperatures.
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.7,0.7])
    plt.plot(thrm[:,0], thrm[:,2], "-ob", lw=0.5, ms=3.0, label=r"$300~{\rm K}$")
    plt.plot(thrm[:,0], thrm[:,3], "-or", lw=0.5, ms=3.0, label=r"$2000~{\rm K}$")
    plt.xlabel(r"$\varepsilon_{\rm F}~[{\rm eV}]$")
    plt.ylabel(r"$\mu~[{\rm eV}]$")
    plt.legend(loc="lower right", frameon=False, fontsize="small")
    plt.savefig("%s/plot-thrm-mu-vs-eF.png" % runDir, dpi=300)


def check_equilibrium_optical_conductivity():
    '''Equilibrium optical conductivity.'''

    # Equilibrium optical conductivity, as a function of frequency.
    ww = np.loadtxt("%s/out-test-cond-freq.csv" % runDir)
    dd = np.loadtxt("%s/out-test-cond-eta.csv" % runDir)
    ss = np.loadtxt("%s/out-test-cond-unscr.csv" % runDir, delimiter=",")

    cmap_name = "gnuplot_r"
    mycolors = [mpl.colormaps[cmap_name](u) for u in np.linspace(0.2, 0.8, len(dd))]
    plt.figure(figsize=(3.5,3.0), frameon=True)
    plt.axes([0.2,0.2,0.53,0.7])
    for i,(d,c) in enumerate(zip(dd,mycolors)):
        plt.plot(ww, ss[:,i], "-", c=c, lw=0.5, label=r"$%.3f$" % d)
    plt.xlabel(r"$\hbar \omega_{\rm pr}~[{\rm eV}]$")
    plt.ylabel(r"$\sigma_{\rm eq}(\omega_{\rm pr})/\sigma_{0}$")
    plt.figtext(0.80, 0.18, r"$\eta~[{\rm eV}]$")
    plt.legend(loc="lower left", bbox_to_anchor=(1.0,0.0), frameon=False, fontsize="small", labelspacing=0.3)
    plt.savefig("%s/plot-eq-cond-unscr.png" % runDir, dpi=300)

    # Optionally plot the dielectric function at q=0 as a function of frequency.
    if os.path.exists("%s/out-test-cond-eps.csv" % runDir):
        ee = np.loadtxt("%s/out-test-cond-eps.csv" % runDir)
        plt.figure(figsize=(3.5,3.0), frameon=True)
        plt.axes([0.2,0.2,0.7,0.7])
        plt.plot(ww, ee, "-", lw=0.5)
        plt.xlabel(r"$\hbar \omega_{\rm pr}~[{\rm eV}]$")
        plt.ylabel(r"${\rm Re}[\varepsilon(q=0,\omega_{\rm pr})]$")
        plt.savefig("%s/plot-eq-cond-eps.png" % runDir, dpi=300)


def check_spectral_density_optical():
    '''Spectral density for optical transitions.'''

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Load the probe frequencies.
    wp = np.loadtxt("%s/out-opt-diff-transm-freqs.csv" % runDir, delimiter=",")

    # Load the spectral density for several probe frequencies.
    spdens = np.loadtxt("%s/out-diff-transm-dens.csv" % runDir, delimiter=",")

    # Distance between k points.
    dk = np.sqrt(np.vdot(kk[1,:] - kk[0,:],kk[1,:] - kk[0,:]))

    # Plot the spectral density for the probe frequencies.
    cmap_name = "gnuplot_r"
    for iw,w in enumerate(wp):
        plt.figure(figsize=(3.5,3.5), frameon=True)
        ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
        # plt.scatter(kk[:,0], kk[:,1], s=4.0, c=ff[:,1], cmap=mpl.colormaps[cmap_name])
        p = mpl.collections.PatchCollection([plt.Circle((k[0],k[1]), radius=0.4*dk, linewidth=0.0) for k in kk], cmap=mpl.colormaps[cmap_name])
        p.set_array(spdens[:,iw])
        plt.gca().add_collection(p)
        plt.xlim([-1.0, 30.0])
        plt.ylim([-30.0, 30.0])
        plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
        plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
        plt.title(r"\scriptsize $\lambda = %.0f~{\rm nm}$" % (np.round(1240.0/w*0.01)*100.0))
        ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
        mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
        plt.savefig("%s/plot-diff-transm-dens-%1d.png" % (runDir, iw), dpi=300)


def check_lindhard_azimuthal():
    '''Lindhard function, averaged over the angle.'''

    # Wave vectors.
    qq = np.loadtxt("%s/out-geo-lindhard-wv-mesh.csv" % runDir, delimiter=",")

    # Transition frequencies.
    ww = np.loadtxt("%s/out-geo-frequency-interp-mesh.csv" % runDir, delimiter=",")

    # Save times.
    times = np.loadtxt("%s/out-snap-times.csv" % runDir, delimiter=",")

    # Load data at all save times.
    chi_re_t = []
    chi_im_t = []
    for it,t in enumerate(times):
        chi_re_t.append(np.loadtxt("%s/out-snap-%04d-lindhard-re.csv" % (runDir, it+1), delimiter=","))
        chi_im_t.append(np.loadtxt("%s/out-snap-%04d-lindhard-im.csv" % (runDir, it+1), delimiter=","))

    # Light cone.
    # In the code, we fix the energy at the M point, so the slope at the K point
    # is different from the standard value.  Fix this value fitting the DOS.
    hbarvf = 0.50  # eV nm
    lc_q = np.array([0.0,10.0])
    lc_w = hbarvf * lc_q

    # Plot real and imaginary part separately because they have different
    # axes limits and labels.
    cmap_name = "seismic"
    for it,(t,chi_re,chi_im) in enumerate(zip(times,chi_re_t,chi_im_t)):
        # Real part.
        valmax = 4.0
        plt.figure(figsize=(3.5,3.5), frameon=True)
        ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
        plt.imshow(chi_re, origin="lower", aspect="auto", interpolation="none", extent=(qq[0],qq[-1],ww[0],ww[-1]), vmin=-valmax, vmax=valmax, cmap=cmap_name)
        plt.plot(lc_q, lc_w, "--", lw=0.5, c="black")
        plt.plot(lc_q, -lc_w, "--", lw=0.5, c="black")
        plt.axhline(0.0, lw=0.5, c="black")
        plt.xlim([0.0, 5.0])
        plt.ylim([-2.0, 2.0])
        plt.xlabel(r"$q~[{\rm nm}^{-1}]$")
        plt.ylabel(r"$\hbar\omega~[{\rm eV}]$")
        plt.title(r"\scriptsize $t = %.0f~{\rm fs}$" % t)
        ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
        mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
        yticks = np.arange(int(-valmax),int(valmax+1))
        yticklabels = np.array(["%.1f" % yt for yt in yticks])
        ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
        ax2.set_yticklabels(yticklabels)
        plt.savefig("%s/plot-lindhard-re-%04d.png" % (runDir, it), dpi=300)
        # Imaginary part.
        valmax = 4.0
        plt.figure(figsize=(3.5,3.5), frameon=True)
        ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
        plt.imshow(chi_im, origin="lower", aspect="auto", interpolation="none", extent=(qq[0],qq[-1],ww[0],ww[-1]), vmin=-valmax, vmax=valmax, cmap=cmap_name)
        plt.plot(lc_q, lc_w, "--", lw=0.5, c="black")
        plt.plot(lc_q, -lc_w, "--", lw=0.5, c="black")
        plt.axhline(0.0, lw=0.5, c="black")
        plt.xlim([0.0, 5.0])
        plt.ylim([-2.0, 2.0])
        plt.xlabel(r"$q~[{\rm nm}^{-1}]$")
        plt.ylabel(r"$\hbar\omega~[{\rm eV}]$")
        plt.title(r"\scriptsize $t = %.0f~{\rm fs}$" % t)
        ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
        mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
        yticks = np.arange(int(-valmax),int(valmax+1))
        yticklabels = np.array(["%.1f" % yt for yt in yticks])
        ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
        ax2.set_yticklabels(yticklabels)
        plt.savefig("%s/plot-lindhard-im-%04d.png" % (runDir, it), dpi=300)


def check_spectral_density_path():
    '''Spectral density on a high-symmetry path in the Brillouin Zone.'''

    # Wave vectors.
    qq = np.loadtxt("%s/out-geo-lindhard-wv-mesh.csv" % runDir, delimiter=",")

    # Transition frequencies.
    ww = np.loadtxt("%s/out-geo-frequency-interp-mesh.csv" % runDir, delimiter=",")

    # Save times.
    times = np.loadtxt("%s/out-snap-times.csv" % runDir, delimiter=",")

    # Load data at all save times.
    chi_re_t = []
    chi_im_t = []
    for it,t in enumerate(times):
        chi_re_t.append(np.loadtxt("%s/out-snap-%04d-lindhard-re.csv" % (runDir, it+1), delimiter=","))
        chi_im_t.append(np.loadtxt("%s/out-snap-%04d-lindhard-im.csv" % (runDir, it+1), delimiter=","))

    # Light cone.
    # In the code, we fix the energy at the M point, so the slope at the K point
    # is different from the standard value.  Fix this value fitting the DOS.
    hbarvf = 0.50  # eV nm
    lc_q = np.array([0.0,10.0])
    lc_w = hbarvf * lc_q

    # Plot real and imaginary part separately because they have different
    # axes limits and labels.
    cmap_name = "seismic"
    for it,(t,chi_re,chi_im) in enumerate(zip(times,chi_re_t,chi_im_t)):
        # Real part.
        valmax = 4.0
        plt.figure(figsize=(3.5,3.5), frameon=True)
        ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
        plt.imshow(chi_re, origin="lower", aspect="auto", interpolation="none", extent=(qq[0],qq[-1],ww[0],ww[-1]), vmin=-valmax, vmax=valmax, cmap=cmap_name)
        plt.plot(lc_q, lc_w, "--", lw=0.5, c="black")
        plt.plot(lc_q, -lc_w, "--", lw=0.5, c="black")
        plt.axhline(0.0, lw=0.5, c="black")
        plt.xlim([0.0, 5.0])
        plt.ylim([-2.0, 2.0])
        plt.xlabel(r"$q~[{\rm nm}^{-1}]$")
        plt.ylabel(r"$\hbar\omega~[{\rm eV}]$")
        plt.title(r"\scriptsize $t = %.0f~{\rm fs}$" % t)
        ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
        mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
        yticks = np.arange(int(-valmax),int(valmax+1))
        yticklabels = np.array(["%.1f" % yt for yt in yticks])
        ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
        ax2.set_yticklabels(yticklabels)
        plt.savefig("%s/plot-lindhard-re-%04d.png" % (runDir, it), dpi=300)
        # Imaginary part.
        valmax = 4.0
        plt.figure(figsize=(3.5,3.5), frameon=True)
        ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
        plt.imshow(chi_im, origin="lower", aspect="auto", interpolation="none", extent=(qq[0],qq[-1],ww[0],ww[-1]), vmin=-valmax, vmax=valmax, cmap=cmap_name)
        plt.plot(lc_q, lc_w, "--", lw=0.5, c="black")
        plt.plot(lc_q, -lc_w, "--", lw=0.5, c="black")
        plt.axhline(0.0, lw=0.5, c="black")
        plt.xlim([0.0, 5.0])
        plt.ylim([-2.0, 2.0])
        plt.xlabel(r"$q~[{\rm nm}^{-1}]$")
        plt.ylabel(r"$\hbar\omega~[{\rm eV}]$")
        plt.title(r"\scriptsize $t = %.0f~{\rm fs}$" % t)
        ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
        mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
        yticks = np.arange(int(-valmax),int(valmax+1))
        yticklabels = np.array(["%.1f" % yt for yt in yticks])
        ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
        ax2.set_yticklabels(yticklabels)
        plt.savefig("%s/plot-lindhard-im-%04d.png" % (runDir, it), dpi=300)


def check_gamma_matrices_brillouin():
    '''Gamma matrices (scattering rates) in the Brillouin zone.'''

    # Wave vectors in the primitive cell.
    kk = np.loadtxt("%s/out-geo-primitive-cell.csv" % runDir, delimiter=",")

    # Distance between k points.
    dk = np.sqrt(np.vdot(kk[1,:] - kk[0,:],kk[1,:] - kk[0,:]))

    # Load the gamma matrices.
    gamma_names = ["ebp-gammas-in", "ebp-gammas-out"]
    gamma_list = []
    for gamma_name in gamma_names:
        gamma_list.append(np.loadtxt("%s/out-stop-%s.csv" % (runDir, gamma_name), delimiter=","))

    # Logarithmic scale for the colorbar.
    yticks = [-7, -6, -5, -4, -3, -2, -1]
    yticklabels = [r"$<\!10^{-7}$", r"$10^{-6}$", r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$>\!10^{-1}$"]
    fmin = 10**(yticks[0])
    fmax = 10**(yticks[-1])

    # Plot the distribution in each band in time.
    cmap_name = "gnuplot_r"
    for gamma, gamma_name in zip(gamma_list, gamma_names):
        for ib,bn in zip([0,1],["va","co"]):  # iterate over bands
            plt.figure(figsize=(3.5,3.5), frameon=True)
            ax1 = plt.axes([0.25,0.2,0.7*(31.0/60.0),0.7])
            p = mpl.collections.PatchCollection([plt.Circle((k[0],k[1]), radius=0.4*dk, linewidth=0.0) for k in kk], cmap=mpl.colormaps[cmap_name])
            p.set_array(trunc_g(gamma[:,ib],fmin,fmax))
            plt.gca().add_collection(p)
            plt.xlim([-1.0, 30.0])
            plt.ylim([-30.0, 30.0])
            plt.xlabel(r"$k_{x}~[{\rm nm}^{-1}]$")
            plt.ylabel(r"$k_{y}~[{\rm nm}^{-1}]$")
            ax2 = plt.axes([0.3+0.7*(31.0/60.0),0.2,0.07,0.7])
            mpl.colorbar.ColorbarBase(ax2, orientation="vertical", cmap=mpl.colormaps[cmap_name])
            ax2.set_yticks(np.linspace(0.0,1.0,len(yticks)))
            ax2.set_yticklabels(yticklabels)
            plt.savefig("%s/plot-%s-%s.png" % (runDir, gamma_name, bn), dpi=300)


def check_spectral_density_support():
    '''Distribution of exchanged energies at fixed exchanged wave vector in electron-electron scattering'''

    



################################################################################
## CROSS-RUN ANALYSES, ORGANIZATION

def overview_html():
    '''Generate an overview HTML file with useful plots for several runs.'''

    # HTML templates.

    html_template_a = '''<!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <meta http-equiv="X-UA-Compatible" content="ie=edge">
            <title>Results overview</title>
            <style>
            body {
                margin: 0;
                padding: 20px;
                font-family: Arial, sans-serif;
                background-color: #f4f4f4;
            }
            .image-section {
                width: 840px;
                margin: auto;
                margin-bottom: 30px;
            }
            .image-row {
                background-color: #bac;
                padding: 10px;
            }
            .image-row img {
                width: 250px;
                margin: 10px;
            }
            .description {
                text-align: left;
                margin-top: 5px;
                padding: 10px;
                color: #ccc;
                background-color: #3aa;
            }
            </style>
        </head>
        <body>
    '''

    html_template_z = '''    </body>
    </html>
    '''

    html_template_img = '''
        <div class="image-section">
            <div class="image-row">
                <img src="{rundir}/plot-diff-transm.png">
                <img src="{rundir}/plot-diff-transm-freq.png">
                <img src="{rundir}/plot-density.png">
            </div>
            <p class="description">{description}</p>
        </div>
    '''

    # Ignore runs before this date.

    day0 = datetime.datetime(2024,11,1)

    dataDir = "../data"

    # Build run list from data on disk.
    # Remember that glob does not sort!
    runDirList = glob.glob("%s/????-??-??_*" % dataDir)

    # For each run, define its date.
    runDateList = []
    for runDir in runDirList:
        m = re.search(r'\.\.\/data\/(....)\-(..)\-(..)\_*', runDir)
        runDateList.append(datetime.datetime(int(m.group(1)), int(m.group(2)), int(m.group(3))))

    # Prepare the text for the HTML file.
    html_text = html_template_a

    # Indices to sort by folder name, reverse order.
    ii = np.argsort(runDirList)[::-1]

    # For the runs in the desired time interval, create the HTML entry.
    for i in ii:
        runDir = runDirList[i]
        runDate = runDateList[i]
        if (runDate > day0):
            data = {
                "rundir": runDir[8:],
                "description": runDir[8:].replace("_"," ")
            }
            html_img = html_template_img.format(**data)
            html_text = html_text + html_img

    html_text = html_text + html_template_z

    # Save file to disk.
    f = open("%s/results_summary.html" % dataDir, "w")
    f.write(html_text)
    f.close()


################################################################################
## RUN A SET OF ANALYSES

if __name__ == "__main__":

    # Folder with the data.
    runDir = "../run"
    #runDir = "../data/2025-01-20_1712_DYN_f10_x09_p30/"
    
    electron_distro_primitive_cell()
    electron_density_in_time()
    energy_density_in_time()
    differential_transmission_in_time()


