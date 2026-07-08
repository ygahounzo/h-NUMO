#!/usr/bin/env python3
"""
Analysis for h-NUMO's 'gravity_wave' test case (Chen 2025, QJRMS,
10.1002/qj.4994, Sec 5.2 / Figure 5).

Estimates BOTH gravity-wave speeds from radial thickness-deviation profiles:
  - barotropic mode (~140 m/s): needs a SHORT snapshot pair, well inside one
    circumnavigation time (circumference / speed ~ 79h), since past that the
    front wraps around, converges/focuses at the antipode, and the simple
    "leading edge from source" tracking below becomes meaningless (a
    negative or >20015km "distance" in the printed result is the tell).
  - baroclinic mode (~10-14 m/s): needs a LONG snapshot pair (Chen's paper
    used ~150-225h) for the front to separate from the initial hump's own
    footprint (radius ~1112km).
Because these two modes need snapshot pairs from very different times --
often from entirely different runs, if a single run didn't save early
snapshots -- the barotropic and baroclinic pairs each have their own
RUNDIR/FPREFIX/indices/timing settings below (defaulting to the shared
RUNDIR/FPREFIX if not overridden).

Also renders flattened lon-lat maps of the barotropic-pair snapshots
(Chen's Fig 5 style).

Usage:
    python3 analyze_gravity_wave.py [rundir]

Edit the SETTINGS block below for run configuration.
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

from numo_vtk import (read_vtk_ascii, great_circle_km, radial_bin_mean, plot_latlon_map,
                       snapshot_time_hours, CIRCUMFERENCE_KM)

# ---------------- SETTINGS ----------------
RUNDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
FPREFIX = 'mlswe_sphere_300.0000_dg_rk35_l001_'
TIME_RESTART_HOURS = 0.5  # time_restart in numo3d.in, shared default

# Short pair for the barotropic mode -- must stay well inside one
# circumnavigation time (see the warning printed below if it doesn't).
# Override RUNDIR_BTP/FPREFIX_BTP/TIME_RESTART_HOURS_BTP if the run with
# early snapshots is a *different* run/directory than the one used for the
# baroclinic pair below (common: a long run kept only sparse late snapshots
# and has nothing early enough for a clean barotropic measurement).
RUNDIR_BTP = RUNDIR
FPREFIX_BTP = FPREFIX
TIME_RESTART_HOURS_BTP = TIME_RESTART_HOURS
IDX1_BTP, IDX2_BTP = 16, 32

# Long pair for the baroclinic mode -- needs to span far enough for the
# front to separate from the source (Chen's paper: ~150-225h of simulated
# time).
RUNDIR_BCL = RUNDIR
FPREFIX_BCL = FPREFIX
TIME_RESTART_HOURS_BCL = TIME_RESTART_HOURS
IDX1_BCL, IDX2_BCL = IDX1_BTP, IDX2_BTP

BACKGROUND = 1000.0     # resting layer thickness (m)
LON0, LAT0 = 0.0, 0.0   # hump center (radians)
BINWIDTH_KM = 500.0
BTP_THRESHOLD = 0.01    # m; deviation magnitude marking the leading edge of
                        # the fast (barotropic) wave, well above the ~0.001-
                        # 0.002m interpolation noise floor but well below
                        # real signal (tens of m near source, ~0.05m at the
                        # barotropic front itself)
BCL_THRESHOLD = 2.0     # m; deviation magnitude marking the edge of the
                        # baroclinic disturbance -- much larger than
                        # BTP_THRESHOLD since the baroclinic signal itself is
                        # O(10s of m).
RHO1, RHO2 = 1000.0, 1020.0  # layer densities (kg/m^3), from
                              # initial_conditions_sphere.F90's 'gravity_wave'
                              # case (paper's other tested config: 1040)
GRAVITY = 9.81
PHI1, PHI2 = 1000.0, 1000.0   # background layer thicknesses (m)
PREDICTED_BTP_SPEED = np.sqrt(GRAVITY * (PHI1 + PHI2))
G_REDUCED = GRAVITY * (RHO2 - RHO1) / RHO1
PREDICTED_BCL_SPEED = np.sqrt(G_REDUCED * PHI1 * PHI2 / (PHI1 + PHI2))
# -------------------------------------------


def _dt_seconds(fprefix):
    # dt (bcl time step, seconds) is embedded literally in h-NUMO's filename,
    # e.g. "..._400.0000_...". Needed because the *actual* spacing between
    # VTK snapshots is irestart*dt (irestart = NINT(time_restart_seconds/dt)),
    # which only equals the nominal time_restart when dt evenly divides
    # 1800s -- see numo_vtk.snapshot_time_hours. A wrong/hardcoded
    # time-per-snapshot here silently mislabels snapshot times and biases
    # the wave-speed estimate below (this bit us already with a dt=400 run
    # -- don't reintroduce it).
    return float(re.search(r'_(\d+\.\d+)_', fprefix).group(1))


def load(rundir, fprefix, idx, time_restart_hours):
    dt_seconds = _dt_seconds(fprefix)
    t = snapshot_time_hours(idx, dt_seconds, time_restart_hours)
    f = os.path.join(rundir, f'{fprefix}{idx:04d}.vtk')
    d = read_vtk_ascii(f, fields=['h'])
    print(f"Loaded {len(d['h'])} points from {f} at t={t:.1f}h")
    return t, d


def radial_profile(d):
    dist = great_circle_km(LON0, LAT0, d['lon'], d['lat'])
    dev = d['h'] - BACKGROUND
    return radial_bin_mean(dist, dev, BINWIDTH_KM)


def leading_edge(centers, mean, threshold):
    """Outermost radius where |mean| still exceeds `threshold`, i.e. the
    leading edge of whichever mode that threshold is tuned for -- NOT
    wherever |mean| happens to be locally largest (that can instead pick
    out a different mode's larger-amplitude residual nearer the source)."""
    absmean = np.nan_to_num(np.abs(mean))
    above = np.where(absmean > threshold)[0]
    if len(above) == 0:
        return np.nan
    i = above[-1]
    # Linear-interpolate the threshold crossing between bin i (above
    # threshold) and bin i+1 (below) for sub-bin precision -- snapping to
    # the raw bin edge otherwise makes the speed estimate sensitive to bin
    # width (e.g. a single 500km-bin shift over an 8h window is already a
    # ~17 m/s error).
    if i == len(centers) - 1 or np.isnan(mean[i + 1]):
        return centers[i]
    y0, y1 = absmean[i], absmean[i + 1]
    if y0 == y1:
        return centers[i]
    frac = (y0 - threshold) / (y0 - y1)
    return centers[i] + frac * (centers[i + 1] - centers[i])


def check_circumnavigation(label, t_hours, predicted_speed_ms):
    """Warn if the wave (at its predicted speed) would already be past the
    antipode by this snapshot's time -- beyond that point great-circle
    distance starts reporting the *shorter* path back toward the source
    (wrapping), and the simple outermost-leading-edge tracking below no
    longer measures a meaningful "distance travelled"."""
    expected_r_km = predicted_speed_ms * t_hours * 3600.0 / 1000.0
    half_circumference = CIRCUMFERENCE_KM / 2.0
    if expected_r_km > half_circumference:
        n_wraps = expected_r_km / CIRCUMFERENCE_KM
        print(f"*** WARNING: at t={t_hours:.1f}h the {label} wave (predicted "
              f"{predicted_speed_ms:.1f} m/s) would already have travelled "
              f"~{expected_r_km:.0f} km, past the antipode ({half_circumference:.0f} km) "
              f"-- roughly {n_wraps:.2f} circumnavigations. The leading-edge "
              f"measurement below is likely contaminated by wraparound/antipodal "
              f"focusing; use a shorter snapshot pair (from an earlier-snapshot run "
              f"if this run only saved sparse late-time output). ***")


def main():
    check_circumnavigation("barotropic",
                            snapshot_time_hours(IDX2_BTP, _dt_seconds(FPREFIX_BTP), TIME_RESTART_HOURS_BTP),
                            PREDICTED_BTP_SPEED)

    t1_btp, d1_btp = load(RUNDIR_BTP, FPREFIX_BTP, IDX1_BTP, TIME_RESTART_HOURS_BTP)
    t2_btp, d2_btp = load(RUNDIR_BTP, FPREFIX_BTP, IDX2_BTP, TIME_RESTART_HOURS_BTP)
    centers_btp, mean1_btp = radial_profile(d1_btp)
    _, mean2_btp = radial_profile(d2_btp)

    same_pair = ((RUNDIR_BCL, FPREFIX_BCL, IDX1_BCL, IDX2_BCL, TIME_RESTART_HOURS_BCL) ==
                 (RUNDIR_BTP, FPREFIX_BTP, IDX1_BTP, IDX2_BTP, TIME_RESTART_HOURS_BTP))
    if same_pair:
        t1_bcl, t2_bcl = t1_btp, t2_btp
        centers_bcl, mean1_bcl, mean2_bcl = centers_btp, mean1_btp, mean2_btp
        d1_bcl, d2_bcl = d1_btp, d2_btp
    else:
        t1_bcl, d1_bcl = load(RUNDIR_BCL, FPREFIX_BCL, IDX1_BCL, TIME_RESTART_HOURS_BCL)
        t2_bcl, d2_bcl = load(RUNDIR_BCL, FPREFIX_BCL, IDX2_BCL, TIME_RESTART_HOURS_BCL)
        centers_bcl, mean1_bcl = radial_profile(d1_bcl)
        _, mean2_bcl = radial_profile(d2_bcl)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(centers_btp, mean1_btp, '-o', label=f't={t1_btp:.1f}h (btp pair)')
    ax.plot(centers_btp, mean2_btp, '-s', label=f't={t2_btp:.1f}h (btp pair)')
    if not same_pair:
        ax.plot(centers_bcl, mean1_bcl, '--o', label=f't={t1_bcl:.1f}h (bcl pair)')
        ax.plot(centers_bcl, mean2_bcl, '--s', label=f't={t2_bcl:.1f}h (bcl pair)')
    ax.set_xlabel('great-circle distance from hump center (km)')
    ax.set_ylabel('mean layer thickness deviation (m)')
    ax.set_title('Gravity-wave radial profile')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    r1_btp = leading_edge(centers_btp, mean1_btp, BTP_THRESHOLD)
    r2_btp = leading_edge(centers_btp, mean2_btp, BTP_THRESHOLD)
    speed_btp = (r2_btp - r1_btp) * 1000.0 / ((t2_btp - t1_btp) * 3600.0)

    r1_bcl = leading_edge(centers_bcl, mean1_bcl, BCL_THRESHOLD)
    r2_bcl = leading_edge(centers_bcl, mean2_bcl, BCL_THRESHOLD)
    speed_bcl = (r2_bcl - r1_bcl) * 1000.0 / ((t2_bcl - t1_bcl) * 3600.0)

    print(f"\nBarotropic front: t={t1_btp:.1f}h -> r={r1_btp:.0f} km,  t={t2_btp:.1f}h -> r={r2_btp:.0f} km")
    print(f"Estimated barotropic speed: {speed_btp:.1f} m/s  (predicted ~{PREDICTED_BTP_SPEED:.1f} m/s)")

    print(f"\nBaroclinic front: t={t1_bcl:.1f}h -> r={r1_bcl:.0f} km,  t={t2_bcl:.1f}h -> r={r2_bcl:.0f} km")
    print(f"Estimated baroclinic speed: {speed_bcl:.1f} m/s  (predicted ~{PREDICTED_BCL_SPEED:.1f} m/s)")
    if same_pair:
        print("(Baroclinic estimate reuses the barotropic pair's short snapshot times -- unreliable; "
              "set IDX1_BCL/IDX2_BCL (and RUNDIR_BCL/FPREFIX_BCL if it's a different run) to a long "
              "run's indices, e.g. ~150-225h apart, for a clean measurement.)")

    # Lon-lat maps (Chen Fig 5 style) -- barotropic-pair snapshots
    fig2, axes = plt.subplots(1, 2, figsize=(13, 5))
    plot_latlon_map(axes[0], d1_btp['lon'], d1_btp['lat'], d1_btp['h'],
                     title=f'layer thickness, t={t1_btp:.1f}h',
                     mark_point=(np.degrees(LON0), np.degrees(LAT0)))
    plot_latlon_map(axes[1], d2_btp['lon'], d2_btp['lat'], d2_btp['h'],
                     title=f'layer thickness, t={t2_btp:.1f}h',
                     mark_point=(np.degrees(LON0), np.degrees(LAT0)))
    fig2.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
