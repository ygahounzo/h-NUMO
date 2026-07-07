#!/usr/bin/env python3
"""
Analysis for h-NUMO's 'gravity_wave' test case (Chen 2025, QJRMS,
10.1002/qj.4994, Sec 5.2 / Figure 5).

For two layer-1 snapshots, this:
  1. Bins the thickness deviation from background by great-circle distance
     from the hump center, and estimates the wavefront propagation speed
     from the shift in the far-field anomaly peak between the two times.
  2. Renders flattened lon-lat maps of both snapshots (Chen's Fig 5 style).

Usage:
    python3 analyze_gravity_wave.py [rundir] [idx1] [idx2]

Defaults assume h-NUMO's 'l001_XXXX.vtk' naming and time_restart=0.5h.
Edit the SETTINGS block below for a different run configuration.
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

from numo_vtk import read_vtk_ascii, great_circle_km, radial_bin_mean, plot_latlon_map, snapshot_time_hours

# ---------------- SETTINGS ----------------
RUNDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
FPREFIX = 'mlswe_sphere_300.0000_dg_rk35_l001_'
IDX1 = int(sys.argv[2]) if len(sys.argv) > 2 else 16
IDX2 = int(sys.argv[3]) if len(sys.argv) > 3 else 32
TIME_RESTART_HOURS = 0.5  # time_restart in numo3d.in
# dt (bcl time step, seconds) -- auto-detected from FPREFIX (h-NUMO embeds
# it literally, e.g. "..._400.0000_..."); override manually if that fails.
# This matters because the *actual* spacing between VTK snapshots is
# irestart*dt (irestart = NINT(time_restart_seconds/dt)), which only equals
# TIME_RESTART_HOURS exactly when dt evenly divides 1800s -- see
# numo_vtk.snapshot_time_hours.
DT_SECONDS = float(re.search(r'_(\d+\.\d+)_', FPREFIX).group(1))
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
                        # O(10s of m). CAVEAT: at short run times the
                        # baroclinic front has barely moved beyond the
                        # initial hump's own footprint (radius ~1112km), so
                        # this measures the edge of the still-adjusting
                        # source region more than a cleanly separated
                        # propagating front -- Chen's own paper only gets a
                        # clean baroclinic speed by running ~150-225h, far
                        # longer than typical short test runs here.
RHO1, RHO2 = 1000.0, 1020.0  # layer densities (kg/m^3), from
                              # initial_conditions_sphere.F90's 'gravity_wave'
                              # case (paper's other tested config: 1040)
GRAVITY = 9.81
PHI1, PHI2 = 1000.0, 1000.0   # background layer thicknesses (m)
PREDICTED_BTP_SPEED = np.sqrt(GRAVITY * (PHI1 + PHI2))
G_REDUCED = GRAVITY * (RHO2 - RHO1) / RHO1
PREDICTED_BCL_SPEED = np.sqrt(G_REDUCED * PHI1 * PHI2 / (PHI1 + PHI2))
# -------------------------------------------


def main():
    t1 = snapshot_time_hours(IDX1, DT_SECONDS, TIME_RESTART_HOURS)
    t2 = snapshot_time_hours(IDX2, DT_SECONDS, TIME_RESTART_HOURS)
    f1 = os.path.join(RUNDIR, f'{FPREFIX}{IDX1:04d}.vtk')
    f2 = os.path.join(RUNDIR, f'{FPREFIX}{IDX2:04d}.vtk')

    d1 = read_vtk_ascii(f1, fields=['h'])
    d2 = read_vtk_ascii(f2, fields=['h'])
    print(f"Loaded {len(d1['h'])} points from {f1} at t={t1:.1f}h")
    print(f"Loaded {len(d2['h'])} points from {f2} at t={t2:.1f}h")

    dist1 = great_circle_km(LON0, LAT0, d1['lon'], d1['lat'])
    dist2 = great_circle_km(LON0, LAT0, d2['lon'], d2['lat'])
    dev1 = d1['h'] - BACKGROUND
    dev2 = d2['h'] - BACKGROUND

    centers, mean1 = radial_bin_mean(dist1, dev1, BINWIDTH_KM)
    _, mean2 = radial_bin_mean(dist2, dev2, BINWIDTH_KM)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(centers, mean1, '-o', label=f't={t1:.1f}h')
    ax.plot(centers, mean2, '-s', label=f't={t2:.1f}h')
    ax.set_xlabel('great-circle distance from hump center (km)')
    ax.set_ylabel('mean layer-1 thickness deviation (m)')
    ax.set_title('Gravity-wave radial profile')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    # Wave-speed estimate: a mode's front position is the outermost radius
    # where the signal is still distinguishable from the background at that
    # mode's characteristic amplitude -- i.e. the leading edge -- not
    # wherever |deviation| happens to be locally largest. The barotropic
    # mode uses a small threshold (its own signal is only ~0.05m) and the
    # baroclinic mode a much larger one (its signal is O(10s of m)), so each
    # picks out its own front rather than the other's.
    def leading_edge(centers, mean, threshold):
        absmean = np.nan_to_num(np.abs(mean))
        above = np.where(absmean > threshold)[0]
        if len(above) == 0:
            return np.nan
        i = above[-1]
        # Linear-interpolate the threshold crossing between bin i (above
        # threshold) and bin i+1 (below) for sub-bin precision -- snapping
        # to the raw bin edge otherwise makes the speed estimate sensitive
        # to bin width (e.g. a single 500km-bin shift over an 8h window is
        # already a ~17 m/s error).
        if i == len(centers) - 1 or np.isnan(mean[i + 1]):
            return centers[i]
        y0, y1 = absmean[i], absmean[i + 1]
        if y0 == y1:
            return centers[i]
        frac = (y0 - threshold) / (y0 - y1)
        return centers[i] + frac * (centers[i + 1] - centers[i])

    r1_btp = leading_edge(centers, mean1, BTP_THRESHOLD)
    r2_btp = leading_edge(centers, mean2, BTP_THRESHOLD)
    speed_btp = (r2_btp - r1_btp) * 1000.0 / ((t2 - t1) * 3600.0)

    r1_bcl = leading_edge(centers, mean1, BCL_THRESHOLD)
    r2_bcl = leading_edge(centers, mean2, BCL_THRESHOLD)
    speed_bcl = (r2_bcl - r1_bcl) * 1000.0 / ((t2 - t1) * 3600.0)

    print(f"\nBarotropic front: t={t1:.1f}h -> r={r1_btp:.0f} km,  t={t2:.1f}h -> r={r2_btp:.0f} km")
    print(f"Estimated barotropic speed: {speed_btp:.1f} m/s  (predicted ~{PREDICTED_BTP_SPEED:.1f} m/s)")

    print(f"\nBaroclinic front: t={t1:.1f}h -> r={r1_bcl:.0f} km,  t={t2:.1f}h -> r={r2_bcl:.0f} km")
    print(f"Estimated baroclinic speed: {speed_bcl:.1f} m/s  (predicted ~{PREDICTED_BCL_SPEED:.1f} m/s)")
    print("(Baroclinic estimate is unreliable at short run times -- see BCL_THRESHOLD comment above; "
          "Chen's paper needed ~150-225h of simulated time for a clean baroclinic-front measurement.)")

    # Lon-lat maps (Chen Fig 5 style)
    fig2, axes = plt.subplots(1, 2, figsize=(13, 5))
    plot_latlon_map(axes[0], d1['lon'], d1['lat'], d1['h'],
                     title=f'layer 1 thickness, t={t1:.1f}h',
                     mark_point=(np.degrees(LON0), np.degrees(LAT0)))
    plot_latlon_map(axes[1], d2['lon'], d2['lat'], d2['h'],
                     title=f'layer 1 thickness, t={t2:.1f}h',
                     mark_point=(np.degrees(LON0), np.degrees(LAT0)))
    fig2.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
