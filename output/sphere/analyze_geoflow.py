#!/usr/bin/env python3
"""
Analysis for h-NUMO's 'geoFlow' test case (adaptation of Williamson et al.
1992 test #2 / Chen 2025 Sec 5.1: steady zonal geostrophic flow).

Chen's own diagnostic for this test (his Figs 2-4) is L2/Linf error in
thickness, vorticity, and divergence relative to a *closed-form analytic
solution* -- but that solution only holds exactly under his Boussinesq-type
approximation (uniform reference density across layers). h-NUMO's geoFlow
uses each layer's REAL density (no Boussinesq substitution -- see project
history), so there is no matching closed-form target to diff against.

The appropriate analog here is a steadiness/drift diagnostic: since the
flow is initialized in (approximate) geostrophic balance, a well-behaved
model should stay close to its initial state over time. This script tracks,
per layer, the RMS drift of thickness and horizontal speed away from their
t=0 values as a function of time -- growth here indicates loss of balance
(expected to leading order from real density stratification -- see
project notes on the layer-2 drift investigation -- but should saturate,
not diverge).

Usage:
    python3 analyze_geoflow.py [rundir] [nlayers]
"""

import sys
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

from numo_vtk import read_vtk_ascii, plot_latlon_map, snapshot_time_hours

# ---------------- SETTINGS ----------------
RUNDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
NLAYERS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
FPREFIX_TEMPLATE = 'mlswe_sphere_*_dg_rk35_l{layer:03d}_'
TIME_RESTART_HOURS = 0.5   # time_restart in numo3d.in
# -------------------------------------------


def find_snapshots(rundir, layer):
    pattern = os.path.join(rundir, FPREFIX_TEMPLATE.format(layer=layer) + '*.vtk')
    files = sorted(glob.glob(pattern))
    idx = [int(re.search(r'_(\d+)\.vtk$', f).group(1)) for f in files]
    order = np.argsort(idx)
    files, idx = [files[i] for i in order], [idx[i] for i in order]
    # dt (bcl time step) is embedded in the filename itself, e.g.
    # "mlswe_sphere_400.0000_dg_rk35_...". Snapshot spacing is
    # irestart*dt (irestart = NINT(time_restart_seconds/dt)), which only
    # equals TIME_RESTART_HOURS exactly when dt evenly divides 1800s -- so
    # dt must be known to label snapshot times correctly (see
    # numo_vtk.snapshot_time_hours).
    dt_seconds = float(re.search(r'_(\d+\.\d+)_dg', files[0]).group(1)) if files else None
    return files, idx, dt_seconds


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for layer in range(1, NLAYERS + 1):
        files, indices, dt_seconds = find_snapshots(RUNDIR, layer)
        if not files:
            print(f"layer {layer}: no snapshots found, skipping")
            continue

        d0 = read_vtk_ascii(files[0], fields=['h', 'u', 'v'])
        h0, u0, v0 = d0['h'], d0['u'], d0['v']
        speed0 = np.hypot(u0, v0)

        times, h_drift, speed_drift = [], [], []
        for f, idx in zip(files, indices):
            d = read_vtk_ascii(f, fields=['h', 'u', 'v'])
            speed = np.hypot(d['u'], d['v'])
            times.append(snapshot_time_hours(idx, dt_seconds, TIME_RESTART_HOURS))
            h_drift.append(np.sqrt(np.mean((d['h'] - h0) ** 2)))
            speed_drift.append(np.sqrt(np.mean((speed - speed0) ** 2)))

        print(f"layer {layer}: {len(files)} snapshots, "
              f"final RMS thickness drift = {h_drift[-1]:.3f} m, "
              f"final RMS speed drift = {speed_drift[-1]:.3f} m/s")

        axes[0].plot(times, h_drift, '-o', label=f'layer {layer}')
        axes[1].plot(times, speed_drift, '-o', label=f'layer {layer}')

    axes[0].set_xlabel('time (h)')
    axes[0].set_ylabel('RMS thickness drift from t=0 (m)')
    axes[0].set_title('geoFlow: thickness drift')
    axes[0].legend()
    axes[0].grid(True)

    axes[1].set_xlabel('time (h)')
    axes[1].set_ylabel('RMS speed drift from t=0 (m/s)')
    axes[1].set_title('geoFlow: horizontal-speed drift')
    axes[1].legend()
    axes[1].grid(True)

    fig.tight_layout()

    # Final-time lon-lat thickness maps, one panel per layer
    fig2, axes2 = plt.subplots(1, NLAYERS, figsize=(6 * NLAYERS, 5), squeeze=False)
    for layer in range(1, NLAYERS + 1):
        files, indices, dt_seconds = find_snapshots(RUNDIR, layer)
        if not files:
            continue
        d_final = read_vtk_ascii(files[-1], fields=['h'])
        t_final = snapshot_time_hours(indices[-1], dt_seconds, TIME_RESTART_HOURS)
        plot_latlon_map(axes2[0][layer - 1], d_final['lon'], d_final['lat'], d_final['h'],
                         title=f'layer {layer} thickness, t={t_final:.1f}h')
    fig2.tight_layout()

    plt.show()


if __name__ == '__main__':
    main()
