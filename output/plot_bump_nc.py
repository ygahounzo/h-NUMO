
"""
Requirements:
  pip install numpy scipy matplotlib netCDF4

FFmpeg is required for MP4 writing:
  - macOS (brew): brew install ffmpeg
  - Linux (apt):  sudo apt-get install ffmpeg
"""

# Python port of the MATLAB animation script
# deps: pip install numpy scipy matplotlib netCDF4

import os
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib import cm

# -------------------- user params --------------------
file_num_start = 0
file_num_end   = 108
delta_file_num = 1

name_root   = "./" # root folder for input NetCDF files
save_folder = "./Animation"
out_root    = os.path.join(save_folder, "bump_2layers")  # output file (no extension)
fps         = 10
nxx, nyy    = 201, 201
# -----------------------------------------------------

os.makedirs(save_folder, exist_ok=True)

# Figure ~ 1600x1200 px
fig = plt.figure(figsize=(16, 12), dpi=100)
ax = fig.add_subplot(111, projection="3d")

plt.rcParams.update({"font.size": 16, "font.family": "Times New Roman"})

writer = FFMpegWriter(fps=fps, metadata={"title": "Layer interfaces", "artist": "python"})

def read_nc(path):
    with Dataset(path, "r") as ds:
        dt = ds.variables["dt"][:]
        x  = ds.variables["x"][:].astype(float)
        y  = ds.variables["y"][:].astype(float)
        z  = ds.variables["eta"][:].astype(float)  # shape (N, 3)
    dt = float(np.array(dt).squeeze())
    return dt, x, y, z

# Build list of frames
frames = list(range(file_num_start, file_num_end + 1, delta_file_num))

# Grid (determined from the first file)
first_path = os.path.join(name_root, f"mlswe{frames[0]:04d}.nc")
dt0, xe0, ye0, _ = read_nc(first_path)

xmin, xmax = float(np.min(xe0)), float(np.max(xe0))
ymin, ymax = float(np.min(ye0)), float(np.max(ye0))
dx = (xmax - xmin) / nxx
dy = (ymax - ymin) / nyy
xg = np.arange(xmin, xmax + 0.5*dx, dx)
yg = np.arange(ymin, ymax + 0.5*dy, dy)
XI, YI = np.meshgrid(xg, yg, indexing="xy")

with writer.saving(fig, out_root + ".mp4", dpi=100):
    for ifile in frames:
        print('ifile = ',ifile)
        nc_path = os.path.join(name_root, f"mlswe{ifile:04d}.nc")
        dt, xe, ye, z = read_nc(nc_path)

        z1 = z[0,:,]
        z2 = (z[1,:,] + 20.0) * 20.0 - 20.0
        z3 = z[2,:]

        pts = np.column_stack([xe, ye])
        Q1 = griddata(pts, z1, (XI, YI), method="cubic")
        Q2 = griddata(pts, z2, (XI, YI), method="cubic")
        Q3 = griddata(pts, z3, (XI, YI), method="linear")

        ax.cla()

        # Layer 1: custom RGB like MATLAB C0
        # G = 0.5 + 1e4*Q1 (clamped to [0,1] so facecolors are valid)
        if Q1 is not None:
            G = np.clip(0.5 + 1e4 * np.nan_to_num(Q1, nan=0.0), 0.0, 1.0)
            R = np.zeros_like(G)
            B = np.ones_like(G)
            facecolors1 = np.stack([R, G, B], axis=-1)
            ax.plot_surface(XI, YI, Q1, rstride=1, cstride=1,
                            facecolors=facecolors1, linewidth=0,
                            antialiased=False, shade=False, alpha=0.7)

        # Layer 2: use colormap, apply clim([-20 -15])
        if Q2 is not None:
            s2 = ax.plot_surface(XI, YI, Q2, rstride=1, cstride=1,
                                 cmap=cm.viridis, linewidth=0,
                                 antialiased=False, shade=False)
            s2.set_clim(vmin=-20, vmax=-15)

        # Bottom topography: same cmap
        if Q3 is not None:
            ax.plot_surface(XI, YI, Q3, rstride=1, cstride=1,
                            cmap=cm.viridis, linewidth=0,
                            antialiased=False, shade=False)

        ax.set_zlim(-40.0, 0.0)
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")

        minutes = (ifile * dt) / 60.0
        ax.set_title(f"Layer interfaces at {minutes:.0f} min", pad=18)

        ax.view_init(elev=15, azim=-15)

        plt.tight_layout()
        writer.grab_frame()

plt.close(fig)
