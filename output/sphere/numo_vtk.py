"""
Shared VTK reading / geometry / plotting helpers for h-NUMO sphere test-case
analysis (geoFlow, gravity_wave, mountain_2layer), in the spirit of the
diagnostics in Chen (2025, QJRMS, 10.1002/qj.4994).

No third-party VTK library needed -- h-NUMO writes legacy ASCII VTK
(POLYDATA-style POINTS + POINT_DATA SCALARS/VECTORS blocks), which this
module parses directly with numpy for speed on the ~150k-point files.
"""

import re
import numpy as np

EARTH_R_KM = 6371.0
CIRCUMFERENCE_KM = 2.0 * np.pi * EARTH_R_KM

# h-NUMO's fixed-width Fortran ASCII output sometimes glues a negative
# number directly onto the previous one with no separating space (the '-'
# sign occupies the field position a separating space would otherwise take,
# e.g. "...E+007-0.514013997789297E+007"). A plain whitespace split (or
# np.fromstring(sep=' ')) mis-parses this, so floats must be extracted with
# an explicit regex instead.
_FLOAT_RE = re.compile(r'-?\d+\.\d+E[+-]\d+')


def _parse_floats(s):
    return np.array(_FLOAT_RE.findall(s), dtype=float)


def read_vtk_ascii(path, fields=None):
    """Parse a legacy-ASCII VTK file written by h-NUMO.

    Returns a dict with:
      'points' : (N,3) array of xyz coordinates (meters)
      'lon'    : (N,) longitude in radians, [-pi, pi]
      'lat'    : (N,) latitude in radians, [-pi/2, pi/2]
      one entry per requested scalar field name -> (N,) array

    If `fields` is None, all SCALARS fields present in the file are read.
    """
    with open(path, 'r') as f:
        text = f.read()

    m = re.search(r'POINTS\s+(\d+)\s+\S+\s*\n', text)
    npts = int(m.group(1))
    start = m.end()
    stop = text.find('CELLS', start)
    if stop == -1:
        stop = text.find('POLYGONS', start)
    if stop == -1:
        stop = text.find('POINT_DATA', start)
    pts = _parse_floats(text[start:stop])
    pts = pts[:3 * npts].reshape(npts, 3)

    pd = text.find('POINT_DATA')
    tail = text[pd:]

    available = re.findall(r'SCALARS\s+(\S+)\s+\S+\s+\d+\s*\nLOOKUP_TABLE\s+\S+\s*\n', tail)
    if fields is None:
        fields = available

    out = {'points': pts}
    for name in fields:
        pat = re.compile(r'SCALARS\s+' + re.escape(name) +
                          r'\s+\S+\s+\d+\s*\nLOOKUP_TABLE\s+\S+\s*\n')
        mm = pat.search(tail)
        if mm is None:
            raise KeyError(f"field '{name}' not found in {path} (available: {available})")
        sstart = mm.end()
        nxt = tail.find('SCALARS', sstart)
        nxtv = tail.find('VECTORS', sstart)
        candidates = [x for x in (nxt, nxtv) if x != -1]
        send = min(candidates) if candidates else len(tail)
        vals = _parse_floats(tail[sstart:send])
        out[name] = vals[:npts]

    r = np.linalg.norm(pts, axis=1)
    out['lon'] = np.arctan2(pts[:, 1], pts[:, 0])
    out['lat'] = np.arcsin(np.clip(pts[:, 2] / r, -1.0, 1.0))
    return out


def _fortran_nint(x):
    """Fortran NINT: round to nearest integer, ties away from zero (Python's
    round() uses banker's rounding/ties-to-even instead, which disagrees with
    Fortran exactly at .5 ties -- e.g. round(4.5)=4 in Python vs NINT(4.5)=5
    in Fortran)."""
    return int(np.floor(x + 0.5)) if x >= 0 else int(np.ceil(x - 0.5))


def snapshot_time_hours(idx, dt_seconds, time_restart_hours):
    """Actual elapsed time (hours) for VTK output snapshot `idx`.

    h-NUMO writes a snapshot every `irestart = NINT(time_restart_seconds /
    dt)` steps (mod_time_loop.F90), NOT every `time_restart_hours` of
    physical time exactly -- those only coincide when dt evenly divides
    time_restart_seconds (1800s for the usual time_restart=0.5h). For any
    other dt the actual snapshot spacing drifts from the nominal value (e.g.
    dt=400s with time_restart=0.5h writes every irestart=5 steps = 2000s =
    0.5556h, not 0.5h) -- using `idx * time_restart_hours` directly silently
    gives the wrong time and, e.g., biases any wave-speed-from-snapshots
    calculation.
    """
    time_restart_seconds = time_restart_hours * 3600.0
    irestart = _fortran_nint(time_restart_seconds / dt_seconds)
    if irestart == 0:
        irestart = 1
    actual_interval_hours = irestart * dt_seconds / 3600.0
    return idx * actual_interval_hours


def great_circle_km(lon0, lat0, lon, lat, R_km=EARTH_R_KM):
    """Great-circle distance (km) from a fixed point (lon0,lat0) [radians]
    to arrays of (lon,lat) [radians]."""
    cosang = np.sin(lat0) * np.sin(lat) + np.cos(lat0) * np.cos(lat) * np.cos(lon - lon0)
    return R_km * np.arccos(np.clip(cosang, -1.0, 1.0))


def radial_bin_mean(dist_km, values, binwidth_km=500.0, max_km=None):
    """Bin `values` by `dist_km` into fixed-width bins; return bin centers
    and the mean value in each bin (nan where empty)."""
    if max_km is None:
        max_km = np.pi * EARTH_R_KM
    edges = np.arange(0.0, max_km + binwidth_km, binwidth_km)
    centers = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(dist_km, edges) - 1
    means = np.full(len(centers), np.nan)
    for b in range(len(centers)):
        sel = idx == b
        if np.any(sel):
            means[b] = values[sel].mean()
    return centers, means


def plot_latlon_map(ax, lon, lat, values, title=None, cmap='RdBu_r',
                     levels=41, mark_point=None, **kwargs):
    """Flattened lon-lat (Plate-Carree-style) filled contour map of a
    scattered field on the sphere, using matplotlib's unstructured
    tricontourf (no gridding / cartopy needed).

    mark_point: optional (lon_deg, lat_deg) tuple to mark with a black 'x'
    (e.g. the mountain peak or the gravity-wave hump center).
    """
    lon_deg = np.degrees(lon)
    lat_deg = np.degrees(lat)
    cf = ax.tricontourf(lon_deg, lat_deg, values, levels=levels, cmap=cmap, **kwargs)
    ax.set_xlabel('longitude (deg)')
    ax.set_ylabel('latitude (deg)')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    if mark_point is not None:
        ax.plot(mark_point[0], mark_point[1], 'kx', markersize=10, markeredgewidth=2)
    if title:
        ax.set_title(title)
    return cf


def read_mass_cons(path, dt_seconds, nlayers=2):
    """Parse h-NUMO's mass_mlswe.cons diagnostic file.

    Format: one row per diagnostic write, columns = [itime, mass_layer1,
    mass_layer2, ...]. itime is the outer (bcl) time-step counter, so
    physical time (hours) = itime * dt_seconds / 3600.

    Returns (time_hours, mass) where mass has shape (n_rows, nlayers).
    """
    data = np.loadtxt(path)
    itime = data[:, 0]
    mass = data[:, 1:1 + nlayers]
    time_hours = itime * dt_seconds / 3600.0
    return time_hours, mass
