#
# pick a random “slice” along the direct ray (excluding ±5° endcaps)
# pick a random angle phi in the cross-section
#
# pick a random radius rho in the cross-section
# get (lat, lon, depth) for that scatterer point
# compute the two legs + compare dt/dbaz/dp constraints

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple, Dict
import sys

from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

EARTH_R_KM = 6371.0  # spherical approx for geometry sampling

def wrap180(a):
    return (a + 180.0) % 360.0 - 180.0

def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v

def latlon_to_unit(lat_deg, lon_deg):
    lat = np.deg2rad(lat_deg); lon = np.deg2rad(lon_deg)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.array([x, y, z], float)

def unit_to_latlon(u):
    x, y, z = u
    lat = np.rad2deg(np.arcsin(np.clip(z, -1, 1)))
    lon = np.rad2deg(np.arctan2(y, x))
    lon = (lon + 540.0) % 360.0 - 180.0  # to [-180,180)
    return float(lat), float(lon)

def rotate_about_axis(v, axis, ang_rad):
    # Rodrigues’ rotation
    axis = unit(axis)
    return (v*np.cos(ang_rad) +
            np.cross(axis, v)*np.sin(ang_rad) +
            axis*np.dot(axis, v)*(1 - np.cos(ang_rad)))

def first_P_time_p_sdeg(model: TauPyModel, delta_deg: float, z_src: float, z_rcv: float) -> Tuple[float, float]:
    arr = model.get_travel_times(
        source_depth_in_km=float(z_src),
        receiver_depth_in_km=float(z_rcv),
        distance_in_degree=float(delta_deg),
        phase_list=["P"]
    )
    if not arr:
        raise ValueError("No P")
    a0 = arr[0]
    return float(a0.time), float(a0.ray_param_sec_degree)  # s, s/deg  :contentReference[oaicite:2]{index=2}

def in_range_and_bin(x, lo, hi, step):
    if x < lo or x > hi:
        return False, np.nan
    k = round((x - lo)/step)
    xb = lo + k*step
    ok = abs(x - xb) <= 0.5*step + 1e-12
    return ok, float(xb)

@dataclass
class Candidate:
    scat_lat: float
    scat_lon: float
    scat_depth_km: float
    dt_obs_s: float
    dp_obs_sdeg: float
    dbaz_obs_deg: float
    dt_bin_s: float
    dp_bin_sdeg: float
    dbaz_bin_deg: float
    t_total_s: float
    p2_sdeg: float
    baz_in_deg: float
    delta1_deg: float
    delta2_deg: float

def sample_tube_point_on_direct_ray(
    src_lat, src_lon, rcv_lat, rcv_lon,
    direct_path_dist_deg: np.ndarray,
    direct_path_depth_km: np.ndarray,
    *,
    tube_radius_km: float,
    min_endcap_deg: float,
    delta_sr_deg: float,
    rng: np.random.Generator
) -> Tuple[float, float, float]:
    """
    Sample a 3D point within a tube of radius tube_radius_km around the direct-ray curve.
    Returns (lat, lon, depth_km) of sampled interior point.
    """
    # Great-circle plane axis
    s_hat = latlon_to_unit(src_lat, src_lon)
    r_hat = latlon_to_unit(rcv_lat, rcv_lon)
    axis = unit(np.cross(s_hat, r_hat))
    if np.linalg.norm(axis) == 0:
        raise ValueError("Source and receiver are antipodal or identical; axis ill-defined.")

    # Choose a ray-path index away from the ends in terms of surface distance (deg from source)
    # direct_path_dist_deg is "deg from source" along ray.
    valid = np.where((direct_path_dist_deg >= min_endcap_deg) &
                     (direct_path_dist_deg <= (delta_sr_deg - min_endcap_deg)))[0]
    if len(valid) < 3:
        raise ValueError("Not enough ray-path points after applying endcaps.")
    i = int(rng.choice(valid[1:-1]))

    # 3D point on the *central* direct ray curve at that dist/depth:
    ang = np.deg2rad(direct_path_dist_deg[i])
    u_surf = rotate_about_axis(s_hat, axis, ang)          # unit vector on great circle
    rad = EARTH_R_KM - direct_path_depth_km[i]
    p0 = rad * u_surf                                      # 3D position (km)

    # Tangent direction along the ray curve (approx from neighbors)
    ang_prev = np.deg2rad(direct_path_dist_deg[i-1])
    ang_next = np.deg2rad(direct_path_dist_deg[i+1])
    u_prev = rotate_about_axis(s_hat, axis, ang_prev)
    u_next = rotate_about_axis(s_hat, axis, ang_next)
    p_prev = (EARTH_R_KM - direct_path_depth_km[i-1]) * u_prev
    p_next = (EARTH_R_KM - direct_path_depth_km[i+1]) * u_next
    t_hat = unit(p_next - p_prev)

    # Build orthonormal basis (n1, n2) perpendicular to tangent
    # pick a helper not parallel to t_hat
    helper = axis if abs(np.dot(axis, t_hat)) < 0.9 else u_surf
    n1 = unit(np.cross(t_hat, helper))
    n2 = unit(np.cross(t_hat, n1))

    # Sample inside tube cross-section (uniform area)
    rho = tube_radius_km * np.sqrt(rng.uniform(0, 1))
    phi = rng.uniform(0, 2*np.pi)
    offset = rho * (np.cos(phi)*n1 + np.sin(phi)*n2)

    p = p0 + offset
    r = np.linalg.norm(p)
    depth = EARTH_R_KM - r
    if depth < 0:  # above surface -> reject by resampling in caller
        return np.nan, np.nan, np.nan

    lat, lon = unit_to_latlon(p / r)
    return lat, lon, float(depth)

def find_scatterers_tube(
    src_lat, src_lon, src_depth_km,
    rcv_lat, rcv_lon,
    *,
    model_name="iasp91",
    tube_radius_km=300.0,
    min_endcap_deg=5.0,
    n_samples=500_000,

    # Ranges + steps (your request)
    dt_range_s=(50.0, 100.0), dt_step_s=2.0,
    dbaz_range_deg=(-5.0, 5.0), dbaz_step_deg=1.0,
    dp_range_sdeg=(0.0, 2.0), dp_step_sdeg=0.5,
) -> Dict[str, object]:

    model = TauPyModel(model=model_name)

    # Source-receiver geometry (for direct baz and Δ)
    dist_m, az_sr, baz_rs = gps2dist_azimuth(src_lat, src_lon, rcv_lat, rcv_lon)
    delta_sr_deg = kilometers2degrees(dist_m / 1000.0)
    baz_dir_in = float(baz_rs)  # receiver->source, 0..360

    # Direct P (time/slowness reference)
    tP_dir, pP_dir = first_P_time_p_sdeg(model, delta_sr_deg, src_depth_km, 0.0)

    # Direct ray path (dist, depth) to define the tube axis curve
    rps = model.get_ray_paths(source_depth_in_km=src_depth_km,
                              distance_in_degree=delta_sr_deg,
                              phase_list=["P"])
    if not rps:
        raise ValueError("No direct P ray path")
    rp = rps[0]

    path_dist_rad = np.asarray(rp.path["dist"], float)
    path_dist_deg = np.rad2deg(path_dist_rad)
    path_depth_km = np.asarray(rp.path["depth"], float)

    rng = np.random.default_rng(0)
    cands: List[Candidate] = []

    dt_lo, dt_hi = dt_range_s
    db_lo, db_hi = dbaz_range_deg
    dp_lo, dp_hi = dp_range_sdeg

    for _ in range(n_samples):
        scat_lat, scat_lon, scat_z = sample_tube_point_on_direct_ray(
            src_lat, src_lon, rcv_lat, rcv_lon,
            path_dist_deg, path_depth_km,
            tube_radius_km=tube_radius_km,
            min_endcap_deg=min_endcap_deg,
            delta_sr_deg=delta_sr_deg,
            rng=rng)
        if not np.isfinite(scat_z):
            continue

        # Epicentral distances for each leg (use angular separation of radial vectors)
        # (TauP's "distance_in_degree" is epicentral distance; receiver_depth handles depth.) :contentReference[oaicite:3]{index=3}
        d1 = kilometers2degrees(gps2dist_azimuth(src_lat, src_lon, scat_lat, scat_lon)[0] / 1000.0)
        d2 = kilometers2degrees(gps2dist_azimuth(scat_lat, scat_lon, rcv_lat, rcv_lon)[0] / 1000.0)

        if d1 < min_endcap_deg or d2 < min_endcap_deg:
            continue

        # Two-leg TauP
        try:
            t1, _p1 = first_P_time_p_sdeg(model, d1, src_depth_km, scat_z)
            t2, p2 = first_P_time_p_sdeg(model, d2, scat_z, 0.0)
        except ValueError:
            continue

        ttot = t1 + t2
        dt_obs = ttot - tP_dir
        dp_obs = p2 - pP_dir
        baz_in = gps2dist_azimuth(rcv_lat, rcv_lon, scat_lat, scat_lon)[1]  # receiver->scatterer azimuth
        dbaz_obs = wrap180(baz_in - baz_dir_in)

        ok_t, dt_bin = in_range_and_bin(dt_obs, dt_lo, dt_hi, dt_step_s)
        if not ok_t: continue
        ok_p, dp_bin = in_range_and_bin(dp_obs, dp_lo, dp_hi, dp_step_sdeg)
        if not ok_p: continue
        ok_b, db_bin = in_range_and_bin(dbaz_obs, db_lo, db_hi, dbaz_step_deg)
        if not ok_b: continue

        cands.append(Candidate(
            scat_lat=scat_lat, scat_lon=scat_lon, scat_depth_km=scat_z,
            dt_obs_s=dt_obs, dp_obs_sdeg=dp_obs, dbaz_obs_deg=dbaz_obs,
            dt_bin_s=dt_bin, dp_bin_sdeg=dp_bin, dbaz_bin_deg=db_bin,
            t_total_s=ttot, p2_sdeg=p2, baz_in_deg=baz_in,
            delta1_deg=d1, delta2_deg=d2
        ))

    return {
        "direct": {
            "delta_sr_deg": float(delta_sr_deg),
            "tP_s": float(tP_dir),
            "pP_sdeg": float(pP_dir),
            "baz_in_deg": float(baz_dir_in),
        },
        "tube": {
            "tube_radius_km": float(tube_radius_km),
            "min_endcap_deg": float(min_endcap_deg),
            "n_samples": int(n_samples),
        },
        "ranges": {
            "dt_range_s": dt_range_s, "dt_step_s": float(dt_step_s),
            "dbaz_range_deg": dbaz_range_deg, "dbaz_step_deg": float(dbaz_step_deg),
            "dp_range_sdeg": dp_range_sdeg, "dp_step_sdeg": float(dp_step_sdeg),
        },
        "candidates": cands,
    }
######
# sys.exit()
if __name__ == "__main__":

    src_lat, src_lon, src_depth_km = 10.0, 20.0, 300.0
    rcv_lat, rcv_lon = 35.0, 90.0

    out = find_scatterers_tube(
        src_lat, src_lon, src_depth_km,
        rcv_lat, rcv_lon,
        model_name="iasp91",

        # tube geometry
        tube_radius_km=300.0,      # radius of the tube around the ray (km)
        min_endcap_deg=5.0,        # don't sample within 5° of src/rcv
        n_samples=200_000,         # number of random tube points to test

        # constraints (your ranges/steps)
        dt_range_s=(50.0, 100.0),  dt_step_s=2.0,
        dbaz_range_deg=(-5.0, 5.0), dbaz_step_deg=1.0,
        dp_range_sdeg=(0.0, 3.0),  dp_step_sdeg=0.5,
    )

    print("Direct (reference):", out["direct"])
    print("Tube:", out["tube"])
    print("Ranges:", out["ranges"])
    print("Matches:", len(out["candidates"]))

    # Print a few matches
    for c in out["candidates"][:10]:
        print(
            f"scat(lat,lon,z)=({c.scat_lat:+.3f},{c.scat_lon:+.3f},{c.scat_depth_km:.1f} km) "
            f"Δ1={c.delta1_deg:.2f} Δ2={c.delta2_deg:.2f} "
            f"dt={c.dt_obs_s:.2f}s(dpbin {c.dt_bin_s:.0f}) "
            f"dp={c.dp_obs_sdeg:.3f}s/deg(dpbin {c.dp_bin_sdeg:.1f}) "
            f"dbaz={c.dbaz_obs_deg:+.2f}°(bin {c.dbaz_bin_deg:+.0f}) "
            f"t_total={c.t_total_s:.2f}s p2={c.p2_sdeg:.3f}s/deg baz_in={c.baz_in_deg:.2f}°"
        )
