#!/usr/bin/env python3
"""
Deterministic tube/cylinder sampling around the direct P ray path using:
  - deterministic slices along ray in degrees (ddeg)
  - deterministic Fibonacci (sunflower) disk points in the plane perpendicular to local ray tangent

For each sampled 3D scatterer point (lat, lon, depth_km), compute two-leg P travel times:
  leg1: source(z_src) -> scatter(z_scat) at distance Δ1
  leg2: scatter(z_scat) -> receiver(0) at distance Δ2

And check constraints in ranges with step/binning:
  dt in [dt_min, dt_max] step dt_step
  dbaz in [dbaz_min, dbaz_max] step dbaz_step
  dp in [dp_min, dp_max] step dp_step

Requires:
  pip install obspy numpy geographiclib
"""

from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import sys
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

EARTH_R_KM = 6371.0  # spherical geometry approximation for sampling


# ----------------------------
# Small utilities
# ----------------------------
def wrap180(a_deg: float) -> float:
    """Wrap an angle to [-180, 180)."""
    return (a_deg + 180.0) % 360.0 - 180.0


def unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n == 0.0:
        return v
    return v / n


def latlon_to_unit(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return np.array([math.cos(lat) * math.cos(lon),
                     math.cos(lat) * math.sin(lon),
                     math.sin(lat)], dtype=float)


def unit_to_latlon(u: np.ndarray) -> Tuple[float, float]:
    x, y, z = float(u[0]), float(u[1]), float(u[2])
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon = math.degrees(math.atan2(y, x))
    # normalize lon to [-180, 180)
    lon = (lon + 540.0) % 360.0 - 180.0
    return lat, lon


def rotate_about_axis(v: np.ndarray, axis: np.ndarray, ang_rad: float) -> np.ndarray:
    """
    Rodrigues rotation formula: rotate vector v about unit axis by ang_rad.
    """
    axis = unit(axis)
    c = math.cos(ang_rad)
    s = math.sin(ang_rad)
    return v * c + np.cross(axis, v) * s + axis * float(np.dot(axis, v)) * (1.0 - c)


def in_range_and_bin(x: float, lo: float, hi: float, step: float) -> Tuple[bool, float]:
    """
    Deterministic binning: accept if x is within [lo, hi] AND within half-step of the nearest bin.
    Returns (ok, x_bin).
    """
    if x < lo or x > hi:
        return False, float("nan")
    k = round((x - lo) / step)
    xb = lo + k * step
    ok = abs(x - xb) <= 0.5 * step + 1e-12
    return ok, float(xb)


# def first_P_time_p_sdeg(model: TauPyModel, delta_deg: float, z_src: float, z_rcv: float) -> Tuple[float, float]:
#     """
#     First-arrival P travel time and ray parameter in s/deg for a path between depths.
#     """
#     arr = model.get_travel_times(
#         source_depth_in_km=float(z_src),
#         receiver_depth_in_km=float(z_rcv),
#         distance_in_degree=float(delta_deg),
#         phase_list=["P","p"],
#     )
#     if not arr:
#         raise ValueError("No P")
#     a0 = arr[0]
#     return float(a0.time), float(a0.ray_param_sec_degree)  # s, s/deg

def first_P_time_p_sdeg(
    model: TauPyModel,
    delta_deg: float,
    z_src: float,
    z_rcv: float,
    *,
    phase_list=("P", "p"),
    try_swap: bool = True,) -> Tuple[float, float]:
    """
    Robust first-arrival compressional time + ray parameter (s/deg) between two depths.

    TauP can return empty for phase_list=["P","p"] when z_src < z_rcv (deep "receiver").
    This function tries (z_src -> z_rcv); if empty and try_swap=True, tries (z_rcv -> z_src).
    Returns the earliest-time arrival among the attempted orderings.

    Note: ray_param is path-invariant; for leg-1 you usually only need time anyway.
    """
    best_time: Optional[float] = None
    best_p: Optional[float] = None

    def _try(zA: float, zB: float):
        arr = model.get_travel_times(
            source_depth_in_km=float(zA),
            receiver_depth_in_km=float(zB),
            distance_in_degree=float(delta_deg),
            phase_list=list(phase_list),
        )
        if not arr:
            return None
        a0 = arr[0]  # already earliest time
        return float(a0.time), float(a0.ray_param_sec_degree)

    out = _try(z_src, z_rcv)
    if out is not None:
        best_time, best_p = out

    if try_swap:
        out2 = _try(z_rcv, z_src)
        if out2 is not None:
            t2, p2 = out2
            if (best_time is None) or (t2 < best_time):
                best_time, best_p = t2, p2

    if best_time is None or best_p is None:
        raise ValueError(
            f"No P/p found for Δ={delta_deg:.2f}°, z_src={z_src:.1f} km, z_rcv={z_rcv:.1f} km "
            f"(tried swap={try_swap})"
        )

    return best_time, best_p

def central_angle_deg(lat1, lon1, lat2, lon2):
    u1 = latlon_to_unit(lat1, lon1)
    u2 = latlon_to_unit(lat2, lon2)
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(u1, u2))))))

# ----------------------------
# Data structures
# ----------------------------
@dataclass
class Candidate:
    scat_lat: float
    scat_lon: float
    scat_depth_km: float

    delta1_deg: float
    delta2_deg: float

    t_total_s: float
    p2_sdeg: float
    baz_in_deg: float

    dt_obs_s: float
    dp_obs_sdeg: float
    dbaz_obs_deg: float

    dt_bin_s: float
    dp_bin_sdeg: float
    dbaz_bin_deg: float

    # for debugging
    slice_dist_deg: float
    rho_km: float
    phi_deg: float


# ----------------------------
# Fibonacci disk in cross-section
# ----------------------------
def fibonacci_disk_points(N: int, R_km: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns arrays (rho_km, phi_rad, phi_deg) of length N.

    rho_k = R * sqrt((k+0.5)/N)  -> roughly uniform in area
    phi_k = k * golden_angle
    """
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))  # ~2.39996323 rad
    k = np.arange(N, dtype=float)
    rho = R_km * np.sqrt((k + 0.5) / N)
    phi = (k * golden_angle) % (2.0 * math.pi)
    phi_deg = np.degrees(phi)
    return rho, phi, phi_deg

# ----------------------------
# Core deterministic sampler
# ----------------------------
def find_scatterers_tube_deterministic(
    src_lat: float, src_lon: float, src_depth_km: float,
    rcv_lat: float, rcv_lon: float,
    *,
    model_name: str = "iasp91",

    # tube sampling controls
    tube_radius_km: float = 300.0,
    min_endcap_deg: float = 5.0,
    ddeg: float = 0.2,          # along-ray slice step in degrees
    N_disk: int = 200,          # Fibonacci points per slice

    # constraint ranges + steps
    dt_range_s: Tuple[float, float] = (50.0, 100.0),
    dt_step_s: float = 2.0,
    dbaz_range_deg: Tuple[float, float] = (-5.0, 5.0),
    dbaz_step_deg: float = 1.0,
    dp_range_sdeg: Tuple[float, float] = (0.0, 2.0),
    dp_step_sdeg: float = 0.5,

    # safety
    max_candidates: Optional[int] = None,  # stop after collecting this many matches (None = no cap)
) -> Dict[str, object]:
    """
    Deterministically sample a tube around the direct P ray and return matching scatterers.
    """

    model = TauPyModel(model=model_name)

    # Source-receiver geometry
    dist_m, az_sr, baz_rs = gps2dist_azimuth(src_lat, src_lon, rcv_lat, rcv_lon)
    delta_sr_deg = kilometers2degrees(dist_m / 1000.0)
    baz_dir_in = float(baz_rs)  #

    if delta_sr_deg <= 2.0 * min_endcap_deg:
        raise ValueError(f"Δ_sr={delta_sr_deg:.2f}° too small for endcaps ±{min_endcap_deg}°")

    # Direct P reference
    tP_dir, pP_dir = first_P_time_p_sdeg(model, delta_sr_deg, src_depth_km, 0.0)

    # Direct ray path
    rps = model.get_ray_paths(
        source_depth_in_km=float(src_depth_km),
        distance_in_degree=float(delta_sr_deg),
        phase_list=["P"])
    if not rps:
        raise ValueError("No direct P ray path")
    rp = rps[0]

    # ObsPy TauP ray-path distances
    path_dist_rad = np.asarray(rp.path["dist"], dtype=float)
    path_dist_deg = np.degrees(path_dist_rad)
    path_depth_km = np.asarray(rp.path["depth"], dtype=float)

    # Ensure monotonic dist for interpolation (usually is)
    # If there are repeats, keep last occurrence
    _, idx_rev = np.unique(path_dist_deg[::-1], return_index=True)
    keep = (len(path_dist_deg) - 1 - idx_rev)[::-1]
    keep = np.sort(keep)
    path_dist_deg = path_dist_deg[keep]
    path_depth_km = path_depth_km[keep]

    # Define slice distances deterministically
    start = min_endcap_deg
    stop = delta_sr_deg - min_endcap_deg
    # Clamp stop to available path max (rare edge case)
    stop = min(stop, float(path_dist_deg[-1]))

    slice_dists = np.arange(start, stop + 0.5 * ddeg, ddeg)
    if slice_dists.size < 3:
        raise ValueError("Not enough slices; increase Δ or reduce endcap/step.")

    # Fibonacci disk points (deterministic)
    rho_km, phi_rad, phi_deg = fibonacci_disk_points(N_disk, tube_radius_km)

    # Great-circle axis (defines the source->receiver plane)
    s_hat = latlon_to_unit(src_lat, src_lon)
    r_hat = latlon_to_unit(rcv_lat, rcv_lon)
    axis = unit(np.cross(s_hat, r_hat))
    if float(np.linalg.norm(axis)) < 1e-12:
        raise ValueError("Source and receiver nearly colinear/antipodal; great-circle axis ill-defined.")

    # Ranges
    dt_lo, dt_hi = dt_range_s
    db_lo, db_hi = dbaz_range_deg
    dp_lo, dp_hi = dp_range_sdeg

    cands: List[Candidate] = []

    # Helper for along-ray centerline 3D point
    def centerline_point(d_deg: float) -> np.ndarray:
        z = float(np.interp(d_deg, path_dist_deg, path_depth_km))
        u_surf = rotate_about_axis(s_hat, axis, math.radians(d_deg))
        return (EARTH_R_KM - z) * u_surf

    # Main loops: slices (deterministic) x Fibonacci points (deterministic)
    for si, d_deg in enumerate(slice_dists):
        # local tangent via neighboring slice points
        d_prev = slice_dists[max(si - 1, 0)]
        d_next = slice_dists[min(si + 1, slice_dists.size - 1)]
        if d_next == d_prev:
            continue

        p0 = centerline_point(float(d_deg))
        p_prev = centerline_point(float(d_prev))
        p_next = centerline_point(float(d_next))
        t_hat = unit(p_next - p_prev)

        # choose a helper not parallel to tangent
        helper = axis if abs(float(np.dot(axis, t_hat))) < 0.9 else unit(p0)
        n1 = unit(np.cross(t_hat, helper))
        if float(np.linalg.norm(n1)) < 1e-12:
            # fallback helper
            helper2 = unit(p0)
            n1 = unit(np.cross(t_hat, helper2))
        n2 = unit(np.cross(t_hat, n1))

        for k in range(N_disk):
            rho = float(rho_km[k])
            phi = float(phi_rad[k])

            offset = rho * (math.cos(phi) * n1 + math.sin(phi) * n2)
            p = p0 + offset
            r0 = float(np.linalg.norm(p0))
            u0 = p0 / r0
            p0_lat, p0_lon = unit_to_latlon(u0)
            p0_depth = EARTH_R_KM - r0

            print(f"[slice] d_deg={d_deg:.2f}  p0(lat,lon,depth)=({p0_lat:+.3f},{p0_lon:+.3f},{p0_depth:.1f} km)")

            rmag = float(np.linalg.norm(p))
            scat_depth = EARTH_R_KM - rmag
            CMB_KM = 2891.0
            if scat_depth > CMB_KM:
                continue
            if scat_depth < 0.0:
                # outside Earth (above surface); skip
                continue

            u = p / rmag
            scat_lat, scat_lon = unit_to_latlon(u)
            lat_ray, lon_ray = unit_to_latlon(p0)


            print('scat_lat, scat_lon,scat_depth =', scat_lat, scat_lon,scat_depth)
            # print('lat_ray, lon_ray =', lat_ray, lon_ray)

            # sys.exit()

            # Epicentral distances for legs (use surface lat/lon; depth handled in TauP)
            # d1 = kilometers2degrees(gps2dist_azimuth(src_lat, src_lon, scat_lat, scat_lon)[0] / 1000.0)
            # d2 = kilometers2degrees(gps2dist_azimuth(scat_lat, scat_lon, rcv_lat, rcv_lon)[0] / 1000.0)
            d1 = central_angle_deg(src_lat, src_lon, scat_lat, scat_lon)
            d2 = central_angle_deg(scat_lat, scat_lon, rcv_lat, rcv_lon)

            # Enforce endcaps in terms of epicentral distance too (optional but usually desired)
            if d1 < min_endcap_deg or d2 < min_endcap_deg:
                print(f"d1={d1:.2f}")
                continue

            # Two-leg TauP
            try:
                t1, _p1 = first_P_time_p_sdeg(model, d1, src_depth_km, scat_depth)
                t2, p2 = first_P_time_p_sdeg(model, d2, scat_depth, 0.0)
                print('P found \n')
            except ValueError:
                print(f"P not found; d1={d1:.2f}, src_depth_km={src_depth_km:.2f}, scat_depth={scat_depth:.2f} ")
                break
                # continue

            ttot = t1 + t2
            dt_obs = ttot - tP_dir
            dp_obs = p2 - pP_dir

            # Incoming direction at receiver: azimuth receiver->scatterer
            _, baz_in, _ = gps2dist_azimuth(rcv_lat, rcv_lon, scat_lat, scat_lon)
            baz_in = float(baz_in)
            dbaz_obs = wrap180(baz_in - baz_dir_in)

            ok_t, dt_bin = in_range_and_bin(dt_obs, dt_lo, dt_hi, dt_step_s)
            if not ok_t:
                continue
            ok_p, dp_bin = in_range_and_bin(dp_obs, dp_lo, dp_hi, dp_step_sdeg)
            if not ok_p:
                continue
            ok_b, db_bin = in_range_and_bin(dbaz_obs, db_lo, db_hi, dbaz_step_deg)
            if not ok_b:
                continue

            cands.append(Candidate(
                scat_lat=scat_lat,
                scat_lon=scat_lon,
                scat_depth_km=float(scat_depth),
                delta1_deg=float(d1),
                delta2_deg=float(d2),
                t_total_s=float(ttot),
                p2_sdeg=float(p2),
                baz_in_deg=float(baz_in),
                dt_obs_s=float(dt_obs),
                dp_obs_sdeg=float(dp_obs),
                dbaz_obs_deg=float(dbaz_obs),
                dt_bin_s=float(dt_bin),
                dp_bin_sdeg=float(dp_bin),
                dbaz_bin_deg=float(db_bin),
                slice_dist_deg=float(d_deg),
                rho_km=float(rho),
                phi_deg=float(phi_deg[k]),))

            if max_candidates is not None and len(cands) >= max_candidates:
                break
        # sys.exit()
        if max_candidates is not None and len(cands) >= max_candidates:
            break

    return {
        "direct": {
            "delta_sr_deg": float(delta_sr_deg),
            "tP_s": float(tP_dir),
            "pP_sdeg": float(pP_dir),
            "baz_in_deg": float(baz_dir_in),
        },
        "sampling": {
            "tube_radius_km": float(tube_radius_km),
            "min_endcap_deg": float(min_endcap_deg),
            "ddeg": float(ddeg),
            "N_disk": int(N_disk),
            "n_slices": int(slice_dists.size),
            "n_total_points": int(slice_dists.size * N_disk),
        },
        "ranges": {
            "dt_range_s": dt_range_s, "dt_step_s": float(dt_step_s),
            "dbaz_range_deg": dbaz_range_deg, "dbaz_step_deg": float(dbaz_step_deg),
            "dp_range_sdeg": dp_range_sdeg, "dp_step_sdeg": float(dp_step_sdeg),
        },
        "candidates": cands,
    }

# sys.exit()
if __name__ == "__main__":

    src_lat, src_lon, src_depth_km = 10.0, 20.0, 300.0
    rcv_lat, rcv_lon = 35.0, 70.0
    model = TauPyModel(model="iasp91")

    out = find_scatterers_tube_deterministic(
        src_lat, src_lon, src_depth_km,
        rcv_lat, rcv_lon,
        model_name="iasp91",

        tube_radius_km=100.0,
        min_endcap_deg=10.0,

        # resolution
        ddeg=0.2,
        N_disk=50,

        #  constraint ranges/steps
        dt_range_s=(50.0, 200.0), dt_step_s=1.0,
        dbaz_range_deg=(-5.0, 5.0), dbaz_step_deg=0.5,
        dp_range_sdeg=(0.0, 2.0), dp_step_sdeg=0.25,

        # optional: stop early possibility
        max_candidates=None,
    )

    print("Direct:", out["direct"])
    print("Sampling:", out["sampling"])
    print("Ranges:", out["ranges"])
    print("Matches:", len(out["candidates"]))

    # Print
    for c in out["candidates"][:15]:
        print(
            f"scat(lat,lon,z)=({c.scat_lat:+.3f},{c.scat_lon:+.3f},{c.scat_depth_km:7.1f} km) "
            f"Δ1={c.delta1_deg:6.2f} Δ2={c.delta2_deg:6.2f} "
            f"dt={c.dt_obs_s:7.2f}s(bin {c.dt_bin_s:5.0f}) "
            f"dp={c.dp_obs_sdeg:6.3f}s/deg(bin {c.dp_bin_sdeg:3.1f}) "
            f"dbaz={c.dbaz_obs_deg:+6.2f}°(bin {c.dbaz_bin_deg:+3.0f}) "
            f"sliceΔ={c.slice_dist_deg:6.2f}° rho={c.rho_km:6.1f}km phi={c.phi_deg:6.1f}° "
            f"t_total={c.t_total_s:8.2f}s p2={c.p2_sdeg:7.3f}s/deg baz_in={c.baz_in_deg:7.2f}°"
        )

# sys.exit()
