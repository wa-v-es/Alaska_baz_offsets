# feb 2
#!/usr/bin/env python3
"""
Depth-consistent single-scatterer search using TauP (ObsPy) in 1-D Earth.

Key constraint implemented:
  Scatterer lies on the *downgoing* branch of the *direct P* ray from source to receiver.
  => scatterer depth z_scat and first-leg time t1 come from interpolating the direct-ray path.

Scattered arrival modeled as:
  t_scatt = t1(direct-ray to depth point at Δ1) + t2(P from depth z_scat to receiver over Δ2)

Observed at receiver:
  slowness proxy = ray_param of leg2 (p2, s/deg)
  incoming baz   = azimuth from receiver to scatterer (receiver->scat)

Dependencies:
  pip install obspy geographiclib numpy
"""
import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from geographiclib.geodesic import Geodesic

def wrap180(angle_deg: float) -> float:
    return (angle_deg + 180.0) % 360.0 - 180.0

def gc_distance_deg(lat1, lon1, lat2, lon2) -> float:
    dist_m, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return kilometers2degrees(dist_m / 1000.0)

def azimuth_deg(lat1, lon1, lat2, lon2) -> float:
    _, az, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
    return float(az)

def first_arrival_P_time_p_sdeg(model: TauPyModel, delta_deg: float, source_depth_km: float) -> Tuple[float, float]:
    arrs = model.get_travel_times(
        source_depth_in_km=source_depth_km,
        distance_in_degree=delta_deg,
        phase_list=["P"],
    )
    if not arrs:
        raise ValueError("No P")
    a0 = arrs[0]
    return float(a0.time), float(a0.ray_param_sec_degree)  # seconds, s/deg

def direct_P_downgoing_profile(model: TauPyModel, delta_sr_deg: float, src_depth_km: float) -> Dict[str, np.ndarray]:
    rps = model.get_ray_paths(
        source_depth_in_km=src_depth_km,
        distance_in_degree=delta_sr_deg,
        phase_list=["P"],
    )
    if not rps:
        raise ValueError("No direct P ray path")
    rp = rps[0]

    dist = np.asarray(rp.path["dist"], float)   # deg from source
    depth = np.asarray(rp.path["depth"], float) # km
    time = np.asarray(rp.path["time"], float)   # s

    i_turn = int(np.argmax(depth))
    dist_dn = dist[: i_turn + 1]
    depth_dn = depth[: i_turn + 1]
    time_dn = time[: i_turn + 1]

    # De-dup any repeated distances (keep last)
    _, idx_rev = np.unique(dist_dn[::-1], return_index=True)
    keep = (len(dist_dn) - 1 - idx_rev)[::-1]
    keep = np.sort(keep)
    dist_dn, depth_dn, time_dn = dist_dn[keep], depth_dn[keep], time_dn[keep]

    return {
        "dist_dn": dist_dn,
        "depth_dn": depth_dn,
        "time_dn": time_dn,
        "delta_turn_deg": float(dist_dn[-1]),
        "z_turn_km": float(depth_dn[-1]),
        "t_dir_s": float(rp.time),
        "p_dir_sdeg": float(rp.ray_param_sec_degree),
    }

def point_on_corridor(src_lat, src_lon, f, xtrack_km, az_sr, dist_sr_m) -> Tuple[float, float]:
    geod = Geodesic.WGS84
    g1 = geod.Direct(src_lat, src_lon, az_sr, f * dist_sr_m)
    lat_mid, lon_mid = g1["lat2"], g1["lon2"]
    az_perp = az_sr + (90.0 if xtrack_km >= 0 else -90.0)
    g2 = geod.Direct(lat_mid, lon_mid, az_perp, abs(xtrack_km) * 1000.0)
    return float(g2["lat2"]), float(g2["lon2"])

def in_range_and_bin(x: float, lo: float, hi: float, step: float) -> Tuple[bool, float]:
    """
    Check if x is within [lo, hi] (inclusive) and return nearest bin center on that grid.
    Accept if within half-step of that bin (so "every step" makes sense).
    """
    if x < lo or x > hi:
        return (False, np.nan)
    k = round((x - lo) / step)
    x_bin = lo + k * step
    ok = abs(x - x_bin) <= 0.5 * step + 1e-12
    return (ok, float(x_bin))

@dataclass
class Candidate:
    scat_lat: float
    scat_lon: float
    delta1_deg: float
    delta2_deg: float
    z_scat_km: float
    t_total_s: float
    p2_sdeg: float
    baz_in_deg: float

    dt_obs_s: float
    dp_obs_sdeg: float
    dbaz_obs_deg: float

    dt_bin_s: float
    dp_bin_sdeg: float
    dbaz_bin_deg: float

def find_scatterers_ranges(
    src_lat: float, src_lon: float, src_depth_km: float,
    rcv_lat: float, rcv_lon: float,
    *,
    model_name: str = "iasp91",

    # RANGES & STEPS (your request)
    dt_range_s: Tuple[float, float] = (50.0, 100.0),
    dt_step_s: float = 2.0,

    dbaz_range_deg: Tuple[float, float] = (-5.0, 5.0),
    dbaz_step_deg: float = 1.0,

    dp_range_sdeg: Tuple[float, float] = (0.0, 2.0),
    dp_step_sdeg: float = 0.5,

    # geometry constraints
    min_endcap_deg: float = 10.0,
    xtrack_max_km: float = 1000.0,
    f_min: float = 0.05,
    f_max: float = 0.95,

    # sampling
    n_random: int = 300_000,
) -> Dict[str, object]:

    model = TauPyModel(model=model_name)

    dist_sr_m, az_sr, baz_rs = gps2dist_azimuth(src_lat, src_lon, rcv_lat, rcv_lon)
    delta_sr_deg = kilometers2degrees(dist_sr_m / 1000.0)

    # direct P
    prof = direct_P_downgoing_profile(model, delta_sr_deg, src_depth_km)
    t_dir = prof["t_dir_s"]
    p_dir = prof["p_dir_sdeg"]
    baz_dir_in = float(baz_rs)  # receiver->source, 0..360

    dist_dn = prof["dist_dn"]
    depth_dn = prof["depth_dn"]
    time_dn = prof["time_dn"]
    delta_turn = prof["delta_turn_deg"]

    rng = np.random.default_rng(0)
    f_samps = rng.uniform(f_min, f_max, size=n_random)
    x_samps = rng.uniform(-xtrack_max_km, xtrack_max_km, size=n_random)

    cands: List[Candidate] = []

    dt_lo, dt_hi = dt_range_s
    db_lo, db_hi = dbaz_range_deg
    dp_lo, dp_hi = dp_range_sdeg

    for f, xkm in zip(f_samps, x_samps):
        scat_lat, scat_lon = point_on_corridor(src_lat, src_lon, float(f), float(xkm), float(az_sr), float(dist_sr_m))

        d1 = gc_distance_deg(src_lat, src_lon, scat_lat, scat_lon)
        d2 = gc_distance_deg(scat_lat, scat_lon, rcv_lat, rcv_lon)

        if d1 < min_endcap_deg or d2 < min_endcap_deg:
            continue
        if d1 > delta_turn:  # must be on downgoing branch of direct P
            continue

        z_scat = float(np.interp(d1, dist_dn, depth_dn))
        t1 = float(np.interp(d1, dist_dn, time_dn))

        try:
            t2, p2 = first_arrival_P_time_p_sdeg(model, d2, z_scat)
        except ValueError:
            continue

        ttot = t1 + t2
        baz_in = azimuth_deg(rcv_lat, rcv_lon, scat_lat, scat_lon)  # receiver->scatterer

        # Observed offsets relative to DIRECT P at receiver
        dt_obs = ttot - t_dir
        dp_obs = p2 - p_dir
        dbaz_obs = wrap180(baz_in - baz_dir_in)  # in [-180,180)

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
            scat_lat=scat_lat, scat_lon=scat_lon,
            delta1_deg=float(d1), delta2_deg=float(d2),
            z_scat_km=float(z_scat),
            t_total_s=float(ttot),
            p2_sdeg=float(p2),
            baz_in_deg=float(baz_in),
            dt_obs_s=float(dt_obs),
            dp_obs_sdeg=float(dp_obs),
            dbaz_obs_deg=float(dbaz_obs),
            dt_bin_s=float(dt_bin),
            dp_bin_sdeg=float(dp_bin),
            dbaz_bin_deg=float(db_bin),
        ))

    return {
        "direct": {
            "delta_sr_deg": float(delta_sr_deg),
            "tP_s": float(t_dir),
            "pP_sdeg": float(p_dir),
            "baz_in_deg": float(baz_dir_in),
            "delta_turn_deg": float(delta_turn),
            "z_turn_km": float(prof["z_turn_km"]),
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
    # Example
    SRC = (10.0, 20.0, 300.0)    # lat, lon, depth_km
    RCV = (35.0, 90.0)         # lat, lon
    out = find_scatterers_ranges(
        SRC[0], SRC[1], SRC[2],
        RCV[0], RCV[1],
        dt_range_s=(50.0, 150.0), dt_step_s=2.0,
        dbaz_range_deg=(-5.0, 5.0), dbaz_step_deg=1.0,
        dp_range_sdeg=(0.0, 2.0), dp_step_sdeg=0.5,
        min_endcap_deg=5.0,
        xtrack_max_km=1000.0,
        n_random=500_000,)

    # out = find_depth_scatterers(
    #     SRC[0], SRC[1], SRC[2],
    #     RCV[0], RCV[1],
    #     model_name="iasp91",
    #     dt_target_s=50.0,
    #     dp_target_s_per_deg=1.0,
    #     dbaz_target_deg=2.0,
    #     tol_t_s=2.0,
    #     tol_p_s_per_deg=0.25,
    #     tol_baz_deg=0.5,
    #     min_endcap_deg=3.0,
    #     xtrack_max_km=1000.0,
    #     mode="random",
    #     n_random=100_000,
    # )

    print("Direct:", out["direct"])
    # print("Target:", out["target"])
    print("Matches:", len(out["candidates"]))

    for c in out["candidates"][:10]:
        print(
            f"scat=({c.scat_lat:+.3f},{c.scat_lon:+.3f}) "
            f"Δ1={c.delta1_deg:.2f} Δ2={c.delta2_deg:.2f} z={c.z_scat_km:.1f}km "
            f"t={c.t_total_s:.2f}s p2={c.p2_s_per_deg:.3f}s/deg baz={c.baz_in_deg:.2f} "
            f"misfit(dt,dp,dbaz)=({c.misfit_t_s:+.2f},{c.misfit_p:+.3f},{c.misfit_baz_deg:+.2f})"
        )
