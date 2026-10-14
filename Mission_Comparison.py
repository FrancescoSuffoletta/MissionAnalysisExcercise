# -*- coding: utf-8 -*-
import sys
sys.stdout.reconfigure(encoding="utf-8")

"""
Mission_Comparison.py
=====================
Side-by-side parametric comparison of two LEO/VLEO mission archetypes:
  1. Direct-to-Device Telecommunications Constellation  (Telecom_Constellation.py)
  2. Thermal IR Wildfire Detection Constellation        (FireDetect.py)

Both missions study VLEO (200 km) vs LEO (500 km) design points over Spain/Galicia.
This script reproduces the key figures of merit from each mission and explicitly
highlights where they agree and where they diverge — and explains why.

Structure
---------
Section 1 : Shared orbital mechanics (applies to both missions)
Section 2 : Telecom-specific sizing  (Walker Delta, link budget summary)
Section 3 : Fire detection-specific sizing (SSO, pushbroom resolution)
Section 4 : Head-to-head comparison table
"""

import math

# =============================================================================
# SECTION 1 — SHARED ORBITAL MECHANICS
# Both missions use these constants and functions.
# =============================================================================

# Physical constants
R_E    = 6_371.0          # Earth mean radius, km
mu     = 398_600.4418     # Earth gravitational parameter, km³/s²
J2     = 1.08263e-3       # Earth oblateness coefficient
omega_E = 7.2921150e-5    # Earth rotation rate, rad/s  (≈ 360°/day)
omega_sun = 0.9856        # Solar mean motion, deg/day  (used for SSO)

# Mission study altitudes (both missions compare these two)
ALTITUDES = {"VLEO": 200.0, "LEO": 500.0}   # km


def orbital_period(h_km):
    """Keplerian orbital period [min] for circular orbit at altitude h_km."""
    a = R_E + h_km                    # semi-major axis, km
    T_sec = 2 * math.pi * math.sqrt(a**3 / mu)
    return T_sec / 60.0               # convert to minutes


def orbital_velocity(h_km):
    """Circular orbital velocity [km/s]."""
    return math.sqrt(mu / (R_E + h_km))


def footprint_half_angle(h_km, elev_min_deg):
    """
    Earth central half-angle λ [deg] — the angular radius of the coverage
    footprint as seen from Earth's centre, for a satellite at altitude h_km
    and a minimum elevation angle elev_min_deg.

    Geometry: λ = 90° - elev_min - arcsin(R_E / (R_E + h) * cos(elev_min))
    """
    e  = math.radians(elev_min_deg)
    rho = R_E / (R_E + h_km)          # Earth angular radius from satellite
    lam = math.pi/2 - e - math.asin(rho * math.cos(e))
    return math.degrees(lam)


def slant_range(h_km, elev_min_deg):
    """Slant range [km] from satellite to ground user at elevation ε_min."""
    e   = math.radians(elev_min_deg)
    lam = math.radians(footprint_half_angle(h_km, elev_min_deg))
    # Sine rule in the ground-satellite-centre triangle
    return (R_E + h_km) * math.sin(lam) / math.cos(e)


def fspl_dB(range_km, freq_GHz):
    """Free-space path loss [dB].  FSPL = 20*log10(4π R f / c)."""
    c = 3e5  # km/s
    return 20 * math.log10(4 * math.pi * range_km * freq_GHz * 1e9 / (c * 1e3))


def drag_thrust_mN(h_km, mass_kg=10.0, Cd=2.2, Am_ratio=0.01):
    """
    Required thrust [mN] to compensate atmospheric drag at altitude h_km.
    Uses Vallado (2013) piecewise exponential atmosphere.
    Am_ratio = A/m [m²/kg].
    """
    # Piecewise exponential density table (reference altitude, density, scale height)
    atm_table = [
        (0,    1.225,      8.44e3),
        (25e3, 3.899e-2,   6.45e3),
        (50e3, 1.027e-3,   7.71e3),
        (100e3, 5.297e-7,  5.877e3),
        (150e3, 2.076e-9,  26.3e3),
        (200e3, 2.789e-10, 37.1e3),
        (300e3, 1.916e-11, 45.5e3),
        (400e3, 2.803e-12, 63.3e3),
        (500e3, 5.215e-13, 88.7e3),
        (600e3, 1.137e-13, 124.64e3),
    ]
    h_m = h_km * 1e3
    # Find bracket
    rho0, href, H = atm_table[5][1], atm_table[5][0], atm_table[5][2]  # default 200 km
    for i in range(len(atm_table) - 1):
        if atm_table[i][0] <= h_m < atm_table[i+1][0]:
            href, rho0, H = atm_table[i][0], atm_table[i][1], atm_table[i][2]
            break
    rho = rho0 * math.exp(-(h_m - href) / H)

    v = orbital_velocity(h_km) * 1e3   # m/s
    F_N = 0.5 * rho * Cd * Am_ratio * mass_kg * v**2  # Drag force, N
    return F_N * 1e3                                    # convert to mN


# --- Print shared parameters ---
print("=" * 70)
print("SECTION 1 — SHARED ORBITAL MECHANICS")
print("=" * 70)
print(f"{'Parameter':<40} {'VLEO (200 km)':>14} {'LEO (500 km)':>14}")
print("-" * 70)
for label, h in ALTITUDES.items():
    pass  # we print together below

h_vleo, h_leo = ALTITUDES["VLEO"], ALTITUDES["LEO"]

rows = [
    ("Altitude [km]",           h_vleo,                    h_leo),
    ("Orbital period [min]",    orbital_period(h_vleo),    orbital_period(h_leo)),
    ("Orbital velocity [km/s]", orbital_velocity(h_vleo),  orbital_velocity(h_leo)),
    ("Drag thrust @ 10 kg [mN]",drag_thrust_mN(h_vleo),    drag_thrust_mN(h_leo)),
]
for name, v, l in rows:
    print(f"  {name:<38} {v:>14.2f} {l:>14.2f}")

print()


# =============================================================================
# SECTION 2 — TELECOM MISSION (Direct-to-Device, Walker Delta)
# =============================================================================
# Key parameters
ELEV_MIN_TELECOM = 20.0   # deg — minimum elevation for link closure (suburban)
FREQ_GHz         = 2.0    # S-band, compatible with unmodified smartphones

print("=" * 70)
print("SECTION 2 — TELECOM MISSION (Direct-to-Device, Walker Delta)")
print("=" * 70)
print(f"  Min. elevation angle: {ELEV_MIN_TELECOM}°   Frequency: {FREQ_GHz} GHz (S-band)")
print()

# --- Footprint and link budget ---
print(f"  {'Parameter':<42} {'VLEO 200 km':>12} {'LEO 500 km':>12}")
print("  " + "-" * 68)

for label, h in ALTITUDES.items():
    pass  # printed together below

lam_v   = footprint_half_angle(h_vleo, ELEV_MIN_TELECOM)
lam_l   = footprint_half_angle(h_leo,  ELEV_MIN_TELECOM)
Rs_v    = slant_range(h_vleo, ELEV_MIN_TELECOM)
Rs_l    = slant_range(h_leo,  ELEV_MIN_TELECOM)
fspl_v  = fspl_dB(Rs_v, FREQ_GHz)
fspl_l  = fspl_dB(Rs_l, FREQ_GHz)
lat_v   = Rs_v / 3e5 * 1e3   # one-way latency, ms
lat_l   = Rs_l / 3e5 * 1e3

telecom_rows = [
    ("Footprint half-angle λ [deg]",       lam_v,   lam_l),
    ("Footprint diameter [km]",            2*lam_v*111.2, 2*lam_l*111.2),
    ("Slant range @ ε_min [km]",           Rs_v,    Rs_l),
    ("Free-space path loss [dB]",          fspl_v,  fspl_l),
    ("FSPL advantage of VLEO [dB]",        fspl_l - fspl_v,  0.0),
    ("One-way latency [ms]",               lat_v,   lat_l),
]
for name, v, l in telecom_rows:
    if "advantage" in name.lower():
        print(f"  {name:<42} {v:>+12.1f} {l:>+12.1f}")
    else:
        print(f"  {name:<42} {v:>12.1f} {l:>12.1f}")

# --- Walker Delta constellation sizing for continuous coverage of Spain ---
# Spain: 36–44°N.  Coverage latitude = 44°N (worst case).
phi_max = 44.0   # deg

print()
print("  Walker Delta sizing — continuous coverage of Spain (36-44°N):")
print(f"  {'Tier':<30} {'N (VLEO)':>10} {'N (LEO)':>10}  Description")
print("  " + "-" * 70)

# Walker: P = ceil(360 cos(phi_max) / (2 lambda))
#         S = ceil(180 / lambda)
#         N = P * S
for tier_name, elev_req in [("Tier 1 IoT/SOS (ε≥5°)",   5),
                             ("Tier 2 Messaging (ε≥15°)", 15),
                             ("Tier 3 Voice (ε≥20°)",     20)]:
    lam_v_t = footprint_half_angle(h_vleo, elev_req)
    lam_l_t = footprint_half_angle(h_leo,  elev_req)
    P_v = math.ceil(360 * math.cos(math.radians(phi_max)) / (2 * lam_v_t))
    S_v = math.ceil(180 / lam_v_t)
    P_l = math.ceil(360 * math.cos(math.radians(phi_max)) / (2 * lam_l_t))
    S_l = math.ceil(180 / lam_l_t)
    N_v = P_v * S_v
    N_l = P_l * S_l
    print(f"  {tier_name:<30} {N_v:>10,} {N_l:>10,}  P={P_v}×S={S_v} / P={P_l}×S={S_l}")

# The reference continuous-voice case (ε=20°) is the headline number
lam_v20 = footprint_half_angle(h_vleo, 20)
lam_l20 = footprint_half_angle(h_leo,  20)
P_v20 = math.ceil(360 * math.cos(math.radians(phi_max)) / (2 * lam_v20))
S_v20 = math.ceil(180 / lam_v20)
P_l20 = math.ceil(360 * math.cos(math.radians(phi_max)) / (2 * lam_l20))
S_l20 = math.ceil(180 / lam_l20)
N_v20 = P_v20 * S_v20
N_l20 = P_l20 * S_l20
print()
print(f"  → Headline continuous-voice constellation:")
print(f"    VLEO: {N_v20:,} sats   LEO: {N_l20:,} sats   ratio: {N_v20/N_l20:.1f}×")
print()


# =============================================================================
# SECTION 3 — FIRE DETECTION MISSION (SSO, Pushbroom Thermal IR)
# =============================================================================
# Key parameters
PIXEL_SIZE_um  = 15.0    # µm — detector pixel pitch
FOCAL_LENGTH_mm = 150.0  # mm — camera focal length
FOV_deg        = 25.0    # deg full cross-track field of view (half = 12.5°)
TARGET_LAT     = 42.8    # deg N — Galicia
REVISIT_REQ_h  = 2.0     # hours — required mean revisit time

print("=" * 70)
print("SECTION 3 — FIRE DETECTION MISSION (SSO, Pushbroom Thermal IR)")
print("=" * 70)
print(f"  Focal length: {FOCAL_LENGTH_mm} mm   Pixel: {PIXEL_SIZE_um} µm   FOV: {FOV_deg}°")
print(f"  Target: Galicia {TARGET_LAT}°N   Revisit target: {REVISIT_REQ_h} h")
print()

# --- SSO inclination (solving J2 nodal regression = ω_sun) ---
def sso_inclination(h_km):
    """
    SSO inclination [deg] such that nodal regression dΩ/dt = ω_sun.
    dΩ/dt = -(3/2) n J2 (R_E/a)^2 cos(i)  [rad/s]
    Solve for i.
    """
    a  = (R_E + h_km) * 1e3               # m
    n  = math.sqrt(mu * 1e9 / a**3)        # mean motion, rad/s  (mu in m³/s²)
    # Note: mu must be in m³/s² here
    mu_m = 3.986004418e14                  # m³/s²
    a_m  = (R_E + h_km) * 1e3
    n_m  = math.sqrt(mu_m / a_m**3)
    target = math.radians(omega_sun / 365.25 * (365.25 / 365.25))  # rad/day → rad/s
    # ω_sun = 0.9856 deg/day
    omega_sun_rads = math.radians(omega_sun) / 86400.0
    cos_i = -omega_sun_rads / ((3/2) * n_m * J2 * (6_371_000 / a_m)**2)
    return math.degrees(math.acos(cos_i))

print(f"  {'Parameter':<42} {'VLEO 200 km':>12} {'LEO 500 km':>12}")
print("  " + "-" * 68)

# Ground resolution (GSD at nadir)
def gsd_nadir(h_km, pixel_um, focal_mm):
    """Ground sampling distance at nadir [m]."""
    return h_km * 1e3 * (pixel_um * 1e-6) / (focal_mm * 1e-3)

# Swath width at nadir (full FOV)
def swath_nadir(h_km, fov_deg):
    """Ground swath width [km] for full-FOV pushbroom at nadir."""
    return 2 * h_km * math.tan(math.radians(fov_deg / 2))

# Single-satellite revisit (all-day, both ascending and descending passes)
def single_sat_revisit_h(h_km, swath_km, lat_deg):
    """
    Mean revisit time [h] for a single SSO satellite over a point at lat_deg.
    Uses both ascending and descending passes (thermal IR works day and night).
    Coverage fraction f = W / (π R_E cos(φ))
    T_rev = T_orbit / f
    """
    W  = swath_km
    f  = W / (math.pi * R_E * math.cos(math.radians(lat_deg)))
    T  = orbital_period(h_km) / 60.0     # hours
    return T / f

i_v  = sso_inclination(h_vleo)
i_l  = sso_inclination(h_leo)
gsd_v = gsd_nadir(h_vleo, PIXEL_SIZE_um, FOCAL_LENGTH_mm)
gsd_l = gsd_nadir(h_leo,  PIXEL_SIZE_um, FOCAL_LENGTH_mm)
sw_v  = swath_nadir(h_vleo, FOV_deg)
sw_l  = swath_nadir(h_leo,  FOV_deg)
rev_v = single_sat_revisit_h(h_vleo, sw_v, TARGET_LAT)
rev_l = single_sat_revisit_h(h_leo,  sw_l, TARGET_LAT)

fire_rows = [
    ("SSO inclination [deg]",             i_v,    i_l),
    ("GSD at nadir [m]",                  gsd_v,  gsd_l),
    ("Swath width at nadir [km]",         sw_v,   sw_l),
    ("Single-sat revisit (all-day) [h]",  rev_v,  rev_l),
]
for name, v, l in fire_rows:
    print(f"  {name:<42} {v:>12.1f} {l:>12.1f}")

# --- Constellation sizing for T_r = 2 h ---
# N = ceil(T_rev_single / T_r),  one satellite per plane (S=1, optimal phasing)
N_fire_v = math.ceil(rev_v / REVISIT_REQ_h)
N_fire_l = math.ceil(rev_l / REVISIT_REQ_h)

print()
print(f"  Walker SSO constellation for ≤{REVISIT_REQ_h} h revisit (S=1 sat/plane):")
print(f"    VLEO: {N_fire_v} planes × 1 = {N_fire_v} sats")
print(f"    LEO : {N_fire_l} planes × 1 = {N_fire_l} sats")
print(f"    Ratio: {N_fire_v/N_fire_l:.1f}×")
print()

# --- Off-nadir agility bonus (VLEO only) ---
# With ±α_max roll, effective swath grows; revisit shrinks accordingly.
alpha_max_deg = 52.0   # deg maximum off-nadir roll
gsd_limit_m   = 50.0   # m  maximum acceptable GSD (LEO-parity)

# Maximum allowed nadir angle given GSD limit
# GSD(α) ≈ GSD_nadir / cos(α)   (first-order approximation)
alpha_gsd_limit = math.degrees(math.acos(gsd_v / gsd_limit_m))

# Effective swath with agile pointing: W_eff = 2 * h * tan(alpha_eff)
# Use the GSD-limited angle, capped by mechanical limit
alpha_eff = min(alpha_max_deg, alpha_gsd_limit)
W_eff_km  = 2 * h_vleo * math.tan(math.radians(alpha_eff))
rev_agile = single_sat_revisit_h(h_vleo, W_eff_km, TARGET_LAT)
N_agile   = math.ceil(rev_agile / REVISIT_REQ_h)

print(f"  VLEO off-nadir agility (±{alpha_max_deg}°, GSD ≤ {gsd_limit_m} m):")
print(f"    Effective α [deg]  : {alpha_eff:.1f}°")
print(f"    Effective swath    : {W_eff_km:.1f} km  (vs {sw_v:.1f} km nadir-only)")
print(f"    Single-sat revisit : {rev_agile:.1f} h  (vs {rev_v:.1f} h nadir-only)")
print(f"    Constellation size : {N_agile} sats  (vs {N_fire_v} nadir-only)")
print()


# =============================================================================
# SECTION 4 — HEAD-TO-HEAD COMPARISON TABLE
# =============================================================================
print("=" * 70)
print("SECTION 4 — HEAD-TO-HEAD COMPARISON")
print("=" * 70)

comparison = [
    # (Feature, Telecom value/note, Fire Detection value/note)
    ("Orbit type",
     "Walker Delta (prograde ~49°)",
     "Sun-Synchronous (retrograde ~96-97°)"),

    ("Reason for orbit choice",
     "Cover Spain latitude band efficiently",
     "Fixed local solar time → uniform thermal IR conditions"),

    ("VLEO payload",
     "Large RF phased array (64 m² class)",
     "Compact pushbroom IR camera (150 mm)"),

    ("LEO payload",
     "Smaller RF array, still large",
     "Same camera, wider swath, lower resolution"),

    (f"VLEO constellation N",
     f"{N_v20:,} sats (continuous voice)",
     f"{N_fire_v} sats (2 h revisit, nadir-only)"),

    (f"LEO constellation N",
     f"{N_l20:,} sats (continuous voice)",
     f"{N_fire_l} sats (2 h revisit)"),

    ("VLEO/LEO ratio",
     f"{N_v20/N_l20:.1f}× more sats at VLEO",
     f"{N_fire_v/N_fire_l:.1f}× more sats at VLEO"),

    ("Coverage requirement",
     "Simultaneous (hard — no gap anywhere)",
     "Statistical (mean revisit ≤ 2 h)"),

    ("Sats per plane (S)",
     "20–42 (fill orbit arc continuously)",
     "1  (one pass per plane suffices)"),

    (f"VLEO path-loss advantage",
     f"+{fspl_l-fspl_v:.1f} dB  BUT wasted in urban blockage",
     "N/A (passive IR, no link budget)"),

    ("VLEO resolution advantage",
     "N/A (resolution not a driver)",
     f"{gsd_l/gsd_v:.1f}× better GSD ({gsd_v:.0f} m vs {gsd_l:.0f} m)"),

    ("Off-nadir agility benefit",
     "Minor (blockage limits usable ε anyway)",
     f"Major: {N_fire_v}→{N_agile} sats ({N_fire_v//max(N_agile,1):.0f}× reduction)"),

    ("ABEP requirement",
     f"Yes — {drag_thrust_mN(h_vleo):.2f} mN per sat × {N_v20:,} units",
     f"Yes — {drag_thrust_mN(h_vleo):.2f} mN per sat × {N_fire_v} units"),

    ("ABEP supply-chain risk",
     "CRITICAL (1000+ thrusters, unproven at fleet scale)",
     "Manageable (100-unit order, viable supply chain)"),

    ("Preferred altitude",
     "LEO (risk and cost favour LEO strongly)",
     "VLEO viable (resolution benefit justifies penalty)"),
]

for feature, telecom_val, fire_val in comparison:
    print(f"\n  {feature}")
    print(f"    Telecom   : {telecom_val}")
    print(f"    Fire Det. : {fire_val}")

print()
print("=" * 70)
print("KEY INSIGHT")
print("=" * 70)
print("""
  Both missions use the same orbital-mechanics toolkit (Walker geometry,
  J2, footprint/swath sizing, ABEP drag compensation) but reach opposite
  conclusions on VLEO viability:

  TELECOM  → LEO preferred
    • Continuous zonal coverage forces N ∝ 1/λ² → 4.5× satellite penalty
    • VLEO path-loss advantage neutralised by urban building blockage
    • 1,260 ABEP thrusters is a logistical and financial blocker

  FIRE DETECTION  → VLEO competitive
    • Point/regional coverage: N ∝ 1/T_rev → only 2.3× penalty
    • 2.5× GSD improvement (20 m vs 50 m) directly improves hotspot detection
    • Off-nadir agility can compress N further (e.g. 122 → {N_agile} sats)
    • 122 ABEP units is a manageable programme risk

  The deciding factor is MISSION ARCHITECTURE:
    Zonal continuous coverage  → LEO wins
    Regional periodic coverage → VLEO is worth the trade
""".format(N_agile=N_agile))
