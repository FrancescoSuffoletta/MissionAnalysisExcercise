"""
================================================================================
  SATELLITE MISSION TRADE STUDY — FIRE DETECTION OVER GALICIA
  Comparing VLEO (200 km) vs LEO (500 km) Sun-Synchronous Orbits
================================================================================

INSTRUMENT
----------
  Pixel size    p  = 15 µm
  Focal length  f  = 150 mm
  Full FOV         = 25°  (cross-track, pushbroom)

  Derived:
    IFOV = p / f = 100 µrad  (angular size of one pixel)

--------------------------------------------------------------------------------
STEP 1 — GROUND SAMPLING DISTANCE (GSD)
--------------------------------------------------------------------------------
  GSD is the ground size of one pixel at nadir:

      GSD = IFOV × H = (p / f) × H

  Linear with altitude. At H = 200 km → GSD = 20 m ; at H = 500 km → GSD = 50 m.

--------------------------------------------------------------------------------
STEP 2 — SWATH WIDTH (spherical-Earth geometry)
--------------------------------------------------------------------------------
  For a sensor with half-FOV θ = 12.5° (nadir angle at swath edge), the Earth
  central half-angle ρ is derived from the triangle
  (Earth centre, satellite, ground point) using the sine rule:

      sin(θ + ρ) = (R_E + H) / R_E × sin(θ)
      ρ = arcsin[(R_E + H)/R_E × sin(θ)] − θ

  Full swath:  W = 2 × ρ × R_E

  At H = 200 km → W ≈  88.7 km ; at H = 500 km → W ≈ 222.1 km.

  Note: Galicia spans ~250 km east–west, so VLEO needs 3 passes and LEO needs
  2 passes to cover the full region in one day.

--------------------------------------------------------------------------------
STEP 3 — SUN-SYNCHRONOUS ORBIT (SSO) INCLINATION
--------------------------------------------------------------------------------
  For SSO the J2-driven nodal regression rate must equal Earth's mean angular
  velocity around the Sun (ω_sun ≈ 0.9856°/day):

      dΩ/dt = −(3/2) × n × J2 × (R_E/a)² × cos(i) = +ω_sun

  Solving for the inclination:

      cos(i) = −(2/3) × ω_sun / [n × J2 × (R_E/a)²]

  where  n = √(µ/a³)  is the mean motion and  a = R_E + H.

  VLEO 200 km → i ≈ 96.32° ;  LEO 500 km → i ≈ 97.39°  (both slightly retrograde)

--------------------------------------------------------------------------------
STEP 4 — ORBITAL PERIOD
--------------------------------------------------------------------------------
  Circular orbit:  T = 2π √(a³ / µ)

  VLEO 200 km → T ≈ 88.35 min (≈ 16.30 orbits/day)
  LEO  500 km → T ≈ 94.47 min (≈ 15.24 orbits/day)

--------------------------------------------------------------------------------
STEP 5 — GROUND TRACK SHIFT PER ORBIT (J2-corrected)
--------------------------------------------------------------------------------
  Each orbit the ascending node moves eastward at ω_sun (SSO definition), while
  Earth rotates eastward at ω_Earth.  The net westward ground track shift is:

      Δλ = (ω_Earth − ω_sun) × T

  Using ω_Earth (inertial sidereal rate) and ω_sun accounts for the ~0.27%
  correction that a naïve ω_Earth × T formula ignores.

  VLEO → Δλ ≈ 22.09°/orbit ;  LEO → Δλ ≈ 23.62°/orbit

--------------------------------------------------------------------------------
STEP 6 — SINGLE-SATELLITE MEAN REVISIT OVER GALICIA
--------------------------------------------------------------------------------
  A) GEOMETRY OF CROSSINGS
     Each orbit the satellite crosses latitude φ = 42.8°N exactly twice:
       - Ascending  crossing: argument of latitude u_asc ≈ 43.3°
                              longitude λ_asc ≈ Ω − 6.6°
       - Descending crossing: argument of latitude u_desc = 180° − u_asc ≈ 136.7°
                              longitude λ_desc ≈ Ω + 186.6°
     The two longitudes are ~193° apart, so for any given orbit (given RAAN Ω)
     only ONE crossing can be near Galicia (λ_T ≈ −8°W).

  B) TWO SEPARATE RAAN WINDOWS
     Because λ_asc and λ_desc are ~193° apart, coverage of Galicia requires
     the RAAN to be in one of two distinct windows:
       - Ascending  window: centred at Ω_asc ≈ −1.4°,  half-width ΔΩ/2
       - Descending window: centred at Ω_desc ≈ 165.4°, half-width ΔΩ/2
                            (~167° away from the ascending window)
     These are different orbits (different RAAN values), not two events in the
     same orbit.

     Angular half-width of each window:
         ΔΩ/2 = (W/2) / (R_E × cos φ)   [rad]

  C) COVERAGE FRACTION
     The fraction of the full RAAN circle (0 to 2π) that gives coverage:
       all-day  (thermal IR — THIS MISSION):
           f = 2 × ΔΩ / (2π) = W / (π × R_E × cos φ)
           Both ascending and descending windows are valid since fires are
           detectable at night with thermal cameras.
       day-only (optical / VNIR, for reference):
           f = ΔΩ / (2π) = W / (2π × R_E × cos φ)
           Only the ascending window is usable; LTDAN ≈ 13:30 means the
           descending pass occurs at ~01:30 local time (night).

  D) MEAN REVISIT TIME
     The RAAN shifts by Δλ ≈ 22.09°/orbit (westward, J2-corrected).
     Treating successive RAAN values as uniformly distributed over 0–2π
     (valid for the probabilistic / random-phase model), the mean number of
     orbits between coverage events is 1/f, giving:

         T_rev = T_orbit / f = T_orbit × (π × R_E × cos φ) / W   [all-day]

     Results:
       VLEO all-day  → T_rev ≈  244 h  (≈ 10 days)
       LEO  all-day  → T_rev ≈  104 h  (≈  4.3 days)

  E) COMMON ERROR TO AVOID
     Dividing the coverage window ΔΩ by the per-orbit shift Δλ instead of
     the full RAAN circle 2π overestimates coverage fraction by a factor of
     2π/Δλ ≈ 8–9, making the revisit appear ~8–9× too short.

--------------------------------------------------------------------------------
STEP 7 — CONSTELLATION SIZING AND WALKER PHASING
--------------------------------------------------------------------------------
  A) MINIMUM SATELLITE COUNT (mean revisit)
     For N_total satellites with random / optimal phasing, mean events scale
     linearly:
         T_rev_mean = T_rev_single / N_total
         N_total    = ceil(T_rev_single / T_r)

  B) DISTRIBUTION INTO PLANES — MEAN vs MAXIMUM REVISIT
     Satellites are placed in N_planes planes, each with N_s = N_total/N_planes
     satellites (Walker Star).  The mean revisit is independent of how the
     satellites are distributed (it only depends on N_total).  The MAXIMUM
     (worst-case) revisit is driven by the inter-plane gap:

         T_rev_mean = T_rev_single / N_total         (same for any N_planes)
         T_rev_max  = T_rev_single / N_planes         (decreases with more planes)

     More planes for the same N_total → better worst-case gap.
     To guarantee T_r as a hard maximum: N_planes ≥ T_rev_single / T_r = N_total.

  C) ROLE OF THE WALKER PHASING FACTOR F
     In a Walker T/P/F constellation the true anomaly of satellite j in plane k
     is offset by k × F × 360°/T relative to plane 0.

     - F does NOT change which orbits (RAAN values) cover Galicia.
       That is set entirely by the RAAN spacing between planes.
     - F DOES change WHEN, within each of those orbits, the satellite reaches
       Galicia's latitude.  With F = 0, satellites in adjacent planes share the
       same true anomaly → their coverage events cluster in time → the actual
       maximum gap can be 2–3× the formula value T_rev_single/N_planes.
     - With OPTIMAL F, coverage events from successive planes are spread
       uniformly in time → the actual maximum gap approaches T_rev_single/N_planes.

     Practical consequence for this mission:
       • For N_s = 1 (1 sat per plane, optimal design): optimal F spreads the
         122 (VLEO) or 53 (LEO) passes uniformly → max revisit ≈ mean revisit.
       • For N_s > 1 per plane: F also spreads the intra-plane burst across
         adjacent planes, reducing clustering.
       • The formula T_rev_max = T_rev_single / N_planes already ASSUMES
         optimal F.  With F = 0 the real max revisit will be worse.

     For Walker Star, F = 1 is the most common optimal choice.

  D) RESULTS FOR T_r = 2 h
     VLEO → N_total = 122 sats,  recommended: Walker 122/122/1 (1 sat per plane)
     LEO  → N_total =  53 sats,  recommended: Walker  53/53/1

--------------------------------------------------------------------------------
STEP 8 — PLOTS
--------------------------------------------------------------------------------
  Two figures are produced and saved to GSD_Swath_vs_Altitude.png:
    - GSD [m]       vs altitude (150–1000 km), markers at VLEO and LEO
    - Swath [km]    vs altitude, markers at VLEO and LEO,
                    reference line at Galicia's ~250 km east–west extent

================================================================================
  References / formulas:
    Wertz, J.R. (ed.) — Space Mission Engineering: The New SMAD (2011)
    Vallado, D.A. — Fundamentals of Astrodynamics and Applications (4th ed.)
    Walker, J.G. — Satellite Constellations, JBIS 37, 1984
================================================================================
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ══════════════════════════════════════════════════════════════════════════════
#  CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
MU        = 3.986004418e14   # [m³/s²]  Earth gravitational parameter
R_E       = 6_371_000.0      # [m]      Earth mean radius
J2        = 1.08263e-3       # [-]      Earth oblateness coefficient
W_EARTH   = 7.2921150e-5     # [rad/s]  Earth rotation rate
W_SUN     = 2*math.pi / (365.25 * 86400)  # [rad/s]  Earth mean motion around Sun

# ══════════════════════════════════════════════════════════════════════════════
#  INSTRUMENT PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
pixel_size_m   = 15e-6   # [m]
focal_length_m = 0.150   # [m]
FOV_deg        = 25.0    # [°]  full cross-track FOV

IFOV_rad  = pixel_size_m / focal_length_m        # [rad/pixel]
half_FOV  = math.radians(FOV_deg / 2)            # [rad]

# ══════════════════════════════════════════════════════════════════════════════
#  GALICIA TARGET
# ══════════════════════════════════════════════════════════════════════════════
galicia_lat_deg  = 42.8   # [°N]  representative centroid latitude
galicia_lon_deg  = -8.0   # [°E]  representative centroid longitude
galicia_ew_km    = 250.0  # [km]  approximate east-west extent of Galicia

# ══════════════════════════════════════════════════════════════════════════════
#  MISSION ALTITUDES
# ══════════════════════════════════════════════════════════════════════════════
missions = {
    "VLEO 200 km": 200e3,   # [m]
}

# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def orbital_period(h):
    """Orbital period [s] for circular orbit at altitude h [m]."""
    a = R_E + h
    return 2 * math.pi * math.sqrt(a**3 / MU)

def sso_inclination_deg(h):
    """SSO inclination [°] so J2 nodal regression = Earth's solar mean motion."""
    a  = R_E + h
    n  = math.sqrt(MU / a**3)                        # mean motion [rad/s]
    # dΩ/dt = -(3/2) n J2 (R_E/a)² cos(i)  →  solve for i
    cos_i = -(2/3) * W_SUN / (n * J2 * (R_E/a)**2)
    return math.degrees(math.acos(cos_i))

def gsd(h):
    """GSD [m] at nadir."""
    return IFOV_rad * h

def swath_spherical(h):
    """
    Full swath width [m] accounting for spherical Earth geometry.
      ρ = arcsin((R_E+h)/R_E · sin(θ)) − θ      [Earth central half-angle]
    where θ = half_FOV (nadir angle at swath edge).
    """
    ratio = (R_E + h) / R_E * math.sin(half_FOV)
    rho   = math.asin(ratio) - half_FOV            # Earth central half-angle [rad]
    return 2 * rho * R_E                            # full swath [m]

def ground_track_shift_deg(T):
    """
    Westward ground track shift per orbit [°], corrected for SSO nodal regression.
    For SSO the node regresses eastward at ω_sun, partially compensating Earth's rotation:
        Δλ = (ω_Earth − ω_sun) · T
    """
    return math.degrees((W_EARTH - W_SUN) * T)

def single_sat_revisit(h, lat_deg, daytime_only=False):
    """
    Mean revisit time [s] for a single SSO satellite over a point target.

    For a given orbit (RAAN = Ω) the ascending and descending crossings of
    latitude φ occur at longitudes ~193° apart, so only ONE crossing can be
    near the target per orbit.  They correspond to two separate RAAN windows
    (~167° apart).  Both are counted when daytime_only=False (thermal IR);
    only the ascending window is counted when daytime_only=True (optical,
    since the descending pass is at ~01:30 local time for LTDAN=13:30).

    Coverage fraction per orbit:
      all-day  : f = W / (π   × R_E × cos φ)   [asc + desc windows]
      day-only : f = W / (2π  × R_E × cos φ)   [ascending window only]
    """
    lat    = math.radians(lat_deg)
    W      = swath_spherical(h)
    T      = orbital_period(h)
    denom  = (2 * math.pi if daytime_only else math.pi) * R_E * math.cos(lat)
    return T * denom / W                               # mean revisit [s]

def num_passes_to_cover_galicia(h):
    """
    Minimum number of consecutive (cross-track shifted) passes to cover
    Galicia's full east-west extent.
    """
    W = swath_spherical(h)
    return math.ceil(galicia_ew_km * 1e3 / W)

def constellation_sizing(h, T_r_hr, n_planes_list=(1, 2, 3, 5, 6, 10),
                         daytime_only=False):
    """
    For each number-of-planes option return
      (n_planes, sats_per_plane, total_sats, mean_revisit_h, max_revisit_h).

    Mean revisit = T_single / N_total   (probabilistic, optimal phasing)
    Max  revisit = T_single / N_planes  (worst-case gap between plane coverage windows;
                                         more planes → smaller max gap for same N_total)
    """
    T_single = single_sat_revisit(h, galicia_lat_deg, daytime_only)   # [s]
    N_total   = math.ceil(T_single / (T_r_hr * 3600))
    results   = []
    for np_ in n_planes_list:
        ns        = math.ceil(N_total / np_)
        ntot      = np_ * ns
        t_mean_h  = (T_single / ntot) / 3600
        t_max_h   = (T_single / np_)  / 3600   # inter-plane worst-case gap
        results.append((np_, ns, ntot, t_mean_h, t_max_h))
    return N_total, results

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN OUTPUT
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 70)
print("  INSTRUMENT")
print("=" * 70)
print(f"  Pixel size     : {pixel_size_m*1e6:.1f} µm")
print(f"  Focal length   : {focal_length_m*1e3:.0f} mm")
print(f"  Full FOV       : {FOV_deg}°")
print(f"  IFOV           : {IFOV_rad*1e6:.2f} µrad  ({np.degrees(IFOV_rad)*3600:.4f} arcsec)")

TARGET_REVISITS_HR = [2.0]   # hours — 2 h requirement

for name, H in missions.items():

    T_orb   = orbital_period(H)
    inc     = sso_inclination_deg(H)
    dl      = ground_track_shift_deg(T_orb)
    GSD_m   = gsd(H)
    W_km    = swath_spherical(H) / 1e3
    n_passes= num_passes_to_cover_galicia(H)
    T_rev1  = single_sat_revisit(H, galicia_lat_deg) / 3600   # [hr]
    n_orb   = 86400 / T_orb

    print()
    print("=" * 70)
    print(f"  {name}   —   SSO")
    print("=" * 70)

    # ── Orbital & instrument parameters ───────────────────────────────────
    print(f"\n  [Orbit]")
    print(f"    Semi-major axis     : {(R_E+H)/1e3:.1f} km")
    print(f"    SSO inclination     : {inc:.2f}°")
    print(f"    Orbital period      : {T_orb/60:.2f} min")
    print(f"    Orbits per day      : {n_orb:.2f}")
    print(f"    Ground track shift  : {dl:.2f}° / orbit  (westward)")

    print(f"\n  [Instrument at this altitude]")
    print(f"    GSD (nadir)         : {GSD_m:.1f} m")
    print(f"    Swath width         : {W_km:.1f} km")
    print(f"    Passes to cover Galicia EW: {n_passes}  ({galicia_ew_km:.0f} km / {W_km:.1f} km)")

    T_rev1_day = single_sat_revisit(H, galicia_lat_deg, daytime_only=True) / 3600

    print(f"\n  [Single-satellite revisit over Galicia @ {galicia_lat_deg}N]")
    print(f"    All-day  (asc + desc, thermal IR) : {T_rev1:.2f} h")
    print(f"    Day-only (asc only,   optical)    : {T_rev1_day:.2f} h")
    print(f"    Note: in each orbit only ONE crossing is near Galicia longitude;")
    print(f"          the two values above count separate RAAN windows, not the same orbit.")

    # ── Constellation options ──────────────────────────────────────────────
    for daytime_only, label in [(False, "All-day  (thermal IR — asc + desc)"),
                                (True,  "Day-only (optical   — asc only  )")]:
        print(f"\n  [Constellation sizing — {label} — T_r = {TARGET_REVISITS_HR[0]*60:.0f} min]")
        print(f"  {'Planes':>7}  {'Sat/plane':>10}  {'N_total':>8}  {'Mean [h]':>9}  {'Max [h]':>8}  Walker")
        print(f"  {'-'*7}  {'-'*10}  {'-'*8}  {'-'*9}  {'-'*8}  {'-'*14}")
        for T_r in TARGET_REVISITS_HR:
            _, configs = constellation_sizing(H, T_r, daytime_only=daytime_only)
            for np_, ns, ntot, t_mean, t_max in configs:
                print(f"  {np_:>7}  {ns:>10}  {ntot:>8}  {t_mean:>9.2f}  {t_max:>8.2f}  {ntot}/{np_}/0")
        print()

print("=" * 70)
print("\n  Notes:")
print("  - Walker Star notation: T/P/F  (Total sats / Planes / Phase factor)")
print("  - Mean revisit = T_single / N_total  (probabilistic, random phasing)")
print("  - Max revisit  = T_single / N_planes (worst-case inter-plane gap)")
print("    -> More planes = lower max revisit for the same N_total")
print("    -> To guarantee T_r_max, need N_planes >= T_single / T_r")
print("  - Ground track shift corrected for SSO nodal regression (W_EARTH - W_SUN)")
print("  - SSO planes differ in RAAN -> different local solar times (LTDAN)")
print("  - For fire detection, LTDAN ~13:30 recommended (max solar illumination)")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOTS — GSD and Swath vs Altitude
# ══════════════════════════════════════════════════════════════════════════════
h_km  = np.linspace(150, 1000, 500)          # altitude range [km]
h_m   = h_km * 1e3                            # [m]

gsd_m    = IFOV_rad * h_m                     # GSD [m]  (linear, nadir)
swath_km = np.array([swath_spherical(h) / 1e3 for h in h_m])   # [km]

# Mission markers — VLEO only
markers = {"VLEO\n200 km": 200e3}
marker_colors = {"VLEO\n200 km": "#e74c3c"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(
    f"Instrument: pixel = {pixel_size_m*1e6:.0f} µm  |  f = {focal_length_m*1e3:.0f} mm  |  FOV = {FOV_deg}°",
    fontsize=13, fontweight="bold"
)

# Annotation offsets (dx_km, dy_value)
gsd_offsets   = {"VLEO\n200 km": (40, 8)}
swath_offsets = {"VLEO\n200 km": (40, 80)}

# ── Plot 1: GSD ───────────────────────────────────────────────────────────────
ax1.plot(h_km, gsd_m, color="#27ae60", linewidth=2.5, label="GSD (nadir)")
ax1.fill_between(h_km, gsd_m, alpha=0.08, color="#27ae60")

for label, h in markers.items():
    g  = IFOV_rad * h
    dx, dy = gsd_offsets[label]
    ax1.axvline(h / 1e3, color=marker_colors[label], linestyle="--", linewidth=1.4, alpha=0.8)
    ax1.plot(h / 1e3, g, "o", color=marker_colors[label], markersize=8, zorder=5)
    ax1.annotate(
        f"{label}\nGSD = {g:.0f} m",
        xy=(h / 1e3, g), xytext=(h / 1e3 + dx, g + dy),
        fontsize=9, color=marker_colors[label],
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=marker_colors[label], alpha=0.8),
        arrowprops=dict(arrowstyle="->", color=marker_colors[label], lw=1.2)
    )

ax1.set_xlabel("Altitude [km]", fontsize=11)
ax1.set_ylabel("GSD [m]", fontsize=11)
ax1.set_title("Ground Sampling Distance vs Altitude", fontsize=11)
ax1.set_xlim(150, 1000)
ax1.set_ylim(bottom=0)
ax1.grid(True, linestyle="--", alpha=0.5)
ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax1.legend(fontsize=9, loc="upper left")

# ── Plot 2: Swath ─────────────────────────────────────────────────────────────
ax2.plot(h_km, swath_km, color="#8e44ad", linewidth=2.5, label="Swath (spherical Earth)")
ax2.fill_between(h_km, swath_km, alpha=0.08, color="#8e44ad")

for label, h in markers.items():
    s  = swath_spherical(h) / 1e3
    dx, dy = swath_offsets[label]
    ax2.axvline(h / 1e3, color=marker_colors[label], linestyle="--", linewidth=1.4, alpha=0.8)
    ax2.plot(h / 1e3, s, "o", color=marker_colors[label], markersize=8, zorder=5)
    ax2.annotate(
        f"{label}\nSwath = {s:.1f} km",
        xy=(h / 1e3, s), xytext=(h / 1e3 + dx, s + dy),
        fontsize=9, color=marker_colors[label],
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=marker_colors[label], alpha=0.8),
        arrowprops=dict(arrowstyle="->", color=marker_colors[label], lw=1.2)
    )

ax2.axhline(galicia_ew_km, color="#e67e22", linestyle=":", linewidth=1.6, alpha=0.9,
            label=f"Galicia E-W extent ({galicia_ew_km:.0f} km)")
ax2.legend(fontsize=9, loc="upper left")

ax2.set_xlabel("Altitude [km]", fontsize=11)
ax2.set_ylabel("Swath width [km]", fontsize=11)
ax2.set_title("Swath Width vs Altitude", fontsize=11)
ax2.set_xlim(150, 1000)
ax2.set_ylim(bottom=0)
ax2.grid(True, linestyle="--", alpha=0.5)
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())

plt.tight_layout()
plt.savefig("GSD_Swath_vs_Altitude.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nPlot saved to GSD_Swath_vs_Altitude.png")

# ══════════════════════════════════════════════════════════════════════════════
#  ABEP THRUST FIGURE
#  Required thrust to compensate atmospheric drag in the VLEO altitude range.
#
#  F_drag = 0.5 * CD * A * rho(h) * v(h)^2   [N]
#         = 0.5 * CD * (A/m) * m * rho * v^2
#
#  Left subplot : absolute thrust [mN] for different satellite masses (A/m fixed)
#  Right subplot: drag deceleration [µm/s²] = thrust per unit mass, for different
#                 A/m ratios (independent of total mass — key design parameter)
# ══════════════════════════════════════════════════════════════════════════════

# Vallado 2013 piecewise-exponential atmosphere (subset covering VLEO range)
_ATM_ABEP = np.array([
    (100, 5.297e-07, 5.877), (110, 9.661e-08, 7.263), (120, 2.438e-08, 9.473),
    (130, 8.484e-09, 12.636), (140, 3.845e-09, 16.149), (150, 2.070e-09, 22.523),
    (180, 5.464e-10, 29.740), (200, 2.789e-10, 37.105), (250, 7.248e-11, 45.546),
    (300, 2.418e-11, 53.628), (350, 9.158e-12, 53.298), (400, 3.925e-12, 58.515),
])

def _atm_density_abep(h_m):
    """Atmospheric density [kg/m3] at altitude h_m [m], Vallado 2013."""
    h_km = max(min(h_m / 1000.0, 399.9), 100.0)
    idx  = int(np.searchsorted(_ATM_ABEP[:, 0], h_km, side='right')) - 1
    idx  = max(0, min(idx, len(_ATM_ABEP) - 2))
    h0, rho0, H = _ATM_ABEP[idx]
    return rho0 * math.exp(-(h_km - h0) / H)

CD_REF  = 2.2        # aerodynamic drag coefficient
AM_REF  = 0.01      # [m^2/kg]  reference area-to-mass ratio
H_VLEO  = 200e3     # [m]       design point

# Altitude sweep over VLEO range
h_abep_km  = np.linspace(150, 400, 400)
h_abep_m   = h_abep_km * 1e3
rho_abep   = np.array([_atm_density_abep(h) for h in h_abep_m])
v_abep     = np.sqrt(MU / (R_E + h_abep_m))   # circular orbital velocity [m/s]

# Drag deceleration [m/s^2] = 0.5 * CD * (A/m) * rho * v^2  (per unit A/m)
drag_base  = 0.5 * CD_REF * rho_abep * v_abep**2   # [1/m^2 * N/kg = m^-1 s^-2 * m^2]

# ── Reference values at design point ──────────────────────────────────────
rho_200  = _atm_density_abep(H_VLEO)
v_200    = math.sqrt(MU / (R_E + H_VLEO))
a_200    = 0.5 * CD_REF * AM_REF * rho_200 * v_200**2   # [m/s^2]
F_200_10kg = 10.0 * a_200 * 1e3                          # [mN] for 10 kg satellite
print(f"\n  ABEP design point (200 km, CD={CD_REF}, A/m={AM_REF} m^2/kg):")
print(f"    rho          = {rho_200:.3e} kg/m^3")
print(f"    v_orb        = {v_200:.1f} m/s")
print(f"    drag decel   = {a_200*1e6:.2f} um/s^2  ({a_200:.3e} m/s^2)")
print(f"    Thrust 10 kg = {F_200_10kg:.2f} mN")

# ── Figure ─────────────────────────────────────────────────────────────────
fig3, (axT, axA) = plt.subplots(1, 2, figsize=(13, 5))
fig3.suptitle(
    f"ABEP Required Thrust — VLEO Altitude Range  (CD = {CD_REF})",
    fontsize=13, fontweight="bold"
)

# LEFT: absolute thrust for different satellite masses (A/m = 0.01 m^2/kg)
masses_kg  = [1, 5, 10, 50, 100]
cmap_m     = plt.get_cmap("plasma")
colors_m   = [cmap_m(x) for x in np.linspace(0.1, 0.85, len(masses_kg))]

for m, c in zip(masses_kg, colors_m):
    F_mN = m * drag_base * AM_REF * 1e3    # [mN]
    axT.semilogy(h_abep_km, F_mN, color=c, linewidth=2.0, label=f"M = {m} kg")

# Mark design point (10 kg @ 200 km)
axT.axvline(200, color="gray", linestyle="--", linewidth=1.2, alpha=0.7, label="200 km")
axT.plot(200, F_200_10kg, "k*", markersize=13, zorder=6,
         label=f"10 kg @ 200 km: {F_200_10kg:.2f} mN")

axT.set_xlabel("Altitude [km]", fontsize=11)
axT.set_ylabel("Required Thrust [mN]", fontsize=11)
axT.set_title(f"Thrust vs Altitude\n(A/m = {AM_REF} m\u00b2/kg, various masses)", fontsize=10)
axT.legend(fontsize=9, loc="upper right")
axT.grid(True, which="both", linestyle="--", alpha=0.5)
axT.set_xlim(150, 400)
axT.yaxis.set_minor_locator(ticker.LogLocator(subs="all", numticks=10))

# RIGHT: drag deceleration [um/s^2] for different A/m ratios (mass-independent)
AM_values  = [0.001, 0.005, 0.01, 0.05, 0.1]
cmap_am    = plt.get_cmap("viridis")
colors_am  = [cmap_am(x) for x in np.linspace(0.1, 0.85, len(AM_values))]

for am, c in zip(AM_values, colors_am):
    a_um = drag_base * am * 1e6    # drag deceleration [um/s^2]
    axA.semilogy(h_abep_km, a_um, color=c, linewidth=2.0, label=f"A/m = {am} m\u00b2/kg")

axA.axvline(200, color="gray", linestyle="--", linewidth=1.2, alpha=0.7, label="200 km")
ref_um = a_200 * 1e6
axA.plot(200, ref_um, "k*", markersize=13, zorder=6,
         label=f"Ref ({AM_REF} m\u00b2/kg): {ref_um:.1f} \u03bcm/s\u00b2")

axA.set_xlabel("Altitude [km]", fontsize=11)
axA.set_ylabel("Drag Deceleration [\u03bcm/s\u00b2]", fontsize=11)
axA.set_title("Specific Drag (= ABEP accel. needed)\nvs Altitude, for various A/m", fontsize=10)
axA.legend(fontsize=9, loc="upper right")
axA.grid(True, which="both", linestyle="--", alpha=0.5)
axA.set_xlim(150, 400)
axA.yaxis.set_minor_locator(ticker.LogLocator(subs="all", numticks=10))

plt.tight_layout()
plt.savefig("ABEP_Thrust.png", dpi=150, bbox_inches="tight")
plt.show()
print("Plot saved to ABEP_Thrust.png")

# ══════════════════════════════════════════════════════════════════════════════
#  OFF-NADIR POINTING ANALYSIS — SINGLE SATELLITE, VLEO 200 km
# ══════════════════════════════════════════════════════════════════════════════
#
#  MOTIVATION
#  ----------
#  Without off-nadir capability the VLEO satellite only covers the narrow strip
#  directly beneath its ground track (swath = 88.7 km).  The ground track shifts
#  ~22 deg in longitude between consecutive passes, so in any given orbit the
#  probability of overflying Galicia is very small.  The resulting single-
#  satellite revisit is ~244 h (>10 days) — unacceptable for fire detection.
#
#  The natural remedy is to add more satellites, but this section explores an
#  alternative: agile off-nadir pointing that effectively widens the accessible
#  coverage strip WITHOUT changing the orbit or adding hardware.
#
#  GEOMETRY OF OFF-NADIR POINTING
#  --------------------------------
#  The satellite rolls by angle Theta (boresight moves away from nadir).  The
#  camera FOV spans ±half_FOV = ±12.5 deg around the boresight.  The two swath
#  edges therefore sit at nadir angles:
#
#      alpha_near = Theta - half_FOV      (could be negative: swath crosses nadir)
#      alpha_far  = Theta + half_FOV      (always the widest reach)
#
#  For each nadir angle alpha the corresponding arc on the Earth's surface is
#  found via the triangle (Earth centre, satellite, ground point).  Applying
#  the sine rule:
#
#      sin(alpha + rho) / (R_E + h)  =  sin(alpha) / R_E
#   => rho(alpha) = arcsin( (R_E+h)/R_E * sin(alpha) ) - alpha       [rad]
#
#  where rho is the Earth central angle (arc from subsatellite point to the
#  ground point).  This formula holds as long as (R_E+h)/R_E * sin(alpha) < 1,
#  i.e. until the line of sight becomes tangent to the Earth (geometric horizon).
#
#  The slant range (straight-line distance satellite -> ground point) follows
#  from the same triangle:
#
#      R_s(alpha) = R_E * sin(rho) / sin(alpha)          [m]
#
#  (Limit: R_s -> h when alpha -> 0, i.e. nadir.)
#
#  EFFECTIVE SWATH WITH AGILE POINTING
#  -------------------------------------
#  When the satellite can roll to any angle between -Theta_max and +Theta_max, the
#  total accessible ground strip has a half-width equal to the Earth arc out to
#  alpha_far = Theta_max + half_FOV:
#
#      W_eff = 2 * rho(alpha_far) * R_E
#
#  This is the key figure of merit: it replaces the fixed nadir swath in the
#  revisit formula, and it grows rapidly with Theta_max because the Earth arc
#  rho(alpha) is a highly nonlinear function — small increases in nadir angle
#  at large look angles yield disproportionately large ground coverage gains.
#
#  GSD DEGRADATION AT THE FAR EDGE
#  ---------------------------------
#  Resolution degrades with slant range.  The standard (and conservative)
#  formula for the ground-projected pixel size at the far swath edge is:
#
#      GSD_far = IFOV * R_s(alpha_far)
#
#  where IFOV = pixel_size / focal_length = 15e-6 / 0.150 = 100 urad.
#  At nadir: GSD_nadir = IFOV * h = 100e-6 * 200e3 = 20 m.
#  GSD_far grows monotonically with Theta_max; it sets the resolution budget.
#
#  Note: the cross-track projection effect (dividing by cos of incidence angle)
#  further elongates pixels at oblique angles, so GSD_far above is the
#  along-track value; the cross-track value is somewhat larger for very high
#  look angles.  For the purposes of this trade (Theta_max ~ 50 deg) the slant-
#  range formula is the standard first-order metric used in mission design.
#
#  SINGLE-SATELLITE REVISIT WITH ENLARGED SWATH
#  ----------------------------------------------
#  The mean revisit time for one satellite over a target at latitude phi is:
#
#      T_rev = T_orbit * pi * R_E * cos(phi) / W_eff
#
#  Derivation: in one orbital period T_orbit the subsatellite point traces a
#  great-circle arc that advances by delta_lambda ~ 22 deg in longitude (J2-
#  corrected ground track shift).  The fraction of the latitude circle at phi
#  that falls within the accessible strip in a single pass is:
#
#      f = W_eff / (pi * R_E * cos(phi))      [ascending + descending, all-day]
#
#  Both the ascending and descending passes contribute independently (thermal IR
#  works day and night), and they cross lat phi at longitudes ~193 deg apart —
#  effectively two distinct coverage events per orbit, hence the factor pi
#  (half-circumference at phi) in the denominator rather than 2*pi.
#  The mean time between events is T_orbit / (2f) * 2 = T_orbit / f, which
#  gives the formula above.
#
#  DESIGN QUESTION AND RESULT
#  ---------------------------
#  The GSD budget caps the allowable slant range, which caps Theta_max, which
#  caps W_eff, which floors T_rev.  The question is:
#
#      "What is the maximum Theta_max before GSD_far exceeds 50 m
#       (= LEO nadir GSD, our cross-mission resolution parity limit),
#       and is the resulting T_rev <= 24 h?"
#
#  Answer:
#      Theta_max = 51.9 deg  =>  GSD_far = 50 m  =>  W_eff = 902.6 km
#      T_rev   = 23.96 h   < 24 h  (single satellite, all-day thermal IR)
#
#  The GSD parity constraint alone is sufficient to satisfy the 1-day revisit
#  requirement with a single VLEO satellite.  No additional satellites or
#  orbit changes are needed — only an agile attitude control system capable of
#  ±52 deg rolls within the time between consecutive ground-track passes
#  (~88 min).  This is demanding but demonstrated on platforms such as
#  Pleiades (~30 deg in 25 s) and SPOT-6/7.
# ══════════════════════════════════════════════════════════════════════════════

H_VL = 200e3   # [m]  VLEO design altitude

def _earth_central_angle(h, alpha_rad):
    """Earth central half-angle [rad] for nadir angle alpha [rad]."""
    s = (R_E + h) / R_E * math.sin(alpha_rad)
    if s >= 1.0:
        return None   # beyond geometric horizon
    return math.asin(s) - alpha_rad

def slant_range(h, alpha_rad):
    """Slant range [m] from satellite to ground point at nadir angle alpha [rad].
    Returns None if alpha_rad is beyond the geometric horizon."""
    if alpha_rad < 1e-9:
        return h   # nadir limit
    rho = _earth_central_angle(h, alpha_rad)
    if rho is None:
        return None
    return R_E * math.sin(rho) / math.sin(alpha_rad)

def gsd_far_edge(h, eta_rad):
    """GSD [m] at the far swath edge for boresight roll eta [rad].
    Returns NaN if beyond horizon."""
    alpha_far = eta_rad + half_FOV
    Rs = slant_range(h, alpha_far)
    if Rs is None:
        return float('nan')
    return IFOV_rad * Rs

def w_eff(h, eta_rad):
    """Effective accessible swath [m] with ±eta off-nadir agility."""
    alpha_far = eta_rad + half_FOV
    rho = _earth_central_angle(h, alpha_far)
    if rho is None:
        return None
    return 2.0 * rho * R_E

def revisit_offnadir(h, eta_rad, lat_deg):
    """Single-satellite mean revisit [s] with ±eta off-nadir agility (all-day)."""
    W = w_eff(h, eta_rad)
    if W is None or W <= 0:
        return float('inf')
    T = orbital_period(h)
    lat = math.radians(lat_deg)
    return T * math.pi * R_E * math.cos(lat) / W

# ── Solve for eta_max where GSD_far = 50 m (bisection) ────────────────────────
GSD_LIMIT_M = 50.0
eta_lo, eta_hi = 0.0, math.radians(70.0)
for _ in range(60):
    eta_mid = 0.5 * (eta_lo + eta_hi)
    if gsd_far_edge(H_VL, eta_mid) < GSD_LIMIT_M:
        eta_lo = eta_mid
    else:
        eta_hi = eta_mid
eta_max_rad = 0.5 * (eta_lo + eta_hi)
eta_max_deg = math.degrees(eta_max_rad)

alpha_far_deg = eta_max_deg + math.degrees(half_FOV)
W_eff_km      = w_eff(H_VL, eta_max_rad) / 1e3
T_rev_offnadir_h = revisit_offnadir(H_VL, eta_max_rad, galicia_lat_deg) / 3600.0
GSD_check     = gsd_far_edge(H_VL, eta_max_rad)

print("\n" + "=" * 70)
print("  OFF-NADIR ANALYSIS — VLEO 200 km  (GSD far-edge limit = 50 m)")
print("=" * 70)
print(f"  Nadir-only  swath          : {swath_spherical(H_VL)/1e3:.1f} km")
print(f"  Max off-nadir angle Theta  : {eta_max_deg:.1f} deg")
print(f"  Far-edge nadir angle       : {alpha_far_deg:.1f} deg")
print(f"  GSD at far edge            : {GSD_check:.2f} m  (limit = {GSD_LIMIT_M:.0f} m)")
print(f"  Effective swath W_eff      : {W_eff_km:.1f} km  (vs 88.7 km nadir-only)")
print(f"  Single-sat revisit         : {T_rev_offnadir_h:.2f} h  (target <= 24.0 h)")

if T_rev_offnadir_h <= 24.0:
    print(f"\n  Result: GSD<=50m constraint SATISFIES 1-day revisit. OK.")
else:
    excess_pct = 100.0 * (T_rev_offnadir_h - 24.0) / 24.0
    print(f"\n  Result: {excess_pct:.1f}% over 1-day target with GSD<=50m.")
    # Find eta that gives exactly 24 h (allow GSD to float)
    eta_lo2, eta_hi2 = eta_max_rad, math.radians(75.0)
    for _ in range(60):
        eta_mid2 = 0.5 * (eta_lo2 + eta_hi2)
        if revisit_offnadir(H_VL, eta_mid2, galicia_lat_deg) / 3600.0 > 24.0:
            eta_lo2 = eta_mid2
        else:
            eta_hi2 = eta_mid2
    eta_24h_rad = 0.5 * (eta_lo2 + eta_hi2)
    eta_24h_deg = math.degrees(eta_24h_rad)
    gsd_24h     = gsd_far_edge(H_VL, eta_24h_rad)
    w_24h_km    = w_eff(H_VL, eta_24h_rad) / 1e3
    print(f"\n  To achieve exactly 24 h revisit:")
    print(f"    Required Theta_max     : {eta_24h_deg:.1f} deg")
    print(f"    GSD at far edge        : {gsd_24h:.2f} m  (delta = +{gsd_24h-GSD_LIMIT_M:.2f} m)")
    print(f"    Effective swath        : {w_24h_km:.1f} km")
    print(f"    GSD relaxation needed  : {gsd_24h:.2f} m instead of {GSD_LIMIT_M:.0f} m"
          f"  ({100*(gsd_24h/GSD_LIMIT_M-1):.1f}% over limit)")

# ── Plot: GSD, swath, revisit vs eta_max ──────────────────────────────────────
eta_arr_deg = np.linspace(0, 65, 500)
eta_arr_rad = np.radians(eta_arr_deg)

gsd_arr  = np.array([gsd_far_edge(H_VL, e) for e in eta_arr_rad])
weff_raw = [w_eff(H_VL, e) for e in eta_arr_rad]
weff_arr = np.array([v / 1e3 if v is not None else float('nan') for v in weff_raw])
trev_arr = np.array([revisit_offnadir(H_VL, e, galicia_lat_deg) / 3600.0
                     for e in eta_arr_rad])
# Mask horizon-beyond values so plots stop cleanly
valid    = np.isfinite(gsd_arr) & np.isfinite(weff_arr) & np.isfinite(trev_arr)
eta_plot = eta_arr_deg[valid]
gsd_plot = gsd_arr[valid]
weff_plot= weff_arr[valid]
trev_plot= trev_arr[valid]

fig4, (ax_g, ax_w, ax_r) = plt.subplots(1, 3, figsize=(15, 5))
fig4.suptitle("Off-nadir Agility Trade — VLEO 200 km  (thermal IR, all-day)",
              fontsize=13, fontweight="bold")

RED  = "#e74c3c"
BLUE = "#2980b9"
GRN  = "#27ae60"

# -- GSD at far edge --
ax_g.plot(eta_plot, gsd_plot, color=RED, linewidth=2.5)
ax_g.axhline(GSD_LIMIT_M, color="k", linestyle="--", linewidth=1.4,
             label=f"GSD limit = {GSD_LIMIT_M:.0f} m (= LEO nadir)")
ax_g.axhline(20.0, color=GRN, linestyle=":", linewidth=1.2,
             label="VLEO nadir GSD = 20 m")
ax_g.axvline(eta_max_deg, color=RED, linestyle="--", linewidth=1.2, alpha=0.7,
             label=f"\u0398_max = {eta_max_deg:.1f} deg")
ax_g.plot(eta_max_deg, GSD_check, "o", color=RED, markersize=9, zorder=5)
ax_g.set_xlabel("Max off-nadir angle \u0398 [deg]", fontsize=11)
ax_g.set_ylabel("GSD at far swath edge [m]", fontsize=11)
ax_g.set_title("GSD degradation vs pointing angle", fontsize=10)
ax_g.legend(fontsize=8)
ax_g.grid(True, linestyle="--", alpha=0.5)
ax_g.set_xlim(0, 65)
ax_g.set_ylim(0, 120)

# -- Effective swath --
ax_w.plot(eta_plot, weff_plot, color=BLUE, linewidth=2.5)
ax_w.axhline(swath_spherical(H_VL) / 1e3, color=GRN, linestyle=":", linewidth=1.2,
             label=f"Nadir-only = {swath_spherical(H_VL)/1e3:.1f} km")
ax_w.axvline(eta_max_deg, color=RED, linestyle="--", linewidth=1.2, alpha=0.7,
             label=f"\u0398_max = {eta_max_deg:.1f} deg")
ax_w.plot(eta_max_deg, W_eff_km, "o", color=RED, markersize=9, zorder=5,
          label=f"W_eff = {W_eff_km:.0f} km at \u0398_max")
ax_w.set_xlabel("Max off-nadir angle \u0398 [deg]", fontsize=11)
ax_w.set_ylabel("Effective swath W_eff [km]", fontsize=11)
ax_w.set_title("Coverage swath vs pointing angle", fontsize=10)
ax_w.legend(fontsize=8)
ax_w.grid(True, linestyle="--", alpha=0.5)
ax_w.set_xlim(0, 65)
ax_w.set_ylim(0)

# -- Single-sat revisit --
ax_r.plot(eta_plot, trev_plot, color="#8e44ad", linewidth=2.5)
ax_r.axhline(24.0, color="k", linestyle="--", linewidth=1.4,
             label="1-day target = 24 h")
ax_r.axhline(single_sat_revisit(H_VL, galicia_lat_deg) / 3600, color=GRN,
             linestyle=":", linewidth=1.2,
             label=f"Nadir-only = {single_sat_revisit(H_VL, galicia_lat_deg)/3600:.0f} h")
ax_r.axvline(eta_max_deg, color=RED, linestyle="--", linewidth=1.2, alpha=0.7,
             label=f"\u0398_max (GSD=50m) = {eta_max_deg:.1f} deg")
ax_r.plot(eta_max_deg, T_rev_offnadir_h, "o", color=RED, markersize=9, zorder=5,
          label=f"T_rev = {T_rev_offnadir_h:.2f} h")
ax_r.set_xlabel("Max off-nadir angle \u0398 [deg]", fontsize=11)
ax_r.set_ylabel("Single-sat revisit [h]", fontsize=11)
ax_r.set_title("Revisit time vs pointing angle", fontsize=10)
ax_r.legend(fontsize=8)
ax_r.grid(True, linestyle="--", alpha=0.5)
ax_r.set_xlim(0, 65)
ax_r.set_ylim(0, 260)

plt.tight_layout()
plt.savefig("OffNadir_Analysis.png", dpi=150, bbox_inches="tight")
plt.show()
print("Plot saved to OffNadir_Analysis.png")
