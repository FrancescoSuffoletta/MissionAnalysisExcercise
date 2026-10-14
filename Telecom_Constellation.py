"""
================================================================================
  DIRECT-TO-DEVICE TELECOM CONSTELLATION — INITIAL MISSION ANALYSIS
  Focus design points: VLEO 200 km  vs  LEO 500 km (Starlink-class)
================================================================================

MISSION CONTEXT
---------------
The mission must provide DIRECT-TO-DEVICE (D2D) telecommunications service
over a defined geographic region.  "Direct-to-device" means the satellite
communicates directly with a commercial user handset (smartphone, tablet, IoT
sensor) WITHOUT any intermediate ground relay, gateway, or specialised terminal.

This is the architecture pioneered by AST SpaceMobile (BlueWalker 3 / BlueBird),
Lynk Global, SpaceX Starlink (direct-to-cell), and similar programmes.

WHY D2D CHANGES THE MISSION ANALYSIS COMPARED TO CLASSIC TELECOM
-----------------------------------------------------------------
In a classic VSAT or GEO broadcast system the user terminal has a dedicated
dish antenna (30–90 cm) that provides 30–50 dB of antenna gain, and the user
can tolerate 250 ms round-trip latency.  In a D2D system:

  1. USER TERMINAL IS A SMARTPHONE:
       - Antenna gain ≈ 0 dBi (near-isotropic, no tracking)
       - Transmit power ≈ 0.2–2 W (23–33 dBm)
       - Receiver noise figure ≈ 5–9 dB
       - Frequency band: must match existing cellular allocations
         (e.g. LTE Band 25 / Band 14 at S-band ~2 GHz, or B3/B7 at 1.7/2.6 GHz)
     → The satellite must compensate for the low terminal gain with a very
       large onboard antenna and high-sensitivity receiver.
     → Typical D2D EIRP deficit vs VSAT: 30–40 dB, requiring large phased arrays
       or deployable reflectors (AST SpaceMobile uses 64 m² aperture).

  2. SINGLE-SATELLITE REVISIT TIME SETS THE SERVICE GAP:
       - With only one satellite in the constellation, a ground user must wait
         for the next pass: T_rev = T_orb · π · R_E · cos(φ) / (2 λ · R_E)
       - At 44° N (northern Spain), ε_min = 20°:
           VLEO 200 km: T_rev ≈ 21.9 h (single sat)
           LEO  500 km: T_rev ≈ 10.9 h (single sat)
       - LEO has a dramatically shorter revisit because its wider footprint (λ)
         sweeps a larger swath, catching the target more frequently.
       - This is the KEY ORBIT-LEVEL DRIVER of constellation size:
           N_tier = ceil(T_rev / T_gap_allowed)
       - Doubling T_rev roughly doubles the minimum constellation needed for
         each service tier — a factor of ~2× penalty for VLEO vs LEO.

  3. MINIMUM ELEVATION ANGLE IS LINK-BUDGET DRIVEN, NOT GEOMETRY DRIVEN:
       - The link budget closes at a specific minimum elevation ε_min that
         depends on the satellite EIRP and the phone's G/T.
       - Urban shadowing and building blockage worsen the effective ε_min
         further: even if the geometry allows ε = 10°, a phone in a city
         may need ε ≥ 30–40° for reliable reception.
       - The mission analyst must input ε_min as a TRADE VARIABLE and iterate
         with the payload team until link budget and geometry converge.
       - For this study: ε_min = 20° is the reference (suburban D2D assumption, Spain).

  4. SERVICE CONTINUITY IS APPLICATION-DEPENDENT:
       - Voice call: needs uninterrupted link for the call duration → continuous
         coverage OR enough satellites for seamless handover (< 100 ms gap).
       - SMS / messaging: tolerates short outages → intermittent OK.
       - IoT / SOS: store-and-forward acceptable → single-satellite revisit
         of minutes to hours may be sufficient.
       - D2D constellations must define their service tier early, as it drives
         N by a factor of 10–100× between intermittent and continuous.

  5. PASS DURATION IS A HARD LIMIT ON SESSION LENGTH:
       - Even with perfect link quality, the satellite is only visible for
         T_pass minutes per orbit.
       - At VLEO 200 km, T_pass ≈ 1.4 min at ε = 30°: a phone call is barely
         possible within a single pass.
       - At LEO 500 km, T_pass ≈ 3.5 min: a short call is feasible.
       - For longer sessions, inter-satellite handover or a dense constellation
         (always another satellite above ε_min) is mandatory.

NOTATION
--------
  h         altitude above Earth surface                [m]
  R_E       Earth mean radius = 6371 km                [m]
  epsilon   minimum elevation angle                    [rad / deg]
  lambda    Earth central half-angle (footprint)       [rad / deg]
  R_s       slant range at min elevation               [m]
  T_orb     orbital period                             [s]
  v_orb     circular orbital velocity                  [m/s]
  T_pass    pass duration above epsilon_min            [s]
  T_rev     single-satellite mean revisit time         [s]
  N, P, S   Walker T/P/F constellation parameters      [-]
  FSPL      free-space path loss                       [dB]
================================================================================
"""

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.stdout.reconfigure(encoding='utf-8')

# ══════════════════════════════════════════════════════════════════════════════
#  PHYSICAL CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
R_E  = 6_371_000.0      # [m]       Earth mean radius
MU   = 3.986004418e14   # [m³/s²]  gravitational parameter
C    = 299_792_458.0    # [m/s]    speed of light
W_E  = 7.2921150e-5     # [rad/s]  Earth sidereal rotation rate

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 1 — MISSION INPUTS
#
#  WHAT MUST BE LOCKED BEFORE ANY ANALYSIS
#  -----------------------------------------
#  (a) SERVICE AREA:  latitude/longitude bounds.
#      Here: Spain 36–44° N (Iberian Peninsula + Balearic Islands).
#      → Sets required orbital inclination (i ≥ 49°) and cross-plane spacing.
#
#  (b) MINIMUM ELEVATION ANGLE ε_min [deg]:
#      For D2D the link budget closes only above a certain elevation.
#      Below ε_min the path loss is too high or shadowing blocks the link.
#      Reference values for D2D:
#        ε_min = 10°  → optimistic, flat open terrain, large sat antenna
#        ε_min = 20°  → realistic average, suburban environment          ← used here
#        ε_min = 30°  → conservative, urban canyon, high-rise density
#      Raising ε_min → shorter slant range (better link) but smaller footprint
#      → more satellites needed.  The payload team must confirm which value
#      the link budget supports.
#
#  (c) FREQUENCY BAND:
#      D2D must use licensed cellular spectrum so unmodified phones work.
#      Options:
#        L-band  (~1.5 GHz) — Satellite mobile (Inmarsat, Iridium heritage)
#        S-band  (~2.0 GHz) — MSS allocations, AST SpaceMobile, Lynk
#        AWS/PCS (~1.7/2.1 GHz) — Starlink D2D, T-Mobile partnership
#      Here: f = 2.0 GHz (S-band reference).
#
#  (d) LATENCY BUDGET [ms]:
#      Voice over IP (VoIP) target: one-way ≤ 150 ms (ITU-T G.114).
#      For D2D with LEO: easily met. GEO (~119 ms) is borderline.
#      Here: budget = 150 ms one-way.
#
#  (e) SERVICE CONTINUITY:
#      D2D requires at least INTERMITTENT service (short message / emergency).
#      Continuous voice service requires more satellites.
#      The reference design aims at continuous service.
#
#  (f) TARGET ALTITUDES FOR FOCUSED COMPARISON:
#      VLEO: 200 km  — very low Earth orbit, high drag, needs ABEP
#      LEO:  500 km  — Starlink-class, manageable drag, proven architecture
# ══════════════════════════════════════════════════════════════════════════════

EPSILON_MIN_DEG   = 20.0    # [deg]  reference min elevation (D2D, suburban average)
TARGET_LAT_MAX    = 44.0    # [deg]  northernmost latitude of service area (northern Spain)
TARGET_LAT_MIN    = 36.0    # [deg]  southernmost latitude (southern Spain)
FREQ_GHZ          = 2.0     # [GHz]  S-band D2D (unmodified phone compatible)
LATENCY_BUDGET_MS = 150.0   # [ms]   one-way latency budget (VoIP standard)

H_VLEO = 200e3   # [m]  VLEO design point
H_LEO  = 500e3   # [m]  LEO design point (Starlink-class)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 2 — FOOTPRINT GEOMETRY
#
#  EARTH CENTRAL HALF-ANGLE  λ(h, ε_min)
#  ----------------------------------------
#  The triangle Earth-centre (C) / satellite (S) / user (P) gives:
#
#      Angle at P = 90° + ε_min   (elevation is from horizontal, interior
#                                   angle toward Earth centre is 90° + ε)
#      Sine rule: (R_E+h)/cos(ε) = R_E/sin(angle_S)
#      → cos(ε + λ) = R_E/(R_E+h) · cos(ε)
#      → λ = arccos( R_E/(R_E+h) · cos(ε) ) − ε                          (1)
#
#  λ is the angular radius of the coverage disk (spherical cap) visible
#  above ε_min, measured at Earth's centre.
#
#  SLANT RANGE AT MIN ELEVATION
#  ------------------------------
#  Projecting along the line of sight from the ground:
#      R_s = sqrt( (R_E+h)² − R_E²cos²(ε) ) − R_E·sin(ε)               (2)
#
#  This is the MAXIMUM slant range within a pass (worst-case link range).
#  Minimum range = h (directly overhead).
#
#  PASS DURATION (overhead pass — best case)
#  -------------------------------------------
#  The satellite covers an arc of 2λ above ε_min in the orbital plane.
#  As a fraction of the full orbit:
#      T_pass = T_orb · λ / π                                             (3)
#  This is the maximum duration; an off-plane pass is shorter.
#  For D2D this is a KEY metric: T_pass sets the maximum uninterrupted
#  session length per satellite without handover.
#
#  SINGLE-SATELLITE REVISIT TIME
#  --------------------------------
#  For a satellite at altitude h over a target at latitude φ (all-day):
#      T_rev = T_orb · π · R_E · cos(φ) / W_eff                         (4)
#  where W_eff = 2·λ·R_E is the footprint diameter on the ground.
#  This is the mean time between consecutive passes of ONE satellite.
#  With N satellites: mean gap = T_rev / N.
# ══════════════════════════════════════════════════════════════════════════════

def earth_central_halfangle(h_m, eps_deg):
    """lambda(h, eps) per eq. (1). Returns NaN beyond geometric horizon."""
    eps = math.radians(eps_deg)
    arg = R_E / (R_E + h_m) * math.cos(eps)
    return float('nan') if arg > 1.0 else math.acos(arg) - eps

def slant_range_min(h_m, eps_deg):
    """Max slant range [m] at elevation eps_deg, eq. (2)."""
    eps = math.radians(eps_deg)
    return math.sqrt((R_E + h_m)**2 - (R_E * math.cos(eps))**2) - R_E * math.sin(eps)

def orbital_period(h_m):
    """Circular orbital period [s]."""
    return 2.0 * math.pi * math.sqrt((R_E + h_m)**3 / MU)

def orbital_velocity(h_m):
    """Circular orbital speed [m/s]."""
    return math.sqrt(MU / (R_E + h_m))

def pass_duration(h_m, eps_deg):
    """Max pass duration [s] for an overhead pass, eq. (3)."""
    lam = earth_central_halfangle(h_m, eps_deg)
    return 0.0 if math.isnan(lam) else orbital_period(h_m) * lam / math.pi

def single_sat_revisit(h_m, eps_deg, lat_deg):
    """
    Mean single-satellite revisit [s] over target at lat_deg, eq. (4).
    Uses W_eff = 2*lambda*R_E (nadir-only, no off-nadir agility assumed here).
    """
    lam = earth_central_halfangle(h_m, eps_deg)
    if math.isnan(lam) or lam <= 0:
        return float('inf')
    W = 2.0 * lam * R_E
    T = orbital_period(h_m)
    return T * math.pi * R_E * math.cos(math.radians(lat_deg)) / W

def fspl_dB(R_s_m, freq_GHz):
    """Free-space path loss [dB] = 20 log10(4 pi R_s f / c)."""
    return 20.0 * math.log10(4.0 * math.pi * R_s_m * freq_GHz * 1e9 / C)

def latency_ms(R_s_m):
    """One-way propagation latency [ms]."""
    return R_s_m / C * 1e3


def n_zonal_estimate(h_m, eps_deg, lat_max_deg):
    """
    Estimated N for continuous zonal coverage (0..lat_max_deg) using
    the Walker coverage conditions:
      P planes: 360° cos(lat_max) / (2 lambda) planes needed cross-track
      S per plane: 180° / lambda  (along-track continuous coverage)
    N = P * S.
    """
    lam = earth_central_halfangle(h_m, eps_deg)
    if math.isnan(lam) or lam <= 0:
        return float('inf')
    lam_d = math.degrees(lam)
    P = math.ceil(360.0 * math.cos(math.radians(lat_max_deg)) / (2.0 * lam_d))
    S = math.ceil(180.0 / lam_d)
    return int(P * S)

def walker_sizing(h_m, eps_deg, lat_max_deg):
    """Returns (P, S, N, F) for a Walker Delta covering 0..lat_max_deg."""
    lam_d = math.degrees(earth_central_halfangle(h_m, eps_deg))
    P = math.ceil(360.0 * math.cos(math.radians(lat_max_deg)) / (2.0 * lam_d))
    S = math.ceil(180.0 / lam_d)
    N = P * S
    F = round(P / 2)
    return P, S, N, F

def eclipse_fraction(h_m):
    """
    Fraction of orbit spent in Earth's shadow (penumbra ignored).
    Derived from the geometry: the shadow cone half-angle = arcsin(R_E/(R_E+h)).
    The shadow arc spans 2 * arcsin(R_E/(R_E+h)) out of 2π of the orbit.
    eclipse_frac = (1/π) * arcsin(R_E / (R_E+h))
    """
    return (1.0 / math.pi) * math.asin(R_E / (R_E + h_m))

def safe_arr(func, h_arr, *args):
    """Evaluate func(h, *args) over h_arr, returning NaN for failures."""
    out = []
    for h in h_arr:
        try:
            v = func(h, *args)
        except Exception:
            v = float('nan')
        out.append(v if not (isinstance(v, float) and v == float('inf')) else float('nan'))
    return np.array(out, dtype=float)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 3 — BROAD ORBIT SCREENING (LEO to GEO)
#
#  Even though the focused design is VLEO vs LEO, a first pass should always
#  sweep the full altitude range to understand where the design space sits:
#
#  KEY TRADE-OFF AXES:
#  - Altitude ↑ → λ ↑ (less satellites) but R_s ↑ (more path loss, more latency)
#  - For D2D the latency constraint generally rules out GEO and pushes
#    strongly toward LEO / MEO.
#  - The Van Allen belts (1000–6000 km and 13000–20000 km) make MEO hostile
#    for commercial satellites without expensive radiation hardening.
#  - GEO latency (~119 ms) is borderline for VoIP but possible; however GEO
#    cannot serve polar users and the large satellite is expensive.
#  → For D2D broadband the rational first shortlist is LEO 400–1200 km.
#  → VLEO (< 400 km) is unconventional but attractive if ABEP solves the drag.
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 76)
print("  DIRECT-TO-DEVICE TELECOM CONSTELLATION — INITIAL MISSION ANALYSIS")
print("=" * 76)
print(f"\n  Service area     : Spain  {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f} deg N  (Iberian Peninsula)")
print(f"  Min elevation    : {EPSILON_MIN_DEG:.0f} deg  (D2D suburban reference)")
print(f"  Frequency        : {FREQ_GHZ:.1f} GHz  (S-band, mobile cellular allocation)")
print(f"  Latency budget   : {LATENCY_BUDGET_MS:.0f} ms one-way  (ITU-T G.114 VoIP)")
print(f"  Service target   : continuous coverage, voice + data capable")

ORBIT_TABLE = [
    ("VLEO  200 km",  200e3),
    ("VLEO  300 km",  300e3),
    ("LEO   500 km",  500e3),
    ("LEO   800 km",  800e3),
    ("LEO  1200 km", 1200e3),
    ("MEO  5000 km", 5000e3),
    ("MEO 10000 km",10000e3),
    ("GEO 35786 km", 35786e3),
]

print("\n" + "-" * 76)
print(f"  {'Orbit':<17} {'lam[deg]':>8} {'Rs[km]':>7} {'FSPL[dB]':>9}"
      f" {'Lat[ms]':>8} {'Torb[min]':>10} {'Tpass[min]':>11}"
      f" {'Ecl[%]':>7} {'N_cont':>7}")
print("-" * 76)
for name, h in ORBIT_TABLE:
    lam  = earth_central_halfangle(h, EPSILON_MIN_DEG)
    Rs   = slant_range_min(h, EPSILON_MIN_DEG)
    fspl  = fspl_dB(Rs, FREQ_GHZ)
    tau   = latency_ms(Rs)
    Torb  = orbital_period(h) / 60.0
    Tpass = pass_duration(h, EPSILON_MIN_DEG) / 60.0
    ecl   = eclipse_fraction(h) * 100.0
    ncont = n_zonal_estimate(h, EPSILON_MIN_DEG, TARGET_LAT_MAX)
    lam_d = math.degrees(lam)
    flag  = "!LATE" if tau > LATENCY_BUDGET_MS else ""
    ncont_s = str(ncont) if ncont < 9999 else ">9999"
    print(f"  {name:<17} {lam_d:>8.1f} {Rs/1e3:>7.0f} {fspl:>9.1f}"
          f" {tau:>7.1f}{flag:<2} {Torb:>10.1f} {Tpass:>11.1f}"
          f" {ecl:>7.1f} {ncont_s:>7}")
print("-" * 76)
print("  !LATE = latency exceeds budget. Ecl = eclipse fraction per orbit.")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4 — FOCUSED COMPARISON: VLEO 200 km  vs  LEO 500 km
#
#  Both orbits are in the LEO family but their characteristics differ
#  significantly.  The following comparison is structured around the FOMs
#  that matter most for a D2D mission.
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 76)
print("  FOCUSED D2D COMPARISON — VLEO 200 km  vs  LEO 500 km")
print("=" * 76)

COMPARE = [
    ("VLEO 200 km", H_VLEO, "#e74c3c"),
    ("LEO  500 km", H_LEO,  "#2980b9"),
]

for name, h, _ in COMPARE:
    lam     = earth_central_halfangle(h, EPSILON_MIN_DEG)
    lam_d   = math.degrees(lam)
    Rs      = slant_range_min(h, EPSILON_MIN_DEG)
    fspl    = fspl_dB(Rs, FREQ_GHZ)
    tau     = latency_ms(Rs)
    Torb    = orbital_period(h) / 60.0
    Tpass   = pass_duration(h, EPSILON_MIN_DEG) / 60.0
    ecl     = eclipse_fraction(h) * 100.0
    v_orb   = orbital_velocity(h) / 1e3
    Trev1   = single_sat_revisit(h, EPSILON_MIN_DEG, TARGET_LAT_MAX) / 3600.0
    P, S, N, F = walker_sizing(h, EPSILON_MIN_DEG, TARGET_LAT_MAX)

    # Number of satellites for 2h mean revisit (intermittent D2D / SOS tier)
    N_2h = math.ceil(Trev1 / 2.0)

    print(f"\n  ── {name} ──────────────────────────────────────────────")
    print(f"    Footprint half-angle  λ    = {lam_d:.2f} deg")
    print(f"    Footprint diameter         = {2*lam_d*111.1:.0f} km  "
          f"  [note: 1 deg arc = 111.1 km at Earth surface]")
    print(f"    Slant range at ε = {EPSILON_MIN_DEG:.0f} deg = {Rs/1e3:.0f} km  (worst-case link range)")
    print(f"    Free-space path loss       = {fspl:.1f} dB  at {FREQ_GHZ:.1f} GHz")
    print(f"    One-way latency            = {tau:.2f} ms")
    print(f"    Orbital period             = {Torb:.1f} min")
    print(f"    Pass duration (overhead)   = {Tpass:.1f} min  → max uninterrupted session")
    print(f"    Orbital velocity           = {v_orb:.3f} km/s")
    print(f"    Eclipse fraction           = {ecl:.1f}%  (battery must cover this fraction)")
    print(f"    Single-sat revisit (44 N)  = {Trev1:.1f} h  (northern Spain)")
    print(f"    N for 2 h revisit          = {N_2h}  (intermittent SOS tier)")
    print(f"    Walker continuous coverage : {P} planes x {S} sats/plane = {N} total  (F={F})")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 4 COMMENTARY — WHY THE NUMBERS LOOK THIS WAY
#
#  FOOTPRINT SIZE:
#    At 200 km, λ ≈ 2.9°; at 500 km, λ ≈ 6.7°.  The footprint area scales
#    roughly as λ², so LEO 500 km has about (9.4/4.3)² ≈ 4.8× larger coverage
#    disk.  This directly multiplies the required number of satellites:
#    ~1260 needed at VLEO vs ~280 at LEO for Spain (36–44° N) coverage.
#
#  PATH LOSS ADVANTAGE OF VLEO:
#    FSPL = 20 log10(R_s · f / const).  At 200 km (R_s ≈ 383 km) vs 500 km
#    (R_s ≈ 909 km): ΔR_s = factor 2.37 → ΔFSPL = 20 log10(2.37) ≈ 7.5 dB.
#    Seven decibels is enormous in a D2D link budget: it can be traded for
#    lower satellite EIRP, smaller antenna, or higher user data rate.
#    This is the primary argument FOR VLEO in a D2D application.
#
#  PASS DURATION — THE D2D SESSION LIMIT:
#    VLEO 200 km: T_pass ≈ 2.1 min (overhead pass at ε = 20°).
#    LEO 500 km:  T_pass ≈ 4.9 min.
#    For a voice call:
#      - 2.1 min: tight — short message exchange feasible; voice requires handover.
#      - 4.9 min: marginal for a voice call; setup signalling leaves ~4 min.
#      - Either case requires seamless handover to the next satellite for
#        uninterrupted service → drives dense constellation or complex
#        handover protocol.
#    For IoT / SOS text message: both are more than sufficient.
#
#  SINGLE-SATELLITE REVISIT TIME:
#    T_rev(h, ε, φ) = T_orb · π · R_E · cos(φ) / (2 λ(h,ε) · R_E)
#    At 44° N (northern Spain), ε = 20°:  VLEO 200 km → 21.9 h;  LEO 500 km → 10.9 h.
#    LEO has a 2.0× shorter revisit because its footprint is 2.2× wider.
#    This is the most consequential orbit-level difference between the two
#    design points: it doubles the constellation size needed for every
#    non-continuous service tier (IoT SOS, messaging) at VLEO vs LEO.
#    Revisit scales with latitude (worse at high latitudes due to cos(φ))
#    and degrades rapidly as ε_min increases (smaller footprint).
#
#  SUSTAINABILITY OF VLEO — ABEP:
#    At 200 km, atmospheric drag is severe: without compensation, a
#    representative satellite (10 kg, A/m = 0.01 m²/kg) re-enters in < 1 day.
#    An Air-Breathing Electric Propulsion (ABEP) system is required, collecting
#    N₂/O molecules from the atmosphere to feed the thruster.
#    Required thrust ≈ 0.19 mN for a 10 kg satellite (see Propagate_Orbits.py).
#    This is a SYSTEM-LEVEL CONSTRAINT with no equivalent at 500 km where
#    chemical or electric thrusters with propellant are sufficient.
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 5 — CONSTELLATION ARCHITECTURE: HOW TO REASON ABOUT N, P, S, F
#
#  CONTINUOUS vs INTERMITTENT
#  ---------------------------
#  For D2D we define three service tiers requiring very different N:
#
#   Tier 1 – Emergency / IoT (intermittent, delay-tolerant):
#     A single pass per few hours is enough for an SOS beacon or IoT report.
#     N_tier1 = ceil(T_rev_single / T_gap_allowed)
#     At 200 km (Spain, 44°N): T_rev_single ≈ 21.9 h; for 2h gap → N ≈ 11 sats.
#     At 500 km (Spain, 44°N): T_rev_single ≈ 10.9 h; for 2h gap → N ≈ 6 sats.
#     Both very achievable with a small constellation.
#
#   Tier 2 – Messaging / data bursts (short gaps tolerated, ~30 min):
#     N_tier2 = ceil(T_rev / 0.5h)
#     At 200 km (Spain) → N ≈ 44; at 500 km → N ≈ 22.
#
#   Tier 3 – Voice / continuous broadband (continuous coverage):
#     N_tier3 = P * S (Walker continuous coverage condition)
#     At 200 km (Spain) → N ≈ 1260; at 500 km → N ≈ 280.
#
#  ONE PLANE vs MULTIPLE PLANES
#  ------------------------------
#  A single-plane constellation covers only one orbital corridor at a time.
#  As the plane precesses or the Earth rotates, users at other longitudes
#  get no service.  For D2D this is generally unacceptable unless the service
#  area is very narrow (e.g. equatorial belt) or the service is heavily
#  intermittent.
#
#  Rule: P ≥ ceil(360° cos(φ_max) / (2λ_deg))
#  At 200 km (ε=20°), λ = 4.3°: P ≥ ceil(360 × cos44° / (2 × 4.3)) = ceil(29.8) = 30 planes.
#  At 500 km (ε=20°), λ = 9.4°: P ≥ ceil(360 × cos44° / (2 × 9.4)) = ceil(13.8) = 14 planes.
#  Multiple planes are unavoidable for continuous coverage.
#
#  PHASING (F PARAMETER)
#  -----------------------
#  In a Walker T/P/F constellation the true anomaly offset between adjacent
#  planes is F × 360°/T.  For CONTINUOUS service F barely matters (there is
#  always a satellite overhead regardless of F).  For INTERMITTENT service
#  (Tier 1/2), F = P/2 distributes coverage events uniformly in time,
#  minimising the maximum gap.  F = 0 clusters events and is worst.
#
#  P vs S TRADE
#  -------------
#  Given a required N = P × S, choosing more planes (large P, small S)
#  improves geographic redundancy: losing one plane removes a narrow RAAN
#  corridor, not a large fraction of the Earth.  Choosing fewer planes
#  (small P, large S) reduces launch complexity but creates single-plane
#  failure vulnerability.  For D2D:
#    - P ≈ S ≈ sqrt(N) is the balanced starting point.
#    - Then adjust P downward to the minimum that closes the cross-track
#      coverage condition (eq. 6).
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 76)
print("  SERVICE TIER ANALYSIS — VLEO 200 km  vs  LEO 500 km")
print("=" * 76)
print(f"\n  {'Tier':<30} {'Condition':<25} {'N (VLEO)':>10} {'N (LEO)':>10}")
print("-" * 76)

TIERS = [
    ("Tier 1: Emergency/IoT",   "mean gap <= 2 h",  2.0),
    ("Tier 2: Messaging",       "mean gap <= 0.5 h", 0.5),
    ("Tier 3: Continuous voice","continuous, Spain", None),
]

for tier_name, cond, gap_h in TIERS:
    if gap_h is not None:
        Trev_v = single_sat_revisit(H_VLEO, EPSILON_MIN_DEG, TARGET_LAT_MAX) / 3600.0
        Trev_l = single_sat_revisit(H_LEO,  EPSILON_MIN_DEG, TARGET_LAT_MAX) / 3600.0
        Nv = math.ceil(Trev_v / gap_h)
        Nl = math.ceil(Trev_l / gap_h)
    else:
        Nv = n_zonal_estimate(H_VLEO, EPSILON_MIN_DEG, TARGET_LAT_MAX)
        Nl = n_zonal_estimate(H_LEO,  EPSILON_MIN_DEG, TARGET_LAT_MAX)
    print(f"  {tier_name:<30} {cond:<25} {Nv:>10} {Nl:>10}")
print("-" * 76)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 6 — ELEVATION ANGLE SENSITIVITY FOR D2D
#
#  The minimum elevation angle ε_min is the single most important parameter
#  the mission analyst does NOT own: it is set by the link budget.
#  But the mission analyst MUST communicate how sensitively N depends on it.
#
#  The logic is: ε_min drives λ (eq.1) → λ drives W_eff → W_eff drives N.
#  Increasing ε_min reduces λ (smaller footprint) → more satellites needed.
#
#  Sweep ε_min from 5° to 40° to show the sensitivity.
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 76)
print("  ELEVATION ANGLE SENSITIVITY  (continuous coverage, Spain 36–44 deg N)")
print("=" * 76)
print(f"\n  {'eps_min':>8} | {'VLEO 200km: lambda':>18} {'N_cont':>7}"
      f" | {'LEO 500km: lambda':>17} {'N_cont':>7}")
print("-" * 76)
for eps in [5, 10, 15, 20, 25, 30, 35, 40]:
    lv = math.degrees(earth_central_halfangle(H_VLEO, eps))
    ll = math.degrees(earth_central_halfangle(H_LEO,  eps))
    nv = n_zonal_estimate(H_VLEO, eps, TARGET_LAT_MAX)
    nl = n_zonal_estimate(H_LEO,  eps, TARGET_LAT_MAX)
    nv_s = str(nv) if nv < 99999 else ">99999"
    nl_s = str(nl) if nl < 99999 else ">99999"
    print(f"  {eps:>7} deg | {lv:>16.2f} deg {nv_s:>7}"
          f" | {ll:>15.2f} deg {nl_s:>7}")
print("-" * 76)
print("  Key insight: every 5 deg increase in eps_min roughly doubles N at VLEO.")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 7 — TELECOM CONSTRAINTS THE MISSION ANALYST MUST TRACK FOR D2D
#
#  The payload is defined by the telecom engineer.  However, the mission
#  analyst provides the GEOMETRIC INPUTS that constrain the link budget and
#  must flag any constraints that couple back to orbit/constellation design.
#
#  1. SLANT RANGE → PATH LOSS → MINIMUM SATELLITE EIRP:
#       FSPL = 20 log10(4π R_s f / c)  grows with altitude.
#       At 200 km (ε=20°): FSPL ≈ 152.9 dB; at 500 km: FSPL ≈ 160.0 dB.
#       The 7.5 dB difference at VLEO means the satellite can use a smaller
#       antenna or support higher user data rate for the same transmit power.
#       Mission analyst input: provide R_s(h, ε_min) for the link budget.
#
#  2. SINGLE-SATELLITE REVISIT TIME → SERVICE GAP:
#       T_rev(h, ε, φ) = T_orb · π · R_E · cos(φ) / (2 λ · R_E)
#       At 44° N (northern Spain), ε = 20°:
#         VLEO 200 km: T_rev ≈ 21.9 h  (single-sat)
#         LEO  500 km: T_rev ≈ 10.9 h  (single-sat)
#       N for a target max gap T_gap: N = ceil(T_rev / T_gap).
#       Revisit is WORSE at VLEO despite lower altitude, because the
#       footprint shrinks faster than the orbit period decreases.
#       Mission analyst output: T_rev(h, ε, φ) table for service tiering.
#
#  3. PASS DURATION → SESSION MANAGEMENT:
#       VLEO 200 km: T_pass ≈ 2.1 min.  The signalling overhead for LTE
#       attach, authentication, and call setup is 1–3 seconds → acceptable,
#       but handover must begin with >20 s margin → effective data time ~1.5 min.
#       LEO 500 km: T_pass ≈ 4.9 min → feasible for a short VoIP session.
#       Mission analyst output: T_pass as a function of h and ε_min must be
#       delivered to the network architecture team for handover protocol design.
#
#  4. INCLINATION → COVERAGE LATITUDE:
#       For Spain 36–44° N: inclination i ≈ 49–54° (Walker Delta, efficient).
#       This covers all of Spain and also southern France, Italy, and Greece.
#       Polar orbits (i = 90°) also cover Spain but waste capacity over poles.
#       Mission analyst must confirm with the constellation team that the chosen
#       inclination satisfies i ≥ lat_max + λ (Walker coverage condition).
#
#  5. ECLIPSE FRACTION → POWER SYSTEM:
#       At LEO the satellite enters Earth's shadow for a fraction of each orbit.
#       f_ecl = (1/π) arcsin(R_E / (R_E + h)):
#         200 km: ≈ 42%;   500 km: ≈ 38%.
#       Mission analyst provides this to the power team for battery sizing.
#       During eclipse the satellite must maintain the telecom payload on battery,
#       which affects the power/mass budget.
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 76)
print("  D2D TELECOM INTERFACE PARAMETERS")
print("=" * 76)
for name, h in [("VLEO 200 km", H_VLEO), ("LEO 500 km", H_LEO)]:
    Rs   = slant_range_min(h, EPSILON_MIN_DEG)
    fspl = fspl_dB(Rs, FREQ_GHZ)
    tpass= pass_duration(h, EPSILON_MIN_DEG) / 60.0
    ecl  = eclipse_fraction(h) * 100.0
    v    = orbital_velocity(h) / 1e3
    print(f"\n  {name}:")
    print(f"    FSPL at ε = {EPSILON_MIN_DEG:.0f} deg          = {fspl:.1f} dB")
    print(f"    Max pass duration          = {tpass:.1f} min  (hard session length limit)")
    print(f"    Orbital velocity           = {v:.3f} km/s")
    print(f"    Eclipse fraction           = {ecl:.1f}% of orbit")
    print(f"    Required inclination       >= {TARGET_LAT_MAX:.0f} + {math.degrees(earth_central_halfangle(h, EPSILON_MIN_DEG)):.1f}"
          f" = {TARGET_LAT_MAX + math.degrees(earth_central_halfangle(h, EPSILON_MIN_DEG)):.1f} deg  (Walker Delta condition)")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 7b — PATH-LOSS EQUIVALENCE: CAN VLEO AFFORD A LOWER ELEVATION ANGLE?
#
#  THE CORE QUESTION
#  ------------------
#  At equal ε_min both missions see the same sky geometry, and VLEO has 7 dB
#  less path loss.  But those 7 dB are "free" and could instead be spent on a
#  LONGER SLANT RANGE — that is, the VLEO satellite could still close the link
#  at a LOWER minimum elevation angle than LEO, achieving a LARGER footprint.
#
#  If the footprint grows enough, VLEO might require FEWER satellites than LEO
#  despite being at a much lower altitude.  This is the key trade-off.
#
#  MATHEMATICAL FORMULATION
#  -------------------------
#  The link closes when:
#      EIRP_sat + G/T_user - FSPL(R_s, f) >= C/N0_required       (link budget)
#  Since the payload is FIXED (same EIRP, same G/T for a given mission),
#  the maximum tolerable path loss is the same for both orbits: FSPL_max.
#
#  For LEO at 500 km and ε_min = ε_LEO:
#      FSPL_LEO = 20 log10(4π · R_s(h_LEO, ε_LEO) · f / c) = FSPL_max
#
#  For VLEO at 200 km, the equivalent ε_min is the elevation at which VLEO
#  reaches the SAME path loss as LEO:
#      FSPL_VLEO(ε_eq) = FSPL_max
#      => R_s(h_VLEO, ε_eq) = R_s(h_LEO, ε_LEO)     [same slant range!]
#
#  In other words: the maximum usable slant range is dictated by the link
#  budget (one value for each payload design), and the question is simply
#  what elevation angle that maximum slant range corresponds to at each altitude.
#  VLEO, being lower, must look further out to reach the same slant range,
#  meaning it looks at a lower elevation angle — which means its footprint is wider.
#
#  SLANT RANGE INVERSION
#  ----------------------
#  R_s(h, ε) = sqrt((R_E+h)² - R_E²cos²(ε)) - R_E sin(ε)
#  Given R_s_target and h, solve for ε:
#      Let X = R_E + h,  q = R_s_target + R_E sin(ε)
#      q² = X² - R_E² cos²(ε) = X² - R_E²(1-sin²(ε))
#      (R_s + R_E sin(ε))² = X² - R_E² + R_E² sin²(ε)
#      R_s² + 2 R_s R_E sin(ε) + R_E² sin²(ε) = X² - R_E² + R_E² sin²(ε)
#      R_s² + 2 R_s R_E sin(ε) = X² - R_E²
#      sin(ε) = (X² - R_E² - R_s²) / (2 R_s R_E)
#      ε_eq = arcsin( ((R_E+h)² - R_E² - R_s²) / (2 R_s R_E) )         (5)
#
#  This closed-form formula gives the minimum elevation for which the slant
#  range equals R_s_target exactly.
# ══════════════════════════════════════════════════════════════════════════════

VLEO_col = "#e74c3c"
LEO_col  = "#2980b9"

def eps_from_slant_range(h_m, Rs_target_m):
    """
    Invert the slant-range formula: given altitude h and a target slant range,
    return the elevation angle [deg] at which R_s(h, eps) = Rs_target.
    Formula (5): sin(eps) = ((R_E+h)^2 - R_E^2 - R_s^2) / (2 * R_s * R_E)
    Returns None if the target range is unreachable from altitude h.
    """
    X   = R_E + h_m
    num = X**2 - R_E**2 - Rs_target_m**2
    den = 2.0 * Rs_target_m * R_E
    s   = num / den
    if abs(s) > 1.0:
        return None     # target range is beyond horizon or below satellite
    return math.degrees(math.asin(s))

# ── Reference: LEO 500 km at ε = EPSILON_MIN_DEG (the "payload design point") ──
Rs_LEO_ref   = slant_range_min(H_LEO, EPSILON_MIN_DEG)
FSPL_LEO_ref = fspl_dB(Rs_LEO_ref, FREQ_GHZ)

# ── Equivalent ε_min for VLEO that gives the same slant range (same FSPL) ──
eps_eq_vleo = eps_from_slant_range(H_VLEO, Rs_LEO_ref)

# ── FOMs at the equivalent elevation ──
lam_eq_vleo  = math.degrees(earth_central_halfangle(H_VLEO, eps_eq_vleo))
Tpass_eq_min = pass_duration(H_VLEO, eps_eq_vleo) / 60.0
N_eq_vleo    = n_zonal_estimate(H_VLEO, eps_eq_vleo, TARGET_LAT_MAX)
P_eq, S_eq, N_eq_w, F_eq = walker_sizing(H_VLEO, eps_eq_vleo, TARGET_LAT_MAX)
N_LEO_ref    = n_zonal_estimate(H_LEO, EPSILON_MIN_DEG, TARGET_LAT_MAX)

print("\n" + "=" * 76)
print("  PATH-LOSS EQUIVALENCE ANALYSIS")
print("  Question: if VLEO uses its 7 dB advantage for a lower ε_min, does")
print("  it need FEWER satellites than LEO at its standard ε_min?")
print("=" * 76)
print(f"\n  REFERENCE (payload design point):")
print(f"    LEO  500 km at ε = {EPSILON_MIN_DEG:.0f} deg:")
print(f"      Slant range        = {Rs_LEO_ref/1e3:.0f} km")
print(f"      FSPL               = {FSPL_LEO_ref:.1f} dB   ← FSPL_max (link closes here)")
print(f"      N continuous       = {N_LEO_ref}")

print(f"\n  EQUIVALENT VLEO 200 km OPERATING POINT (same FSPL_max = {FSPL_LEO_ref:.1f} dB):")
print(f"    Equivalent ε_min    = {eps_eq_vleo:.1f} deg  (vs {EPSILON_MIN_DEG:.0f} deg for LEO)")
print(f"    Slant range         = {Rs_LEO_ref/1e3:.0f} km  (identical — same path loss by construction)")
print(f"    Footprint half-angle= {lam_eq_vleo:.1f} deg  (vs {math.degrees(earth_central_halfangle(H_LEO, EPSILON_MIN_DEG)):.1f} deg for LEO)")
print(f"    Pass duration       = {Tpass_eq_min:.1f} min  (vs {pass_duration(H_LEO, EPSILON_MIN_DEG)/60:.1f} min for LEO)")
print(f"    N continuous        = {N_eq_vleo}  (vs {N_LEO_ref} for LEO)")
print(f"    Walker sizing       = {P_eq} planes x {S_eq} sats/plane = {N_eq_w} total (F={F_eq})")

verdict = "VLEO WINS on constellation size" if N_eq_vleo < N_LEO_ref else "LEO still wins on constellation size"
print(f"\n  Verdict: {verdict}")
delta_N = N_LEO_ref - N_eq_vleo
delta_pct = 100.0 * delta_N / N_LEO_ref
if N_eq_vleo < N_LEO_ref:
    print(f"  VLEO needs {abs(delta_N)} fewer satellites ({abs(delta_pct):.0f}% reduction)")
else:
    print(f"  LEO needs {abs(delta_N)} fewer satellites ({abs(delta_pct):.0f}% reduction)")

# ══════════════════════════════════════════════════════════════════════════════
#  BUT WAIT — THE PRACTICAL LIMIT ON ε_min FOR D2D
#  -------------------------------------------------
#  The analysis above assumes the link budget is the only constraint on ε_min.
#  In reality, TWO independent constraints limit how low ε_min can go:
#
#  1. LINK BUDGET LIMIT (ε_link):
#     The minimum elevation at which FSPL ≤ FSPL_max.  This is what we just
#     computed.  For VLEO it is ε_link ≈ 4.4° given the same payload as LEO.
#
#  2. ENVIRONMENTAL / BLOCKAGE LIMIT (ε_block):
#     A phone at 4.4° elevation in a city sees a satellite near the horizon,
#     easily blocked by any building, hill, or even a parked truck.
#     Empirical data from terrestrial mobile networks:
#       Open rural:     effective ε_min ≈  5–10°  (mostly unobstructed)
#       Suburban:       effective ε_min ≈ 15–25°
#       Urban / dense:  effective ε_min ≈ 30–45°
#     The link cannot close if the satellite is physically blocked, regardless
#     of how strong the transmitted signal is.
#     For D2D, the OPERATIVE ε_min is:  ε_op = max(ε_link, ε_block)
#
#  THE CONCLUSION:
#  - In RURAL / OPEN terrain: ε_block is small → VLEO can exploit its link
#    advantage to operate near ε_link ≈ 4.4°.  Its footprint expands and it
#    NEEDS FEWER SATELLITES than LEO.
#  - In URBAN terrain: ε_block dominates and is LARGE (25–40°) for BOTH orbits.
#    Here VLEO's link advantage is partially wasted: both must point high to
#    clear buildings.  VLEO still has better FSPL margin (which can be used
#    for higher data rate or smaller antenna), but the footprint is no longer
#    the limiting factor — building blockage is.
#  - The CROSSOVER ELEVATION (ε_cross) is where N_VLEO = N_LEO.
#    Below ε_cross: VLEO needs fewer sats.  Above ε_cross: LEO is more efficient.
# ══════════════════════════════════════════════════════════════════════════════

# ── Find the crossover elevation (where N_VLEO = N_LEO at the SAME ε_min) ──
# Note: this is a DIFFERENT question — it asks at what common ε_min VLEO
# and LEO need the same number of satellites.  Below that, VLEO's larger
# footprint wins; above it, LEO is more efficient per satellite.
eps_cross = None
for eps_test in np.linspace(1, 45, 5000):
    nv = n_zonal_estimate(H_VLEO, eps_test, TARGET_LAT_MAX)
    nl = n_zonal_estimate(H_LEO,  eps_test, TARGET_LAT_MAX)
    if nv <= nl:
        eps_cross = eps_test
        break

print("\n" + "-" * 76)
print("  CROSSOVER ANALYSIS (where N_VLEO = N_LEO at the SAME ε_min):")
if eps_cross:
    nv_cross = n_zonal_estimate(H_VLEO, eps_cross, TARGET_LAT_MAX)
    nl_cross = n_zonal_estimate(H_LEO,  eps_cross, TARGET_LAT_MAX)
    print(f"  At ε_min < {eps_cross:.1f} deg: VLEO footprint large enough → N_VLEO <= N_LEO")
    print(f"  At ε_min > {eps_cross:.1f} deg: LEO is more efficient → N_VLEO > N_LEO")
    print(f"  At ε_cross = {eps_cross:.1f} deg: N_VLEO = {nv_cross}, N_LEO = {nl_cross}")
else:
    print("  No crossover found in 1–45 deg range: LEO always more efficient at equal ε_min.")

print(f"\n  D2D practical ε_min ranges:")
print(f"    Rural open terrain : ε ~ 5–10 deg   → VLEO advantage USABLE")
print(f"    Suburban           : ε ~ 15–25 deg  → mixed, depends on exact link budget")
print(f"    Urban dense        : ε ~ 30–45 deg  → blockage dominates, VLEO advantage REDUCED")

# ══════════════════════════════════════════════════════════════════════════════
#  REVERSE QUESTION: IF VLEO IS FIXED AT ε = 20°, WHAT MUST LEO DO TO MATCH?
#
#  The previous section asked: "what ε can VLEO afford if LEO operates at 20°?"
#  This section asks the REVERSE: "if both are forced to ε = 20° by suburban
#  blockage, what elevation would LEO need to match VLEO's path loss?"
#
#  The answer is direct from the same formula:
#    FSPL_VLEO(h=200km, ε=20°) = 152.9 dB  →  R_s_VLEO = 530 km
#    For LEO to also have R_s = 530 km:
#      ε_eq_LEO = arcsin( ((R_E+h_LEO)² - R_E² - R_s²) / (2 R_s R_E) )
#
#  Physically: a satellite at 500 km altitude has a slant range of 530 km
#  only when the user is looking nearly overhead — the satellite subtends a
#  very high elevation angle.  This reveals the asymmetry:
#
#    VLEO at 530 km slant range → ε = 20°  (suburban average, viable)
#    LEO  at 530 km slant range → ε ≈ 70°  (nearly directly overhead)
#
#  A minimum elevation of 77° means you only serve users when the satellite
#  is almost at zenith.  This is NOT a viable system design:
#    - The footprint at 77° elevation is tiny (λ ≈ 0.8°, ~180 km diameter)
#    - N_continuous would be enormous (> 5000 satellites)
#    - A pass at ε_min = 77° lasts only ~2 minutes at most elevations
#      but you are DISCARDING the entire low-elevation part of the pass
#
#  CONCLUSION (both directions of the argument combined):
#  ───────────────────────────────────────────────────────
#  The 7.5 dB path loss advantage of VLEO over LEO at equal elevation is
#  STRUCTURAL — it comes from the ~2.37× shorter slant range and cannot
#  be recovered by LEO through a different pointing geometry.
#
#  The advantage can be used in three mutually non-exclusive ways:
#    (1) Lower ε_min → wider footprint → fewer satellites   (rural D2D)
#    (2) Higher data rate for the same satellite EIRP        (capacity gain)
#    (3) Smaller satellite antenna for the same link margin  (mass/cost saving)
#
#  In urban D2D, where blockage forces both to ε ≈ 30°+, option (1) is
#  unavailable.  Options (2) and (3) remain valid and still justify VLEO —
#  but now the decisive disadvantage of VLEO comes back: ~4.5× more satellites
#  needed (1260 vs 280 for Spain) plus the ABEP propulsion requirement.  Whether the
#  path loss gain (→ smaller/cheaper satellite antenna) outweighs the higher
#  N is the core techno-economic trade for a D2D VLEO mission.
# ══════════════════════════════════════════════════════════════════════════════

# Reference: VLEO at ε = 30° (the blocked-environment constraint for both)
Rs_VLEO_20   = slant_range_min(H_VLEO, EPSILON_MIN_DEG)
FSPL_VLEO_20 = fspl_dB(Rs_VLEO_20, FREQ_GHZ)

# What elevation must LEO reach to hit the same slant range?
eps_eq_leo = eps_from_slant_range(H_LEO, Rs_VLEO_20)

print("\n" + "=" * 76)
print("  REVERSE EQUIVALENCE: VLEO fixed at ε = 20 deg — what must LEO match?")
print("=" * 76)
print(f"\n  VLEO 200 km at ε = {EPSILON_MIN_DEG:.0f} deg (blockage-constrained):")
print(f"    Slant range                = {Rs_VLEO_20/1e3:.0f} km")
print(f"    FSPL                       = {FSPL_VLEO_20:.1f} dB  ← VLEO path loss ceiling")

if eps_eq_leo is not None and eps_eq_leo <= 90:
    lam_eq_leo   = math.degrees(earth_central_halfangle(H_LEO, eps_eq_leo))
    Tpass_eq_leo = pass_duration(H_LEO, eps_eq_leo) / 60.0
    N_eq_leo     = n_zonal_estimate(H_LEO, eps_eq_leo, TARGET_LAT_MAX)
    Rs_check     = slant_range_min(H_LEO, eps_eq_leo)
    print(f"\n  For LEO 500 km to achieve FSPL = {FSPL_VLEO_20:.1f} dB:")
    print(f"    Required elevation         = {eps_eq_leo:.1f} deg  (vs {EPSILON_MIN_DEG:.0f} deg for VLEO)")
    print(f"    Slant range check          = {Rs_check/1e3:.0f} km  (= VLEO slant range ✓)")
    print(f"    Footprint half-angle       = {lam_eq_leo:.2f} deg  (vs {math.degrees(earth_central_halfangle(H_VLEO, EPSILON_MIN_DEG)):.2f} deg for VLEO)")
    print(f"    Footprint diameter         = {2*lam_eq_leo*111.1:.0f} km  (vs {2*math.degrees(earth_central_halfangle(H_VLEO, EPSILON_MIN_DEG))*111.1:.0f} km for VLEO)")
    print(f"    Pass duration at this ε    = {Tpass_eq_leo:.1f} min  (satellite usable ONLY above {eps_eq_leo:.0f} deg)")
    print(f"    N continuous (Spain)       = {N_eq_leo}  (vs {n_zonal_estimate(H_VLEO, EPSILON_MIN_DEG, TARGET_LAT_MAX)} for VLEO at ε={EPSILON_MIN_DEG:.0f} deg)")
    print(f"\n  Interpretation:")
    print(f"    LEO must be viewed at ε >= {eps_eq_leo:.0f} deg — nearly OVERHEAD — to match")
    print(f"    VLEO's path loss at ε = {EPSILON_MIN_DEG:.0f} deg.  This is operationally absurd:")
    print(f"    setting ε_min = {eps_eq_leo:.0f} deg discards the entire low-elevation portion")
    print(f"    of every pass, shrinks the footprint to {2*lam_eq_leo*111.1:.0f} km diameter, and")
    print(f"    requires {N_eq_leo} satellites for continuous coverage — far MORE than")
    print(f"    the {n_zonal_estimate(H_VLEO, EPSILON_MIN_DEG, TARGET_LAT_MAX)} VLEO satellites at normal ε = {EPSILON_MIN_DEG:.0f} deg.")
    print(f"\n  The 7 dB FSPL gap is STRUCTURAL and cannot be closed by LEO geometry.")
    print(f"  At equal ε_min LEO always has {FSPL_LEO_ref - FSPL_VLEO_20:.1f} dB more path loss than VLEO.")
else:
    print(f"\n  No valid LEO elevation exists: the target slant range ({Rs_VLEO_20/1e3:.0f} km)")
    print(f"  is shorter than LEO altitude ({H_LEO/1e3:.0f} km) — geometrically impossible.")

# ── Summary table: path-loss comparison at equal ε and at equivalent ε ───────
print("\n" + "-" * 76)
print("  SUMMARY TABLE — Four scenarios compared")
print("-" * 76)
print(f"  {'Scenario':<42} {'ε_VLEO':>7} {'ε_LEO':>7} {'N_VLEO':>7} {'N_LEO':>7}")
print("-" * 76)

def row(label, ev, el):
    nv = n_zonal_estimate(H_VLEO, ev, TARGET_LAT_MAX)
    nl = n_zonal_estimate(H_LEO,  el, TARGET_LAT_MAX)
    winner = "<" if nv < nl else (">" if nv > nl else "=")
    return f"  {label:<42} {ev:>6.1f}° {el:>6.1f}° {nv:>7} {winner} {nl:<7}"

print(row(f"Both at reference ε = {EPSILON_MIN_DEG:.0f}°  (suburban D2D)", EPSILON_MIN_DEG, EPSILON_MIN_DEG))
print(row(f"VLEO equiv-ε (LEO ref {EPSILON_MIN_DEG:.0f}°, same FSPL)",
          eps_eq_vleo, EPSILON_MIN_DEG))
print(row(f"LEO equiv-ε (VLEO ref {EPSILON_MIN_DEG:.0f}°, same FSPL)",
          EPSILON_MIN_DEG, eps_eq_leo if eps_eq_leo else 0))
print(row("Both at open-terrain ε = 5°", 5, 5))
print("-" * 76)
print("  < means VLEO needs fewer satellites; > means LEO needs fewer.")

# ══════════════════════════════════════════════════════════════════════════════
#  HOW P AND S WERE SELECTED — STEP-BY-STEP REASONING
#
#  After computing the footprint half-angle λ, the entire constellation sizing
#  reduces to answering two independent geometric questions:
#
#    Question 1 (P): How many orbital PLANES do I need so no user is ever in a
#                    longitude gap between adjacent planes?
#    Question 2 (S): How many satellites per PLANE do I need so the coverage
#                    disk never has a gap along the direction of flight?
#
#  These two questions are independent because planes handle the east-west
#  (RAAN) dimension while satellites-per-plane handle the along-track dimension.
#
#  ── SELECTING P (NUMBER OF PLANES) ──────────────────────────────────────────
#
#  Imagine looking down at the Earth from the North Pole.  Each orbital plane
#  is a vertical "slice" through the Earth, equally spaced in RAAN.  The
#  angular gap between adjacent planes is 360°/P.
#
#  At latitude φ the east-west ground distance between adjacent ascending nodes
#  is:  d_cross = (360°/P) × cos(φ) × R_E × π/180  [m]
#
#  For continuous coverage, adjacent footprint disks must just touch at the
#  worst-case latitude (φ_max).  The disk radius on the ground is λ × R_E.
#  Two adjacent disks just touch when the centre-to-centre distance ≤ 2λ R_E:
#
#      (360°/P) × cos(φ_max) ≤ 2λ_deg
#      → P ≥ ceil( 360° × cos(φ_max) / (2 × λ_deg) )               (6)
#
#  KEY INSIGHT: the cos(φ_max) factor is crucial.  At the equator (φ=0),
#  adjacent orbital planes are far apart in longitude, so you need many planes.
#  At high latitude, the planes converge, so fewer are needed.  We size for
#  the hardest case, which is the lowest latitude in the service area (φ = 0
#  for a mid-latitude constellation covering down to the equator).
#
#  Wait — for a zone 0..60°N, the WIDEST east-west spacing at the boundary is
#  at φ = 0° (equator).  But for a Walker Delta inclined at i ≥ φ_max + λ,
#  the ground track reaches φ_max not at the ascending node but near the top
#  of the arc.  The critical latitude for cross-track closure is φ_max because
#  that is the edge of the service area — if coverage closes there it closes
#  everywhere below.  Hence we use cos(φ_max) in eq.(6), not cos(0).
#
#  ── SELECTING S (SATELLITES PER PLANE) ──────────────────────────────────────
#
#  Within one orbital plane the satellites are equally spaced by 360°/S in
#  true anomaly (i.e., along the orbit circle).  As each satellite moves, it
#  sweeps the ground.  The along-track gap between consecutive satellite
#  footprints must be zero at all times.
#
#  The angular separation between consecutive satellites at Earth's centre is
#  360°/S.  For the two footprint disks to just touch:
#
#      360°/S ≤ 2λ_deg
#      → S ≥ ceil( 180° / λ_deg )                                    (7)
#
#  This is simpler than the plane condition because there is no cos(φ) factor
#  — the along-track chain must be unbroken everywhere along the orbital arc,
#  and the worst-case gap occurs at the equator crossing where the apparent
#  speed on the ground is highest.
#
#  ── WHY VLEO NEEDS MORE OF BOTH ─────────────────────────────────────────────
#
#  Both P and S are inversely proportional to λ.  Since λ_VLEO < λ_LEO
#  (smaller footprint), VLEO requires more planes AND more satellites per plane.
#  Specifically:  N = P × S ∝ 1/λ²  (scales as the SQUARE of the footprint).
#  Halving λ quadruples N — this is why the difference between 200 km and
#  500 km is so dramatic (1260 vs 280 for Spain at equal ε = 20°).
#
#  ── THE F PARAMETER — WHY IT DOES NOT CHANGE N ──────────────────────────────
#
#  Walker's phasing parameter F (0 ≤ F < P) shifts adjacent planes' satellites
#  in true anomaly by F × 360°/N relative to each other.  This controls WHEN
#  each satellite crosses a given latitude — it distributes coverage events
#  in time.  It does NOT change the footprint size or the number of gaps in
#  space — coverage closure is guaranteed by P and S alone.
#  For continuous coverage F is irrelevant (there is always a satellite up
#  regardless).  For intermittent service, F = P/2 minimises the maximum gap.
# ══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 76)
print("  CONSTELLATION SIZING — WHY THESE P AND S VALUES")
print("=" * 76)

for name, h, eps in [(f"VLEO 200 km  (suburban, ε={EPSILON_MIN_DEG:.0f}°)", H_VLEO, EPSILON_MIN_DEG),
                     (f"LEO  500 km  (suburban, ε={EPSILON_MIN_DEG:.0f}°)", H_LEO,  EPSILON_MIN_DEG),
                     (f"VLEO 200 km  (rural, equiv-ε={eps_eq_vleo:.1f}°)", H_VLEO, eps_eq_vleo)]:
    lam_d = math.degrees(earth_central_halfangle(h, eps))
    P = math.ceil(360.0 * math.cos(math.radians(TARGET_LAT_MAX)) / (2.0 * lam_d))
    S = math.ceil(180.0 / lam_d)
    N = P * S
    F = round(P / 2)
    # Cross-track gap at phi_max (should be ~= 2*lambda or just under)
    gap_cross_deg = 360.0 * math.cos(math.radians(TARGET_LAT_MAX)) / P
    # Along-track gap
    gap_along_deg = 360.0 / S

    print(f"\n  {name}:")
    print(f"    λ = {lam_d:.2f} deg")
    print(f"")
    print(f"    P = ceil( 360 x cos({TARGET_LAT_MAX:.0f}) / (2 x {lam_d:.2f}) )")
    print(f"      = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX)):.1f} / {2*lam_d:.2f} )")
    print(f"      = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX))/(2*lam_d):.2f} ) = {P} planes")
    print(f"      Cross-track arc between planes at {TARGET_LAT_MAX:.0f} deg N: "
          f"{gap_cross_deg:.2f} deg  (= {gap_cross_deg/lam_d:.2f} x lambda  ← must be <= 2)")

    print(f"")
    print(f"    S = ceil( 180 / {lam_d:.2f} )")
    print(f"      = ceil( {180/lam_d:.2f} ) = {S} sats/plane")
    print(f"      Along-track arc between sats: "
          f"{gap_along_deg:.2f} deg  (= {gap_along_deg/lam_d:.2f} x lambda  ← must be <= 2)")

    print(f"")
    print(f"    N = {P} x {S} = {N}  (F = {F})")
# This plot goes into a separate figure for clarity.
eps_sweep   = np.linspace(1, 45, 300)
N_vleo_arr  = np.array([n_zonal_estimate(H_VLEO, e, TARGET_LAT_MAX) for e in eps_sweep], dtype=float)
N_leo_arr   = np.array([n_zonal_estimate(H_LEO,  e, TARGET_LAT_MAX) for e in eps_sweep], dtype=float)

# Equivalent-ε curve: VLEO operates at eps_eq(ε_LEO) where slant ranges match
eps_leo_arr = np.linspace(1, 45, 300)
N_vleo_equiv = []
for el in eps_leo_arr:
    Rs_l = slant_range_min(H_LEO, el)
    eq   = eps_from_slant_range(H_VLEO, Rs_l)
    if eq is None or eq < 0:
        N_vleo_equiv.append(float('nan'))
    else:
        N_vleo_equiv.append(n_zonal_estimate(H_VLEO, eq, TARGET_LAT_MAX))
N_vleo_equiv = np.array(N_vleo_equiv, dtype=float)
N_leo_ref_arr = np.array([n_zonal_estimate(H_LEO, el, TARGET_LAT_MAX)
                           for el in eps_leo_arr], dtype=float)

fig2, (axL, axR) = plt.subplots(1, 2, figsize=(14, 6))
fig2.suptitle("VLEO 200 km vs LEO 500 km — Path-Loss Equivalence Trade",
              fontsize=12, fontweight="bold")

# LEFT: N vs ε_min at the SAME elevation (equal playing field)
axL.semilogy(eps_sweep, N_vleo_arr, color=VLEO_col, linewidth=2.5,
             label="VLEO 200 km (same ε_min as LEO)")
axL.semilogy(eps_sweep, N_leo_arr,  color=LEO_col,  linewidth=2.5,
             label="LEO  500 km")
if eps_cross:
    axL.axvline(eps_cross, color="green", linestyle="--", linewidth=1.4,
                label=f"Crossover ε = {eps_cross:.1f}°  (VLEO wins below)")
axL.axvline(EPSILON_MIN_DEG, color="gray", linestyle=":", linewidth=1.2,
            label=f"Reference ε = {EPSILON_MIN_DEG:.0f}°")
axL.fill_betweenx([1, 10000], 0, eps_cross if eps_cross else 0,
                  alpha=0.08, color="green", label="VLEO footprint advantage")
axL.fill_betweenx([1, 10000], eps_cross if eps_cross else 0, 45,
                  alpha=0.08, color=LEO_col,  label="LEO efficiency advantage")
axL.set_xlabel("Min elevation angle ε [deg]  (same for both)", fontsize=10)
axL.set_ylabel("N satellites (continuous, Spain 36–44° N)", fontsize=10)
axL.set_title("(A) Equal-ε comparison\nVLEO's larger footprint wins at low elevation",
              fontsize=10)
axL.legend(fontsize=8)
axL.grid(True, which="both", linestyle="--", alpha=0.4)
axL.set_xlim(1, 45)

# RIGHT: N vs LEO's ε_min — VLEO uses equivalent (lower) ε to match path loss
axR.semilogy(eps_leo_arr, N_leo_ref_arr, color=LEO_col, linewidth=2.5,
             label=f"LEO 500 km at ε_LEO")
axR.semilogy(eps_leo_arr, N_vleo_equiv,  color=VLEO_col, linewidth=2.5,
             linestyle="--",
             label="VLEO 200 km at equiv. ε (same FSPL as LEO)")
axR.axvline(EPSILON_MIN_DEG, color="gray", linestyle=":", linewidth=1.2,
            label=f"Reference ε_LEO = {EPSILON_MIN_DEG:.0f}°")
# Mark the design point
Rs_ref = slant_range_min(H_LEO, EPSILON_MIN_DEG)
eq_ref = eps_from_slant_range(H_VLEO, Rs_ref)
nv_ref_eq = n_zonal_estimate(H_VLEO, eq_ref, TARGET_LAT_MAX)
axR.plot(EPSILON_MIN_DEG, N_LEO_ref,  "o", color=LEO_col,  markersize=10, zorder=6,
         label=f"LEO design pt: N={N_LEO_ref}")
axR.plot(EPSILON_MIN_DEG, nv_ref_eq,  "o", color=VLEO_col, markersize=10, zorder=6,
         label=f"VLEO equiv pt: N={nv_ref_eq} (ε={eq_ref:.1f}°)")
axR.set_xlabel("LEO reference ε_min [deg]  (VLEO ε auto-adjusted for equal FSPL)", fontsize=10)
axR.set_ylabel("N satellites (continuous, Spain 36–44° N)", fontsize=10)
axR.set_title("(B) Equal-FSPL comparison\n'VLEO spends its 7 dB on a lower elevation'",
              fontsize=10)
axR.legend(fontsize=8)
axR.grid(True, which="both", linestyle="--", alpha=0.4)
axR.set_xlim(1, 45)

plt.tight_layout()
plt.savefig("PathLoss_Equivalence.png", dpi=150, bbox_inches="tight")
plt.show()
print("\n  Plot saved to PathLoss_Equivalence.png")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 8 — RECOMMENDED STUDY FLOW
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 76)
print("  RECOMMENDED STUDY FLOW — Phase 1: Requirements to Constellation Candidates")
print("=" * 76)
FLOW = [
    ("1", "Lock D2D requirements",
     "Service area, user terminal type (phone/IoT), service tier (emergency/\n"
     "        voice/broadband), frequency band (must match handset licensing)."),
    ("2", "Define FOMs",
     "lambda, R_s, FSPL, latency, T_pass, T_rev (single-sat revisit), N per service tier."),
    ("3", "Payload interface",
     "Obtain ε_min from link budget (do NOT assume 10 deg — get it from payload\n"
     "        team). This single number drives all downstream geometry."),
    ("4", "Altitude trade sweep",
     "Plot all FOMs vs h (300–1200 km). Highlight Van Allen belt (skip 1000-\n"
     "        6000 km). Mark VLEO 200 km and LEO 500 km as design points."),
    ("5", "Service tier N estimate",
     "For each design altitude compute N for Tier 1 (SOS), Tier 2 (messaging),\n"
     "        Tier 3 (continuous). Present cost/complexity vs capability."),
    ("6", "Inclination selection",
     "Walker Delta at i = lat_max + lambda is minimum energy. Compare with\n"
     "        polar (i=90) for global coverage if required."),
    ("7", "Walker T/P/F pre-sizing",
     "Apply coverage conditions to get P, S, F for each shortlisted orbit.\n"
     "        Confirm N x satellite cost is within programme budget."),
    ("8", "Revisit time vs service tier",
     "Compute T_rev(h, ε, φ) at target latitude. Verify N = T_rev/T_gap\n"
     "        is consistent with the mission cost/launch budget."),
    ("9", "ABEP / propulsion check",
     "For VLEO: confirm ABEP thrust budget closes (0.19 mN per 10 kg sat).\n"
     "        For LEO 500 km: confirm electric prop propellant budget for lifetime."),
    ("10","Handover protocol",
     "Derive T_pass curve. Confirm network team can design handover within\n"
     "        T_pass - margin (target: handover starts 20s before loss of signal)."),
    ("11","Down-select to 2 candidates",
     "Score VLEO vs LEO on: N, FSPL gain, T_pass, T_rev, propulsion complexity,\n"
     "        regulatory risk. Present to system PDR."),
]
for step, title, desc in FLOW:
    print(f"\n  {step:>2}. {title}")
    print(f"      {desc}")

# ══════════════════════════════════════════════════════════════════════════════
#  STEP 9 — MAIN RISKS AND OPEN ITEMS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 76)
print("  MAIN RISKS AND UNKNOWNS  (geometry alone cannot close these)")
print("=" * 76)
RISKS = [
    ("ε_min not yet from link budget",
     "All N estimates now assume ε_min = 20 deg (suburban D2D baseline, Spain).\n"
     "     N at VLEO is ~1260 vs ~2835 at ε_min = 30 deg. If the link budget\n"
     "     forces ε_min = 30 deg, the constellation nearly triples.\n"
     "     This is the highest-priority iteration with the payload team."),
    ("Revisit time gap at VLEO",
     "Single-sat T_rev ≈ 21.9 h at VLEO vs 10.9 h at LEO (44° N Spain, ε=20°).\n"
     "     Every additional service tier requires proportionally more satellites.\n"
     "     The 2× revisit penalty at VLEO is a hard orbital geometry consequence."),
    ("ABEP readiness for VLEO (TRL)",
     "Air-Breathing Electric Propulsion is still TRL 4–5 (2025). If ABEP is\n"
     "     not available by launch date, VLEO is not sustainable and falls back to\n"
     "     LEO, eliminating the path-loss advantage."),
    ("Handset software compatibility",
     "D2D requires over-the-air firmware updates to the phone or coordination\n"
     "     with MNO/chipset vendors. Regulatory approval (FCC, ETSI) takes 1–3\n"
     "     years and must start early. Failure = entire market access blocked."),
    ("Interference with terrestrial networks",
     "The satellite downlink overlaps terrestrial base station frequencies.\n"
     "     The ITU/national regulators may impose power flux density (PFD) limits\n"
     "     that reduce the satellite EIRP → increases required ε_min → more sats."),
    ("Constellation deployment phasing",
     "Revenue requires a minimum functional constellation. For Tier 3 at LEO\n"
     "     (~280 sats for Spain) this may take 2–4 launch campaigns. Must define the minimum\n"
     "     operational constellation for Tier 1/2 service early revenue."),
    ("Radiation and lifetime at VLEO",
     "At 200 km the atmosphere provides partial shielding from cosmic rays,\n"
     "     but solar particle events can still damage electronics. Lower altitude\n"
     "     reduces design lifetime marginally; more frequent replacements drive cost."),
    ("Ground segment architecture",
     "LEO handover rate (one pass every ~90 min per satellite) requires either\n"
     "     many ground stations or inter-satellite links (ISL). ISL add mass/cost.\n"
     "     Mission analyst must size the gateway network alongside the space segment."),
]
for risk, desc in RISKS:
    print(f"\n  [!] {risk}")
    print(f"      {desc}")

print("\n" + "=" * 76)

# ══════════════════════════════════════════════════════════════════════════════
#  PLOTS
# ══════════════════════════════════════════════════════════════════════════════

h_sweep = np.logspace(np.log10(150e3), np.log10(15000e3), 500)
eps_list   = [5.0, 10.0, 20.0, 30.0]
cols_eps   = ["#27ae60", "#2980b9", "#e67e22", "#e74c3c"]
VLEO_col   = "#e74c3c"
LEO_col    = "#2980b9"

fig, axes = plt.subplots(2, 3, figsize=(17, 10))
fig.suptitle(
    f"D2D Telecom Constellation — FOMs vs Altitude  "
    f"(f = {FREQ_GHZ:.1f} GHz S-band, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N)",
    fontsize=12, fontweight="bold"
)

def mark_design_pts(ax, y_func, eps=EPSILON_MIN_DEG, y_scale=1.0):
    """Add VLEO/LEO markers to an axis."""
    for h, label, col in [(H_VLEO, "VLEO\n200 km", VLEO_col),
                          (H_LEO,  "LEO\n500 km",  LEO_col)]:
        yv = y_func(h, eps) * y_scale
        ax.axvline(h / 1e3, color=col, linestyle="--", linewidth=1.2, alpha=0.6)
        ax.plot(h / 1e3, yv, "o", color=col, markersize=9, zorder=6)

# ── A: Footprint half-angle λ ──────────────────────────────────────────────
ax = axes[0, 0]
for eps, col in zip(eps_list, cols_eps):
    lam_arr = safe_arr(lambda h, e=eps: math.degrees(earth_central_halfangle(h, e)), h_sweep)
    ax.semilogx(h_sweep / 1e3, lam_arr, color=col, linewidth=2.0, label=f"ε={eps:.0f}°")
mark_design_pts(ax, lambda h, e: math.degrees(earth_central_halfangle(h, e)))
ax.set_xlabel("Altitude [km]", fontsize=10)
ax.set_ylabel("Footprint half-angle λ [deg]", fontsize=10)
ax.set_title("(A) Coverage footprint radius", fontsize=10)
ax.legend(fontsize=8); ax.grid(True, which="both", linestyle="--", alpha=0.4)

# ── B: Continuous coverage N ───────────────────────────────────────────────
ax = axes[0, 1]
for eps, col in zip(eps_list, cols_eps):
    N_arr = safe_arr(lambda h, e=eps: n_zonal_estimate(h, e, TARGET_LAT_MAX), h_sweep)
    ax.semilogy(h_sweep / 1e3, N_arr, color=col, linewidth=2.0, label=f"ε={eps:.0f}°")
for h, col in [(H_VLEO, VLEO_col), (H_LEO, LEO_col)]:
    ax.axvline(h / 1e3, color=col, linestyle="--", linewidth=1.2, alpha=0.6)
ax.set_xlabel("Altitude [km]", fontsize=10)
ax.set_ylabel("N satellites (continuous coverage)", fontsize=10)
ax.set_title(f"(B) Constellation size (continuous, 0–{TARGET_LAT_MAX:.0f}° N)", fontsize=10)
ax.legend(fontsize=8); ax.grid(True, which="both", linestyle="--", alpha=0.4)
ax.set_xlim(150, 15000)

# ── C: Free-space path loss ────────────────────────────────────────────────
ax = axes[0, 2]
Rs_arr   = safe_arr(slant_range_min, h_sweep, EPSILON_MIN_DEG)
fspl_arr = np.array([fspl_dB(r, FREQ_GHZ) for r in Rs_arr])
ax.semilogx(h_sweep / 1e3, fspl_arr, color="#8e44ad", linewidth=2.5)
for h, label, col in [(H_VLEO, "VLEO", VLEO_col), (H_LEO, "LEO", LEO_col)]:
    Rs  = slant_range_min(h, EPSILON_MIN_DEG)
    yv  = fspl_dB(Rs, FREQ_GHZ)
    ax.axvline(h / 1e3, color=col, linestyle="--", linewidth=1.2, alpha=0.6)
    ax.plot(h / 1e3, yv, "o", color=col, markersize=9, zorder=6,
            label=f"{label}: {yv:.1f} dB")
ax.set_xlabel("Altitude [km]", fontsize=10)
ax.set_ylabel("FSPL [dB]", fontsize=10)
ax.set_title(f"(C) Path loss  (ε={EPSILON_MIN_DEG:.0f}°, f={FREQ_GHZ:.1f} GHz)", fontsize=10)
ax.legend(fontsize=8); ax.grid(True, which="both", linestyle="--", alpha=0.4)

# ── D: Pass duration ───────────────────────────────────────────────────────
ax = axes[1, 0]
for eps, col in zip(eps_list, cols_eps):
    tp_arr = safe_arr(lambda h, e=eps: pass_duration(h, e) / 60.0, h_sweep)
    ax.semilogx(h_sweep / 1e3, tp_arr, color=col, linewidth=2.0, label=f"ε={eps:.0f}°")
mark_design_pts(ax, lambda h, e: pass_duration(h, e) / 60.0)
ax.axhline(5.0, color="gray", linestyle=":", linewidth=1.2,
           label="5 min reference (short call)")
ax.set_xlabel("Altitude [km]", fontsize=10)
ax.set_ylabel("Max pass duration [min]", fontsize=10)
ax.set_title("(D) Session length limit per satellite", fontsize=10)
ax.legend(fontsize=8); ax.grid(True, which="both", linestyle="--", alpha=0.4)

# ── E: Single-satellite mean revisit time vs altitude ──────────────────────
ax = axes[1, 1]
lat_list   = [0.0, 30.0, 45.0, 60.0]
lat_colors = ["#16a085", "#2980b9", "#e67e22", "#e74c3c"]
for lat, col in zip(lat_list, lat_colors):
    trev_arr = safe_arr(
        lambda h, e=EPSILON_MIN_DEG, phi=lat: single_sat_revisit(h, e, phi) / 3600.0,
        h_sweep)
    ax.semilogy(h_sweep / 1e3, trev_arr, color=col, linewidth=2.0,
                label=f"φ = {lat:.0f}° N")
for h, label, col in [(H_VLEO, "VLEO", VLEO_col), (H_LEO, "LEO", LEO_col)]:
    trev_ref = single_sat_revisit(h, EPSILON_MIN_DEG, TARGET_LAT_MAX) / 3600.0
    ax.axvline(h / 1e3, color=col, linestyle="--", linewidth=1.2, alpha=0.6)
    ax.plot(h / 1e3, trev_ref, "o", color=col, markersize=9, zorder=6,
            label=f"{label}: {trev_ref:.1f} h  (60° N)")
ax.set_xlabel("Altitude [km]", fontsize=10)
ax.set_ylabel("Single-sat revisit time [h]", fontsize=10)
ax.set_title(f"(E) Single-sat revisit  (ε = {EPSILON_MIN_DEG:.0f}°, service 0–60° N)", fontsize=10)
ax.legend(fontsize=8); ax.grid(True, which="both", linestyle="--", alpha=0.4)

# ── F: VLEO vs LEO head-to-head N per service tier ──────────────────────────
ax = axes[1, 2]
eps_range = np.linspace(5, 45, 100)

N_cont_vleo = [n_zonal_estimate(H_VLEO, e, TARGET_LAT_MAX) for e in eps_range]
N_cont_leo  = [n_zonal_estimate(H_LEO,  e, TARGET_LAT_MAX) for e in eps_range]

# Tier 2: 30 min gap
N_t2_vleo = [math.ceil(single_sat_revisit(H_VLEO, e, TARGET_LAT_MAX)/3600 / 0.5)
             for e in eps_range]
N_t2_leo  = [math.ceil(single_sat_revisit(H_LEO,  e, TARGET_LAT_MAX)/3600 / 0.5)
             for e in eps_range]

ax.semilogy(eps_range, N_cont_vleo, color=VLEO_col, linewidth=2.5,
            label="VLEO 200 km — continuous")
ax.semilogy(eps_range, N_cont_leo,  color=LEO_col,  linewidth=2.5,
            label="LEO  500 km — continuous")
ax.semilogy(eps_range, N_t2_vleo,   color=VLEO_col, linewidth=1.5, linestyle="--",
            label="VLEO 200 km — 30 min gap")
ax.semilogy(eps_range, N_t2_leo,    color=LEO_col,  linewidth=1.5, linestyle="--",
            label="LEO  500 km — 30 min gap")
ax.axvline(EPSILON_MIN_DEG, color="gray", linestyle=":", linewidth=1.2,
           label=f"Reference ε = {EPSILON_MIN_DEG:.0f} deg")
ax.set_xlabel("Min elevation angle ε [deg]", fontsize=10)
ax.set_ylabel("N satellites required", fontsize=10)
ax.set_title("(F) VLEO vs LEO: N vs ε_min  (0–60° N)", fontsize=10)
ax.legend(fontsize=7.5); ax.grid(True, which="both", linestyle="--", alpha=0.4)
ax.set_xlim(5, 45)

plt.tight_layout()
plt.savefig("Telecom_Constellation_FOMs.png", dpi=150, bbox_inches="tight")
plt.show()
print("\n  Plot saved to Telecom_Constellation_FOMs.png")
