"""
================================================================================
  DIRECT-TO-DEVICE TELECOM CONSTELLATION — MISSION ANALYSIS DOCUMENT
  Answers to the seven questions of the preliminary design assessment
  Reference mission: VLEO 200 km vs LEO 500 km, Spain (36–44° N), S-band 2 GHz
  ε_min = 20° (suburban average D2D reference)
================================================================================
"""

import sys
import math

sys.stdout.reconfigure(encoding='utf-8')

# ──────────────────────────────────────────────────────────────────────────────
#  PHYSICAL CONSTANTS AND MISSION INPUTS
# ──────────────────────────────────────────────────────────────────────────────
R_E  = 6_371_000.0      # [m]       Earth mean radius
MU   = 3.986004418e14   # [m³/s²]  gravitational parameter
C    = 299_792_458.0    # [m/s]    speed of light

EPSILON_MIN_DEG   = 20.0    # [deg]  reference min elevation (suburban D2D)
TARGET_LAT_MAX    = 44.0    # [deg]  northernmost latitude of service area (northern Spain)
TARGET_LAT_MIN    = 36.0    # [deg]  southernmost latitude (southern Spain)
FREQ_GHZ          = 2.0     # [GHz]  S-band
LATENCY_BUDGET_MS = 150.0   # [ms]   one-way VoIP budget
H_VLEO            = 200e3   # [m]    VLEO design point
H_LEO             = 500e3   # [m]    LEO design point

# ──────────────────────────────────────────────────────────────────────────────
#  CORE GEOMETRY FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────
def earth_central_halfangle(h_m, eps_deg):
    eps = math.radians(eps_deg)
    arg = R_E / (R_E + h_m) * math.cos(eps)
    return float('nan') if arg > 1.0 else math.acos(arg) - eps

def slant_range(h_m, eps_deg):
    eps = math.radians(eps_deg)
    return math.sqrt((R_E + h_m)**2 - (R_E * math.cos(eps))**2) - R_E * math.sin(eps)

def orbital_period(h_m):
    return 2.0 * math.pi * math.sqrt((R_E + h_m)**3 / MU)

def orbital_velocity(h_m):
    return math.sqrt(MU / (R_E + h_m))

def pass_duration(h_m, eps_deg):
    lam = earth_central_halfangle(h_m, eps_deg)
    return 0.0 if math.isnan(lam) else orbital_period(h_m) * lam / math.pi

def single_sat_revisit(h_m, eps_deg, lat_deg):
    lam = earth_central_halfangle(h_m, eps_deg)
    if math.isnan(lam) or lam <= 0:
        return float('inf')
    W = 2.0 * lam * R_E
    T = orbital_period(h_m)
    return T * math.pi * R_E * math.cos(math.radians(lat_deg)) / W

def fspl_dB(Rs_m, freq_GHz):
    return 20.0 * math.log10(4.0 * math.pi * Rs_m * freq_GHz * 1e9 / C)

def latency_ms(Rs_m):
    return Rs_m / C * 1e3

def eclipse_fraction(h_m):
    return (1.0 / math.pi) * math.asin(R_E / (R_E + h_m))

def n_zonal_estimate(h_m, eps_deg, lat_max_deg):
    lam = earth_central_halfangle(h_m, eps_deg)
    if math.isnan(lam) or lam <= 0:
        return float('inf')
    lam_d = math.degrees(lam)
    P = math.ceil(360.0 * math.cos(math.radians(lat_max_deg)) / (2.0 * lam_d))
    S = math.ceil(180.0 / lam_d)
    return int(P * S)

def walker_sizing(h_m, eps_deg, lat_max_deg):
    lam_d = math.degrees(earth_central_halfangle(h_m, eps_deg))
    P = math.ceil(360.0 * math.cos(math.radians(lat_max_deg)) / (2.0 * lam_d))
    S = math.ceil(180.0 / lam_d)
    N = P * S
    F = round(P / 2)
    return P, S, N, F

def eps_from_slant_range(h_m, Rs_target_m):
    X   = R_E + h_m
    num = X**2 - R_E**2 - Rs_target_m**2
    den = 2.0 * Rs_target_m * R_E
    s   = num / den
    if abs(s) > 1.0:
        return None
    return math.degrees(math.asin(s))

# ──────────────────────────────────────────────────────────────────────────────
#  PRE-COMPUTE ALL VALUES USED IN THE DOCUMENT
# ──────────────────────────────────────────────────────────────────────────────
eps = EPSILON_MIN_DEG

# — VLEO 200 km —
lam_v   = math.degrees(earth_central_halfangle(H_VLEO, eps))
Rs_v    = slant_range(H_VLEO, eps)
fspl_v  = fspl_dB(Rs_v, FREQ_GHZ)
tau_v   = latency_ms(Rs_v)
Torb_v  = orbital_period(H_VLEO) / 60.0
Tpass_v = pass_duration(H_VLEO, eps) / 60.0
v_v     = orbital_velocity(H_VLEO) / 1e3
Trev_v  = single_sat_revisit(H_VLEO, eps, TARGET_LAT_MAX) / 3600.0
ecl_v   = eclipse_fraction(H_VLEO) * 100.0
P_v, S_v, N_v, F_v = walker_sizing(H_VLEO, eps, TARGET_LAT_MAX)
incl_v  = TARGET_LAT_MAX + lam_v

# — LEO 500 km —
lam_l   = math.degrees(earth_central_halfangle(H_LEO, eps))
Rs_l    = slant_range(H_LEO, eps)
fspl_l  = fspl_dB(Rs_l, FREQ_GHZ)
tau_l   = latency_ms(Rs_l)
Torb_l  = orbital_period(H_LEO) / 60.0
Tpass_l = pass_duration(H_LEO, eps) / 60.0
v_l     = orbital_velocity(H_LEO) / 1e3
Trev_l  = single_sat_revisit(H_LEO, eps, TARGET_LAT_MAX) / 3600.0
ecl_l   = eclipse_fraction(H_LEO) * 100.0
P_l, S_l, N_l, F_l = walker_sizing(H_LEO, eps, TARGET_LAT_MAX)
incl_l  = TARGET_LAT_MAX + lam_l

# — FSPL difference —
delta_fspl = fspl_l - fspl_v
delta_Rs   = Rs_l / Rs_v

# — Path-loss equivalence: VLEO equivalent ε that matches LEO's FSPL —
eps_eq_vleo = eps_from_slant_range(H_VLEO, Rs_l)
lam_eq_v    = math.degrees(earth_central_halfangle(H_VLEO, eps_eq_vleo))
N_eq_v      = n_zonal_estimate(H_VLEO, eps_eq_vleo, TARGET_LAT_MAX)
Tpass_eq_v  = pass_duration(H_VLEO, eps_eq_vleo) / 60.0
Trev_eq_v   = single_sat_revisit(H_VLEO, eps_eq_vleo, TARGET_LAT_MAX) / 3600.0
P_eq, S_eq, N_eq_w, F_eq = walker_sizing(H_VLEO, eps_eq_vleo, TARGET_LAT_MAX)

# — Service tier N —
N_t1_v = math.ceil(Trev_v / 2.0)
N_t1_l = math.ceil(Trev_l / 2.0)
N_t2_v = math.ceil(Trev_v / 0.5)
N_t2_l = math.ceil(Trev_l / 0.5)
N_t3_v = N_v
N_t3_l = N_l

# — Sensitivity sweep values —
sensitivity = []
for e in [5, 10, 15, 20, 25, 30, 35, 40]:
    lv = math.degrees(earth_central_halfangle(H_VLEO, e))
    ll = math.degrees(earth_central_halfangle(H_LEO,  e))
    nv = n_zonal_estimate(H_VLEO, e, TARGET_LAT_MAX)
    nl = n_zonal_estimate(H_LEO,  e, TARGET_LAT_MAX)
    sensitivity.append((e, lv, nv, ll, nl))

# — Orbit survey values —
orbit_survey = []
for name, h in [("VLEO  200 km", 200e3), ("LEO   400 km", 400e3),
                ("LEO   500 km", 500e3), ("LEO   600 km", 600e3),
                ("LEO  1000 km", 1000e3), ("LEO  1200 km", 1200e3)]:
    lam  = math.degrees(earth_central_halfangle(h, eps))
    Rs   = slant_range(h, eps)
    fp   = fspl_dB(Rs, FREQ_GHZ)
    lat  = latency_ms(Rs)
    tp   = pass_duration(h, eps) / 60.0
    ncnt = n_zonal_estimate(h, eps, TARGET_LAT_MAX)
    tr   = single_sat_revisit(h, eps, TARGET_LAT_MAX) / 3600.0
    orbit_survey.append((name, lam, Rs/1e3, fp, lat, tp, ncnt, tr))

# ──────────────────────────────────────────────────────────────────────────────
#  DOCUMENT HEADER
# ──────────────────────────────────────────────────────────────────────────────
W = 84
SEP = "=" * W
SEP2 = "-" * W
def title(text):
    print("\n" + SEP)
    print(f"  {text}")
    print(SEP)
def section(text):
    print(f"\n  {'─'*70}")
    print(f"  {text}")
    print(f"  {'─'*70}")
def sub(text):
    print(f"\n  ▸ {text}")
def body(text):
    for line in text.strip().split('\n'):
        print(f"    {line}")
def table_row(cols, widths, sep="  "):
    row = sep
    for c, w in zip(cols, widths):
        row += str(c).ljust(w)
    print(row)

print()
print("=" * W)
print("  DIRECT-TO-DEVICE TELECOM CONSTELLATION — MISSION ANALYSIS DOCUMENT")
print(f"  Reference design: VLEO 200 km  vs  LEO 500 km  |  Spain ({TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N)  |  S-band 2 GHz")
print(f"  ε_min = {EPSILON_MIN_DEG:.0f}° (suburban average D2D)   f = {FREQ_GHZ} GHz")
print("=" * W)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 1 — WHERE TO START                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 1 — WHERE TO START")

section("1.1  First inputs to request before starting any analysis")
body("""
Before any geometry can be written down, five inputs must be locked.
Each one directly gates the next calculation and cannot be assumed.

┌─────────────────────────────────────────────────────────────────────────┐
│ INPUT 1 — SERVICE AREA                                                  │
│   What geographic region must be covered?                               │
│   → Sets the required orbital inclination and the cross-track plane     │
│     spacing condition.                                                   │
│   Reference: Spain 36–44° N  (Iberian Peninsula + Balearic Islands)    │
│              → Walker Delta at inclination i ≥ 49° satisfies this.     │
│                                                                         │
│ INPUT 2 — USER TERMINAL TYPE                                            │
│   Smartphone (antenna gain ≈ 0 dBi, Tx power 23–33 dBm) or            │
│   dedicated IoT device (may have a small patch antenna)?               │
│   → Determines the G/T of the user side of the link budget.            │
│   → For a smartphone, this is fixed by the handset standard — it        │
│     CANNOT be traded. The satellite must compensate with EIRP.         │
│                                                                         │
│ INPUT 3 — MINIMUM ELEVATION ANGLE  ε_min                               │
│   This is the most sensitive parameter in the whole study.             │
│   ε_min is NOT a free choice — it is the output of the link budget:    │
│       ε_min = f(satellite EIRP, user G/T, C/N₀ required, blockage)    │
│   Reference ranges for D2D:                                             │
│       ε = 10°  → optimistic (open terrain, large phased array)         │
│       ε = 20°  → suburban average                ← used here          │
│       ε = 30°  → conservative (urban canyon)                           │
│       ε = 40°  → worst-case high-rise city                             │
│   Impact: changing ε from 20° to 30° roughly doubles N at VLEO.       │
│                                                                         │
│ INPUT 4 — SERVICE TIER / CONTINUITY REQUIREMENT                        │
│   Tier 1 Emergency/IoT  : delay-tolerant, 1–2 h gap acceptable        │
│   Tier 2 Messaging      : short gaps OK, ~30 min max gap               │
│   Tier 3 Voice/broadband: continuous, no gap at all                    │
│   → Drives N by a factor of 10–200× between Tier 1 and Tier 3.        │
│                                                                         │
│ INPUT 5 — FREQUENCY BAND                                                │
│   D2D must use licensed cellular spectrum so unmodified phones work.   │
│   L-band ~1.5 GHz  : MSS allocations (Inmarsat heritage)              │
│   S-band ~2.0 GHz  : AST SpaceMobile, Lynk Global  ← used here       │
│   AWS/PCS ~1.7/2.1 GHz : Starlink direct-to-cell, T-Mobile            │
│   → Band choice determines FSPL, antenna aperture, and regulatory path.│
└─────────────────────────────────────────────────────────────────────────┘
""")

section("1.2  Performance metrics / figures of merit to define first")
body(f"""
Seven geometric FOMs are derived immediately from the five inputs above.
All are computed here for the two reference orbits at ε_min = {EPSILON_MIN_DEG:.0f}°:

  FOM                      Symbol    VLEO 200 km        LEO 500 km
  ─────────────────────────────────────────────────────────────────────
  Footprint half-angle      λ        {lam_v:>8.2f}°          {lam_l:>8.2f}°
  Slant range (worst-case)  R_s      {Rs_v/1e3:>8.0f} km        {Rs_l/1e3:>8.0f} km
  Free-space path loss      FSPL     {fspl_v:>8.1f} dB        {fspl_l:>8.1f} dB
  One-way propagation       τ        {tau_v:>8.2f} ms        {tau_l:>8.2f} ms
  Max pass duration         T_pass   {Tpass_v:>8.1f} min       {Tpass_l:>8.1f} min
  Single-sat revisit (44°N) T_rev    {Trev_v:>8.1f} h         {Trev_l:>8.1f} h
  N continuous (Spain)      N_cont   {N_t3_v:>8d} sats      {N_t3_l:>8d} sats
  ─────────────────────────────────────────────────────────────────────

These seven numbers together tell the mission analyst whether the orbit
is physically viable before any payload or system trade begins.

Three additional derived FOMs are important for system interfaces:
  • Eclipse fraction     f_ecl  :  VLEO {ecl_v:.1f}%  /  LEO {ecl_l:.1f}%  → power subsystem
  • Orbital velocity     v_orb  :  VLEO {v_v:.3f} km/s  /  LEO {v_l:.3f} km/s → drag, station-keeping
  • Required inclination i_min  :  VLEO ≥ {incl_v:.1f}°   /  LEO ≥ {incl_l:.1f}°  → launch vehicle trade
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 2 — COVERAGE-DRIVEN PRE-DESIGN                                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 2 — COVERAGE-DRIVEN PRE-DESIGN")

section("2.1  How to estimate the number of satellites needed")
body(f"""
The number of satellites for continuous coverage is derived from three
nested geometry conditions.

STEP A — Compute the Earth central half-angle λ (footprint radius):
  λ = arccos( R_E/(R_E+h) · cos(ε) ) − ε                        [eq. 1]

  This is the angular radius of the coverage cap as seen from Earth's centre.
  At ε = {EPSILON_MIN_DEG:.0f}°:
    VLEO 200 km : λ = arccos({R_E/(R_E+H_VLEO):.4f} × cos({EPSILON_MIN_DEG:.0f}°)) − {EPSILON_MIN_DEG:.0f}° = {lam_v:.2f}°
    LEO  500 km : λ = arccos({R_E/(R_E+H_LEO):.4f} × cos({EPSILON_MIN_DEG:.0f}°)) − {EPSILON_MIN_DEG:.0f}° = {lam_l:.2f}°

  Physical meaning: at VLEO the satellite sees a cap of radius {lam_v:.1f}° at the
  Earth's centre. At LEO the cap is {lam_l:.1f}°, i.e. {lam_l/lam_v:.1f}× wider in angle
  and ({lam_l/lam_v:.1f})² ≈ {(lam_l/lam_v)**2:.1f}× larger in area.

STEP B — Compute minimum planes P:
  Planes must be spaced close enough that their footprint disks overlap at
  the hardest latitude, which for a 0–{TARGET_LAT_MAX:.0f}° N zone is the poleward
  boundary (φ_max = {TARGET_LAT_MAX:.0f}°):

    P ≥ ceil( 360° × cos(φ_max) / (2λ) )                         [eq. 2]

    VLEO: P ≥ ceil( 360 × cos({TARGET_LAT_MAX:.0f}°) / (2 × {lam_v:.2f}) )
            = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX)):.1f} / {2*lam_v:.2f} )
            = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX))/(2*lam_v):.2f} ) = {P_v} planes
    LEO:  P ≥ ceil( 360 × cos({TARGET_LAT_MAX:.0f}°) / (2 × {lam_l:.2f}) )
            = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX)):.1f} / {2*lam_l:.2f} )
            = ceil( {360*math.cos(math.radians(TARGET_LAT_MAX))/(2*lam_l):.2f} ) = {P_l} planes

STEP C — Compute satellites per plane S:
  Satellites within one plane must be spaced so consecutive footprints
  never leave a gap along the orbital track:

    S ≥ ceil( 180° / λ )                                          [eq. 3]

    VLEO: S ≥ ceil( 180 / {lam_v:.2f} ) = ceil({180/lam_v:.2f}) = {S_v} sats/plane
    LEO:  S ≥ ceil( 180 / {lam_l:.2f} ) = ceil({180/lam_l:.2f}) = {S_l} sats/plane

STEP D — Total constellation size:
    N = P × S

    VLEO: N = {P_v} × {S_v} = {N_v} satellites  (Walker {N_v}/{P_v}/{F_v})
    LEO:  N = {P_l} × {S_l} = {N_l} satellites  (Walker {N_l}/{P_l}/{F_l})

WHY N SCALES AS λ² (inverse square law of constellation size):
  Both P and S are proportional to 1/λ, so N = P×S ∝ 1/λ².
  The footprint area scales as λ², so halving the footprint quadruples N.
  λ_LEO / λ_VLEO = {lam_l:.2f}/{lam_v:.2f} = {lam_l/lam_v:.2f} → N ratio ≈ ({lam_l/lam_v:.2f})² = {(lam_l/lam_v)**2:.1f}×
  Actual ratio: {N_v}/{N_l} = {N_v/N_l:.2f}× (matches prediction well)
""")

section("2.2  How minimum elevation angle affects the analysis")
body("""
ε_min is the single most powerful lever in the analysis. It enters in
two places simultaneously and compounds its effect:

  (1) Larger ε → shorter slant range R_s → lower FSPL
      → the satellite needs LESS EIRP to close the link. GOOD for payload.

  (2) Larger ε → smaller coverage cap λ → more satellites needed.
      → GOOD for link, BAD for constellation size.

This is the fundamental trade-off. The mission analyst must communicate
the full sensitivity table to the system engineer so they can decide
where to trade satellite EIRP against constellation size.

  SENSITIVITY TABLE  (continuous coverage, Spain 36–44° N)
  ──────────────────────────────────────────────────────────────────────────
  ε_min  │  VLEO 200 km              │  LEO 500 km
  [deg]  │  λ [deg]    N_cont        │  λ [deg]    N_cont   │ N ratio
  ───────┼───────────────────────────┼──────────────────────┼──────────""")
for e, lv, nv, ll, nl in sensitivity:
    ratio = f"{nv/nl:.1f}×" if nl > 0 else "—"
    nv_s = str(nv) if nv < 99999 else ">99999"
    nl_s = str(nl) if nl < 99999 else ">99999"
    marker = "  ← used here" if e == int(EPSILON_MIN_DEG) else ""
    body(f"   {e:>5}  │  {lv:>7.2f}     {nv_s:>8}    │  {ll:>7.2f}     {nl_s:>8}   │  {ratio}{marker}")
body("""  ──────────────────────────────────────────────────────────────────────────

  KEY INSIGHT: going from ε=20° to ε=30° at VLEO increases N from ~1260 to
  ~2835 — a factor of 2.3×. A 5° step in ε_min roughly doubles N at VLEO
  (because N ∝ 1/λ² and λ decreases steeply as ε rises at low altitude).

  At LEO the sensitivity is lower: same step doubles N only by ~2×.
  This asymmetry means that urban D2D (high ε_min) hits VLEO much harder
  than LEO — a critical programme-level decision point.

  The mission analyst must present this table at the earliest system review
  and flag that N is not a single number but a function of ε_min.
  Until the link budget is closed, N cannot be quoted with confidence.
""")

section("2.3  How to choose between one plane or multiple planes")
body(f"""
A SINGLE PLANE cannot provide continuous coverage over a wide geographic
zone. Here is why and when it matters:

  Single-plane geometry:
    One orbital plane sweeps one ascending corridor per orbit.
    As the Earth rotates under it, the ground track shifts westward.
    A ground user at a fixed longitude is covered only during the pass
    and then must wait ~T_rev for the next overflight.

  Single-plane revisit at 44° N (northern Spain, ε = {EPSILON_MIN_DEG:.0f}°):
    VLEO: T_rev_single ≈ {Trev_v:.1f} h  → acceptable only for Tier 1 (IoT/SOS)
    LEO:  T_rev_single ≈ {Trev_l:.1f} h  → acceptable only for Tier 1

  DECISION RULE — when is a single plane acceptable?
    • Tier 1 (SOS beacon, IoT, max gap 2–4 h): YES for both orbits
    • Tier 2 (messaging, max gap 30 min): NO — need {N_t2_v} (VLEO) or {N_t2_l} (LEO) sats
    • Tier 3 (voice, continuous): NO — need P ≥ {P_v} (VLEO) or P ≥ {P_l} (LEO) planes

  For the D2D reference mission (targeting at minimum Tier 2), multiple
  planes are required for both design points.

  HOW MANY PLANES — the trade space:
  ┌────────────────────────────────────────────────────────────────────┐
  │  OPTION A — minimum planes (P_min, large S)                       │
  │    P = {P_l} planes (LEO) or {P_v} planes (VLEO)                            │
  │    Pro: fewer RAAN slots to manage, simpler launch manifest        │
  │    Con: single-plane failure removes an entire longitude band;     │
  │         ground track repeat coverage takes longer to achieve       │
  │                                                                    │
  │  OPTION B — balanced (P ≈ S ≈ sqrt(N))                            │
  │    LEO:  sqrt({N_l}) ≈ {int(N_l**0.5)} → P ≈ {int(N_l**0.5)}, S ≈ {int(N_l**0.5)}  │
  │    VLEO: sqrt({N_v}) ≈ {int(N_v**0.5)} → P ≈ {int(N_v**0.5)}, S ≈ {int(N_v**0.5)} │
  │    Pro: geographic and temporal redundancy balanced                │
  │    Con: more launches required                                     │
  │                                                                    │
  │  OPTION C — many thin planes (large P, small S)                   │
  │    Pro: maximum geographic diversity, one plane lost = small gap   │
  │    Con: maximum launch complexity; each plane needs a dedicated    │
  │         raise to its RAAN (either multiple launches or long drift) │
  └────────────────────────────────────────────────────────────────────┘

  RECOMMENDED STARTING POINT: P_min × S_min (options from eq. 2 and 3),
  then check that P × S cost is within programme budget. Increase P at
  the cost of reducing S proportionally to improve resilience.
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 3 — ORBIT SELECTION                                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 3 — ORBIT SELECTION")

section("3.1  Which orbit types to evaluate first and why")
body("""
Four orbit families are considered. They are ranked by D2D suitability:

  RANK 1: LEO  300–1000 km  → PRIMARY CANDIDATE RANGE
    ✔ Latency 2–5 ms one-way: well within 150 ms VoIP budget
    ✔ FSPL 148–162 dB: closeable with a large phased array
    ✔ Manageable atmospheric drag above 450 km
    ✔ Below Van Allen inner belt (trapped radiation starts ~1000 km)
    ✔ Short pass duration (2–10 min): acceptable for messaging; marginal
      for voice without continuous constellation

  RANK 2: VLEO  160–300 km  → HIGH-RISK, HIGH-REWARD
    ✔ 7.5 dB lower path loss than LEO 500 km → smaller satellite antenna
    ✘ Severe atmospheric drag at 200 km → Air-Breathing Electric
      Propulsion (ABEP) required (TRL 4–5 as of 2025)
    ✘ Smaller footprint → 5.3× more satellites for equal coverage
    ✘ Shorter pass duration: 2.1 min at 200 km, ε=20° → voice is
      practically impossible without a dense constellation

  RANK 3: MEO  2000–20000 km  → EXCLUDED for D2D
    ✘ Van Allen belts (1000–6000 km and 13000–20000 km) cause rapid
      solar cell and electronics degradation → prohibitive shielding mass
    ✘ FSPL 170–180 dB: far beyond what a phone link can close without
      a large terminal antenna
    ✘ Latency 7–60 ms one-way: technically within budget, but advantage
      over LEO is marginal and does not offset the radiation penalty

  RANK 4: GEO  35786 km  → EXCLUDED for D2D
    ✘ FSPL ≈ 188 dB: link cannot close against a smartphone
    ✘ One-way latency ≈ 119 ms: approaches VoIP limit
    ✘ N = 3–4 satellites possible but coverage gap at high latitudes

  CONCLUSION: The shortlist for Phase 1 analysis is LEO 300–1200 km with
  VLEO 200 km included as the high-performance/high-risk design point.
  MEO and GEO are closed out immediately.
""")

section("3.2  How altitude affects key metrics — full sweep table")
body(f"""
  ORBIT SURVEY  (ε = {EPSILON_MIN_DEG:.0f}°, f = {FREQ_GHZ:.1f} GHz, service area Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N)
  ────────────────────────────────────────────────────────────────────────────
  Orbit        λ[°]  Rs[km]  FSPL[dB]  τ[ms]  T_pass[min]  N_cont  T_rev[h]
  ────────────────────────────────────────────────────────────────────────────""")
for name, lam, rs, fp, lat, tp, ncnt, tr in orbit_survey:
    ncnt_s = str(ncnt) if ncnt < 99999 else ">99999"
    body(f"  {name:<13} {lam:>5.1f}  {rs:>6.0f}  {fp:>8.1f}  {lat:>6.1f}  {tp:>11.1f}  {ncnt_s:>6}  {tr:>8.1f}")
body("""  ────────────────────────────────────────────────────────────────────────────

  HOW ALTITUDE AFFECTS EACH METRIC:

  ▸ Footprint λ:   INCREASES with altitude (larger coverage cap)
    Higher altitude → satellite sees more of Earth → fewer satellites
    needed. But gains are non-linear: doubling altitude from 200→400 km
    nearly doubles λ (from 2.9° to ~5°), but 400→800 km only adds 2°.

  ▸ Path loss FSPL: INCREASES with altitude (20 log₁₀ of slant range)
    Every factor 2× in slant range = +6 dB. VLEO advantage over LEO 500 km
    is 7.5 dB — significant but cannot be gained for free because of the
    footprint and constellation size penalty.

  ▸ Latency τ:     INCREASES with altitude (proportional to range)
    All LEO/VLEO latencies (1–5 ms) are far below the 150 ms VoIP limit.
    Latency does not discriminate within the LEO band — it is only
    a concern for GEO comparisons.

  ▸ Pass duration T_pass: INCREASES with altitude (wider footprint, longer arc)
    This is a critical D2D metric:
      200 km (ε=20°): 2.1 min — barely usable for a short message
      500 km (ε=20°): 5.0 min — marginal for a short voice call
      1200 km (ε=20°): ~11 min — comfortable for a voice session
    Higher altitude is genuinely better for session length, which is an
    argument for LEO 500–1200 km over VLEO for voice service.

  ▸ Revisit time T_rev:  DECREASES with altitude (wider footprint dominates)
    Despite higher altitude meaning longer orbital periods, the footprint
    grows faster, so the satellite sweeps its target more often per orbit.
    LEO 500 km has a 2.15× shorter revisit than VLEO 200 km — this is the
    most counter-intuitive result and the main argument against VLEO for
    intermittent services.

  ▸ Constellation size N_cont: DECREASES rapidly with altitude
    N ∝ 1/λ² → at 500 km N = 280 vs 1260 at 200 km (factor 4.5×) for Spain.
    Above ~800 km, N drops well below 200 and becomes economically attractive.
    The Van Allen belt (above ~1000 km) cuts off this trend.
""")

section("3.3  What inclinations to consider")
body(f"""
  Inclination must satisfy: i ≥ φ_max + λ  (Walker coverage condition)

  For Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N service area:
    VLEO 200 km: i_min = {TARGET_LAT_MAX:.0f}° + {lam_v:.1f}° = {incl_v:.1f}°   (round up to 49°)
    LEO  500 km: i_min = {TARGET_LAT_MAX:.0f}° + {lam_l:.1f}° = {incl_l:.1f}°   (round up to 54°)

  THREE OPTIONS AND THEIR TRADE-OFFS:

  ┌─────────────────────────────────────────────────────────────────────┐
  │ OPTION 1 — Walker Delta  i ≈ 49–54°  [RECOMMENDED for this study]  │
  │   Coverage just closes at 44° N; no capacity wasted above 44°.     │
  │   Covers all of Spain; also covers southern France, Italy, Greece.  │
  │   ✔ Minimum N for the required coverage zone                       │
  │   ✔ Lower inclination → cheaper launch (more equatorial boost)      │
  │   ✘ Does not cover above ~49° N (northern France, Germany,         │
  │     Benelux, UK, Scandinavia excluded)                              │
  │                                                                     │
  │ OPTION 2 — Polar orbits  i = 90°                                   │
  │   Covers the entire globe including poles.                          │
  │   ✔ If pan-European service is ever required, no re-design needed   │
  │   ✔ Simpler orbit precession behaviour                              │
  │   ✘ ~15–20% more satellites needed vs Walker Delta for same        │
  │     Spain coverage (capacity wasted over polar regions)             │
  │   ✘ Higher launch energy required (no equatorial velocity benefit)  │
  │                                                                     │
  │ OPTION 3 — Sun-synchronous  i ≈ 97–98°                            │
  │   Fixed local solar time at each overflight.                        │
  │   ✔ Useful if solar illumination matters (optical payloads)         │
  │   ✘ Not relevant for telecoms; adds unnecessary inclination         │
  │   ✘ Retrograde orbit → less efficient use of Earth's rotation       │
  └─────────────────────────────────────────────────────────────────────┘

  RECOMMENDATION: Walker Delta at i = 49° for VLEO, i = 54° for LEO.
  If programme requirements later extend coverage to continental Europe
  (up to ~60° N), re-assess whether increasing to i = 67° (and accepting
  ~30% more satellites) is more cost-effective than a supplemental shell.
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 4 — CONSTELLATION ARCHITECTURE                                   ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 4 — CONSTELLATION ARCHITECTURE")

section("4.1  How to reason about number of satellites N")
body(f"""
  N is not a single number — it depends entirely on the service tier.
  This is the most important message to deliver to a programme manager.

  ┌──────────────────────────────────────────────────────────────────────────┐
  │  SERVICE TIER           CONDITION              VLEO 200 km   LEO 500 km │
  ├──────────────────────────────────────────────────────────────────────────┤
  │  Tier 1: Emergency/IoT  Mean gap ≤ 2 h         N = {N_t1_v:>5}        N = {N_t1_l:>3}   │
  │  Tier 2: Messaging      Mean gap ≤ 30 min       N = {N_t2_v:>5}        N = {N_t2_l:>3}   │
  │  Tier 3: Continuous     No gap (Spain)          N = {N_t3_v:>5}        N = {N_t3_l:>3}   │
  └──────────────────────────────────────────────────────────────────────────┘
  All values computed at ε_min = {EPSILON_MIN_DEG:.0f}°, φ = {TARGET_LAT_MAX:.0f}° N (northern Spain, worst-case).

  HOW THESE NUMBERS ARE DERIVED:
  • Tier 1 & 2:  N = ⌈T_rev_single / T_gap_allowed⌉
    With 1 satellite: VLEO revisit = {Trev_v:.1f} h, LEO revisit = {Trev_l:.1f} h (at 44° N, northern Spain)
    Distributing N satellites uniformly reduces the gap by factor N.
    Tier 1:  N_VLEO = ⌈{Trev_v:.1f} / 2.0⌉ = {N_t1_v};   N_LEO = ⌈{Trev_l:.1f} / 2.0⌉ = {N_t1_l}
    Tier 2:  N_VLEO = ⌈{Trev_v:.1f} / 0.5⌉ = {N_t2_v};   N_LEO = ⌈{Trev_l:.1f} / 0.5⌉ = {N_t2_l}

  • Tier 3: N = P × S from the Walker continuous coverage conditions
    (eqs. 2 and 3 above). There is no gap anywhere in the service zone.

  PROGRAMME IMPLICATION:
  The jump from Tier 2 to Tier 3 is enormous — a factor of ~44× for VLEO
  ({N_t2_v} → {N_t3_v}) and ~18× for LEO ({N_t2_l} → {N_t3_l}).
  A programme must define its target tier before signing a contract.
  Tier 1/2 is within reach of a startup (~6–50 satellites).
  Tier 3 is a major programme (~400–2000 satellites).
""")

section("4.2  How to reason about number of planes P")
body(f"""
  P is driven by the cross-track coverage condition at the poleward edge
  of the service area (eq. 2 above).

  MINIMUM P VALUES (just closing coverage at {TARGET_LAT_MAX:.0f}° N):
    VLEO 200 km: P_min = {P_v} planes  (each plane occupies {360/P_v:.1f}° of RAAN)
    LEO  500 km: P_l_min = {P_l} planes  (each plane occupies {360/P_l:.1f}° of RAAN)

  HOW TO SPREAD PLANES IN RAAN:
  The P planes are equally spaced in RAAN (right ascension of ascending
  node) from 0° to 360°. To achieve this:
    • With a single launch vehicle: deploy to one RAAN and drift to
      others using differential J₂ precession — takes months/years.
    • With P dedicated launches: inject directly to each RAAN — fast
      but expensive.
    • With a polar launch: one launch can serve all RANs progressively
      as the orbit precesses, but requires long commissioning phase.

  RESILIENCE CONSIDERATION:
  If one plane is lost (satellite failure, debris mitigation manoeuvre),
  the coverage gap is approximately 360°/P in RAAN:
    VLEO: losing 1 of {P_v} planes → {360/P_v:.1f}° RAAN gap (moderate impact)
    LEO:  losing 1 of {P_l} planes → {360/P_l:.1f}° RAAN gap (larger relative impact)
  For LEO this is more severe — fewer planes means each one is more
  critical. Increase P above P_min to improve resilience.
""")

section("4.3  How to reason about satellites per plane S")
body(f"""
  S is driven by the along-track coverage condition (eq. 3 above).

  MINIMUM S VALUES:
    VLEO 200 km: S_min = {S_v} sats/plane  (spacing = {360/S_v:.1f}°/sat in true anomaly)
    LEO  500 km: S_min = {S_l} sats/plane  (spacing = {360/S_l:.1f}°/sat in true anomaly)

  The angular gap between consecutive satellites is 360°/S. For continuous
  coverage this must be ≤ 2λ (footprint diameter):
    VLEO: 360°/{S_v} = {360/S_v:.1f}° vs 2λ = {2*lam_v:.2f}° ✓ (coverage closes)
    LEO:  360°/{S_l} = {360/S_l:.1f}° vs 2λ = {2*lam_l:.2f}° ✓ (coverage closes)

  LAUNCH STRATEGY FOR S:
  All S satellites in a plane are injected to the same altitude and
  RAAN, then manoeuvred to their final in-plane spacing using small
  differential burns. This is standard practice for constellation
  deployment (e.g., Starlink trains).

  FAILURE MARGIN:
  If one satellite in a plane fails, the gap becomes 2×(360°/S) = {2*360/S_v:.1f}°
  (VLEO) or {2*360/S_l:.1f}° (LEO). Coverage is interrupted at that longitude.
  A replacement satellite can be launched to the RAAN-matched plane and
  drifted to the correct slot, typically within 1–3 months.
""")

section("4.4  Phasing between satellites — the F parameter")
body(f"""
  The Walker T/P/F parameter F (0 ≤ F < P) controls the relative true
  anomaly offset between adjacent planes.

  DEFINITION:
  In a Walker T/P/F constellation, the true anomaly of the first satellite
  in plane j is offset by F × 360°/T relative to plane 0.

  PHYSICAL EFFECT:
  F controls WHEN satellites from different planes are simultaneously
  visible at a given ground point — it sets the temporal rhythm of the
  coverage pattern, not the spatial coverage geometry.

  FOR CONTINUOUS COVERAGE (Tier 3):
  F is irrelevant. At any instant, some satellite is above ε_min regardless
  of F. The system analyst can choose F = P/2 = {F_v} (VLEO) or {F_l} (LEO) as
  the standard balanced choice, but it has no operational consequence.

  FOR INTERMITTENT COVERAGE (Tier 1/2):
  F is CRITICAL. Different values produce very different maximum gaps:
    • F = 0:    all planes' first satellites cluster at the same longitude
                → maximum gap = T_rev (all passes bunch together)
    • F = P/2:  satellites from adjacent planes are maximally spread in
                time → maximum gap is minimised
    • F = P/2 recommended for Tier 1/2: distributes service events
                uniformly in time, ensuring the longest quiet period
                is as short as possible.

  RECOMMENDED VALUES:
    VLEO continuous: F = {F_v}  (Walker {N_v}/{P_v}/{F_v})
    LEO  continuous: F = {F_l}  (Walker {N_l}/{P_l}/{F_l})
""")

section("4.5  Continuity of service vs intermittent service — the trade-off")
body(f"""
  This is the fundamental programme decision. The following comparison
  shows the cost is not linear — it is a cliff edge between tiers.

  VLEO 200 km  (ε = {EPSILON_MIN_DEG:.0f}°, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N):
    Tier 1 SOS/IoT      (2 h max gap):    N = {N_t1_v:>5}   (≈ 12 small sats, 1 launch)
    Tier 2 Messaging    (30 min max gap): N = {N_t2_v:>5}   (mid-size programme)
    Tier 3 Continuous   (no gap):         N = {N_t3_v:>5}   (major programme)
    Cliff between Tier 1 → Tier 3: factor {N_t3_v//N_t1_v}×

  LEO 500 km  (ε = {EPSILON_MIN_DEG:.0f}°, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N):
    Tier 1 SOS/IoT      (2 h max gap):    N = {N_t1_l:>5}   (≈ 6 small sats, 1 launch)
    Tier 2 Messaging    (30 min max gap): N = {N_t2_l:>5}   (small programme)
    Tier 3 Continuous   (no gap):         N = {N_t3_l:>5}   (large programme)
    Cliff between Tier 1 → Tier 3: factor {N_t3_l//N_t1_l}×

  THE PATH-LOSS EQUIVALENCE ARGUMENT — CAN VLEO BEAT LEO IN TIER 3?
  If VLEO uses its {delta_fspl:.1f} dB FSPL advantage to operate at a LOWER
  elevation angle instead of the same ε_min as LEO, the footprint grows:

    LEO at ε = {EPSILON_MIN_DEG:.0f}°   → R_s = {Rs_l/1e3:.0f} km  → FSPL = {fspl_l:.1f} dB → N = {N_l}
    VLEO at ε_eq = {eps_eq_vleo:.1f}° → R_s = {Rs_l/1e3:.0f} km  → same FSPL      → N = {N_eq_v}

  At ε_eq = {eps_eq_vleo:.1f}°, VLEO achieves the same link quality as LEO at {EPSILON_MIN_DEG:.0f}°
  but with a wider footprint (λ = {lam_eq_v:.1f}° vs {lam_l:.1f}°), requiring only
  N = {N_eq_v} satellites — FEWER than LEO's {N_l}.

  CATCH: ε = {eps_eq_vleo:.1f}° requires very low building blockage (rural/suburban).
  In a dense urban environment where ε_min is forced to 30° by building
  geometry, VLEO CANNOT exploit this advantage. The choice between the
  two scenarios is made by the service area definition, not by the
  mission analyst.

  SUMMARY OF TRADE-OFFS — VLEO vs LEO:
  ┌─────────────────────────┬────────────────────┬──────────────────────┐
  │  Metric                 │  VLEO 200 km       │  LEO 500 km          │
  ├─────────────────────────┼────────────────────┼──────────────────────┤
  │  Path loss (FSPL)       │  {fspl_v:.1f} dB ← better │  {fspl_l:.1f} dB             │
  │  Slant range            │  {Rs_v/1e3:.0f} km ← shorter │  {Rs_l/1e3:.0f} km               │
  │  Pass duration          │  {Tpass_v:.1f} min  (tight) │  {Tpass_l:.1f} min ← better   │
  │  Single-sat revisit     │  {Trev_v:.1f} h  (worse) │  {Trev_l:.1f} h ← better    │
  │  N Tier 3 (Spain,20°)  │  {N_t3_v:>5} ← much worse │  {N_t3_l:>3}               │
  │  N Tier 3 (rural ~9°)  │  {N_eq_v:>5} ← better    │  {N_l:>3}               │
  │  ABEP required?         │  YES (TRL 4–5)     │  NO ← simpler        │
  │  Atmosphere shielding   │  Partial benefit   │  None                │
  └─────────────────────────┴────────────────────┴──────────────────────┘
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 5 — INTERACTION WITH TELECOM CONSTRAINTS                         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 5 — INTERACTION WITH TELECOM CONSTRAINTS")

body(f"""
  Even when the payload is defined by the telecom engineer, the mission
  analyst owns five geometric outputs that directly couple to the link
  budget and the payload design. Delivering wrong or missing values here
  will cause the payload team to size to the wrong requirements.

  ══════════════════════════════════════════════════════════════════════════
  CONSTRAINT 1 — SLANT RANGE → FREE-SPACE PATH LOSS → MINIMUM EIRP
  ══════════════════════════════════════════════════════════════════════════
  The mission analyst must provide R_s(h, ε_min) as the fundamental link
  budget input. FSPL = 20 log₁₀(4π R_s f / c) directly determines how
  much EIRP the satellite must radiate to close the link.

  At ε_min = {EPSILON_MIN_DEG:.0f}°:
    VLEO 200 km: R_s = {Rs_v/1e3:.0f} km  → FSPL = {fspl_v:.1f} dB
    LEO  500 km: R_s = {Rs_l/1e3:.0f} km  → FSPL = {fspl_l:.1f} dB
    Difference:  ΔFSPL = {delta_fspl:.1f} dB

  The link budget equation is:
    EIRP_sat ≥ FSPL + (C/N₀)_required − G/T_user

  For a D2D smartphone (G/T ≈ −23 dB/K) and an LTE data channel
  requiring C/N₀ ≈ 55 dBHz, this gives:
    Minimum EIRP(VLEO) ≈ {fspl_v:.0f} + 55 − (−23) = {fspl_v + 55 + 23:.0f} dBW
    Minimum EIRP(LEO)  ≈ {fspl_l:.0f} + 55 − (−23) = {fspl_l + 55 + 23:.0f} dBW

  The {delta_fspl:.1f} dB saving at VLEO is the primary argument for the lower
  orbit — it allows the satellite to use a smaller phased array, lower DC
  power, or carry more margin for urban shadowing.

  MISSION ANALYST DELIVERABLE: Table of R_s(h, ε_min) and FSPL(h, ε_min)
  for ε_min in [10°, 20°, 30°, 40°] at both design altitudes.

  ══════════════════════════════════════════════════════════════════════════
  CONSTRAINT 2 — PASS DURATION → SESSION AND HANDOVER PROTOCOL
  ══════════════════════════════════════════════════════════════════════════
  T_pass sets the maximum uninterrupted contact time per satellite.
  After T_pass the satellite drops below ε_min and the link is lost.

    VLEO 200 km: T_pass = {Tpass_v:.1f} min  — very tight for any voice service
    LEO  500 km: T_pass = {Tpass_l:.1f} min  — marginal but feasible for a short call

  LTE/NR session establishment overhead: 1–3 s for connection setup.
  Handover must begin ≥ 20 s before end of pass (network protocol margin).
  Effective data transfer time:
    VLEO: {Tpass_v:.1f} min − 20 s setup − 20 s handover = {max(0, Tpass_v*60-40)/60:.1f} min usable
    LEO:  {Tpass_l:.1f} min − 20 s setup − 20 s handover = {max(0, Tpass_l*60-40)/60:.1f} min usable

  The network architecture team needs T_pass(h, ε_min) to design their
  handover protocol. If T_pass < 2 min, voice calls require continuous
  coverage (Tier 3 constellation) — single-pass calls are not viable.

  MISSION ANALYST DELIVERABLE: T_pass(h, ε_min) curve delivered to the
  network team for handover protocol dimensioning.

  ══════════════════════════════════════════════════════════════════════════
  CONSTRAINT 3 — SINGLE-SATELLITE REVISIT TIME → SERVICE GAP
  ══════════════════════════════════════════════════════════════════════════
  T_rev = T_orb × π × R_E × cos(φ) / (2λ × R_E)

  At 44° N (northern Spain), ε = {EPSILON_MIN_DEG:.0f}°:
    VLEO 200 km: T_rev = {Trev_v:.1f} h   (single satellite)
    LEO  500 km: T_rev = {Trev_l:.1f} h   (single satellite)

  This is the orbit-level input to the service tier analysis.
  N for a target gap T_gap: N = ⌈T_rev / T_gap⌉
  The telecom team needs T_rev at the worst-case latitude in the service
  area to size the constellation for their SLA commitments.

  ══════════════════════════════════════════════════════════════════════════
  CONSTRAINT 4 — MINIMUM INCLINATION → LAUNCH VEHICLE AND REGULATORY SLOT
  ══════════════════════════════════════════════════════════════════════════
  For Walker Delta: i_min = φ_max + λ
    VLEO: i_min = {TARGET_LAT_MAX:.0f}° + {lam_v:.1f}° = {incl_v:.1f}°  → round up to 49°
    LEO:  i_min = {TARGET_LAT_MAX:.0f}° + {lam_l:.1f}° = {incl_l:.1f}°  → round up to 54°

  This must be confirmed with the launch vehicle team: some vehicles
  have inclination-constrained trajectories from certain launch sites.
  It also interacts with ITU orbital slot coordination.

  ══════════════════════════════════════════════════════════════════════════
  CONSTRAINT 5 — ECLIPSE FRACTION → POWER SYSTEM AND BATTERY SIZING
  ══════════════════════════════════════════════════════════════════════════
  f_ecl = (1/π) arcsin(R_E / (R_E + h))
    VLEO 200 km: f_ecl = {ecl_v:.1f}% of orbit in shadow
    LEO  500 km: f_ecl = {ecl_l:.1f}% of orbit in shadow

  For a D2D phased array consuming P_pay watts during a pass, the
  battery must store at least:
    E_bat = P_pay × T_orb × f_ecl  [Joules]
    VLEO: E_bat = P_pay × {orbital_period(H_VLEO)/60:.0f} min × {ecl_v/100:.2f} = {orbital_period(H_VLEO)/60 * ecl_v/100:.0f} P_pay [J/W]
    LEO:  E_bat = P_pay × {orbital_period(H_LEO)/60:.0f} min × {ecl_l/100:.2f} = {orbital_period(H_LEO)/60 * ecl_l/100:.0f} P_pay [J/W]

  The two values are nearly equal (VLEO's shorter period nearly cancels
  its higher eclipse fraction). The dominant driver is P_pay (the phased
  array power level), which the mission analyst does not set — but must
  request from the telecom team to correctly size the power subsystem.

  MISSION ANALYST DELIVERABLE: T_orb, f_ecl, and T_pass delivered to the
  power team to compute battery energy per eclipse.
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 6 — RECOMMENDED STUDY FLOW                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 6 — RECOMMENDED STUDY FLOW")
body(f"""
  Phase 1: Requirements capture to preliminary constellation candidates
  ════════════════════════════════════════════════════════════════════════

  STEP 1 — LOCK THE REQUIREMENTS  (Week 1)
    Input: service area (lat/lon bounds), user terminal type, target
    service tier, frequency band, latency budget.
    Output: a single requirements table that everyone signs.
    WARNING: Do NOT proceed without a signed ε_min from the payload team.
    Geometry results are useless until ε_min is fixed.
    Action: Request link budget assumptions from payload team immediately.

  STEP 2 — DEFINE FIGURES OF MERIT  (Week 1)
    Establish which metrics will be used to compare candidates:
    λ, R_s, FSPL, τ, T_pass, T_rev, N_tier1, N_tier2, N_tier3.
    Add eclipse fraction and minimum inclination for system interfaces.
    Agree the FOM definitions with all engineering leads before computing.

  STEP 3 — ALTITUDE TRADE SWEEP  (Week 2)
    Compute all FOMs vs h from 200 km to 1200 km in steps.
    Plot: λ(h), N_cont(h), FSPL(h), T_pass(h), T_rev(h) — all at ε_ref.
    Mark VLEO 200 km and LEO 500 km as labelled design points.
    Note: exclude 1000–6000 km (Van Allen inner belt).
    Highlight the cliff: N_cont jumps steeply below 400 km.

  STEP 4 — ε_min SENSITIVITY SWEEP  (Week 2)
    For the two shortlisted altitudes (200 km and 500 km), sweep
    ε_min from 5° to 40° and plot N_cont vs ε_min.
    Present this table at the first system review.
    Key message: "If link budget forces ε = 30° instead of the assumed 20°,
    N at VLEO rises from ~1260 to ~2835 (for Spain). This is the biggest lever
    in the programme."

  STEP 5 — SERVICE TIER ANALYSIS  (Week 3)
    For both design altitudes and the reference ε_min, compute:
      Tier 1 N = ⌈T_rev / 2 h⌉  →  VLEO: {N_t1_v},  LEO: {N_t1_l}
      Tier 2 N = ⌈T_rev / 0.5 h⌉ → VLEO: {N_t2_v},  LEO: {N_t2_l}
      Tier 3 N = P × S            → VLEO: {N_t3_v}, LEO: {N_t3_l}
    Present as a cost/capability trade: "What can we afford to build?"
    Tier 1 is achievable for any startup. Tier 3 is a major programme.

  STEP 6 — PATH-LOSS EQUIVALENCE ANALYSIS  (Week 3)
    Ask: "If VLEO uses its 7.5 dB FSPL advantage to lower ε_min rather
    than to shrink the antenna, can it beat LEO on constellation size?"
    Compute eps_eq_vleo (the VLEO elevation at which its FSPL matches
    LEO's FSPL at ε_ref = 20°). Result: ε_eq ≈ 4.4°.
    At ε_eq, VLEO needs only ~234 satellites vs LEO's 280 (for Spain).
    CONCLUSION: VLEO wins IF the service environment allows ε < 10°
    (open rural terrain). In suburban/urban areas where ε ≥ 20°, LEO wins.

  STEP 7 — WALKER PRE-SIZING  (Week 4)
    For each shortlisted (orbit, ε_min) combination, compute P, S, N, F.
    Check P × satellite_unit_cost is within programme budget envelope.
    Check P is achievable with the planned launch cadence.
    Select F = P/2 as the default phasing for Tier 1/2 service.

  STEP 8 — INCLINATION SELECTION  (Week 4)
    Confirm Walker Delta inclination i = lat_max + λ.
    Check with launch vehicle team: is this inclination reachable from
    the planned launch site with the required payload mass?
    If not, identify the minimum feasible inclination and recheck coverage.

  STEP 9 — REVISIT TIME VALIDATION  (Week 4)
    Compute T_rev(h, ε, φ) at 36°, 40°, 44° N (southern, central, northern Spain).
    Verify that N = ⌈T_rev / T_gap⌉ at every latitude in the service
    zone is consistent with the N computed from Walker coverage.
    The Walker continuous N is always ≥ the intermittent N — confirm this
    is not violated (if it is, the Walker sizing has an error).

  STEP 10 — PROPULSION AND SUSTAINABILITY CHECK  (Week 5)
    For VLEO 200 km: confirm ABEP thrust budget closes.
      Required thrust ≈ 0.19 mN for a 10 kg satellite (from drag calculation).
      Check ABEP TRL and confirm a flight-ready solution exists for the
      target launch date. If not: VLEO is NOT a viable design point.
    For LEO 500 km: compute propellant mass for station-keeping
      and de-orbit at end-of-life. Confirm it fits the mass budget.

  STEP 11 — DOWN-SELECT TO 2 CANDIDATES  (Week 5)
    Score each candidate on a structured trade matrix:
    ┌────────────────────────────────┬────────────┬────────────┬────────┐
    │  Criterion                     │  VLEO 200  │  LEO 500   │ Weight │
    ├────────────────────────────────┼────────────┼────────────┼────────┤
    │  N (Tier 3, ε={EPSILON_MIN_DEG:.0f}°)           │  {N_t3_v:<6d} (−)  │  {N_t3_l:<4d} (+)   │  High  │
    │  Path loss advantage (FSPL)    │  {delta_fspl:.1f} dB (+)│  0 dB  (−) │  Med   │
    │  Pass duration (voice)         │  {Tpass_v:.1f} min(−)│  {Tpass_l:.1f} min(+)│  High  │
    │  Single-sat revisit            │  {Trev_v:.1f} h(−) │  {Trev_l:.1f} h(+) │  High  │
    │  Propulsion complexity (ABEP)  │  HIGH (−)  │  LOW   (+) │  High  │
    │  Technology readiness          │  TRL 4–5(−)│  TRL 8  (+)│  High  │
    │  Urban link margin (ε=30°)     │  EQUAL     │  EQUAL     │  Med   │
    │  Rural link margin (ε=10°)     │  VLEO wins │            │  Med   │
    └────────────────────────────────┴────────────┴────────────┴────────┘
    Deliver trade matrix to system PDR. Recommend LEO 500 km as baseline
    and VLEO 200 km as the high-risk/high-reward alternative for further
    study, conditional on ABEP technology maturation.
""")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  QUESTION 7 — MAIN RISKS AND UNKNOWNS                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
title("QUESTION 7 — MAIN RISKS AND UNKNOWNS")
body(f"""
  These are items that CANNOT be closed by orbital geometry alone.
  Each one requires direct iteration with another engineering discipline.
  Failure to track them will invalidate every number in the analysis.

  ══════════════════════════════════════════════════════════════════════════
  RISK 1 — ε_min IS NOT YET FROM THE LINK BUDGET  [CRITICAL — Priority 1]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    All N estimates in this study assume ε_min = 20° (suburban average).
    This number has NOT been validated by the link budget.

  WHY IT MATTERS (Spain 36–44° N context):
    ε_min = 20° → N_VLEO ≈ 1260,  N_LEO ≈ 280   (manageable programmes)
    ε_min = 30° → N_VLEO ≈ 2835,  N_LEO ≈ 560   (large programmes)
    ε_min = 40° → N_VLEO ≈ 5700+, N_LEO ≈ 1060+ (very large programmes)
    A 10° error in ε_min changes N by a factor of 2–4×. This can flip
    the programme from affordable to unaffordable.

  WHAT IS NEEDED:
    The payload team must supply: satellite EIRP at ε_min, G/T of the
    user terminal, required C/N₀, and an estimate of urban shadowing loss.
    These together give ε_min from the link budget, not from assumption.

  WHO OWNS THE ACTION:
    Mission analyst requests; payload/RF team delivers. Iteration loop
    between geometry (λ, R_s) and payload (EIRP, G/T) must be documented.

  ══════════════════════════════════════════════════════════════════════════
  RISK 2 — ABEP TECHNOLOGY READINESS FOR VLEO  [CRITICAL for VLEO]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    At 200 km altitude, atmospheric drag force requires continuous thrust
    to maintain orbit. The only viable technology is Air-Breathing Electric
    Propulsion (ABEP), which ingests ambient N₂/O molecules as propellant.
    Required thrust: ≈ 0.19 mN for a representative 10 kg satellite.
    ABEP TRL as of 2025: TRL 4–5 (demonstrated in laboratory, not in orbit).

  WHY IT MATTERS:
    Without ABEP, the satellite re-enters within hours at 200 km altitude.
    No conventional propulsion system carries enough propellant for more
    than a few weeks of drag compensation at 200 km. The entire VLEO
    value proposition depends on ABEP becoming flight-ready.

  WHAT IS NEEDED:
    ABEP supplier to confirm TRL roadmap and flight demonstration date.
    If TRL 7 (flight qualified) is not achievable before the target launch
    date, VLEO must be removed from the trade and the analysis collapses
    to LEO 400–600 km only.

  WHO OWNS THE ACTION:
    System engineer, in coordination with propulsion team and ABEP supplier.

  ══════════════════════════════════════════════════════════════════════════
  RISK 3 — SATELLITE EIRP / PHASED ARRAY SIZE  [CRITICAL — Priority 1]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    The D2D link requires the satellite to compensate for the phone's
    0 dBi antenna gain with a very large phased array.
    Estimated minimum EIRP ≈ 231 dBW (VLEO, ε=20°) to 238 dBW (LEO, ε=20°).
    AST SpaceMobile achieves this with a 64 m² deployable array.

  WHY IT MATTERS:
    The array area, mass, and DC power consumption are programme-defining.
    They determine the satellite bus size, launch mass, and cost per unit.
    A 1 dB error in the link budget translates to a ~26% change in
    required array area — a non-trivial mass and cost change.
    If the required EIRP exceeds what fits within the mass/cost envelope,
    ε_min must increase → more satellites → higher constellation cost.
    This is a circular dependency that must be iterated to convergence.

  WHAT IS NEEDED:
    Payload team to provide: array area, gain, efficiency, mass, DC power,
    and the ε_min at which the link closes for the target EIRP.

  WHO OWNS THE ACTION:
    Payload/RF team delivers; mission analyst updates N accordingly.

  ══════════════════════════════════════════════════════════════════════════
  RISK 4 — INTERFERENCE WITH TERRESTRIAL NETWORKS  [REGULATORY]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    The satellite downlink uses the same spectrum as terrestrial base
    stations. ITU and national regulators may impose Power Flux Density
    (PFD) limits on the satellite signal at the Earth's surface.
    These limits could cap the satellite EIRP below the link budget minimum.

  WHY IT MATTERS:
    If PFD limits force EIRP_max < EIRP_required:
      → Link does not close at ε_min → ε_min must increase
      → Larger ε_min → smaller footprint → more satellites
    This can effectively increase N by a factor of 2–4× through regulation
    alone, with no technical solution available.

  WHAT IS NEEDED:
    Legal/regulatory team to assess ITU/FCC PFD limits for the target band
    and geographic service area. Results must be fed back into the link budget
    as an EIRP cap constraint before any satellite sizing begins.

  WHO OWNS THE ACTION:
    Regulatory team, in coordination with payload (EIRP) and mission analyst (N).

  ══════════════════════════════════════════════════════════════════════════
  RISK 5 — HANDSET SOFTWARE COMPATIBILITY  [PROGRAMME RISK]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    D2D requires the smartphone to run modified protocol software to handle
    satellite-specific parameters (different timing advance, pre-compensation
    for geometric delays). This requires cooperation from:
      • Chipset vendors (Qualcomm, MediaTek) to enable satellite mode
      • Mobile Network Operators (MNOs) to provision the service
      • Regulatory bodies (FCC, ETSI) for type approval
    Type approval takes 1–3 years.

  WHY IT MATTERS:
    If chipset or operator agreement is not in place, the service cannot
    launch regardless of how well the satellite constellation is designed.
    This is a market access risk, not a technical risk — it can shut down
    the entire programme at the commercial rather than engineering level.

  WHAT IS NEEDED:
    Business development team to confirm chipset and MNO partnerships
    early in the programme. Target: signed agreement before CDR.

  WHO OWNS THE ACTION:
    Programme director, not the mission analyst — but the analyst must
    flag this risk explicitly in every mission analysis review.

  ══════════════════════════════════════════════════════════════════════════
  RISK 6 — PASS DURATION TOO SHORT FOR VOICE AT VLEO  [ARCHITECTURE RISK]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    At VLEO 200 km and ε = 20°, T_pass = 2.1 min.
    After subtracting connection setup (5 s) and handover margin (20 s),
    effective call time per satellite is ≈ 1.7 min.
    A normal phone call (minimum 3–5 min) requires continuous coverage,
    meaning another satellite must always be above ε_min when the first
    drops below — i.e., Tier 3 coverage is MANDATORY for voice at VLEO.

  WHY IT MATTERS:
    If the mission targets voice service at VLEO over Spain, N ≈ 1260 is non-negotiable.
    There is no architectural workaround: inter-satellite links (ISL) could
    relay the call, but that adds mass, cost, and complexity to every satellite.
    At LEO 500 km (T_pass = 4.9 min), the situation is better
    for voice — handover is still needed every ~5 min, which the network
    must handle seamlessly.

  WHAT IS NEEDED:
    Network architecture team to confirm handover capability at the planned
    T_pass values, and to decide whether ISL are baselined or excluded.

  ══════════════════════════════════════════════════════════════════════════
  RISK 7 — REVISIT TIME ASYMMETRY VLEO vs LEO  [ARCHITECTURE RISK]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    Single-sat revisit at 44° N / northern Spain (ε=20°):
      VLEO 200 km: T_rev = 21.9 h  — worse than LEO despite lower altitude
      LEO  500 km: T_rev = 10.9 h
    This is counter-intuitive. The reason: VLEO's smaller footprint reduces
    the swath width faster than its shorter orbital period compensates.

  WHY IT MATTERS:
    For Tier 1/2 services (IoT, messaging), which are the commercially
    achievable early revenue services, VLEO needs 2× more satellites than
    LEO for the same service quality. This is a permanent geometric penalty
    that cannot be overcome by any payload improvement.
    It means VLEO's path-loss advantage (lower EIRP needed) must be weighed
    against its geometry disadvantage (higher N needed) at every tier.

  WHAT IS NEEDED:
    System engineer to evaluate whether the EIRP savings at VLEO translate
    into cost savings per satellite that exceed the cost of the additional
    satellites. This is a techno-economic analysis, not a geometry analysis.

  ══════════════════════════════════════════════════════════════════════════
  RISK 8 — CONSTELLATION DEPLOYMENT PHASING  [PROGRAMME RISK]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    A full Tier 3 constellation ({N_t3_l} sats for LEO, {N_t3_v} for VLEO) requires
    multiple launch campaigns. Revenue must begin before full deployment.

  WHY IT MATTERS:
    The minimum operational constellation for Tier 1 revenue:
      VLEO: 12 satellites  → achievable in 1 launch
      LEO:  6 satellites   → achievable in 1 launch
    But investors expect a path to Tier 2/3, which requires visible progress.
    The deployment sequence (which planes to launch first) affects the
    geographic coverage ramp-up and the revenue model.

  WHAT IS NEEDED:
    Mission analyst to provide a deployment sequence plan showing coverage
    at N = 6, 12, 50, 100, {N_t3_l} satellites. Programme manager to confirm
    launch cadence is consistent with the business plan.

  ══════════════════════════════════════════════════════════════════════════
  RISK 9 — GROUND SEGMENT ARCHITECTURE  [SYSTEM RISK]
  ══════════════════════════════════════════════════════════════════════════
  WHAT IT IS:
    Each satellite in a LEO/VLEO constellation passes over a ground station
    approximately once per orbital period (every 88–95 min).
    To manage 400–2000 satellites in real time, the ground segment requires
    either a large network of ground stations (every ~1000 km) or ISL.
    ISL add mass, power, and cost to every satellite.

  WHY IT MATTERS:
    Without ISL, time-critical operations (orbit manoeuvres, emergency
    rerouting) can only be commanded during the ~5–10 min window per orbit.
    The network operations centre must be designed for this.
    The mission analyst must size the gateway network alongside the space
    segment, otherwise the total system cost is underestimated.

  WHAT IS NEEDED:
    Ground segment architect to define the number and location of ground
    stations, and the data volume per pass. Feeds into operations cost model.
""")

# ──────────────────────────────────────────────────────────────────────────────
#  CLOSING SUMMARY
# ──────────────────────────────────────────────────────────────────────────────
title("EXECUTIVE SUMMARY AND KEY DESIGN CONCLUSIONS")
body(f"""
  MISSION: Direct-to-Device telecom, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N, S-band 2 GHz, ε_min = {EPSILON_MIN_DEG:.0f}°

  ┌──────────────────────────────────────────────────────────────────────────┐
  │  DESIGN POINT COMPARISON  (ε_min = {EPSILON_MIN_DEG:.0f}°, f = {FREQ_GHZ:.1f} GHz, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N)    │
  │                                                                          │
  │  Parameter              VLEO 200 km      LEO 500 km        Winner       │
  │  ─────────────────────────────────────────────────────────────────────  │
  │  Footprint λ            {lam_v:.2f}°            {lam_l:.2f}°           LEO          │
  │  Path loss FSPL         {fspl_v:.1f} dB        {fspl_l:.1f} dB       VLEO         │
  │  Pass duration T_pass   {Tpass_v:.1f} min          {Tpass_l:.1f} min         LEO          │
  │  Single-sat revisit     {Trev_v:.1f} h           {Trev_l:.1f} h          LEO          │
  │  N Tier 1 (2 h gap)     {N_t1_v:<16d} {N_t1_l:<16d} LEO          │
  │  N Tier 2 (30 min)      {N_t2_v:<16d} {N_t2_l:<16d} LEO          │
  │  N Tier 3 (continuous)  {N_t3_v:<16d} {N_t3_l:<16d} LEO          │
  │  ABEP propulsion        Required         Not required      LEO          │
  │  Technology TRL         4–5              8–9               LEO          │
  └──────────────────────────────────────────────────────────────────────────┘

  RECOMMENDATION:
  1. Baseline: LEO 500 km, Walker {N_l}/{P_l}/{F_l}, i = {incl_l:.0f}°, ε_min = {EPSILON_MIN_DEG:.0f}°
     → N = {N_t3_l} for Tier 3 continuous coverage.
     → N = {N_t1_l} for Tier 1 IoT/SOS (achievable with 1 launch).

  2. VLEO 200 km remains compelling IF:
     (a) Service area is rural/suburban (ε_min ≤ 15°) → N = {N_eq_v}, beating LEO
     (b) ABEP achieves TRL 7 before the target launch date
     (c) Cost savings from smaller satellite antenna (7.5 dB less EIRP)
         outweigh the additional satellites required at urban ε_min

  3. CRITICAL OPEN ITEM: ε_min from the link budget.
     This one number changes N by a factor of 2–4×.
     No programme decision should be made until this is validated.

  All values computed at: ε_min = {EPSILON_MIN_DEG:.0f}°, f = {FREQ_GHZ:.1f} GHz, Spain {TARGET_LAT_MIN:.0f}–{TARGET_LAT_MAX:.0f}° N (φ_max = {TARGET_LAT_MAX:.0f}° N),
  R_E = {R_E/1e3:.0f} km, μ = {MU:.4e} m³/s².
""")

print("\n" + "=" * W)
print("  END OF DOCUMENT")
print("=" * W)
