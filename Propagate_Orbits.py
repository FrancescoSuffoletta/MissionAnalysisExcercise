"""
================================================================================
  ORBIT PROPAGATION  —  J2 + ATMOSPHERIC DRAG
  Ground Tracks and Altitude Decay  |  VLEO (200 km) vs LEO (500 km) SSO
================================================================================

Dynamical model (ECI Cartesian, RK45 adaptive integrator):
  1. Two-body central gravity
  2. J2 oblateness  (dominant secular perturbation → drives SSO nodal precession)
  3. Atmospheric drag with co-rotating exponential atmosphere (Vallado 2013)

Coordinate pipeline:
  ECI (propagation)  →  ECEF (Earth rotation)  →  geodetic lat/lon (ground track)

Two simulations per mission:
  A) Ground track  — 24 h, small A/m to preserve orbit for clean visualisation
  B) Altitude decay — long-term, realistic A/m to quantify drag impact

Outputs:
  Ground_Tracks.png   —  24 h SSO ground track on world map, Galicia starred
  Orbital_Decay.png   —  altitude vs time + atmospheric density profile
================================================================================
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ══════════════════════════════════════════════════════════════════════════════
#  PHYSICAL CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════
MU       = 3.986004418e14   # [m³/s²]  Earth gravitational parameter
R_E      = 6_371_000.0      # [m]      Earth mean radius
J2       = 1.08263e-3       # [-]      Earth oblateness coefficient
W_EARTH  = 7.2921150e-5     # [rad/s]  Earth sidereal rotation rate
H_KARMAN = 100_000.0        # [m]      reentry altitude threshold

# ══════════════════════════════════════════════════════════════════════════════
#  SATELLITE PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
CD       = 2.2     # [-]       aerodynamic drag coefficient (flat-plate approx)
AM_TRACK = 0.001   # [m²/kg]   area/mass for 24 h ground track  (weak drag)
AM_DECAY = 0.01    # [m²/kg]   area/mass for decay study  (e.g. 10 kg, 0.1 m²)

# ── Air-Breathing Electric Propulsion (ABEP) thrust compensation ─────────────
# Required thrust = drag force = 1/2 * CD * (A/m) * rho * v² * m_sat
# For a satellite of mass M [kg]:
#   F_drag = 1/2 * CD * A * rho * v²   [N]  where A = AM_DECAY * M
#   a_drag = F_drag / M = 1/2 * CD * AM_DECAY * rho * v²  [m/s²]
#
# At VLEO 200 km:  rho = 2.789e-10 kg/m³,  v = 7789 m/s
#   a_drag = 0.5 * 2.2 * 0.01 * 2.789e-10 * 7789²  = 1.86e-4 m/s²
#   For M = 10 kg:  F_thrust = M * a_drag = 1.86e-3 N  ≈ 1.86 mN
#   For M = 100 kg: F_thrust = 1.86e-2 N  ≈ 18.6 mN
#
# ABEP concept: intake collects residual atmospheric N₂/O atoms at ~200 km,
# feeds them directly into an ion/Hall thruster as propellant.
# Typical ABEP specific impulse Isp ~ 1500–3000 s (much better than cold gas).
# Required mass flow:  mdot = F / (Isp * g0)
#   For F = 0.19 mN, Isp = 2000 s:  mdot = 1.86e-4 / (2000 * 9.81) ≈ 9.5e-9 kg/s
#   Intake collects:  mdot_atm = rho * v * A_intake
#   For A_intake = A_sat = 0.1 m²:  mdot_atm = 2.789e-10 * 7791 * 0.1 ≈ 2.17e-10 kg/s...
#   → ABEP efficiency must be >> 1 (impossible) unless A_intake >> A_drag
#   In practice ABEP requires large intake (A_intake / A_drag ≈ 10–100),
#   making it a key design driver for VLEO satellites (e.g. JAXA RAM-EP, ESA GOCE).
#
# In the ABEP-compensated simulation below, thrust exactly cancels drag → orbit maintained.
ABEP_ACTIVE = True   # set False to see uncompensated VLEO reentry
#
# Decay rates WITHOUT compensation (analytical):
#   da/dt = -CD * (A/m) * rho * sqrt(mu*a)
#   VLEO 200 km:  rho ~ 2.8e-10 kg/m³  → da/dt ~ -27 km/day  (reentry in ~4 days!)
#   LEO  500 km:  rho ~ 7.0e-13 kg/m³  → da/dt ~ -0.07 km/day (25 km/year)

# ══════════════════════════════════════════════════════════════════════════════
#  EXPONENTIAL ATMOSPHERE MODEL  (Vallado 2013, Table 8-4)
#  Each row: (h_base [km], rho_ref [kg/m³], scale_height [km])
# ══════════════════════════════════════════════════════════════════════════════
_ATM = np.array([
    (  0,  1.225000e+00,  7.249), ( 25,  3.899000e-02,  6.349),
    ( 30,  1.774000e-02,  6.682), ( 40,  3.972000e-03,  7.554),
    ( 50,  1.057000e-03,  8.382), ( 60,  3.206000e-04,  7.714),
    ( 70,  8.770000e-05,  6.549), ( 80,  1.905000e-05,  5.799),
    ( 90,  3.396000e-06,  5.382), (100,  5.297000e-07,  5.877),
    (110,  9.661000e-08,  7.263), (120,  2.438000e-08,  9.473),
    (130,  8.484000e-09, 12.636), (140,  3.845000e-09, 16.149),
    (150,  2.070000e-09, 22.523), (180,  5.464000e-10, 29.740),
    (200,  2.789000e-10, 37.105), (250,  7.248000e-11, 45.546),
    (300,  2.418000e-11, 53.628), (350,  9.158000e-12, 53.298),
    (400,  3.925000e-12, 58.515), (450,  1.585000e-12, 60.828),
    (500,  6.967000e-13, 63.822), (600,  1.454000e-13, 71.835),
    (700,  3.614000e-14, 88.667), (800,  1.170000e-14,124.640),
    (900,  5.245000e-15,181.050), (1000, 3.019000e-15,268.000),
])

def atm_density(h_m):
    """
    Atmospheric density [kg/m³] at altitude h_m [m].
    Piecewise exponential fit: rho = rho_ref * exp(-(h - h_base) / H)
    Valid from 0 to ~1000 km. Representative of moderate solar activity.
    """
    h_km = max(min(h_m / 1000.0, 999.9), 0.0)
    idx  = int(np.searchsorted(_ATM[:, 0], h_km, side='right')) - 1
    idx  = max(0, min(idx, len(_ATM) - 2))
    h0, rho0, H = _ATM[idx]
    return rho0 * math.exp(-(h_km - h0) / H)

# ══════════════════════════════════════════════════════════════════════════════
#  EQUATIONS OF MOTION  (ECI Cartesian, state = [x,y,z, vx,vy,vz])
# ══════════════════════════════════════════════════════════════════════════════
def eom(t, state, CD, AM, compensate_drag=False):
    """
    Equations of motion in ECI.

    Forces:
      1) Two-body:  a = -mu/r^3 * r_vec
      2) J2:        a = -3/2 * mu*J2*R_E^2 / r^5 * [x(1-5(z/r)^2),
                                                      y(1-5(z/r)^2),
                                                      z(3-5(z/r)^2)]
      3) Drag:      a = -1/2 * CD * (A/m) * rho * |v_rel| * v_rel
         where v_rel = v_sat - v_atm,  v_atm = omega_E x r = (-We*y, We*x, 0)
      4) ABEP thrust (optional): cancels drag exactly when compensate_drag=True
         This models an Air-Breathing Electric Propulsion system that collects
         atmospheric particles and expels them to counteract drag.
    """
    x, y, z, vx, vy, vz = state
    r2  = x*x + y*y + z*z
    r   = math.sqrt(r2)

    # 1) Two-body
    c3  = -MU / (r2 * r)
    ax, ay, az = c3*x, c3*y, c3*z

    # 2) J2
    j2c = -1.5 * MU * J2 * R_E**2 / (r2 * r2 * r)
    zr2 = (z / r)**2
    ax += j2c * x * (1.0 - 5.0 * zr2)
    ay += j2c * y * (1.0 - 5.0 * zr2)
    az += j2c * z * (3.0 - 5.0 * zr2)

    # 3) Drag  (co-rotating atmosphere)
    rho = atm_density(r - R_E)
    vrx = vx + W_EARTH * y    # v_rel_x = vx - v_atm_x = vx - (-We*y)
    vry = vy - W_EARTH * x    # v_rel_y = vy - v_atm_y = vy - ( We*x)
    vrz = vz
    vr  = math.sqrt(vrx*vrx + vry*vry + vrz*vrz)
    df  = -0.5 * CD * AM * rho * vr
    ax += df * vrx
    ay += df * vry
    az += df * vrz

    # 4) ABEP thrust: exactly cancels drag (net zero aerodynamic acceleration)
    if compensate_drag:
        ax -= df * vrx
        ay -= df * vry
        az -= df * vrz

    return [vx, vy, vz, ax, ay, az]

# ══════════════════════════════════════════════════════════════════════════════
#  HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════
def initial_state(h_m, inc_deg, raan_deg):
    """
    Initial ECI state for a circular orbit at the ascending node (true anomaly=0).
    Position: on the equatorial plane at the ascending node longitude.
    Velocity: perpendicular to radius in the orbital plane.
    """
    a   = R_E + h_m
    vc  = math.sqrt(MU / a)
    i   = math.radians(inc_deg)
    Om  = math.radians(raan_deg)
    return [
        a * math.cos(Om),
        a * math.sin(Om),
        0.0,
        vc * (-math.sin(Om) * math.cos(i)),
        vc * ( math.cos(Om) * math.cos(i)),
        vc *   math.sin(i),
    ]

def eci_to_latlon(xyz_arr, t_arr, theta0=0.0):
    """
    Vectorised ECI → geodetic lat/lon (spherical Earth).
    xyz_arr : (3, N) position array [m]
    t_arr   : (N,)   times [s]
    theta0  : GAST at t=0 [rad] (default 0 — relative ground track)
    Returns : lat [deg], lon [deg]  each shape (N,)
    """
    theta = theta0 + W_EARTH * t_arr
    xe =  xyz_arr[0] * np.cos(theta) + xyz_arr[1] * np.sin(theta)
    ye = -xyz_arr[0] * np.sin(theta) + xyz_arr[1] * np.cos(theta)
    ze =  xyz_arr[2]
    r  = np.sqrt(xe**2 + ye**2 + ze**2)
    return np.degrees(np.arcsin(ze / r)), np.degrees(np.arctan2(ye, xe))

def make_reentry_event():
    """Terminal event: stops integration when altitude < Kármán line (100 km)."""
    def event(t, y):
        return math.sqrt(y[0]**2 + y[1]**2 + y[2]**2) - (R_E + H_KARMAN)
    event.terminal  = True
    event.direction = -1
    return event

# ══════════════════════════════════════════════════════════════════════════════
#  MISSION DEFINITIONS
#  RAAN = -2° → first ascending crossing at ~42.8°N lands near Galicia (-8°W)
# ══════════════════════════════════════════════════════════════════════════════
MISSIONS = {
    "VLEO 200 km": dict(h=200e3, inc=96.32, raan=-2.0, color="#e74c3c"),
    "LEO  500 km": dict(h=500e3, inc=97.39, raan=-2.0, color="#2980b9"),
}
GALICIA_LAT, GALICIA_LON = 42.8, -8.0

T_TRACK_S = 24 * 3600          # ground track: 1 day
T_DECAY_S = {
    "VLEO 200 km":  7 * 86400, # 7 days   (satellite likely reenters before this)
    "LEO  500 km": 365 * 86400, # 1 year
}

# ══════════════════════════════════════════════════════════════════════════════
#  PROPAGATION
# ══════════════════════════════════════════════════════════════════════════════
results = {}

for name, cfg in MISSIONS.items():
    print(f"\n  Propagating {name} ...")
    s0 = initial_state(cfg["h"], cfg["inc"], cfg["raan"])

    # ── A) Ground track: 24 h, small drag ────────────────────────────────
    t_eval_trk = np.arange(0.0, T_TRACK_S, 10.0)   # 10-s samples
    sol_trk = solve_ivp(
        lambda t, y: eom(t, y, CD, AM_TRACK),
        (0, T_TRACK_S), s0,
        method="RK45", t_eval=t_eval_trk,
        rtol=1e-9, atol=1e-11,
    )
    lat_trk, lon_trk = eci_to_latlon(sol_trk.y[:3], sol_trk.t)
    alt_trk = (np.linalg.norm(sol_trk.y[:3], axis=0) - R_E) / 1e3  # [km]

    # ── B) Altitude decay: long-term, full drag ──────────────────────────
    abep_on = ABEP_ACTIVE and (name == "VLEO 200 km")
    t_eval_dec = np.linspace(0, T_DECAY_S[name], 8000)
    sol_dec = solve_ivp(
        lambda t, y, _n=name: eom(t, y, CD, AM_DECAY,
                                   compensate_drag=(ABEP_ACTIVE and _n == "VLEO 200 km")),
        (0, T_DECAY_S[name]), s0,
        method="RK45", t_eval=t_eval_dec,
        events=make_reentry_event(),
        rtol=1e-8, atol=1e-10,
    )
    alt_dec = (np.linalg.norm(sol_dec.y[:3], axis=0) - R_E) / 1e3  # [km]

    # ── Summary ──────────────────────────────────────────────────────────
    h0 = cfg["h"] / 1e3
    hf = alt_dec[-1]
    tf = sol_dec.t[-1] / 86400
    abep_tag = " [ABEP ON - drag compensated]" if abep_on else ""
    if sol_dec.status == 1:
        print(f"  => REENTRY after {tf:.2f} days  (final h = {hf:.1f} km){abep_tag}")
    else:
        print(f"  => After {tf:.0f} days:  h = {hf:.2f} km  "
              f"(total decay = {h0-hf:.2f} km, "
              f"mean rate = {(h0-hf)/tf:.3f} km/day){abep_tag}")

    results[name] = dict(
        lat=lat_trk, lon=lon_trk,
        t_trk=sol_trk.t / 3600,   alt_trk=alt_trk,
        t_dec=sol_dec.t / 86400,  alt_dec=alt_dec,
        reentry=(sol_dec.status == 1),
        color=cfg["color"], inc=cfg["inc"],
    )

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 1 — GROUND TRACKS  (24 h, cartopy)
# ══════════════════════════════════════════════════════════════════════════════
PC = ccrs.PlateCarree()
fig1, axes = plt.subplots(
    2, 1, figsize=(15, 9),
    subplot_kw=dict(projection=ccrs.Robinson()),
)
fig1.suptitle(
    "SSO Ground Tracks — 24 h propagation  (J2 + drag,  A/m = "
    f"{AM_TRACK*1e3:.1f}×10⁻³ m²/kg)",
    fontsize=13, fontweight="bold",
)

for ax, (name, res) in zip(axes, results.items()):
    # Map features
    ax.add_feature(cfeature.LAND,      facecolor="#e8e4dc", zorder=0)
    ax.add_feature(cfeature.OCEAN,     facecolor="#cce4f5", zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5,        zorder=1)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.3, linestyle=":", zorder=1)
    gl = ax.gridlines(draw_labels=True, linewidth=0.35,
                      color="gray", alpha=0.6, zorder=1)
    gl.top_labels = gl.right_labels = False
    ax.set_global()

    # Colour-code by time
    lats, lons = res["lat"], res["lon"]
    t_h = res["t_trk"]
    cmap = plt.get_cmap("plasma")

    # Break ground track at longitude wraps to avoid horizontal artefacts
    breaks = np.where(np.abs(np.diff(lons)) > 180)[0] + 1
    segs_l = np.split(lons, breaks)
    segs_a = np.split(lats, breaks)
    segs_t = np.split(t_h,  breaks)

    for sl, sa, st in zip(segs_l, segs_a, segs_t):
        sc = ax.scatter(sl, sa, c=st, cmap=cmap, vmin=0, vmax=24,
                        s=0.6, transform=PC, zorder=3)

    # Galicia marker
    ax.plot(GALICIA_LON, GALICIA_LAT, marker="*",
            color="red", markersize=14, transform=PC,
            zorder=5, label="Galicia", linestyle="None")
    ax.annotate("Galicia", xy=(GALICIA_LON, GALICIA_LAT),
                xytext=(GALICIA_LON + 15, GALICIA_LAT + 5),
                fontsize=9, color="red", fontweight="bold",
                transform=PC,
                arrowprops=dict(arrowstyle="->", color="red", lw=1.0))

    cbar = plt.colorbar(sc, ax=ax, orientation="vertical",
                        pad=0.02, shrink=0.7, label="Time [h]")

    ax.set_title(
        f"{name}  —  SSO  i = {res['inc']:.2f}°  |  "
        f"RAAN = -2°  |  start at ascending node",
        fontsize=10,
    )

plt.tight_layout()
plt.savefig("Ground_Tracks.png", dpi=150, bbox_inches="tight")
print("\n  Saved Ground_Tracks.png")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 2 — ALTITUDE DECAY + ATMOSPHERIC DENSITY PROFILE
# ══════════════════════════════════════════════════════════════════════════════
fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))
fig2.suptitle(
    f"Orbital Decay  (J2 + drag,  CD = {CD},  A/m = {AM_DECAY} m²/kg)",
    fontsize=13, fontweight="bold",
)

for ax, (name, res) in zip(axes2[:2], results.items()):
    t   = res["t_dec"]
    h   = res["alt_dec"]
    c   = res["color"]
    h0  = h[0]
    hf  = h[-1]
    tf  = t[-1]

    ax.plot(t, h, color=c, linewidth=1.8, label=name)
    ax.axhline(h0, color="gray",  linestyle="--", linewidth=0.9,
               alpha=0.7, label=f"Initial: {h0:.1f} km")

    if res["reentry"]:
        ax.axhline(H_KARMAN / 1e3, color="red", linestyle=":", linewidth=1.2,
                   label="Karman line (100 km)")
        ax.annotate(f"Reentry\nt = {tf:.2f} d",
                    xy=(tf, hf), xytext=(tf * 0.5, hf + 20),
                    fontsize=9, color="red",
                    arrowprops=dict(arrowstyle="->", color="red"))
    else:
        ax.annotate(
            f"After {tf:.0f} d:\n{hf:.1f} km\n(-{h0-hf:.1f} km)",
            xy=(tf, hf), xytext=(tf * 0.55, hf + (h0 - hf) * 0.4),
            fontsize=9, color=c,
            arrowprops=dict(arrowstyle="->", color=c),
        )

    ax.set_xlabel("Time [days]", fontsize=11)
    ax.set_ylabel("Altitude [km]", fontsize=11)
    ax.set_title(name, fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# ── Atmospheric density profile ───────────────────────────────────────────────
ax3 = axes2[2]
h_arr   = np.linspace(100, 800, 600)
rho_arr = [atm_density(h * 1e3) for h in h_arr]

ax3.semilogy(rho_arr, h_arr, "k-", linewidth=2)
for name, cfg in MISSIONS.items():
    h_km = cfg["h"] / 1e3
    rho  = atm_density(cfg["h"])
    ax3.axhline(h_km, color=results[name]["color"],
                linestyle="--", linewidth=1.3,
                label=f"{name}:  {rho:.2e} kg/m³")
    ax3.plot(rho, h_km, "o", color=results[name]["color"], markersize=8, zorder=5)

ax3.set_xlabel("Density  [kg/m³]", fontsize=11)
ax3.set_ylabel("Altitude [km]",    fontsize=11)
ax3.set_title("Exponential Atmosphere\n(Vallado 2013)", fontsize=11)
ax3.legend(fontsize=9, loc="upper right")
ax3.grid(True, which="both", linestyle="--", alpha=0.5)
ax3.set_ylim(100, 800)
ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())  # linear y-axis: OK
# x-axis is log scale — AutoMinorLocator not used there

plt.tight_layout()
plt.savefig("Orbital_Decay.png", dpi=150, bbox_inches="tight")
print("  Saved Orbital_Decay.png")

plt.show()
