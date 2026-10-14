# Satellite Mission Analysis Scripts

Python scripts for preliminary mission analysis of Low Earth Orbit (LEO) and Very Low Earth Orbit (VLEO) satellites, with a focus on two mission archetypes: **Direct-to-Device (D2D) telecommunications** and **wildfire detection** over Spain.

All scripts compare a **VLEO 200 km** design point against a **LEO 500 km** (Starlink-class) reference.

---

## Scripts

### `Telecom_Constellation.py`
Full Walker constellation sizing for a D2D telecom mission over Spain (36–44° N, S-band 2 GHz).

Covers:
- Footprint geometry and slant range
- Free-space path loss and one-way latency
- Service tier analysis (IoT/SOS, messaging, continuous voice)
- Elevation angle sensitivity (how N scales with ε_min)
- Path-loss equivalence trade: *can VLEO's 7 dB advantage compensate for its smaller footprint?*
- Walker T/P/F sizing with step-by-step reasoning

Outputs: `PathLoss_Equivalence.png`, `Telecom_Constellation_FOMs.png`

---

### `D2D_Mission_Analysis_Document.py`
Structured answers to the seven standard questions of a preliminary design assessment for the same D2D mission. Companion document to `Telecom_Constellation.py`, formatted for reporting.

---

### `FireDetect.py`
Instrument and constellation sizing for a thermal-IR pushbroom fire-detection mission over Galicia.

Covers:
- Ground Sampling Distance (GSD) and swath width vs altitude
- Off-nadir angle analysis
- Revisit time and required number of satellites per service tier

Outputs: `GSD_Swath_vs_Altitude.png`, `OffNadir_Analysis.png`, `ABEP_Thrust.png`

---

### `Propagate_Orbits.py`
Numerical orbit propagation (J2 + atmospheric drag, RK45) for both VLEO and LEO SSO orbits.

Covers:
- 24 h ground tracks (ECI → ECEF → geodetic)
- Long-term altitude decay for realistic ballistic coefficient
- Atmospheric density profile (exponential model, Vallado 2013)

Outputs: `Ground_Tracks.png`, `Orbital_Decay.png`

---

### `Mission_Comparison.py`
Side-by-side parametric comparison of the telecom and fire-detection missions. Reproduces the key figures of merit from both and highlights where the two mission types agree and diverge — and why.

---

## Dependencies

```
numpy
matplotlib
scipy
cartopy   # Propagate_Orbits.py only
```

Install with:
```bash
pip install numpy matplotlib scipy cartopy
```

---

## Quick start

Each script is self-contained and can be run directly:

```bash
python Telecom_Constellation.py
python FireDetect.py
python Propagate_Orbits.py
python Mission_Comparison.py
python D2D_Mission_Analysis_Document.py
```

Results are printed to stdout and plots are saved as PNG files in the same directory.

---

## Context

These scripts were written as part of a mission analysis study for a D2D telecom constellation and an Earth-observation fire-detection constellation, both targeting the Iberian Peninsula. The core trade is whether VLEO (200 km) with Air-Breathing Electric Propulsion (ABEP) is preferable to conventional LEO (500 km) given the path-loss advantage vs the higher satellite count and drag penalty.
