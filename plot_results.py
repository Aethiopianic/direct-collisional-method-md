# -*- coding: utf-8 -*-
"""
analysis.py
-----------
Post-processing and visualization for the Direct Collisional MD Simulation.

Reads three CSV outputs from collision.f90:
  - tracer.csv       : vx, x, y of the tracked particle at each time step
  - initial_v.csv    : initial velocities of all particles
  - final_v.csv      : final velocities of all particles

Author: Getabalew Manegerew
Date  : May 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# =============================================================================
# 1. Parameters
# =============================================================================

path = r"."           # <-- change this to your CSV folder path

m  = 1.0e-26          # particle mass [kg]
kb = 1.38e-23         # Boltzmann constant [J/K]
Np = 150              # number of particles

# T_init is the temperature used to draw the Maxwell-Boltzmann reference curve.
# This should match the temperature you used when initializing velocities in
# the Fortran code (T_init in the PARAMETER block), NOT the computed T_sim.
T_init = 300.0        # [K]  -- initial velocity sampling temperature

# =============================================================================
# 2. Load data
# =============================================================================

tracer = pd.read_csv(f"{path}/tracer.csv", sep=r"\s+")
print("=== Tracer Data Loaded ===")
print(tracer.head())

vel_in = pd.read_csv(
    f"{path}/initial_v.csv", sep=r"\s+", skiprows=1, names=["vx", "vy", "vz"]
)
print("\n=== Initial Velocities Loaded ===")
print(vel_in.head())

vel_f = pd.read_csv(
    f"{path}/final_v.csv", sep=r"\s+", skiprows=1, names=["vfx", "vfy", "vfz"]
)
print("\n=== Final Velocities Loaded ===")
print(vel_f.head())

# =============================================================================
# 3. Extract arrays
# =============================================================================

tracer_vx   = tracer["vx_tracer"]
tracer_posx = tracer["x_tracer"]
tracer_posy = tracer["y_tracer"]

vfx = vel_f["vfx"].to_numpy(dtype=float)
vfy = vel_f["vfy"].to_numpy(dtype=float)
vfz = vel_f["vfz"].to_numpy(dtype=float)

# =============================================================================
# 4. Compute T_sim from the final state (equipartition theorem)
#    T_sim = m * <v^2> / (3 * kb)   -- averaged over all particles and 3 DOF
# =============================================================================

v2_mean = np.mean(vfx**2 + vfy**2 + vfz**2)   # mean squared speed [m^2/s^2]
T_sim   = m * v2_mean / (3 * kb)

print(f"\n=== Temperatures ===")
print(f"  T_init (velocity sampling) : {T_init:.2f} K")
print(f"  T_sim  (from final state)  : {T_sim:.2f} K")

# =============================================================================
# 5. Thermodynamic quantities from the tracer particle
#    Uses T_sim (the actual system temperature) for comparison with theory
# =============================================================================

nsize   = len(tracer_vx)
vx_arr  = tracer_vx.to_numpy(dtype=float)

E_sim   = 0.5 * m * vx_arr**2                  # KE per time step [J]
E_mean  = np.mean(E_sim)
E2_mean = np.mean(E_sim**2)
Cv      = (E2_mean - E_mean**2) / (kb * T_sim**2)

print(f"\n=== Thermodynamic Summary (T_sim = {T_sim:.2f} K) ===")
print(f"  <E>  sim    [1e-21 J]  : {E_mean*1e21:.4f}")
print(f"  <E>  theory [1e-21 J]  : {0.5*kb*T_sim*1e21:.4f}   (= 0.5 kb T_sim)")
print(f"  <Cv> sim    [1e-24 J/K]: {Cv*1e24:.4f}")
print(f"  <Cv> theory [1e-24 J/K]: {0.5*kb*1e24:.4f}   (= 0.5 kb)")

# =============================================================================
# 6. Figure 1 — Velocity distribution vs. Maxwell-Boltzmann theory
#    Histogram : final vx of all particles (full ensemble snapshot)
#    Theory    : drawn at T_init (the temperature velocities were sampled from)
# =============================================================================

a     = m / (2 * kb * T_init)
v     = np.linspace(-2000, 2000, 1000)
pdf   = np.sqrt(a / np.pi) * np.exp(-a * v**2)

plt.figure("hist", figsize=(7, 4.5))
plt.plot(v, pdf, label=f"Theory  (T = {T_init:.0f} K)", color="steelblue", linewidth=2)
plt.hist(vfx, bins=50, density=True, alpha=0.75, color="darkorange", label="Simulation")
plt.xlabel("Velocity in x direction [m/s]")
plt.ylabel("Probability density")
plt.title("Velocity Distribution — x component")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

# =============================================================================
# 7. Figure 2 — Tracer particle vx(t)
# =============================================================================

plt.figure("tracer", figsize=(9, 4))
plt.plot(vx_arr, label="Gas (tracer)", color="red", linewidth=0.8)
plt.xlabel("Time step")
plt.ylabel("vx [m/s]")
plt.title("Tracer Particle — x Velocity vs. Time")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.show()
