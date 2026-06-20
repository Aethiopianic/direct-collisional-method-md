# Direct Collisional Method of Molecular Dynamics Simulation

> **BSc Senior Thesis** — Department of Physics, College of Natural and Computational Sciences, Addis Ababa University (AAU)
> **Author:** Getabalew Manegerew | **Advisor:** Dr. Lemi D. | **April 2026**

---

## Overview

This repository contains the Fortran 90 implementation and supporting analysis scripts for a hybrid gas-dynamics simulation model proposed in my senior thesis. The **Direct Collisional Method** blends key features of two established approaches:

| Method | Strength | Weakness |
|--------|----------|----------|
| Direct Simulation Monte Carlo (DSMC) | Fast — O(N) scaling | No continuous particle tracking |
| Molecular Dynamics (MD) | Full trajectory history preserved | Catastrophic O(N²) scaling |

The proposed model captures the physical accuracy of MD (full phase-space history, explicit particle identities) at a cost closer to DSMC (no continuous force-field evaluation). Collisions are detected geometrically and resolved in the Center-of-Mass Frame (CMF) using isotropic elastic scattering.

---

## Physical Background

### The Knudsen Regime
When the mean free path λ is comparable to the system size L, continuum fluid equations (Navier–Stokes) break down. The **Knudsen number** Kn = λ/L characterizes this:

- **Kn < 0.01** → continuum regime, standard fluid mechanics applies  
- **Kn > 0.01** → rarefied regime, particle-level kinetics must be resolved ← *this model's domain*

### Hard-Sphere Potential
Particles interact via a rigid hard-sphere potential:

```
φ(r) = ∞   for r ≤ dl
φ(r) = 0   for r > dl
```

### Collision Mechanics (CMF)
For two colliding particles with masses m₁, m₂:

- CM velocity **V**_c is conserved (no external forces)
- Elastic scattering: |**v**_r| is conserved, direction randomized uniformly over 4π steradians
- Post-collision velocities recovered via reduced-mass transformation

---

## Repository Structure

```
.
├── collision.f90          # Main Fortran simulation program
├── analysis/
│   └── plot_results.py    # Python script for reading output CSVs and plotting
├── output/                # (generated at runtime)
│   ├── initial_v.csv      # Initial velocity distribution of all particles
│   ├── final_v.csv        # Final velocity distribution of all particles
│   └── tracer.csv         # Phase-space trajectory of the tracked particle
├── thesis/
│   └── Getabalew_Senior_Thesis_Final.pdf
├── defense/
│   └── Getabalew_Thesis_Defense.pdf
└── README.md
```

---

## Simulation Parameters

The default parameters represent a low-density atmospheric gas:

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Number of particles | N_p | 150 |
| Number density | n | 2.34 × 10²⁴ m⁻³ |
| Initial temperature | T | 300 K |
| Molecular diameter | dl | 4.0 × 10⁻¹⁰ m |
| Box volume | V_box | 6.4 × 10⁻²³ m³ |
| Particle mass | M | 10⁻²⁶ kg |
| Time step | dt | 0.1 ps |
| Total simulation time | t_span | 1.7 µs |

Parameters can be modified directly in the `PARAMETER` block at the top of `collision.f90`.

---

## Getting Started

### Requirements
- A Fortran compiler: [gfortran](https://gcc.gnu.org/fortran/) (recommended) or ifort
- Python 3.x with `numpy`, `matplotlib`, and `pandas` (for analysis scripts)

### Compile and Run

```bash
# Compile
gfortran -O2 -o collision collision.f90

# Run
./collision
```

The program prints a simulation summary to stdout and writes three CSV files to the working directory.

### Plot Results

```bash
cd analysis
python plot_results.py
```

This generates:
- Velocity distribution histogram vs. Maxwell-Boltzmann theory curve
- Tracer particle trajectory X(t)

---

## Key Results

At T = 343 K and n = 2.34 × 10²⁴ m⁻³, computed quantities compared to theoretical predictions:

| Quantity | Computed | Theory | Ratio |
|----------|----------|--------|-------|
| Mean energy ⟨E⟩ | 1.927 × 10⁻²¹ J | 2.367 × 10⁻²¹ J | 0.814 |
| Specific heat ⟨c_v⟩ | 5.388 × 10⁻²⁴ J/K | 6.900 × 10⁻²⁴ J/K | 0.781 |
| Mean free path λ | 1.216 × 10⁻⁶ m | 6.001 × 10⁻⁶ m | 0.203 |

⟨E⟩ and ⟨c_v⟩ agree with theory to within ~20%. The mean free path discrepancy (~5×) is attributed to hard-sphere packing effects at the simulated density; this is an active area for model improvement.

The simulated velocity distribution of a single traced particle converges to the Maxwell-Boltzmann distribution over the full simulation window, confirming ergodic behavior.

---

## Application: Radiation Stopping Power

The model was also used to qualitatively reproduce energy-loss behavior of energetic particles traversing the gas — analogous to nuclear radiation passing through a gas detector. The maximum energy transfer per collision is:

```
T_max = [4 m₁ m₂ / (m₁ + m₂)²] × K_in
```

Simulations correctly reproduce: a **massive particle** (m₁ ≫ m₂) traversing nearly undisturbed, and an **equal-mass particle** (m₁ ≈ m₂) stopping rapidly after successive collisions.

---

## Limitations and Future Work

- Hard-sphere model is restricted to **dilute (low-density)** regimes; long-range interactions are neglected
- Single-particle ergodic analysis is sensitive to finite simulation time
- Planned extensions:
  - Longer runs for better statistical convergence of ⟨E⟩ and ⟨c_v⟩
  - External field support (gravity, electromagnetic) via explicit velocity updates
  - Soft potentials (e.g., Lennard-Jones) for dense-gas regimes

---

## Citation

If you use this code or build on this work, please cite:

```
Getabalew Manegerew. "Direct Collisional Method of Molecular Dynamics Simulation."
BSc Senior Thesis, Department of Physics, Addis Ababa University, April 2026.
```

---

## License

This project is released under the [MIT License](LICENSE). You are free to use, modify, and distribute the code with attribution.

---

## Acknowledgements

Advised by **Dr. Lemi D.**, Department of Physics, AAU. Thanks also to Dr. Mesfin (course instructor) and fellow students for discussions throughout the project.
