# Flowforge CFD

An open-source, browser-based learning lab for computational fluid dynamics. Build an experiment, solve a 2D laminar flow field, and inspect velocity, pressure, and thermal transport without installing a solver.

**Live application:** https://raghav1822.github.io/flowforge-cfd/

## What you can explore

- D2Q9 single-relaxation-time (BGK) lattice Boltzmann flow solver with a velocity inlet, open outlet, solid walls, and a movable cylinder or square obstacle.
- Velocity magnitude, relative pressure, and normalized temperature contours, with optional streamline and mesh overlays.
- Reynolds number, pressure-drop, and outlet-temperature readouts.
- A three-resolution, normalized pressure-coefficient mesh-convergence study.
- Adjustable domain, blockage ratio, inlet speed, kinematic viscosity, and thermal diffusivity.

The application runs locally in your browser and does not upload experiment data. It is an educational approximation, not a validated engineering design tool. In particular, results are sensitive to the compact grid and simplified boundary conditions.

## Run locally

Open `site/index.html` in a browser, or serve the `site` directory with any static HTTP server:

```sh
python3 -m http.server 8000 --directory site
```

Visit http://localhost:8000. GitHub Actions deploys the contents of `site/` to GitHub Pages on every push to `main`.
