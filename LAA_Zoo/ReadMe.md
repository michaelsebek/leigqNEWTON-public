# LAA Zoo scripts (paper + Supplement B reproduction)

This folder contains MATLAB **demo scripts** that reproduce the numerical results reported for selected “Zoo” examples in the LAA paper and **Supplement B** (left eigenvalues of quaternion matrices computed by `leigqNEWTON`).

## What is included

### Computation demos (Supplement B examples)
- `LAA_Zoo_Ex_2_1_computation` — computation companion for Supplement B (Section 2 Example 1).
- `LAA_Zoo_Ex_2_2_computation` — computation companion for Supplement B (Section 2 Example 2).
- `LAA_Zoo_Ex_3_3_computation` — computation companion for Supplement B (Section 3 Example 3).
- `LAA_Zoo_Ex_4_2_computation` — computation companion for Supplement B (Section 4 Example 2).
- `LAA_Zoo_Ex_4_5_computation` — computation companion for Supplement B (Section 4 Example 5).
- `LAA_Zoo_Ex_4_6_computation` — computation companion for Supplement B (Section 4 Example 6).

### Evidence demos (optional)
- `LAA_Zoo_Ex_2_1_evidence` — evidence/diagnostics companion for Supplement B (Section 2 Example 1).
- `LAA_Zoo_Ex_2_2_evidence` — evidence/diagnostics companion for Supplement B (Section 2 Example 2).

## Requirements

- MATLAB with a quaternion datatype available (these scripts expect the `quaternion` class).
- The **`leigqNEWTON` toolbox** (and its dependencies) must be on the MATLAB path.
- For the sphere fits used in some demos: `fitSPHEREfromLambdas.m` must be on the MATLAB path.

## How to run

1. Add the `leigqNEWTON` package folder (and this folder) to the MATLAB path.
2. Run a script, for example:
   ```matlab
   run('LAA_Zoo_Ex_4_6_computation.m');
   ```

Each script:
- prints a compact run log,
- produces a set of variables in the base workspace (listed at the end of the script output),
- is deterministic up to the internal randomness controlled by the displayed seeds.

## Notes on “refined spheres”

Some examples contain spherical components in the left spectrum. A post-processing “refinement” step typically improves **residual certificates** for individual candidates. The refined candidates are then optionally **re-fit** to a sphere model (using a geometric sphere fit), which can be viewed as projecting the refined points back onto an idealized spherical component.

## Citation / provenance

These scripts accompany the LAA submission and Supplement B. Please cite the paper when using any of the examples or results.

**Authors:** M. Šebek (CTU Prague) + ChatGPT  
**Last updated:** 2026-01-23
