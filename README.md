## leigqNEWTON v1.0.5 — 2026-10-03

**Numerical core: `LAA-R1-relative-2026-10-03`. Native MATLAB regression validation: 29/29 groups OK (0.98 s, 2026-10-03).**
The validation was run with the stand-alone public installation, with the author's personal quaternion toolbox removed from the MATLAB path. This corrected core uses the scale-invariant relative eigenpair residual for stopping
and acceptance and corrects nullspace and local-polisher defects. See
[`REVISION_NOTES.md`](REVISION_NOTES.md) for the changed tolerance semantics,
compatibility limits and the small regression driver:

```matlab
report = test_leigqNEWTON_revision;
save('leigqNEWTON_revision_report.mat','report');
assert(strcmp(report.status,'OK'),'Revision regression failed.');
```

The corrected core was published as GitHub release `v1.0.5` on 2026-10-03 and archived on Zenodo:
**[10.5281/zenodo.23124184](https://doi.org/10.5281/zenodo.23124184)**.
The release tag and archived snapshot identify the numerical code; subsequent updates
of this README and the citation files only record the assigned DOI and release metadata.
Historical benchmark outputs and prebuilt HTML/PDF documentation have not been regenerated.
Current MATLAB help and `REVISION_NOTES.md` govern the changed functions.
The `private` folder is required; add only the package root to the path.

---

## leigqNEWTON (MATLAB) — public toolbox

**leigqNEWTON** is a stand-alone MATLAB toolbox for computing and refining **left eigenpairs** of quaternion matrices using a Newton-type solver, with **residual certificates** and experimental **sphere** diagnostics.

- **Version:** 1.0.5 (2026-10-03)
- **Zenodo (archival DOI, v1.0.5):** https://doi.org/10.5281/zenodo.23124184
- **Historical first release (v1.0):** https://doi.org/10.5281/zenodo.18410141
- **Requirements:** MATLAB with the `quaternion` class available
- **Funding:** This work was co-funded by the European Union under the project ROBOPROX (reg. no. CZ.02.01.01/00/22_008/0004590).

---

## Install

Use the folder containing `leigqNEWTON.m` and `Contents.m` as the package root.
The optional distribution ZIP uses the folder name `leigqNEWTON_public`; GitHub's
automatically generated source archives may use a repository/tag-based folder name.
No additional nested package folder is needed. Keep the `private`, `examples`,
`docs`, and `LAA_Zoo` subfolders in their supplied locations.
If your MATLAB startup already adds the package root, no path change is required.

1. Download / unzip this toolbox folder (the folder that contains `leigqNEWTON.m` and `Contents.m`).
2. In MATLAB, add the toolbox folder to the path:

```matlab
toolboxRoot = "path/to/leigqNEWTON_public";
addpath(toolboxRoot);                       % recommended (no genpath)
addpath(fullfile(toolboxRoot,"examples"));  % optional: example scripts
rehash toolboxcache
```

**Why not `genpath`?** It can accidentally add unrelated folders and cause shadowing/path chaos.  
This toolbox is intentionally small; adding just the root (and optionally `examples`) is safest.

---

## Quick start (copy/paste)

Copy/paste this whole block after you set `toolboxRoot`:

```matlab
% 1) Release sanity gate (recommended before you start using the toolbox)
R = PACKAGEVerifierNEWTON_fromContents('MetaMode','strict');

% 2) Quick smoke test (solves + certifies + prints a compact report)
A = quaternion([1 2; 3 4], [0 1; 0 0], [0 0; 1 0], [0 0; 0 1]);  % 2x2 test
out = checkNEWTON(A);
```

If everything is installed correctly, you should see a short report with residuals and recommendations.
## Reproduce paper / supplement examples

Run scripts from the `examples/` folder:

```matlab
root = fileparts(which("leigqNEWTON"));   % robust toolbox root
run(fullfile(root,"examples","ExNEWTON_1_HuangSo.m"));
run(fullfile(root,"examples","ExNEWTON_2_MVPS.m"));
run(fullfile(root,"examples","ExNEWTON_3.m"));
```

**Notes about the examples:**
- MATLAB’s `quaternion` does not implement some matrix operations (`mtimes`, `abs`, `rank`).  
  The example scripts therefore use helper functions shipped here (e.g., `qmtimesNEWTON`, `qmldivideNEWTON`) when needed.
- One MVPS subcase is currently skipped because it contains `NaN` placeholders (see the comment inside the script).

---

## Using the solver

Minimal use:

```matlab
[lambda,V,res,info,lambdaU,VU,resU] = leigqNEWTON(A, "Num", 50);
```

- `lambda,V,res` are the accepted “hits” (may contain duplicates).
- `lambdaU,VU,resU` are a **distinct representative set** clustered using a tolerance.
- `info` contains summary details and, optionally, per‑trial logs depending on options.

Start here:
- `help leigqNEWTON`
- `help checkNEWTON`
- `doc leigqNEWTON` (after you publish the docs)

---

## Documentation

Prebuilt documentation is provided in:
- `docs/html/index.html`
- `docs/pdf/`

The folder `docs/source` contains the MATLAB documentation sources used to produce these pages. (Live Script `.mlx` sources and export utilities are maintainer tools and may be kept outside this public distribution.)

---

## Troubleshooting

### “Undefined function or variable 'quaternion'”
You likely do not have the Aerospace Toolbox installed/enabled. Verify:

```matlab
try, quaternion(0,0,0,0); disp("quaternion OK");
catch ME, disp(ME.message);
end
```

### Path shadowing / duplicates
If MATLAB finds an unexpected function version:

```matlab
which -all leigqNEWTON
which -all checkNEWTON
```

If needed (hard reset):

```matlab
restoredefaultpath
rehash toolboxcache
addpath(toolboxRoot);
```

---

## License

See the file **LICENSE** in this folder.

---

## How to cite

If you use this toolbox in academic work, please cite the accompanying paper/preprint.

- `CITATION.cff` (for GitHub/Zenodo)
- `CITATION.bib` (BibTeX)

For the corrected public core, cite **version 1.0.5**, released **2026-10-03**,
with DOI **10.5281/zenodo.23124184**. The older DOI **10.5281/zenodo.18410141**
identifies the first release, not the corrected core. The citation files on `main`
have been updated after archiving to record the newly assigned version DOI;
this metadata-only update does not change the numerical code or the release tag.

---

## Optional: LAA_Zoo companion bundle (paper/supplement demos)

This distribution may include an optional folder `LAA_Zoo` containing reproduction scripts for selected
examples used in the paper and the supplement.

Recommended way to run (no genpath):

```matlab
root = fileparts(which("leigqNEWTON"));
addpath(root); addpath(fullfile(root,"examples")); rehash toolboxcache
ExNEWTON_0_LAA_ZooLauncher
```

Notes:
- `LAA_Zoo_Ex_4_6_computation` uses an optional sphere-refit step via `fitSPHEREfromLambdas`
  (a public, self-contained implementation is shipped inside `LAA_Zoo/`).
