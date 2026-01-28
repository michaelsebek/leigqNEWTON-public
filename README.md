# leigqNEWTON (MATLAB) — public toolbox

**leigqNEWTON** is a stand‑alone MATLAB toolbox for computing and refining **left eigenpairs** of quaternion matrices using a Newton‑type solver, with **residual certificates** and experimental **sphere** diagnostics.

- **Version:** 1 (2026-01-28)
- **Requirements:** MATLAB with the built‑in `quaternion` class (Aerospace Toolbox)

---

## Install

1. Download / unzip this toolbox folder (the folder that contains `leigqNEWTON.m` and `Contents.m`).
2. In MATLAB, add the toolbox folder to the path:

```matlab
toolboxRoot = "path/to/leigqNEWTON-public";
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

Repository URL and DOI can be added later (e.g., when the GitHub repo and/or Zenodo DOI exist).

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

