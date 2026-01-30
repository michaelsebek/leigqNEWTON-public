% leigqNEWTON  Public MATLAB toolbox for left eigenpairs of quaternion matrices.
% Version 1 2026-01-28
%
% Overview
%   leigqNEWTON computes (and refines) left eigenpairs (lambda,v) of quaternion matrices A,
%   i.e., A*v = lambda*v, using a practical Newton-type method with certificate-style
%   residual checks. The toolbox also provides refinement routines and utilities tailored
%   to MATLAB's built-in quaternion class.
%
% Requirements
%   MATLAB with the built-in quaternion class (Aerospace Toolbox).
%
% Getting started
%   1) Add the toolbox to the path (avoid genpath to prevent shadowing):
%        root = "path/to/leigqNEWTON_public";
%        addpath(root); addpath(fullfile(root,"examples")); rehash toolboxcache
%   2) Run a quick smoke test:
%        out = checkNEWTON();   % or: out = checkNEWTON(A)
%   3) Reproduce selected paper examples:
%        root = fileparts(which("leigqNEWTON"));
%        run(fullfile(root,"examples","ExNEWTON_1_HuangSo.m"));
%        run(fullfile(root,"examples","ExNEWTON_2_MVPS.m"));
%
% Documentation
%   docs/source contains documentation source scripts intended for publish/Live Script.
%
% License and citation
%   See LICENSE and README.md in the toolbox root.
%   See CITATION.cff and CITATION.bib for citation metadata.
%
% Author: Michael Sebek (michael.sebek@fel.cvut.cz)
%   Version: v1.0
%
% Note:
%  This function is part of the public MATLAB toolbox leigqNEWTON accompanying the paper:
%     M. Sebek, "Computing Left Eigenvalues of Quaternion Matrices", submitted to
%     Linear Algebra and its Applications, 2026.
%  If you use this software in academic work, please cite the paper and/or the Zenodo archive:
%     https://doi.org/10.5281/zenodo.18410141
%
% Toolbox functions

%   leigqNEWTON                    - Compute left eigenpairs of quaternion matrices using Newton multi-start sampling.
%   checkNEWTON                    - Run a quick smoke test / verifier and print a compact report.
%   leigqNEWTON_quickstart         - Quick-start convenience wrapper to run a standard solver pipeline.
%   leigqNEWTON_init_vec           - Construct initial vectors/guesses used by the Newton iteration.
%   leigqNEWTON_classify_multiplicity - Classify distinctness groups / multiplicity heuristics (diagnostic).
%   leigqNEWTON_cert_resMin        - Compute the eigenvalue-only residual certificate r_min(lambda).
%   leigqNEWTON_cert_resPair       - Compute the eigenpair residual certificate for (lambda,v).
%   leigqNEWTON_refine_auto        - Automatic refinement/polish pipeline for a set of candidates.
%   leigqNEWTON_refine_batch       - Batch refinement for multiple candidate eigenvalues/eigenpairs.
%   leigqNEWTON_refine_lambda      - Refine a single candidate eigenvalue (Newton/LS polishing).
%   leigqNEWTON_refine_polish      - Polish candidates to strengthen certificates and consistency checks.
%   leigqNEWTON_sphere_sample      - Sample candidate eigenvalues for sphere detection (experimental).
%   leigqNEWTON_sphere_detect      - Detect spherical components from sampled candidates (experimental).
%   leigqNEWTON_sphere_refine      - Refine fitted sphere parameters (center/radius/inliers) (experimental).
%   leigqNEWTON_sphere_validate    - Validate a sphere hypothesis (diagnostic checks) (experimental).
%   qcleanNEWTON                   - Clean quaternion arrays (trim tiny components, sanitize numerical noise).
%   qroundNEWTON                   - Round quaternion arrays to a given precision/tolerance.
%   qmtimesNEWTON                  - Quaternion matrix multiplication helper (via embeddings).
%   qmldivideNEWTON                - Left matrix division helper for quaternion matrices (via embeddings).
%   qmrdivideNEWTON                - Right matrix division helper for quaternion matrices (via embeddings).
%
% See also quaternion
