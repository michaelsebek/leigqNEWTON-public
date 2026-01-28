%% leigqNEWTON Toolbox — Documentation index
% This toolbox computes quaternionic **LEFT** eigenpairs
%
% \[
%   A v = \lambda v,\qquad v\neq 0,
% \]
% where \(\lambda\) is a quaternion acting on the **left**.
%
% The package is **stand-alone**: it relies only on MATLAB + the built-in
% |quaternion| class.
%
% Tip: If you prefer Live Scripts, open any file in |docs/source| and convert it via
% *Save As → Live Script (.mlx)*.

%% Quickstart checklist
% 1) Put the toolbox *root folder* on the MATLAB path.
%    (The root is the folder that contains |leigqNEWTON.m|.)
%
%   toolboxRoot = '.../leigqNEWTON_public';
%   addpath(toolboxRoot);
%
% Avoid |genpath| unless you know exactly what is already on your path.
% (Using genpath can accidentally pull in old copies of similarly-named files.)
%
% 2) Run a one-call diagnostic workflow:
%
%   out = checkNEWTON(A);

%% What is where
% * Core solver: |leigqNEWTON|
% * Diagnostic wrapper: |checkNEWTON|
% * Refinement: |leigqNEWTON_refine_*|
% * Certificates: |leigqNEWTON_cert_*|
% * Sphere tooling (experimental): |leigqNEWTON_sphere_*|
% * Quaternion utilities: |qcleanNEWTON|, |qroundNEWTON|, |qmtimesNEWTON|,
%   |qmldivideNEWTON|, |qmrdivideNEWTON|

%% Reproducible paper examples
% The documentation pages intentionally reuse fixed matrices from the shipped scripts:
% * |examples/ExNEWTON_1_HuangSo.m| (Huang–So, explicit 2×2 cases; includes a sphere case)
% * |examples/ExNEWTON_2_MVPS.m|   (MVPS, explicit 3×3 cases)
% * |examples/ExNEWTON_3.m|        (additional small example)

%% Recommended reading order
% # |doc_GettingStarted| — add-to-path, smoke test, and the basic workflow
% # |doc_leigqNEWTON| — solver API + profiles + key options
% # |doc_SolverAndOptions| — option patterns and practical recipes
% # |doc_checkNEWTON| — one-call report workflow
% # |doc_RefinementAndCertificates| — refine + certify (recommended production workflow)
% # |doc_SphereHunting| — sphere sampling/refinement (experimental)
% # |doc_examples_HuangSo|, |doc_examples_MVPS|, |doc_examples_PanNg| — reproducible examples

%   M. Šebek (CTU Prague), 2026-01-24