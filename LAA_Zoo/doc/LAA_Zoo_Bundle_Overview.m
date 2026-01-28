%% Zoo bundle overview and quick-start runnable walkthrough.
% (LAA_Zoo_Bundle_Overview.m)
%
% This script is part of the "Zoo" demo bundle accompanying the paper on
% Newton-based computation of left eigenvalues of quaternion matrices.
% The bundle enables running each of the Zoo examples.
%
% Where are these examples in the paper?
%  - Supplement B, Section 2: more than n distinct isolated eigenvalues (Examples 1 and 2)
%  - Supplement B, Section 2: less than n distinct isolated eigenvalues (Example 3)
%  - Supplement B, Section 4: spherical components (Examples 5, 6, ...)
%
% What do the demo scripts do?
%   Computation scripts (LAA_Zoo_Ex_*_computation.m) typically:
%   1) define the matrix A from the paper,
%   2) run a multi-start Newton solver to collect candidate left eigenvalues,
%   3) filter duplicates (distinctness tolerance) and report K distinct candidates,
%   4) (optional) refine/polish candidates and report improved certificates,
%   5) for spherical cases: detect sphere models and report centers/radii/inliers,
%   6) print a short "workspace guide" listing key output variables.
%
% Evidence scripts (LAA_Zoo_Ex_*_evidence.m) are complementary to the randomly generated exmaples:
%  - they run a larger/batched multi-start search to see whether NEW distinct
%       isolated eigenvalues appear at the stated tolerances,
%  - and/or they probe off-axis points on a hypothesized eigen-sphere to
%       confirm that the reported eigenvalues are isolated (not spherical).
%
% Requirements
%  - The leigqNEWTON toolbox must be on the MATLAB path.
%  - Quaternion support as used by your toolbox (typically MATLAB's quaternion).
%
% Author/citation
%   This demo bundle accompanies the leigqNEWTON reference implementation and the
%   associated manuscript (see the paper's "Data and code availability").
%
% Authors: M. Sebek (CVUT Prague) + ChatGPT, 24-Jan-2026
%
% -------------------------------------------------------------------------

%% Setup: paths and reproducibility
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
addpath(thisDir);

rng(24680,'twister');

if exist('leigqNEWTON','file') ~= 2
    error('Cannot find leigqNEWTON on the MATLAB path.\n');
end

%% Bundle contents (Zoo scripts)
% Computations (Supplement B examples):
%   - LAA_Zoo_Ex_2_1_computation   : 3x3, five distinct isolated eigenvalues
%   - LAA_Zoo_Ex_2_2_computation   : (as in Supplement B)
%   - LAA_Zoo_Ex_3_3_computation   : (as in Supplement B)
%   - LAA_Zoo_Ex_4_2_computation   : 4x4, spherical eigenvalues
%   - LAA_Zoo_Ex_4_5_computation   : 4x4, one sphere + isolated points
%   - LAA_Zoo_Ex_4_6_computation   : 6x6, two spherical components
%
% Evidence (Supplement B):
%   - LAA_Zoo_Ex_2_1_evidence      : evidence for "no extra isolated eigs" + no sphere
%   - LAA_Zoo_Ex_2_2_evidence      : analogous evidence run for Example 2-2

%% What the scripts print: 
% Distinct eigenvalue cases (Supplement B, Section 2 and 3):
%   - distinct eigenvalues as returned by solver (before refinement) 
%     with certificates
%   - distinct eigenvalues after refine with certificates 
%   - tips, where to find results in the workspace
%
% Spherical cases (Supplement B, Section 4)
% Stage 1 (sampling + detection):
%   - collect K distinct candidates,
%   - cluster and fit sphere models in R^4,
%   - report centers/radii/inliers and isolated/outlier candidates.
%
% Stage 2 (refinement + optional sphere re-fit):
%   - refine to improve eigenvalue certificates,
%   - optionally re-fit spheres using refined candidates (diagnostic consistency check).
%   - tips, where to find results in the workspace

%% Quick-start runnable walkthrough: Example 2-1 (fastest)
fprintf('\n=== Quick-start: running LAA_Zoo_Ex_2_1_computation ===\n');
LAA_Zoo_Ex_2_1_computation;

% The computation script should print a short "workspace guide" at the end.
% Below we additionally *show* a few key variables if they exist.

varsWanted = {'A','lamDistinct','resMinDistinct','lamRef','resMinRef'};
for ii = 1:numel(varsWanted)
    if evalin('base', sprintf("exist('%s','var')", varsWanted{ii}))
        fprintf('  (workspace) %s is available.\n', varsWanted{ii});
    end
end

if evalin('base',"exist('A','var')")
    A = evalin('base','A');
    fprintf('\nA (Example 2-1 matrix):\n');
    disp(A);
end

if evalin('base',"exist('lamDistinct','var') && exist('resMinDistinct','var')")
    lamDistinct = evalin('base','lamDistinct');
    resMinDistinct = evalin('base','resMinDistinct');
    fprintf('\nDistinct candidates (pre-polish):\n');
    disp(table(lamDistinct(:), resMinDistinct(:), 'VariableNames', {'lamDistinct','resMinDistinct'}));
end

if evalin('base',"exist('lamRef','var') && exist('resMinRef','var')")
    lamRef = evalin('base','lamRef');
    resMinRef = evalin('base','resMinRef');
    fprintf('\nRefined/polished candidates:\n');
    disp(table(lamRef(:), resMinRef(:), 'VariableNames', {'lambdaRef','resMinRef'}));
end

fprintf('\n=== Quick-start: running LAA_Zoo_Ex_2_1_evidence ===\n');
LAA_Zoo_Ex_2_1_evidence;
