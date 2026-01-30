% LAA Zoo demo scripts
%
% This folder contains small, self-contained MATLAB scripts used to reproduce
% the numerical results reported for selected "Zoo" examples in the paper
% and in Supplement B (left eigenvalues of quaternion matrices).
%
% Computation demos (Supplement B examples):
%   LAA_Zoo_Ex_2_1_computation - computation companion for Supplement B (Section 2 Example 1).
%   LAA_Zoo_Ex_2_2_computation - computation companion for Supplement B (Section 2 Example 2).
%   LAA_Zoo_Ex_3_3_computation - computation companion for Supplement B (Section 3 Example 3).
%   LAA_Zoo_Ex_4_2_computation - computation companion for Supplement B (Section 4 Example 2).
%   LAA_Zoo_Ex_4_5_computation - computation companion for Supplement B (Section 4 Example 5).
%   LAA_Zoo_Ex_4_6_computation - computation companion for Supplement B (Section 4 Example 6).
%
% Evidence demos (optional, for selected examples):
%   LAA_Zoo_Ex_2_1_evidence - evidence/diagnostics companion for Supplement B (Section 2 Example 1).
%   LAA_Zoo_Ex_2_2_evidence - evidence/diagnostics companion for Supplement B (Section 2 Example 2).
%
% Each script prints a short run log and leaves key variables in the workspace.
% See the header of each script for the exact mapping to the paper/Supplement.
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