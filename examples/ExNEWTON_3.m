%% ExNEWTON_3.m  (NEWTON toolbox examples — updated)
% This script collects explicit quaternion-matrix examples and runs the leigqNEWTON toolbox
% (solver + refinement + certificates + optional sphere diagnostics).
%
% IMPORTANT
%   * Matrices are specified explicitly and should not be modified (paper reproducibility).
%   * Requires MATLAB's built-in quaternion class (Aerospace Toolbox).
%   * Requires the leigqNEWTON toolbox on the MATLAB path.
%
% This version replaces any legacy calls (qleig/leigqVANILA/inspectSPHERE/testSPHERE/...)
% by calls to NEWTON toolbox functions only.

clc;

%% 0) Sanity checks
% Quaternion class
try
    quaternion(0,0,0,0);
catch ME
    error('MATLAB quaternion class not available. Install/enable Aerospace Toolbox.\n%s', ME.message);
end

% NEWTON toolbox on path
if exist('leigqNEWTON','file') ~= 2
    error('leigqNEWTON not found on MATLAB path. Add the toolbox, e.g. addpath(genpath(<toolbox_root>)).');
end
if exist('checkNEWTON','file') ~= 2
    error('checkNEWTON not found on MATLAB path. Please update to the current leigqNEWTON toolbox.');
end
if exist('qmtimesNEWTON','file') ~= 2
    error('qmtimesNEWTON not found on MATLAB path. Please update the toolbox (quaternion matrix multiply helper).');
end

%% 1) Pan–Ng (2024) dense 4x4 quaternion circulant matrix (Example for the paper)
% Source: Pan & Ng (2024), SIAM J. Matrix Anal. Appl., Example 1.
%
% A =
%  [a b c d;
%   d a b c;
%   c d a b;
%   b c d a]
%
% where
%   a = -2 + 1*i + 1*j + 4*k
%   b =  2 + 4*i + 1*j + 1*k
%   c =  1 + 3*i + 2*j + 2*k
%   d = -1 + 2*i + 2*j + 3*k

q = @(w,x,y,z) quaternion(w,x,y,z);

a = q(-2, 1, 1, 4);
b = q( 2, 4, 1, 1);
c = q( 1, 3, 2, 2);
d = q(-1, 2, 2, 3);

A = [a b c d;
     d a b c;
     c d a b;
     b c d a];

fprintf('=== Pan--Ng dense 4x4 quaternion circulant matrix ===\n');

%% 2) Solve + report (uses leigqNEWTON + certificates; computes vectors if needed)
% Reliable profile is usually helpful for 4x4:
out = checkNEWTON(A, [], [], ...
    'SolveArgs',  {'SolveProfile','reliable','Seed',1}, ...
    'SphereCheck','auto');

%% 3) Optional: refine/polish the returned eigenvalues, then re-check
if ~isempty(out.lambda)
    fprintf('\nRefinement (batch):\n');
    [lamR, VR, cert] = leigqNEWTON_refine_batch(A, out.lambda, 'DoPolish', true);

    fprintf('  median(resMin)=%.3e, max(resMin)=%.3e\n', median(cert.resMin), max(cert.resMin));

    % Re-check with refined values:
    checkNEWTON(A, lamR, VR, 'SphereCheck','off');
end

%% 4) LaTeX snippet (copy/paste into the paper)
% The following snippet is kept verbatim to match the draft text.
latexLines = [
"\\begin{Ex}[A dense $4\\times4$ quaternion matrix]\\label{ex:PanNg_dense_4x4}"
"Consider the dense quaternion circulant matrix $A\\in\\HH^{4\\times 4}$:"
"\\["
"\\setlength{\\arraycolsep}{2pt}"
"\\begin{bmatrix}"
" -2 + \\i + \\j + 4\\k & 2 + 4\\i + \\j + \\k & 1 + 3\\i + 2\\j + 2\\k & -1 + 2\\i + 2\\j + 3\\k \\"
" -1 + 2\\i + 2\\j + 3\\k & -2 + \\i + \\j + 4\\k & 2 + 4\\i + \\j + \\k & 1 + 3\\i + 2\\j + 2\\k \\"
" 1 + 3\\i + 2\\j + 2\\k & -1 + 2\\i + 2\\j + 3\\k & -2 + \\i + \\j + 4\\k & 2 + 4\\i + \\j + \\k \\"
" 2 + 4\\i + \\j + \\k & 1 + 3\\i + 2\\j + 2\\k & -1 + 2\\i + 2\\j + 3\\k & -2 + \\i + \\j + 4\\k"
"\\end{bmatrix}."
"\\]"
"This matrix is taken from \\cite[Example~1]{PanNg2024}, where it is used to illustrate the"
"block-diagonalization of quaternion circulant matrices; the authors do not compute or report its"
"\\emph{left} eigenvalues."
""
"Using our implementation, we obtained the following (rounded) left eigenvalues:"
"\\["
"\\widehat\\lambda \\approx"
"\\begin{bmatrix}"
" -1.5 - 0.43\\i + 1.5\\j + 4.6\\k \\"
" -2 - 2\\i + 2\\k \\"
" -2.9 + 0.85\\i - 1.2\\j + 5.1\\k \\"
" -5.2 + 1.5\\i - 0.25\\j + 1.9\\k"
"\\end{bmatrix}."
"\\]"
"The corresponding residual norms are around $10^{-16}$, confirming near machine-precision accuracy"
"for the computed eigenpairs."
"\\end{Ex}"
""
"\\bibitem{PanNg2024}"
"J.~Pan and M.~K.~Ng,"
"\\newblock Block-Diagonalization of Quaternion Circulant Matrices with Applications,"
"\\newblock \\emph{SIAM J. Matrix Anal. Appl.} 45(3):1429--1454, 2024."
"\\newblock \\DOI{10.1137/23M1552115}."
];
latex_Ex_PanNg_dense_4x4 = strjoin(latexLines, newline);
disp(latex_Ex_PanNg_dense_4x4);

fprintf('\nDone.\n');

