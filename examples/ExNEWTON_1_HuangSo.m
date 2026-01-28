%% ExNEWTON_1_HuangSo.m (UPDATED for leigqNEWTON toolbox)
% Test cases from:
%   L. Huang and W. So, "On left eigenvalues of a quaternionic matrix,"
%   Linear Algebra and its Applications 323 (2001) 105--116.
%   DOI: 10.1016/S0024-3795(00)00246-9
%
% This script reproduces the paper's explicit 2x2 examples (Examples 2.5--2.7)
% using MATLAB's built-in quaternion class and the leigqNEWTON toolbox.
%
% Requirements
%   - MATLAB class quaternion (Aerospace Toolbox).
%   - leigqNEWTON toolbox on the MATLAB path.
%
% Notes
%   - MATLAB quaternion does NOT implement abs(...) or matrix rank directly.
%     This script relies on checkNEWTON(...) and the sphere sampling utilities
%     from the leigqNEWTON toolbox.
%
% -------------------------------------------------------------------------
clc;

% --- Sanity checks --------------------------------------------------------
try
    quaternion(0,0,0,0);
catch
    error(['This script requires MATLAB''s built-in quaternion class. ' ...
           'Ensure Aerospace Toolbox is installed and quaternion(w,x,y,z) works.']);
end

if exist('leigqNEWTON','file') ~= 2
    error(['leigqNEWTON not found on the MATLAB path. ' ...
           'Add the leigqNEWTON toolbox folder, e.g. addpath(genpath(<toolboxRoot>)).']);
end

% Pretty-print tolerance (component-wise zeroing).
tolShow = 1e-12;

% Quaternion units and constants.
q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);
qj = quaternion(0,0,1,0);
qk = quaternion(0,0,0,1);

fprintf('\nHuang--So (LAA 323, 2001) explicit 2x2 examples\n');
fprintf('DOI: 10.1016/S0024-3795(00)00246-9\n');
fprintf('=============================================================\n');

%% Example 2.5
% A = [ 0, 1+i ; 1-i, 0 ]
% sigma_l(A) = { +sqrt(2), -sqrt(2) }
A25 = [ q0, q1 + qi;
        q1 - qi, q0 ];
lam25_true = [qreal(sqrt(2)); qreal(-sqrt(2))];

fprintf('\nEXAMPLE 2.5\n');
disp('A25 ='); disp(qcleanNEWTON(A25, tolShow));
disp('Expected left eigenvalues (paper) ='); disp(qcleanNEWTON(lam25_true, tolShow));

out25 = checkNEWTON(A25);
local_compare_to_true(out25.lambdac, lam25_true, tolShow);

%% Example 2.6
% A = [ 0, i ; j, 1 ]
% sigma_l(A) = { (1+i+j-k)/2 , (1-i-j-k)/2 }
A26 = [ q0, qi;
        qj, q1 ];
lam26_true = [ (q1 + qi + qj - qk)/2;
               (q1 - qi - qj - qk)/2 ];

fprintf('\nEXAMPLE 2.6\n');
disp('A26 ='); disp(qcleanNEWTON(A26, tolShow));
disp('Expected left eigenvalues (paper) ='); disp(qcleanNEWTON(lam26_true, tolShow));

out26 = checkNEWTON(A26);
local_compare_to_true(out26.lambdac, lam26_true, tolShow);

%% Example 2.7 (spherical spectrum)
% A = [ 2, i ; -i, 2 ]
% The paper states that sigma_l(A) is infinite (a 2-sphere family).
A27 = [ qreal(2),  qi;
       -qi,        qreal(2) ];

fprintf('\nEXAMPLE 2.7\n');
disp('A27 ='); disp(qcleanNEWTON(A27, tolShow));
disp('Paper claim: sigma_l(A27) is infinite (a sphere family).');

% First run the solver as usual (it will return a finite set of samples).
out27 = checkNEWTON(A27, [], [], 'SphereCheck','off');

% Then explicitly run sphere sampling + refinement (NEWTON toolbox).
% Increase Collect if you want stronger evidence.
Collect = 20; % batches (not raw lambda count)
[lamAll, lamS, ~, ~, sph0, infoS] = leigqNEWTON_sphere_sample(A27, ...
    'Collect', Collect, 'Seed0', 1, 'Report', 'progress');
C27 = struct('A',A27,'lamAll',lamAll,'lamSamples',lamS,'sph',sph0,'info',infoS);
%
fprintf('\nPlease note: The refinement may take about ten minutes.\n\n');
%
[A2, lamAll2, resAll, lamS2, resS, sph2, conf, outS] = leigqNEWTON_sphere_validate(C27, 'Verbose', 1);

fprintf('\nSphere report (from leigqNEWTON_sphere_refine):\n');
disp(sph2);
disp(conf);

fprintf('\nDone.\n');

%% ------------------------------------------------------------------------
%% Local helpers

function q = qreal(x)
%QREAL Make a real quaternion scalar/array from a real numeric array.
q = quaternion(double(x), 0*double(x), 0*double(x), 0*double(x));
end

function local_compare_to_true(lamFound, lamTrue, tolShow)
%LOCAL_COMPARE_TO_TRUE Compare a found set to a reference set (2x2 examples).
%
% For these 2x2 examples the reference set size is small.
% We compute a simple nearest-neighbor matching in 4D Euclidean distance.

lamFound = lamFound(:);
lamTrue  = lamTrue(:);

if isempty(lamFound)
    fprintf('  No eigenvalues returned.\n');
    return;
end

D = zeros(numel(lamFound), numel(lamTrue));
for a = 1:numel(lamFound)
    for b = 1:numel(lamTrue)
        D(a,b) = local_qdist(lamFound(a), lamTrue(b));
    end
end

fprintf('  Nearest distances to the reference set (in 4D):\n');
for a = 1:numel(lamFound)
    [dmin, j] = min(D(a,:));
    fprintf('    lambda[%d] ~ ref[%d], dist = %.3e, lambda = ', a, j, dmin);
    disp(qcleanNEWTON(lamFound(a), tolShow));
end
end

function d = local_qdist(q1, q2)
%LOCAL_QDIST Euclidean distance between two quaternions in R^4.
[w1,x1,y1,z1] = parts(q1);
[w2,x2,y2,z2] = parts(q2);
d = sqrt((w1-w2).^2 + (x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2);
end