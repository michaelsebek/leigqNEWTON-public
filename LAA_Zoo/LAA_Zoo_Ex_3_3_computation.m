%LAA_ZOO_EX_3_3_COMPUTATION  Zoo demo — computation companion for Supplement B (Section 3 Example 3).
%
%   - LaTeX block:  Supplement B, Section 3 Example 3
%   - Script:       LAA_Zoo_Ex_3_3_computation.m
%
% This script reproduces (up to numerical tolerances and randomized sampling)
% the numerical results reported for the corresponding Zoo example.
%
% Pipeline
%   Stage 1  Candidate sampling + clustering (leigqNEWTON_sphere_sample)
%   Stage 2  Refinement/certification (leigqNEWTON_refine_batch)
%   Stage 3  Optional sphere refit from refined candidates (fitSPHEREfromLambdas)
%
% How to run
%   >> LAA_Zoo_Ex_3_3_computation
%
% Outputs
%   Creates variables in the base workspace; see the final "Variables created" printout.
%
% Requirements
%   - leigqNEWTON toolbox on the MATLAB path (quaternion class + leigqNEWTON_* routines).
%
% Notes
%   - Option 'Tol' is supported as a shortcut alias for leigqNEWTON_refine_batch (sets TolX/TolFun/TolResPolish).
%     use full names ('TolX','TolFun','TolResPolish',...).
%
%   M. Šebek (CTU Prague) + ChatGPT, 2026-01-23

clear; clc;

%% (0) Basic dependency check
local_require_function('leigqNEWTON');
local_require_function('leigqNEWTON_cert_resMin');

%% (A) Construct the matrix A (component-wise integers)
% A = W + X*i + Y*j + Z*k
W = [  1,   8,  -8,  -7;
       0,   1,   0,   0;
       0,  -1,   2,   2;
       0,   8,  -8,  -6];
X = [  3,   1,   0,  -2;
      -1,   1,   0,   1;
      -1,   0,   1,   1;
       0,   0,   0,   1];
Y = [  0,   1,   0,   4;
      -1,  -2,   0,   1;
      -1,  -6,   4,   1;
       0,   0,   0,   4];
Z = [  1,  -4,   4,   1;
       0,   1,   0,   0;
       0,   3,  -2,  -1;
       0,  -4,   4,   2];

A = quaternion(W, X, Y, Z);
n = size(A,1);

%% (B) User-facing parameters
TOL.UniqueTol     = 1e-5;      % clustering tolerance in R^4 (important for this example)
TOL.TolNewton     = 1e-12;     % leigqNEWTON convergence tolerance
TOL.TargetResMin  = 1e-16;     % target certificate after refine

% Refiner tolerances (must use full names; DO NOT use 'Tol')
TOL.TolX          = 1e-16;
TOL.TolFun        = 1e-16;
TOL.MaxIter       = 600;
TOL.TolResPolish  = 1e-16;
TOL.MaxIterPolish = 60;

CFG.Runs          = 60;
CFG.Restarts      = 1500;
CFG.Krequest      = n;
CFG.baseSeed      = 12345;

fprintf('=== LAA demo: Zoo Example 3-3 (4x4, TWO distinct isolated left eigenvalues) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', TOL.UniqueTol, TOL.TolNewton, TOL.TargetResMin);
fprintf('Runs=%d | Restarts/run=%d | Krequest=%d\n\n', CFG.Runs, CFG.Restarts, CFG.Krequest);

%% (C) Multi-start campaign: collect accepted hits
lamAll = quaternion.empty(0,1);
hitsTotal = 0;

for r = 1:CFG.Runs
    seed = CFG.baseSeed + r;

    [lam, ~, info] = leigqNEWTON(A, ...
        'Num',       CFG.Krequest, ...
        'Restarts',  CFG.Restarts, ...
        'Seed',      seed, ...
        'Tol',       TOL.TolNewton, ...
        'UniqueTol', TOL.UniqueTol, ...
        'Verbose',   0);

    lamAcc = local_getAcceptedLambdas(lam, info);
    hitsTotal = hitsTotal + numel(lamAcc);
    lamAll(end+1:end+numel(lamAcc),1) = lamAcc(:); %#ok<AGROW>
end

fprintf('Collected accepted hits: %d\n', hitsTotal);

%% (D) Cluster to DISTINCT representatives (stable order)
lamAll = local_sort_by_invariants(lamAll);
lamDistinct = local_unique_quat_r4(lamAll, TOL.UniqueTol);

fprintf('Distinct eigenvalues (tol = %.1e): %d\n\n', TOL.UniqueTol, numel(lamDistinct));

% Certificates for the solver output (before refinement)
resMinDistinct = leigqNEWTON_cert_resMin(A, lamDistinct);

% Stable order for display
[lamDistinct, idx0] = local_sort_by_invariants(lamDistinct);
resMinDistinct = resMinDistinct(idx0);

fprintf('--- DISTINCT eigenvalues as returned by solver (before refinement) ---\n');
for k = 1:numel(lamDistinct)
    fprintf('%2d) %s   | r_min = %.3e\n', k, local_qstr(lamDistinct(k), 6), double(resMinDistinct(k)));
end
fprintf('\n');


%% (E) Refine + polish (recommended)

% Refinement stage (improves certificates)
if exist('leigqNEWTON_refine_batch','file') == 2
    refinePlanned = 'refine_batch';
elseif exist('leigqNEWTON_refine_auto','file') == 2
    refinePlanned = 'refine_auto';
else
    refinePlanned = 'no_refine';
end

fprintf('Refinement will now run (%s).\n', refinePlanned);
fprintf('Progress:\n');
tRef__ = tic;
[lamRef, resMinRef, refineNote] = local_refine_list(A, lamDistinct, TOL, CFG.baseSeed+999);
elapsedRef = toc(tRef__);
fprintf('Refinement finished (%s), elapsed %.2fs.\n\n', refineNote, elapsedRef);


% Stable order for display
[lamRef, idxSort] = local_sort_by_invariants(lamRef);
resMinRef = resMinRef(idxSort);

%% (F) Report
fprintf('--- DISTINCT eigenvalues after refine (%s) ---\n', refineNote);
for k = 1:numel(lamRef)
    fprintf('%2d) %s   | r_min = %.3e\n', k, local_qstr(lamRef(k), 6), double(resMinRef(k)));
end

fprintf('\nCertificate range r_min: [%.3e, %.3e]\n', min(resMinRef), max(resMinRef));
fprintf('Minimum pairwise separation (R^4 distance) between listed eigenvalues: %.6g\n', local_minsep_r4(lamRef));
fprintf('  (Diagnostic only: larger values mean the eigenvalues are well separated.)\n');

fprintf('\nWhere to find results in the workspace:\n');
fprintf('  A              : %dx%d quaternion matrix\n', n, n);
fprintf('  lamAll         : raw accepted hits (quaternion, %d-by-1)\n', numel(lamAll));
fprintf('  lamDistinct    : distinct eigenvalues BEFORE refinement (quaternion, %d-by-1)\n', numel(lamDistinct));
fprintf('  resMinDistinct : certificates for lamDistinct (double, %d-by-1)\n', numel(lamDistinct));
fprintf('  lamRef         : refined eigenvalues (quaternion, %d-by-1)\n', numel(lamRef));
fprintf('  resMinRef      : certificates for lamRef (double, %d-by-1)\n', numel(lamRef));
fprintf('  refineNote     : which refiner was used\n');
fprintf('=== End of demo ===\n');

%% ============================================================
% Local helpers (script-local functions)
%% ============================================================

function local_require_function(funName)
    if exist(funName,'file') ~= 2
        error('Required function not found on path: %s', funName);
    end
end

function lamAcc = local_getAcceptedLambdas(lam, info)
    lamAcc = quaternion.empty(0,1);
    if isstruct(info)
        if isfield(info,'lambdaAccepted') && ~isempty(info.lambdaAccepted)
            lamAcc = info.lambdaAccepted;
        elseif isfield(info,'lambda') && ~isempty(info.lambda)
            lamAcc = info.lambda;
        end
    end
    if isempty(lamAcc)
        lamAcc = lam;
    end
    lamAcc = lamAcc(:);
end

function [lamOut, resMinOut, note] = local_refine_list(A, lam0, TOL, seed)
% Refine a list of candidates and return certificates r_min.
%
% This helper prints brief progress so the user sees MATLAB is working.

    if nargin < 4, seed = 0; end
    if isempty(seed), seed = 0; end

    lam0 = lam0(:);
    K = numel(lam0);

    if K == 0
        lamOut = lam0;
        resMinOut = [];
        note = 'empty';
        return;
    end

    % Preferred: batch refiner (used here per-candidate to provide progress)
    if exist('leigqNEWTON_refine_batch','file') == 2
        note = 'refine_batch';
        lamOut = lam0;
        resMinOut = nan(K,1);

        for kk = 1:K
            fprintf('  [%2d/%2d] refining ... ', kk, K);
            drawnow;

            [lamk, ~, certk] = leigqNEWTON_refine_batch(A, lam0(kk), ...
                'Mode',          'auto', ...
                'TargetResMin',  TOL.TargetResMin, ...
                'TolX',          TOL.TolX, ...
                'TolFun',        TOL.TolFun, ...
                'MaxIter',       TOL.MaxIter, ...
                'DoPolish',      true, ...
                'TolResPolish',  TOL.TolResPolish, ...
                'MaxIterPolish', TOL.MaxIterPolish, ...
                'Seed',          seed + kk, ...
                'Verbose',       0);

            lamOut(kk) = lamk;

            if isstruct(certk) && isfield(certk,'resMin') && ~isempty(certk.resMin)
                resMinOut(kk) = double(certk.resMin(1));
            else
                resMinOut(kk) = double(leigqNEWTON_cert_resMin(A, lamk));
            end

            fprintf('done (r_min=%.3e)\n', resMinOut(kk));
        end
        return;
    end

    % Fallback: auto refiner
    if exist('leigqNEWTON_refine_auto','file') == 2
        note = 'refine_auto';
        lamOut = lam0;
        resMinOut = nan(K,1);

        for kk = 1:K
            fprintf('  [%2d/%2d] refining ... ', kk, K);
            drawnow;

            [lamBest, lamPol] = leigqNEWTON_refine_auto(A, lam0(kk), ...
                'Mode',          'auto', ...
                'TargetResMin',  TOL.TargetResMin, ...
                'TolX',          TOL.TolX, ...
                'TolFun',        TOL.TolFun, ...
                'MaxIter',       TOL.MaxIter, ...
                'DoPolish',      true, ...
                'TolResPolish',  TOL.TolResPolish, ...
                'MaxIterPolish', TOL.MaxIterPolish, ...
                'Seed',          seed + kk, ...
                'Verbose',       0);

            if ~isempty(lamPol)
                lamk = lamPol;
            else
                lamk = lamBest;
            end

            lamOut(kk) = lamk;
            resMinOut(kk) = double(leigqNEWTON_cert_resMin(A, lamk));

            fprintf('done (r_min=%.3e)\n', resMinOut(kk));
        end
        return;
    end

    % Last resort: certify the provided values
    lamOut = lam0;
    resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
    note = 'no_refine';
end
function lamU = local_unique_quat_r4(lam, tol)
% Greedy clustering in R^4.
    lamU = quaternion.empty(0,1);
    for i = 1:numel(lam)
        q = lam(i);
        if isempty(lamU)
            lamU(1,1) = q; %#ok<AGROW>
            continue;
        end
        v = local_qvec(q);
        isNew = true;
        for j = 1:numel(lamU)
            if norm(v - local_qvec(lamU(j)), 2) <= tol
                isNew = false;
                break;
            end
        end
        if isNew
            lamU(end+1,1) = q; %#ok<AGROW>
        end
    end
end

function [lamSorted, idx] = local_sort_by_invariants(lam)
    M = local_qmat(lam);
    rePart = M(:,1);
    imNorm = sqrt(sum(M(:,2:4).^2,2));
    [~, idx] = sortrows([rePart(:), imNorm(:), M(:,2:4)], [1 2 3 4 5]);
    lamSorted = lam(idx);
end

function s = local_qstr(q, nd)
    if nargin < 2, nd = 6; end
    v = local_qvec(q);
    fmt = sprintf('%% .%df %+. %df*i %+. %df*j %+. %df*k', nd, nd, nd, nd);
    % MATLAB does not allow spaces inside conversion spec; build without those.
    fmt = sprintf('%% .%df %+. %df*i %+. %df*j %+. %df*k', nd, nd, nd, nd); %#ok<NASGU>
    % Safer explicit:
    s = sprintf(['% .',num2str(nd),'f %+.',num2str(nd),'f*i %+.',num2str(nd),'f*j %+.',num2str(nd),'f*k'], ...
        v(1), v(2), v(3), v(4));
end

function v = local_qvec(q)
    try
        [a,b,c,d] = parts(q);
        v = double([a,b,c,d]);
    catch
        v = double(compact(q));
    end
    v = v(:).';
end

function M = local_qmat(q)
    if isempty(q)
        M = zeros(0,4);
    end
    q = q(:);
    M = zeros(numel(q),4);
    for i = 1:numel(q)
        M(i,:) = local_qvec(q(i));
    end
end

function dmin = local_minsep_r4(lam)
    dmin = inf;
    for i = 1:numel(lam)
        vi = local_qvec(lam(i));
        for j = i+1:numel(lam)
            dmin = min(dmin, norm(vi - local_qvec(lam(j)), 2));
        end
    end
    if isinf(dmin)
        dmin = NaN;
    end
end
