%LAA_ZOO_EX_2_1_COMPUTATION Zoo Example 2--1: 3x3, K=5, computation + refine.
%
%   Reader-facing demo script shipped with the LAA paper and leigqNEWTON toolbox
%
%   Zoo – Supplement B – Section 2 - Example 1: a 3x3 integer quaternion matrix
%     with FIVE distinct isolated left eigenvalues (K=5 > n).
%
%   What this script does
%     1) Construct the 3x3 integer quaternion matrix A.
%     2) Run leigqNEWTON in a reproducible multi-start loop (request K=n per run).
%     3) Cluster accepted hits into DISTINCT representatives (in R^4) using UniqueTol.
%     4) Refine + optionally polish the DISTINCT representatives using
%        leigqNEWTON_refine_batch (fallback: leigqNEWTON_refine_auto).
%     5) Print refined eigenvalues and their certificates r_min.
%
%   Outputs / printing
%     - Prints a compact summary and two eigenvalue lists:
%         (i) DISTINCT candidates returned by the solver, and
%         (ii) the same list after refinement/polishing.
%     - At the end, prints which workspace variables hold the results.
%
%   Reproducibility
%     - Deterministic seeds: Seed0 + runIndex.
%     - Stable ordering of printed results.
%
%   Citation
%     - Please cite the associated LAA article and the released leigqNEWTON code.
%
%   License
%     - Distributed under the same license as the leigqNEWTON toolbox.
%
%   Requirements on MATLAB path
%     - leigqNEWTON
%     - leigqNEWTON_cert_resMin
%     - leigqNEWTON_refine_batch (recommended) or leigqNEWTON_refine_auto
%
%   M. Šebek (CTU Prague) + ChatGPT, 2026-01-24

clear; clc;

%% (0) Basic dependency check
local_require_function('leigqNEWTON');
local_require_function('leigqNEWTON_cert_resMin');

%% (A) Construct the matrix A (component-wise integers)
% A = W + X*i + Y*j + Z*k
W = [
    -7   3  11;
     6   9  -5;
    16  20  -5];

X = [
     6  -7  -9;
     1   7  15;
     6  -3  19];

Y = [
    -6  11   0;
    -5  -6  14;
    14   4   1];

Z = [
    -7  -2   0;
     3 -10  -1;
    11   6  -3];

A = quaternion(W, X, Y, Z);
[n,~] = size(A);

%% (B) User-facing parameters (match Supplement B Ex 2--1)
TOL.UniqueTol     = 1e-8;      % clustering tolerance in R^4
TOL.TolNewton     = 1e-12;     % leigqNEWTON convergence tolerance
TOL.TargetResMin  = 1e-16;     % target certificate after refine

% Refiner tolerances: specify full names (avoid name ambiguity)
TOL.TolX          = 1e-16;
TOL.TolFun        = 1e-16;
TOL.MaxIter       = 600;

TOL.TolResPolish  = 1e-15;
TOL.MaxIterPolish = 50;

CFG.Runs          = 50;        % number of independent solver calls
CFG.Restarts      = 500;       % restarts per solver call
CFG.Krequest      = n;         % request K=n eigenpairs per run
CFG.Seed0         = 24680;     % reproducibility

fprintf('=== LAA demo: Zoo Example 2--1 (3x3, K=5 distinct isolated left eigenvalues) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', TOL.UniqueTol, TOL.TolNewton, TOL.TargetResMin);
fprintf('Runs=%d | Restarts/run=%d | Krequest=%d | Seed0=%d\n\n', CFG.Runs, CFG.Restarts, CFG.Krequest, CFG.Seed0);

%% (C) Multi-start campaign: collect accepted hits
lamAll = quaternion.empty(0,1);
hitsTotal = 0;

for r = 1:CFG.Runs
    seed = CFG.Seed0 + r;

    % Try the most common leigqNEWTON signature (Tol) used in the Zoo scripts.
    try
        [lam, ~, info] = leigqNEWTON(A, ...
            'Num',       CFG.Krequest, ...
            'Restarts',  CFG.Restarts, ...
            'Seed',      seed, ...
            'Tol',       TOL.TolNewton, ...
            'UniqueTol', TOL.UniqueTol, ...
            'Verbose',   0);
    catch ME
        % Fallback for older/newer parameter naming.
        if contains(ME.message, 'TolNewton') || contains(ME.message, 'Unrecognized')
            [lam, ~, info] = leigqNEWTON(A, ...
                'Num',       CFG.Krequest, ...
                'Restarts',  CFG.Restarts, ...
                'Seed',      seed, ...
                'TolNewton', TOL.TolNewton, ...
                'UniqueTol', TOL.UniqueTol, ...
                'Verbose',   0);
        else
            rethrow(ME);
        end
    end

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
[lamRef, resMinRef, refineNote] = local_refine_list(A, lamDistinct, TOL, CFG.Seed0+999);
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

if isempty(lamRef)
    fprintf('\nCertificate range r_min: N/A (no refined eigenvalues)\n');
    fprintf('Minimum pairwise separation (R^4 distance): N/A\n');
else
    fprintf('\nCertificate range r_min: [%.3e, %.3e]\n', min(resMinRef), max(resMinRef));
    fprintf('Minimum pairwise separation (R^4 distance): %.6g\n', local_minsep_r4(lamRef));
end

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

function local_require_function(name)
    if exist(name,'file') ~= 2
        error('Missing required function on MATLAB path: %s', name);
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

function lamU = local_unique_quat_r4(lamAll, tol)
% Greedy clustering in R^4 (keeps first occurrence).
    lamU = quaternion.empty(0,1);
    for i = 1:numel(lamAll)
        q = lamAll(i);
        if isempty(lamU)
            lamU(1,1) = q;
            continue;
        end
        d = arrayfun(@(u) norm(local_qvec(u)-local_qvec(q),2), lamU);
        if all(d > tol)
            lamU(end+1,1) = q; %#ok<AGROW>
        end
    end
end

function [lamS, idx] = local_sort_by_invariants(lam)
% Sort by (real part, imag norm, i, j, k) for stable printing.
    if isempty(lam)
        lamS = lam;
        idx = [];
        return;
    end
    Q = arrayfun(@(q) local_qvec(q), lam, 'UniformOutput', false);
    Q = cat(1, Q{:}); % K-by-4
    imagNorm = sqrt(sum(Q(:,2:4).^2,2));
    key = [Q(:,1), imagNorm, Q(:,2:4)];
    [~, idx] = sortrows(key);
    lamS = lam(idx);
end

function v = local_qvec(q)
% Convert quaternion to a 1x4 real row [w x y z].
% (Works with MATLAB quaternion; for other classes this may need adjustment.)
    a = compact(q);
    v = double(a(:)).';
end

function s = local_qstr(q, nd)
% Short pretty string for a quaternion (nd decimals). Falls back to disp() output.
    fmt = sprintf('%%.%df', nd);
    a = compact(q);
    w = double(a(1)); x = double(a(2)); y = double(a(3)); z = double(a(4));
    s = sprintf('%s%s%s%s', num2str(w,fmt), local_qterm(x,'i',fmt), local_qterm(y,'j',fmt), local_qterm(z,'k',fmt));
    s = strrep(s,'+-','-');
    s = strrep(s,'  ',' ');
    s = strtrim(s);
end

function t = local_qterm(val, symb, fmt)
    if abs(val) < 0.5*10^(-12)
        t = '';
        return;
    end
    if val >= 0
        t = [' +' num2str(val,fmt) symb];
    else
        t = [' ' num2str(val,fmt) symb];
    end
end

function dmin = local_minsep_r4(lam)
% Minimum pairwise separation in R^4.
    K = numel(lam);
    if K <= 1
        dmin = NaN;
        return;
    end
    Q = arrayfun(@(q) local_qvec(q), lam, 'UniformOutput', false);
    Q = cat(1, Q{:});
    dmin = inf;
    for i = 1:K
        for j = i+1:K
            dij = norm(Q(i,:)-Q(j,:), 2);
            if dij < dmin
                dmin = dij;
            end
        end
    end
end

function [lamRef, resMinRef, note] = local_refine_list(A, lam0, TOL, seed)
% Wrapper to use whichever refiner is available.
    lamRef = lam0;
    resMinRef = leigqNEWTON_cert_resMin(A, lam0);

    if isempty(lam0)
        note = 'none';
        return;
    end

    if exist('leigqNEWTON_refine_batch','file') == 2
        note = 'refine_batch';
        [lamRef, ~, cert] = leigqNEWTON_refine_batch(A, lam0, ...
            'Mode',         'auto', ...
            'TargetResMin', TOL.TargetResMin, ...
            'MaxIter',      TOL.MaxIter, ...
            'TolX',         TOL.TolX, ...
            'TolFun',       TOL.TolFun, ...
            'DoPolish',     true, ...
            'TolResPolish', TOL.TolResPolish, ...
            'MaxIterPolish',TOL.MaxIterPolish, ...
            'Seed',         seed, ...
            'Verbose',      0);

        if isstruct(cert) && isfield(cert,'resMin')
            resMinRef = double(cert.resMin(:));
        else
            resMinRef = leigqNEWTON_cert_resMin(A, lamRef);
        end
        return;
    end

    if exist('leigqNEWTON_refine_auto','file') == 2
        note = 'refine_auto';
        try
            [lamRef, vRef] = leigqNEWTON_refine_auto(A, lam0, ...
                'TargetResMin', TOL.TargetResMin, ...
                'Seed', seed, 'Verbose', 0);
            %#ok<NASGU>
        catch
            [lamRef, ~] = leigqNEWTON_refine_auto(A, lam0);
        end
        resMinRef = leigqNEWTON_cert_resMin(A, lamRef);
        return;
    end

    note = 'none';
end
