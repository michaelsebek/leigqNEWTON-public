%LAA_ZOO_EX_2_1_EVIDENCE  Zoo demo — evidence/diagnostics companion for Supplement B (Section 2 Example 1).
%
%   - LaTeX block:  Supplement B, Section 2 Example 1
%   - Script:       LAA_Zoo_Ex_2_1_evidence.m
%
% This script reproduces (up to numerical tolerances and randomized sampling)
% the numerical results reported for the corresponding Zoo example.
%
% What it does
%   - Runs additional solver batches/diagnostics to support the claims in the Example text.
%   - Typical tasks include saturation checks, residual statistics, and sphere/no-sphere probes.
%
% How to run
%   >> LAA_Zoo_Ex_2_1_evidence
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

%% (A) Matrix A (Zoo Example 2–1)
% A = W + X*i + Y*j + Z*k
W = [ -7,   3,  11;
       6,   9,  -5;
      16,  20,  -5 ];
X = [  6,  -7,  -9;
       1,   7,  15;
       6,  -3,  19 ];
Y = [ -6,  11,   0;
      -5,  -6,  14;
      14,   4,   1 ];
Z = [ -7,  -2,   0;
       3, -10,  -1;
      11,   6,  -3 ];

A = quaternion(W, X, Y, Z);
n = size(A,1);

%% (B) Parameters
TOL.UniqueTol     = 1e-8;      % clustering tolerance in R^4
TOL.TolNewton     = 1e-12;     % solver tolerance
TOL.TargetResMin  = 1e-16;     % target certificate after refine

% Refiner tolerances (use full names; DO NOT use 'Tol')
TOL.TolX          = 1e-16;
TOL.TolFun        = 1e-16;
TOL.MaxIter       = 600;
TOL.TolResPolish  = 1e-16;
TOL.MaxIterPolish = 60;

CFG.batchRuns       = 25;      % solver calls per batch
CFG.batchRestarts   = 500;     % restarts per solver call
CFG.maxRuns         = 400;     % hard cap on total solver calls
CFG.patienceBatches = 8;       % stop after so many no-new batches
CFG.baseSeed        = 13579;   % reproducibility

DIAG.invTol      = 1e-8;       % tolerance for invariant grouping
DIAG.probeDirs   = 300;        % directions per candidate sphere
DIAG.probeSeed   = 24680;      % reproducibility

fprintf('=== LAA evidence: Zoo Example 2–1 (3x3, K=5) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', TOL.UniqueTol, TOL.TolNewton, TOL.TargetResMin);
fprintf('Sampling: batchRuns=%d | restarts/run=%d | maxRuns=%d | patience=%d batches\n', ...
    CFG.batchRuns, CFG.batchRestarts, CFG.maxRuns, CFG.patienceBatches);
fprintf('Sphere-probe: %d directions per candidate (seed=%d)\n\n', DIAG.probeDirs, DIAG.probeSeed);

%% (C) Stage 1 — multi-start saturation
lamDistinct = quaternion.empty(0,1);
rawHits = 0;
noNewBatches = 0;

maxBatches = ceil(CFG.maxRuns / CFG.batchRuns);
logRuns = zeros(maxBatches,1);
logK    = zeros(maxBatches,1);

for b = 1:maxBatches
    addedThisBatch = 0;

    for r = 1:CFG.batchRuns
        runIdx = (b-1)*CFG.batchRuns + r;
        if runIdx > CFG.maxRuns, break; end

        seed = CFG.baseSeed + runIdx;
        lamAcc = local_run_one(A, n, CFG.batchRestarts, TOL, seed);
        rawHits = rawHits + numel(lamAcc);

        [lamDistinct, added] = local_unique_append(lamDistinct, lamAcc, TOL.UniqueTol);
        addedThisBatch = addedThisBatch + added;
    end

    runsDone = min(b*CFG.batchRuns, CFG.maxRuns);
    logRuns(b) = runsDone;
    logK(b)    = numel(lamDistinct);

    if addedThisBatch == 0
        noNewBatches = noNewBatches + 1;
    else
        noNewBatches = 0;
    end

    if b == 1
        fprintf('[runs=%4d] DISTINCT increased: 0 -> %d\n', runsDone, logK(b));
    else
        fprintf('[runs=%4d] DISTINCT stays at %d   (no-new batches: %d/%d)\n', ...
            runsDone, logK(b), noNewBatches, CFG.patienceBatches);
    end

    if noNewBatches >= CFG.patienceBatches || runsDone >= CFG.maxRuns
        logRuns = logRuns(1:b);
        logK    = logK(1:b);
        break;
    end
end

fprintf('\nStage 1 done. Total runs=%d, raw hits=%d, distinct=%d\n\n', ...
    logRuns(end), rawHits, numel(lamDistinct));

%% (D) Refine / polish the distinct representatives
[lamRef, resMinRef, refineNote] = local_refine_list(A, lamDistinct, TOL, CFG.baseSeed+999);

% Stable order for reporting
[lamRef, idxSort] = local_sort_by_invariants(lamRef);
resMinRef = resMinRef(idxSort);

fprintf('Distinct eigenvalues found after refine (%s):\n', refineNote);
for k = 1:numel(lamRef)
    fprintf('  %2d) %s   | r_min = %.3e\n', k, local_qstr(lamRef(k),6), double(resMinRef(k)));
end
fprintf('Min pairwise ||lambda_i - lambda_j||_{R^4}: %.6g\n\n', local_minsep_r4(lamRef));

fprintf('Discovery log (runs vs Kdistinct):\n');
fprintf('    Runs    Kdistinct\n');
for t = 1:numel(logRuns)
    fprintf('%8d%11d\n', logRuns(t), logK(t));
end

%% (E) Stage 2 — no-sphere diagnostics
fprintf('\n=== Stage 2: no-sphere diagnostics ===\n');

groups = local_group_by_invariants(lamRef, DIAG.invTol);
fprintf('Invariant groups by (Re, ||Im||) within tol=%.1e:\n', DIAG.invTol);
for g = 1:numel(groups)
    fprintf('  Group %d size %d: indices %s\n', g, numel(groups{g}), mat2str(groups{g}));
end
if all(cellfun(@numel, groups) == 1)
    fprintf('No group contains 2+ points with matching (Re, ||Im||) invariants.\n');
else
    fprintf('WARNING: some invariant groups contain 2+ points (possible spherical structure).\n');
end

fprintf('\nSphere-probe test (r_min over random directions on candidate spheres):\n');
rng(DIAG.probeSeed,'twister');
for k = 1:numel(lamRef)
    v = local_qvec(lamRef(k));
    alpha = v(1);
    rho = norm(v(2:4),2);

    if rho <= 0
        fprintf('  k=%d: alpha=%.6g, rho=0 (real) -> probe skipped\n', k, alpha);
        continue;
    end

    rvals = zeros(DIAG.probeDirs,1);
    for t = 1:DIAG.probeDirs
        u = randn(3,1);
        u = u / norm(u,2);
        lamProbe = quaternion(alpha, rho*u(1), rho*u(2), rho*u(3));
        rvals(t) = leigqNEWTON_cert_resMin(A, lamProbe);
    end

    fprintf('  k=%d: alpha=%.6g, rho=%.6g | min r_min=%.3e | median r_min=%.3e\n', ...
        k, alpha, rho, min(rvals), median(rvals));
end

fprintf('\nInterpretation:\n');
fprintf('  - If an eigen-sphere existed and were reachable, many off-axis probes would yield r_min ~ 0.\n');
fprintf('  - For isolated eigenvalues, off-axis r_min should stay bounded away from 0.\n');

fprintf('\n=== End of evidence script ===\n');

%% ============================================================
% Local helpers (script-local functions)
%% ============================================================

function local_require_function(fname)
    if exist(fname,'file') ~= 2
        error('Missing required function "%s" on the MATLAB path.', fname);
    end
end

function lamAcc = local_run_one(A, n, restarts, TOL, seed)
    [~, ~, info] = leigqNEWTON(A, ...
        'Num',       n, ...
        'Restarts',  restarts, ...
        'Seed',      seed, ...
        'Tol',       TOL.TolNewton, ...
        'UniqueTol', TOL.UniqueTol, ...
        'Verbose',   0);

    lamAcc = [];
    if isstruct(info)
        if isfield(info,'lambdaAccepted') && ~isempty(info.lambdaAccepted)
            lamAcc = info.lambdaAccepted;
        elseif isfield(info,'lambda') && ~isempty(info.lambda)
            lamAcc = info.lambda;
        end
    end
    lamAcc = lamAcc(:);
end

function [lamU, added] = local_unique_append(lamU, lamNew, tol)
    added = 0;
    if isempty(lamNew), return; end

    for i = 1:numel(lamNew)
        q = lamNew(i);
        if isempty(lamU)
            lamU(1,1) = q;
            added = added + 1;
            continue;
        end
        d = arrayfun(@(u) norm(local_qvec(u)-local_qvec(q),2), lamU);
        if all(d > tol)
            lamU(end+1,1) = q; %#ok<AGROW>
            added = added + 1;
        end
    end
end

function [lamOut, resMinOut, note] = local_refine_list(A, lam0, TOL, seed)
    if nargin < 4, seed = []; end

    if isempty(lam0)
        lamOut = lam0;
        resMinOut = [];
        note = 'empty';
        return;
    end

    if exist('leigqNEWTON_refine_batch','file') == 2
        [lamOut, ~, cert] = leigqNEWTON_refine_batch(A, lam0, ...
            'Mode',          'auto', ...
            'TargetResMin',  TOL.TargetResMin, ...
            'TolX',          TOL.TolX, ...
            'TolFun',        TOL.TolFun, ...
            'MaxIter',       TOL.MaxIter, ...
            'DoPolish',      true, ...
            'TolResPolish',  TOL.TolResPolish, ...
            'MaxIterPolish', TOL.MaxIterPolish, ...
            'Seed',          seed, ...
            'Verbose',       0);

        if isstruct(cert) && isfield(cert,'resMin')
            resMinOut = cert.resMin(:);
        else
            resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
        end

        note = 'refine_batch';
        return;
    end

    if exist('leigqNEWTON_refine_auto','file') == 2
        [lamBest, lamPol] = leigqNEWTON_refine_auto(A, lam0, ...
            'Mode',          'auto', ...
            'TargetResMin',  TOL.TargetResMin, ...
            'TolX',          TOL.TolX, ...
            'TolFun',        TOL.TolFun, ...
            'MaxIter',       TOL.MaxIter, ...
            'DoPolish',      true, ...
            'TolResPolish',  TOL.TolResPolish, ...
            'MaxIterPolish', TOL.MaxIterPolish, ...
            'Seed',          seed, ...
            'Verbose',       0);

        lamOut = lamPol;
        if isempty(lamOut), lamOut = lamBest; end
        resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
        note = 'refine_auto';
        return;
    end

    lamOut = lam0(:);
    resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
    note = 'no_refine';
end

function [lamSorted, idx] = local_sort_by_invariants(lam)
    if isempty(lam)
        lamSorted = lam;
        idx = [];
        return;
    end

    M = local_qmat(lam);
    rePart = M(:,1);
    imNorm = sqrt(sum(M(:,2:4).^2,2));
    key = [rePart, imNorm, M(:,2:4)];
    [~, idx] = sortrows(key, [1 2 3 4 5]);
    lamSorted = lam(idx);
end

function groups = local_group_by_invariants(lam, tol)
    M = local_qmat(lam);
    rePart = M(:,1);
    imNorm = sqrt(sum(M(:,2:4).^2,2));

    used = false(numel(lam),1);
    groups = {};

    for i = 1:numel(lam)
        if used(i), continue; end
        idx = find(~used & abs(rePart-rePart(i)) <= tol & abs(imNorm-imNorm(i)) <= tol);
        used(idx) = true;
        groups{end+1} = idx(:).'; %#ok<AGROW>
    end
end

function dmin = local_minsep_r4(lam)
    dmin = inf;
    for i = 1:numel(lam)
        vi = local_qvec(lam(i));
        for j = i+1:numel(lam)
            dmin = min(dmin, norm(vi - local_qvec(lam(j)),2));
        end
    end
end

function s = local_qstr(q, nd)
    if nargin < 2, nd = 6; end
    v = local_qvec(q);
    s = sprintf(['% .',num2str(nd),'f %+.',num2str(nd),'f*i %+.',num2str(nd),'f*j %+.',num2str(nd),'f*k'], ...
        v(1), v(2), v(3), v(4));
end

function v = local_qvec(q)
% Return a 1x4 double row [a b c d] for quaternion q.
    try
        [a,b,c,d] = parts(q);
        v = double([a(:), b(:), c(:), d(:)]);
    catch
        v = double(compact(q));
    end
    if size(v,1) ~= 1
        v = v(1,:);
    end
end

function M = local_qmat(q)
% Return a Kx4 double matrix for quaternion array q.
    q = q(:);
    try
        [a,b,c,d] = parts(q);
        M = double([a(:), b(:), c(:), d(:)]);
    catch
        M = double(compact(q));
    end
end
