%LAA_ZOO_EX_2_2_EVIDENCE  Zoo demo — evidence/diagnostics companion for Supplement B (Section 2 Example 2).
%
%   - LaTeX block:  Supplement B, Section 2 Example 2
%   - Script:       LAA_Zoo_Ex_2_2_evidence.m
%
% This script reproduces (up to numerical tolerances and randomized sampling)
% the numerical results reported for the corresponding Zoo example.
%
% What it does
%   - Runs additional solver batches/diagnostics to support the claims in the Example text.
%   - Typical tasks include saturation checks, residual statistics, and sphere/no-sphere probes.
%
% How to run
%   >> LAA_Zoo_Ex_2_2_evidence
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

%% (A) Matrix A (Zoo Example 2–2)
% A = W + X*i + Y*j + Z*k
W = [
    -6  -2  13   7  -3;
     0   0 -20  10 -17;
     4  -5   7   6  -3;
   -12   5   1  -4  16;
     7 -10  19   6   1];

X = [
     1  -8  -4  -6  16;
    -7 -14  -6  -7   1;
    15 -10  -3   9   6;
    16   8  -3  -2   3;
   -13  -4  -1   8   3];

Y = [
   -10   1   7  -7 -18;
    -7  -4   4   8  -5;
    -2  -8 -19  23  -8;
    -4 -18  -7  -5  -3;
   -11  -2   9   6  -5];

Z = [
    -8  -5  14  14   0;
   -17 -18   4   9 -25;
     3 -18 -17   6  -3;
    -1 -19   5  -7   5;
     1  -8   6   3  -7];

A = quaternion(W,X,Y,Z);
n = size(A,1);

%% (B) Parameters
TOL.UniqueTol     = 1e-8;
TOL.TolNewton     = 1e-12;
TOL.TargetResMin  = 1e-16;

TOL.TolX          = 1e-16;
TOL.TolFun        = 1e-16;
TOL.MaxIter       = 600;
TOL.TolResPolish  = 1e-16;
TOL.MaxIterPolish = 60;

CFG.batchRuns       = 15;      % solver calls per batch
CFG.batchRestarts   = 1200;    % restarts per solver call
CFG.maxRuns         = 250;     % cap on solver calls
CFG.patienceBatches = 8;       % stop after this many no-new batches
CFG.baseSeed        = 24680;

SP.nProbe           = 250;     % directions per candidate
SP.invTol           = 5e-9;    % invariant grouping tolerance
SP.seedProbe        = 777;

fprintf('=== LAA evidence: Zoo Example 2–2 (5x5, K=9 distinct isolated left eigenvalues) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', TOL.UniqueTol, TOL.TolNewton, TOL.TargetResMin);
fprintf('Batches: batchRuns=%d, restarts/run=%d, maxRuns=%d | patience=%d batches\n', ...
    CFG.batchRuns, CFG.batchRestarts, CFG.maxRuns, CFG.patienceBatches);
fprintf('Sphere probe: %d directions/candidate | invTol=%.1e\n\n', SP.nProbe, SP.invTol);

%% (C) Stage 1 — multi-start saturation
lamDistinct = quaternion.empty(0,1);
rawHits = 0;
noNewBatches = 0;

maxBatches = ceil(CFG.maxRuns / CFG.batchRuns);
discLog = zeros(maxBatches,2); % [runsDone, Kdistinct]

for b = 1:maxBatches
    addedThisBatch = 0;

    for r = 1:CFG.batchRuns
        runIdx = (b-1)*CFG.batchRuns + r;
        if runIdx > CFG.maxRuns
            break;
        end

        seed = CFG.baseSeed + runIdx;
        lamAcc = local_run_one(A, n, CFG.batchRestarts, TOL.TolNewton, TOL.UniqueTol, seed);

        rawHits = rawHits + numel(lamAcc);
        [lamDistinct, added] = local_unique_append(lamDistinct, lamAcc, TOL.UniqueTol);
        addedThisBatch = addedThisBatch + added;
    end

    runsDone = min(b*CFG.batchRuns, CFG.maxRuns);
    discLog(b,:) = [runsDone, numel(lamDistinct)];

    if addedThisBatch == 0
        noNewBatches = noNewBatches + 1;
    else
        noNewBatches = 0;
    end

    if b == 1
        fprintf('[runs=%4d] DISTINCT increased: 0 -> %d\n', runsDone, numel(lamDistinct));
    else
        fprintf('[runs=%4d] DISTINCT stays at %d   (no-new batches: %d/%d)\n', ...
            runsDone, numel(lamDistinct), noNewBatches, CFG.patienceBatches);
    end

    if noNewBatches >= CFG.patienceBatches || runsDone >= CFG.maxRuns
        break;
    end
end

discLog = discLog(1:b,:);

fprintf('\nStage 1 done. Total runs=%d, raw hits=%d, distinct=%d\n\n', ...
    discLog(end,1), rawHits, numel(lamDistinct));

%% (D) Refine / polish the distinct representatives
[lamRef, resMinRef, refineNote] = local_refine_list(A, lamDistinct, TOL, CFG.baseSeed+999);

[lamRef, idxSort] = local_sort_by_invariants(lamRef);
resMinRef = resMinRef(idxSort);

fprintf('--- DISTINCT eigenvalues after refine (%s) ---\n', refineNote);
for k = 1:numel(lamRef)
    fprintf('%2d) %s   | r_min = %.3e\n', k, local_qstr(lamRef(k), 6), double(resMinRef(k)));
end
fprintf('Min separation min_{i!=j} ||lambda_i-lambda_j||_{R^4}: %.6g\n\n', local_minsep_r4(lamRef));

fprintf('Discovery log (runsDone vs Kdistinct):\n');
fprintf('    Runs    Kdistinct\n');
for t = 1:size(discLog,1)
    fprintf('%8d%11d\n', discLog(t,1), discLog(t,2));
end

%% (E) Stage 2 — sphere checks
fprintf('\n=== Stage 2: sphere checks (necessary conditions + probe) ===\n');

groups = local_group_by_invariants(lamRef, SP.invTol);
fprintf('Invariant groups by (Re, ||Im||) with tol=%.1e:\n', SP.invTol);
for g = 1:numel(groups)
    fprintf('  Group %d size %d: indices %s\n', g, numel(groups{g}), mat2str(groups{g}));
end

if all(cellfun(@numel, groups) == 1)
    fprintf('No group contains 2+ points with matching invariants -> strong evidence against eigen-spheres.\n\n');
else
    fprintf('WARNING: some invariant groups contain 2+ points -> possible spherical structure.\n\n');
end

fprintf('Sphere-probe test (r_min over random directions on candidate spheres):\n');
rng(SP.seedProbe,'twister');

for k = 1:numel(lamRef)
    v = local_qvec(lamRef(k));
    alpha = v(1);
    rho = sqrt(sum(v(2:4).^2));

    if rho <= 0
        fprintf('  k=%d: alpha=%.6g, rho=0 (real) -> skipped\n', k, alpha);
        continue;
    end

    rvals = zeros(SP.nProbe,1);
    for t = 1:SP.nProbe
        u = randn(3,1);
        u = u / norm(u);
        lamProbe = quaternion(alpha, rho*u(1), rho*u(2), rho*u(3));
        rvals(t) = double(leigqNEWTON_cert_resMin(A, lamProbe));
    end

    fprintf('  k=%d: alpha=%.6g, rho=%.6g | probes=%d | min r_min=%.3e | median r_min=%.3e\n', ...
        k, alpha, rho, SP.nProbe, min(rvals), median(rvals));
end

fprintf('\nInterpretation (rule of thumb):\n');
fprintf('  - If an eigen-sphere existed and were reachable, many off-axis probes would give r_min ≈ 0.\n');
fprintf('  - For isolated eigenvalues, off-axis r_min stays bounded away from 0.\n');
fprintf('=== End of evidence script ===\n');

%% ============================================================
% Local helper functions (kept script-local for portability)
%% ============================================================

function local_require_function(fname)
    if exist(fname,'file') ~= 2
        error('Missing required function on MATLAB path: %s', fname);
    end
end

function lamAcc = local_run_one(A, n, restarts, tolNewton, uniqueTol, seed)
    [lam, ~, info] = leigqNEWTON(A, ...
        'Num',       n, ...
        'Restarts',  restarts, ...
        'Seed',      seed, ...
        'Tol',       tolNewton, ...
        'UniqueTol', uniqueTol, ...
        'Verbose',   0);

    lamAcc = lam;
    if isstruct(info)
        if isfield(info,'lambdaAccepted') && ~isempty(info.lambdaAccepted)
            lamAcc = info.lambdaAccepted;
        elseif isfield(info,'lambda') && ~isempty(info.lambda)
            lamAcc = info.lambda;
        end
    end
    lamAcc = lamAcc(:);
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

        lamOut = lamBest;
        if ~isempty(lamPol)
            lamOut = lamPol;
        end
        resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
        note = 'refine_auto';
        return;
    end

    lamOut = lam0(:);
    resMinOut = leigqNEWTON_cert_resMin(A, lamOut);
    note = 'no_refine';
end

function [lamU, added] = local_unique_append(lamU, lamNew, tol)
    added = 0;
    if isempty(lamNew)
        return;
    end
    lamNew = lamNew(:);
    for i = 1:numel(lamNew)
        if isempty(lamU)
            lamU = lamNew(i);
            added = added + 1;
            continue;
        end
        dmin = inf;
        vi = local_qvec(lamNew(i));
        for k = 1:numel(lamU)
            dmin = min(dmin, norm(vi - local_qvec(lamU(k)), 2));
        end
        if dmin > tol
            lamU(end+1,1) = lamNew(i); %#ok<AGROW>
            added = added + 1;
        end
    end
end

function [lamSorted, idx] = local_sort_by_invariants(lam)
    if isempty(lam)
        lamSorted = lam; idx = [];
        return;
    end
    M = local_qmat(lam);
    rePart = M(:,1);
    imNorm = sqrt(sum(M(:,2:4).^2,2));
    key = [rePart(:), imNorm(:), M(:,2:4)];
    [~, idx] = sortrows(key, [1 2 3 4 5]);
    lamSorted = lam(idx);
end

function groups = local_group_by_invariants(lam, tol)
    if isempty(lam)
        groups = {};
        return;
    end
    M = local_qmat(lam);
    rePart = M(:,1);
    imNorm = sqrt(sum(M(:,2:4).^2,2));

    used = false(numel(lam),1);
    groups = {};
    for i = 1:numel(lam)
        if used(i), continue; end
        idx = find(~used & abs(rePart - rePart(i)) <= tol & abs(imNorm - imNorm(i)) <= tol);
        used(idx) = true;
        groups{end+1} = idx(:).'; %#ok<AGROW>
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
    q = q(:);
    try
        [a,b,c,d] = parts(q);
        M = double([a(:), b(:), c(:), d(:)]);
    catch
        M = double(compact(q));
    end
end

function s = local_qstr(q, nd)
    if nargin < 2, nd = 6; end
    v = local_qvec(q);
    fmt = ['% .',num2str(nd),'f %+.',num2str(nd),'f*i %+.',num2str(nd),'f*j %+.',num2str(nd),'f*k'];
    s = sprintf(fmt, v(1), v(2), v(3), v(4));
end
