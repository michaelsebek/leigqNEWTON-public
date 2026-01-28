%LAA_ZOO_EX_4_5_COMPUTATION  Zoo demo — computation companion for Supplement B (Section 4 Example 5).
%
%   - LaTeX block:  Supplement B, Section 4 Example 5
%   - Script:       LAA_Zoo_Ex_4_5_computation.m
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
%   >> LAA_Zoo_Ex_4_5_computation
%
% Outputs
%   Creates variables in the base workspace; see the final "Variables created" printout.
%
% Requirements
%   - leigqNEWTON toolbox on the MATLAB path (quaternion class + leigqNEWTON_* routines).
%
% Notes
%   - This script temporarily suppresses selected MATLAB warnings and restores them on exit.
%   - Option 'Tol' is supported as a shortcut alias for leigqNEWTON_refine_batch (sets TolX/TolFun/TolResPolish).
%     use full names ('TolX','TolFun','TolResPolish',...).
%
%   M. Šebek (CTU Prague) + ChatGPT, 2026-01-23
clc;
format short g;

% ===================== user-tunable parameters =====================
UniqueTol     = 1.0e-08;  % clustering tol for DISTINCT candidates (R^4)
TolNewton     = 1.0e-12;  % shown for reference (Stage 2 uses refine_batch tolerances)
TargetResMin  = 1.0e-16;  % target residual for refinement/polish

Collect       = 20;       % desired number of DISTINCT candidates (Supplement uses K=20)
RunsMax       = 120;      % max multi-start runs in sampling
Restarts      = 1500;     % restarts per run in sampling
Seed0         = 24680;    % base seed

fprintf('=== LAA demo: Zoo Example 4-5 (4x4, spherical eigenvalues; AAsph2) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', UniqueTol, TolNewton, TargetResMin);
fprintf('Sampling: Collect=%d distinct | RunsMax=%d | Restarts=%d | Seed0=%d\n\n', Collect, RunsMax, Restarts, Seed0);

% ===================== warning handling (silence + restore) =====================
% During sampling/refinement, MATLAB may emit linear-solver conditioning warnings
% (e.g., 'Matrix is close to singular or badly scaled.'). These are expected in
% some Newton steps and would clutter a reader-facing demo. We silence them and
% restore the user's warning state automatically when the script ends.
warnState = warning; %#ok<NASGU>
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:illConditionedMatrix');

% ===================== construct the matrix from the Supplement =====================
A = build_AAsph2();

% (Optional) report a numerical matrix norm used in the Supplement text.
% For quaternion matrices, norm(A,2) may fail; we compute the induced 2-norm
% via a real 4n-by-4n left-multiplication embedding.
nA2 = local_norm2_quat_left(A);

% ===================== Stage 1: sampling + sphere detection =====================
fprintf('Stage 1: sphere sampling/detection ...\n');

[lamAllFull, lamSamplesFull, lambda0Full, clsFull, sphFull, infoFull] = leigqNEWTON_sphere_sample(A, ...
    'Collect', Collect, 'RunsMax', RunsMax, 'Seed0', Seed0, 'Restarts', Restarts, 'UniqueTol', UniqueTol);

KtotFull      = numel(lamAllFull);
KdistinctFull = numel(lamSamplesFull);

fprintf('Stage 1 done.  Ktot=%d, Kdistinct=%d\n\n', KtotFull, KdistinctFull);

% Residuals for Stage 1 lists
resAllFull     = leigqNEWTON_cert_resMin(A, lamAllFull);
resSamplesFull = leigqNEWTON_cert_resMin(A, lamSamplesFull);

fprintf('Stage 1 results (candidates):\n');
print_lam_list('lamAll (before refinement)', lamAllFull, resAllFull);

fprintf('Stage 1 samples (distinct):\n');
print_lam_list('lamSamples (before refinement)', lamSamplesFull, resSamplesFull);

fprintf('Stage 1 classification summary (cls):\n');
print_cls_summary(clsFull);

% Take the first detected sphere (this example is designed to have one)
sph0Full = [];
if ~isempty(sphFull)
    sph0Full = sphFull(1);
end

if isempty(sph0Full)
    fprintf('\nSphere information returned by detection:\n');
    fprintf('  (none)\n');
else
    fprintf('\nSphere information returned by detection:\n');
    print_sphere_info('sph0', sph0Full);
end

% ------------------ downselect to EXACTLY Collect candidates ------------------
% The sampler may overshoot Collect by a small amount (e.g., it may add multiple
% DISTINCT points in the last run). For the Supplement, we present K=20.
% Here we deterministically downselect to K=Collect while preserving:
%   - all isolated/unclassified samples (cls==0), and
%   - the best (Collect-#isolated) inliers of the detected sphere.

[lamAll0, lamSamples0, cls0, sph0, downselectNote] = downselect_to_collect( ...
    lamAllFull, lamSamplesFull, clsFull, sph0Full, resSamplesFull, Collect);

if ~isempty(downselectNote)
    fprintf('\n%s\n', downselectNote);
end

% For this demo we treat "candidates" as the DISTINCT set shown in the Supplement.
% (In this example lamAllFull and lamSamplesFull are typically identical.)
lamAll0 = lamSamples0;
resAll0 = leigqNEWTON_cert_resMin(A, lamAll0);

Ktot      = numel(lamAll0);
Kdistinct = numel(lamSamples0);

fprintf('\nStage 1 (selected for Supplement): K=%d distinct candidates\n', Ktot);
print_lam_list('lamAll0 (selected distinct candidates)', lamAll0, resAll0);

if ~isempty(sph0)
    fprintf('Detected spherical cluster (selected set):\n');
    print_sphere_info('sph0', sph0);
end

idxIso = find(cls0(:) == 0);
idxSphere = [];
if ~isempty(sph0)
    idxSphere = sph0.inliers(:);
end

if ~isempty(idxIso)
    fprintf('Isolated/unclassified DISTINCT samples (before refinement):\n');
    print_lam_list('isolated', lamSamples0(idxIso), leigqNEWTON_cert_resMin(A, lamSamples0(idxIso)));
end

% ===================== Stage 2: refinement/polish =====================
fprintf('\nStage 2: refinement/polish of candidates ...\n');

tStart = tic;
K = numel(lamAll0);
lamRef = lamAll0;
resRef = nan(K,1);
certRef = cell(K,1);

for k = 1:K
    fprintf('  Refinement %d/%d ... ', k, K);
    tOne = tic;

    [lamRef(k), certRef{k}] = local_refine_one(A, lamAll0(k), TargetResMin);
    % Extract numeric certificate r_min (resMin) robustly.
    if isstruct(certRef{k}) && isfield(certRef{k},'resMin') && ~isempty(certRef{k}.resMin)
        resRef(k) = double(certRef{k}.resMin(1));
    else
        % Fallback: compute certificate directly (should not happen in normal use)
        resRef(k) = double(leigqNEWTON_cert_resMin(A, lamRef(k)));
    end

    fprintf('done (dt=%.2fs, r~%.3e)\n', toc(tOne), resRef(k));
end

fprintf('Stage 2 done. elapsed=%.1fs\n\n', toc(tStart));

fprintf('Stage 2 results (refined candidates):\n');
print_lam_list('lamRef (after Stage 2)', lamRef, resRef);

fprintf('Stage 2: residual summary (refined candidates)\n');
fprintf('  med(resMin) ~ %.3e\n', median(resRef));
fprintf('  max(resMin) ~ %.3e\n', max(resRef));
fprintf('  ||A||_2 (left-embedding) ~ %.4g\n\n', nA2);

% In this demo, samples==candidates, so the refined samples are identical.
lamSamplesRef = lamRef;
resSamplesRef = resRef;

fprintf('Refined DISTINCT samples (mapped to Stage 2 refinements):\n');
print_lam_list('lamSamplesRef (after Stage 2)', lamSamplesRef, resSamplesRef);

% ===================== Sphere fit AFTER refinement (no validator) =====================
if ~isempty(sph0) && ~isempty(idxSphere)
    fprintf('Sphere fit AFTER refinement (from refined candidates) ...\n');

    % Primary: use fitSPHEREfromLambdas if available (it also rejects outliers).
    % Fallback: use the detected inliers and a local least-squares fit.
    [sphRef, metaRef, fitNote] = fit_refined_sphere(lamSamplesRef, idxSphere);

    if ~isempty(fitNote)
        fprintf('%s\n', fitNote);
    end

    if ~isempty(sphRef)
        fprintf('Refined sphere parametrization (sphRef):\n');
        print_sphere_info('sphRef', sphRef);

        % Consistency of refined inliers with the REFITTED sphere
        idxIn = sphRef.inliers(:);
        if isempty(idxIn)
            idxIn = idxSphere(:);
        end
        dev = builtin('abs', qabs(lamSamplesRef(idxIn) - sphRef.center) - double(sphRef.radius));
        fprintf('Sphere consistency on refined inliers (using sphRef center/radius):\n');
        fprintf('  median | ||lam-c|| - r | ~ %.3e\n', median(dev));
        fprintf('  max    | ||lam-c|| - r | ~ %.3e\n\n', max(dev));

        % Two representative inliers (as in the Supplement text)
        if numel(idxIn) >= 2
            a = lamSamplesRef(idxIn(1));
            b = lamSamplesRef(idxIn(2));
            fprintf('Two representative refined inliers (from sphRef.inliers):\n');
            fprintf('  lambda_a ~ %s\n', qfmt(a));
            fprintf('  lambda_b ~ %s\n\n', qfmt(b));
        end
    else
        fprintf('Refined sphere parametrization: (not available)\n\n');
    end
else
    sphRef = [];
    metaRef = [];
end

% ===================== Isolated eigenvalues AFTER refinement =====================
if ~isempty(idxIso)
    fprintf('Isolated/unclassified DISTINCT samples (refined):\n');
    print_lam_list('isolatedRef', lamSamplesRef(idxIso), resSamplesRef(idxIso));
end

% ===================== Final hint: where to find results =====================
fprintf('Done.\n\n');
fprintf('Variables created by this demo (inspect in the workspace):\n');
fprintf('  A                 : the 4x4 quaternion matrix (AAsph2)\n');
fprintf('  lamAll0, resAll0   : selected DISTINCT candidates (K=%d) and their residuals\n', K);
fprintf('  cls0               : classification labels (0 = isolated/unclassified, >0 = sphere cluster)\n');
if ~isempty(sph0)
    fprintf('  sph0               : sphere detected from sampled candidates (center/radius/inliers)\n');
end
fprintf('  lamRef, resRef     : refined candidates and refined residuals\n');
if exist('sphRef','var') && ~isempty(sphRef)
    fprintf('  sphRef             : sphere refit from refined candidates (center/radius/inliers)\n');
end
fprintf('  lamSamplesFull      : DISTINCT samples returned by sampler (before downselect)\n');
fprintf('  lamSamples0         : DISTINCT samples used in the Supplement (after downselect)\n');
fprintf('  lamSamplesRef       : refined DISTINCT samples (same ordering as lamSamples0)\n');
if ~isempty(idxIso)
    fprintf('  idxIso              : indices of isolated candidates inside lamSamples0/lamSamplesRef\n');
end
if ~isempty(idxSphere)
    fprintf('  idxSphere           : indices of sphere inliers inside lamSamples0/lamSamplesRef\n');
end

% ====================================================================== %
% Local helper functions (MATLAB R2016b+ allows local functions in scripts)
% ====================================================================== %

function A = build_AAsph2()
% Matrix from the Supplement (Example 4-5).
A = quaternion.zeros(4,4);
A(1,:) = [q(  2,  5, -5,  6), q( 12, -3, -5,  4), q(  8, -1, -1, -2), q( 20, -4,  2,  2)];
A(2,:) = [q(  0,  0,  4,  0), q( 10,  4, -2,  4), q(  0,  0,  4,  0), q(  0,  0,  8,  0)];
A(3,:) = [q(  8, -1,  3, -2), q( 28, -5, -3,  0), q(  2,  5, -1,  6), q( 20, -4, 10,  2)];
A(4,:) = [q(  0,  0, -4,  0), q(-20,  4, -6, -2), q(  0,  0, -4,  0), q(-10,  8,-16,  2)];
end

function qq = q(a,b,c,d)
qq = quaternion(a,b,c,d);
end

function [lam1, cert] = local_refine_one(A, lam0, TargetResMin)
% Refine a single candidate via refine_batch (1-element batch).
% Note: do not pass ambiguous Name,Value pairs.

% Prefer refine_batch if available.
if exist('leigqNEWTON_refine_batch','file') == 2
    % NOTE: leigqNEWTON_refine_batch returns [lambdaRef, vRef, cert, info].
    % Requesting only TWO outputs would return vRef as the 2nd output (quaternion),
    % which is a common source of bugs. We explicitly request the cert struct.
    [lamOut, ~, certOut] = leigqNEWTON_refine_batch(A, lam0, ...
        'TargetResMin', TargetResMin, 'DoPolish', true, 'Verbose', 0);
    lam1 = lamOut(1);
    cert = certOut;
else
    % Fallback to auto-refiner
    [lam1, cert] = leigqNEWTON_refine_auto(A, lam0, 'TargetResMin', TargetResMin, 'DoPolish', true, 'Verbose', 0);
end
end

function print_lam_list(titleStr, lam, res)
% Pretty-print a quaternion list with residuals.

lam = lam(:);
res = double(res(:));

fprintf('--- %s ---\n', titleStr);
for k = 1:numel(lam)
    fprintf('%4d)  %s   |  res=%.3e\n', k, qfmt(lam(k)), res(k));
end
fprintf('\n');
end

function s = qfmt(q)
% ASCII-only quaternion formatter.
[w,x,y,z] = parts(q);
vec = [double(w), double(x), double(y), double(z)];
lab = {'', 'i', 'j', 'k'};

s = sprintf('%.12g', vec(1));
for t = 2:4
    if vec(t) >= 0
        s = [s, sprintf(' +%.12g%s', vec(t), lab{t})]; %#ok<AGROW>
    else
        s = [s, sprintf(' %.12g%s', vec(t), lab{t})]; %#ok<AGROW>
    end
end
end

function print_cls_summary(cls)
%PRINT_CLS_SUMMARY  Print a small summary of integer class labels.
%
% leigqNEWTON_sphere_sample returns cls as an integer vector (same length as
% lamSamples). To make this demo robust against accidental output-misalignment,
% shadowed functions on user paths, or older stored structs, we defensively
% extract numeric labels and avoid calling UNIQUE (which calls SORT internally).

[lbl, ok, why] = local_cls_to_numeric(cls);
if ~ok
    fprintf('  (cls summary unavailable)  cls is %s, size %s.  %s\n\n', ...
        class(cls), mat2str(size(cls)), why);
    return;
end

lbl = lbl(:);
K = numel(lbl);
labels = local_unique_int(lbl);

fprintf('  K=%d, unique labels = {', K);
for i = 1:numel(labels)
    if i > 1, fprintf(', '); end
    fprintf('%d', labels(i));
end
fprintf('}\n');

for i = 1:numel(labels)
    li = labels(i);
    fprintf('    label %d: %d\n', li, sum(lbl==li));
end
fprintf('\n');
end

function [lbl, ok, why] = local_cls_to_numeric(cls)
% Try to obtain a numeric label vector from cls (numeric or struct).
ok = false; why = ''; lbl = [];

if isempty(cls)
    ok = true; lbl = []; return;
end

if isnumeric(cls) || islogical(cls)
    ok = true; lbl = double(cls); return;
end

if isstruct(cls)
    % Common field names used in some meta outputs
    cand = {'cls','class','lbl','label','labels'};
    for k = 1:numel(cand)
        f = cand{k};
        if isfield(cls, f)
            v = cls.(f);
            if isnumeric(v) || islogical(v)
                ok = true; lbl = double(v); return;
            end
        end
    end
    why = 'struct without a usable numeric label field (tried cls/class/lbl/label/labels).';
    return;
end

why = 'unsupported cls type.';
end

function labels = local_unique_int(v)
% Compute unique integers without calling unique/sort.

v = double(v(:));
labels = [];
for k = 1:numel(v)
    if isempty(labels)
        labels = v(k);
    else
        found = false;
        for j = 1:numel(labels)
            if labels(j) == v(k)
                found = true;
                break;
            end
        end
        if ~found
            labels(end+1) = v(k); %#ok<AGROW>
        end
    end
end
% keep stable order of first appearance
labels = labels(:).';
end

function print_sphere_info(tag, sph)
% Print sphere struct info in a compact way.

if isempty(sph)
    fprintf('  %s: (empty)\n', tag);
    return
end

c = sph.center;
r = double(sph.radius);
idx = sph.inliers(:);

fprintf('  %s: center=%s   radius=%.12g\n', tag, qfmt(c), r);
fprintf('  %s: inliers count = %d\n', tag, numel(idx));
if ~isempty(idx)
    fprintf('  %s: inliers idx(1:%d) = [', tag, numel(idx));
    fprintf('%d ', idx);
    fprintf(']\n');
end
fprintf('\n');
end

function nrm = qabs(q)
% Quaternion absolute value (Euclidean norm in R^4), vectorized.
[w,x,y,z] = parts(q);
w = double(w); x = double(x); y = double(y); z = double(z);
nrm = sqrt(w.^2 + x.^2 + y.^2 + z.^2);
end

function [lamAll0, lamSamples0, cls0, sph0, note] = downselect_to_collect(...
    lamAllFull, lamSamplesFull, clsFull, sph0Full, resSamplesFull, Collect)

note = '';
lamAll0 = lamAllFull;
lamSamples0 = lamSamplesFull;
cls0 = clsFull;
sph0 = sph0Full;

K = numel(lamSamplesFull);
if K <= Collect
    % no downselect needed
    if ~isempty(sph0)
        sph0.inliers = sph0.inliers(:).';
    end
    return
end

% If no sphere was detected, keep the first Collect samples.
if isempty(sph0)
    keep = (1:Collect).';
    lamSamples0 = lamSamplesFull(keep);
    cls0 = clsFull(keep);
    lamAll0 = lamAllFull(min(keep, numel(lamAllFull)));
    note = sprintf('NOTE: sampler returned K=%d > Collect=%d. No sphere detected, keeping the first %d samples.', K, Collect, Collect);
    return
end

idxIso = find(double(clsFull(:))==0);
idxIn  = sph0.inliers(:);
idxIn  = idxIn(:);

% Ensure indices are within bounds
idxIso = idxIso(idxIso>=1 & idxIso<=K);
idxIn  = idxIn(idxIn>=1 & idxIn<=K);

kIso = numel(idxIso);
if kIso > Collect
    % Extremely unlikely; keep first Collect isolated.
    keep = idxIso(1:Collect);
else
    kNeed = Collect - kIso;

    % Rank inliers by their deviation from the detected sphere model.
    c0 = sph0.center;
    r0 = double(sph0.radius);
    dev = builtin('abs', qabs(lamSamplesFull(idxIn) - c0) - r0);

    % Tie-break by residual (prefer smaller residuals).
    resIn = double(resSamplesFull(idxIn));
    key = dev + resIn;  % dev dominates in this example

    [~, ord] = builtin('sort', key, 'ascend');
    idxInKeep = idxIn(ord(1:min(kNeed, numel(idxIn))));

    keep = [idxIso; idxInKeep];
end

% Keep ascending order for stable printing (avoid unique/sort shadowing issues)
keep = double(keep(:));
keep = keep(keep>=1 & keep<=K);
[keep, ~] = builtin('sort', keep);
keep = keep([true; diff(keep)~=0]);

% If still too many (should not happen), truncate.
if numel(keep) > Collect
    keep = keep(1:Collect);
end

lamSamples0 = lamSamplesFull(keep);
cls0 = clsFull(keep);

% "Candidates" list: in this example we keep the same subset
if numel(lamAllFull) == K
    lamAll0 = lamAllFull(keep);
else
    lamAll0 = lamSamples0;
end

% Reindex sphere inliers to the kept list
maskKeepIn = ismember(keep, idxIn);
idxInNew = find(maskKeepIn);

sph0 = sph0Full;
sph0.inliers = idxInNew(:).';

note = sprintf('NOTE: sampler returned K=%d > Collect=%d. Downselected to K=%d by keeping all %d isolated samples and the best %d sphere inliers.', ...
    K, Collect, numel(keep), numel(idxIso), numel(idxInNew));
end

function [sphRef, metaRef, note] = fit_refined_sphere(lamRef, idxSphereFromStage1)
% Fit a sphere from refined candidates.
% Primary: fitSPHEREfromLambdas on all refined candidates (robustly rejects outliers).
% Fallback: local LS sphere fit from the Stage-1 inliers.

sphRef = [];
metaRef = [];
note = '';

% Try robust fitter if available.
if exist('fitSPHEREfromLambdas','file') == 2
    try
        [sphRef, metaRef] = fitSPHEREfromLambdas(lamRef, ...
            'InlierMode','auto', ...
            'InlierRelTol', 1e-10, ...
            'InlierAbsTol', 0, ...
            'InlierMadFactor', 3, ...
            'DoRefine', true, ...
            'RefineMaxIter', 8, ...
            'SnapInteger', true, ...
            'SnapTol', 1e-9, ...
            'Verbose', false);

        if isfield(metaRef,'snappedInteger') && metaRef.snappedInteger
            note = '  (fitSPHEREfromLambdas: sphere parameters snapped to nearby integers)';
        else
            note = '  (fitSPHEREfromLambdas: sphere refit completed)';
        end
        return

    catch ME
        % Common user issue: a user-defined abs.m shadows MATLAB builtin abs,
        % which breaks internal computations that call abs on doubles.
        note = sprintf(['  (fitSPHEREfromLambdas failed: %s)\n', ...
                        '  Falling back to a local least-squares fit from Stage-1 inliers.\n', ...
                        '  Tip: if you want the robust fitter, remove any shadowing abs.m from your path or run: clear functions; rehash toolboxcache.'], ME.message);
    end
end

% Fallback: fit from the Stage-1 inliers only
idx = idxSphereFromStage1(:);
idx = idx(idx>=1 & idx<=numel(lamRef));

[sphRef, metaRef] = local_fit_sphere_from_inliers(lamRef(idx));

% The fallback fit returns inliers indices local to the passed subset.
% Re-map to indices into lamRef.
if ~isempty(sphRef) && isfield(sphRef,'inliers')
    sphRef.inliers = idx(sphRef.inliers);
end

note = sprintf('%s\n  (fallback fit: LS sphere fit from Stage-1 inliers)', note);
end

function [sph, meta] = local_fit_sphere_from_inliers(lamIn)
% Minimal sphere fit in R^4 for points known to lie on a 2-sphere.
% Steps: PCA to get the affine 3-space, then LS sphere fit in induced R^3.

meta = struct();
sph = [];

lamIn = lamIn(:);
K = numel(lamIn);
if K < 4
    return
end

% Convert to R^4
[w,x,y,z] = parts(lamIn);
P = [double(w), double(x), double(y), double(z)];

% PCA (affine 3-space)
pbar = mean(P,1);
X = P - pbar;
[U,S,V] = svd(X, 'econ'); %#ok<ASGLU>
basis4x3 = V(:,1:3);
Y = X * basis4x3;

% LS sphere in R^3
[c3, r] = local_fitSphere3(Y);

center4 = pbar + c3*basis4x3.';

sph = struct();
sph.center4  = center4(:).';
sph.center   = quaternion(center4(1),center4(2),center4(3),center4(4));
sph.radius   = r;
sph.basis4x3 = basis4x3;
sph.p0       = sph.center4;
sph.center3  = c3;

% All provided points are treated as inliers in this fallback
sph.inliers  = (1:K);

meta.K = K;
meta.dev = builtin('abs', sqrt(sum((Y - c3).^2,2)) - r);
meta.dev_med = median(meta.dev);
meta.dev_max = max(meta.dev);
end

function [c, r] = local_fitSphere3(Y)
% Least-squares sphere fit in R^3.
Y = double(Y);
b = sum(Y.^2, 2);
A = [2*Y, ones(size(Y,1),1)];
sol = A \ b;
c = sol(1:3).';
d = sol(4);

r2 = d + sum(c.^2);
if r2 < 0
    r2 = 0;
end
r = sqrt(r2);
end

function nrm2 = local_norm2_quat_left(A)
% Induced 2-norm of a quaternion matrix via real left-multiplication embedding.
% For A in H^{n x n}, build B in R^{4n x 4n} so that ||A||_2 = ||B||_2.

[n,m] = size(A);
B = zeros(4*n, 4*m);

for i = 1:n
    for j = 1:m
        qij = A(i,j);
        [a,b,c,d] = parts(qij);
        a = double(a); b = double(b); c = double(c); d = double(d);

        L = [ a, -b, -c, -d; ...
              b,  a, -d,  c; ...
              c,  d,  a, -b; ...
              d, -c,  b,  a ];

        ii = (4*(i-1)+1):(4*i);
        jj = (4*(j-1)+1):(4*j);
        B(ii,jj) = L;
    end
end

nrm2 = builtin('norm', B, 2);
end

function restore_warning_states(ids, states)
% Restore warning states captured at the beginning of the script.
for ii = 1:numel(ids)
    try
        warning(states{ii}.state, ids{ii});
    catch
        % ignore
    end
end
end
