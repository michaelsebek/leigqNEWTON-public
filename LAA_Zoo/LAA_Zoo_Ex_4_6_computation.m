%LAA_ZOO_EX_4_6_COMPUTATION  Zoo demo — computation companion for Supplement B (Section 4 Example 6).
%
%   - LaTeX block:  Supplement B, Section 4 Example 6
%   - Script:       LAA_Zoo_Ex_4_6_computation.m
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
%   >> LAA_Zoo_Ex_4_6_computation
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

%% LAA_Zoo_Ex_4_6_computation
% LAA demo: Zoo Example 4-6 (6x6, two spherical components; AA6s2)
%
% Computational companion of Supplement-B, Example:
%   “A 6x6 example with two spherical components in the left spectrum”.
%
% Stages
%   1) Sample distinct candidate left eigenvalues and detect spherical clusters.
%   2) Refine/polish the selected distinct candidates (certificate r_min).
%   3) Refit spheres from the refined candidates (cluster-wise).
%
% Notes
%   - The sampling stage uses leigqNEWTON_sphere_sample(...), which repeatedly
%     calls leigqNEWTON with random starts and collects distinct candidates.
%   - The detection stage fits up to MaxSpheres affine 2-spheres in H~R^4.
%   - Refinement is done point-wise (each lambda separately), so it does not
%     explicitly enforce sphere membership; we therefore also provide an
%     after-refinement sphere refit (Stage 3).

clear; clc;

% ------------------- parameters (match Supplement style) -------------------
UniqueTol     = 1.0e-08;
TolNewton     = 1.0e-12;   % used in refinement tolerances
TargetResMin  = 1.0e-16;   % refinement target

Collect       = 32;        % target number of DISTINCT candidates
RunsMax       = 250;
Restarts      = 3000;
Seed0         = 24680;

MaxSpheres    = 4;         % we expect 2; allow a little slack

% Printing controls (keep console output compact by default)
PRINT_DEBUG         = false;  % set true to print long candidate lists + label histograms
PRINT_REFINE_EVERY  = 8;      % show progress every N refinements (and on the last one)

fprintf('=== LAA demo: Zoo Example 4-6 (6x6, two spheres; AA6s2) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', UniqueTol, TolNewton, TargetResMin);
fprintf('Sampling: Collect=%d distinct | RunsMax=%d | Restarts=%d | Seed0=%d\n\n', Collect, RunsMax, Restarts, Seed0);

% ------------------- disable warnings (restore automatically) --------------
warnState = warning; %#ok<NASGU>
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:illConditionedMatrix');

% ------------------- build the matrix -------------------------------------
A = build_AA6s2();

% ===========================================================================
% Stage 1: sampling + sphere detection
% ===========================================================================
fprintf('Stage 1: sphere sampling/detection ...\n');

% NOTE: The sampler signature has changed across toolbox revisions.
% Current signature is:
%   [lamAll, lamSamples, lambda0, cls, sph, info] = leigqNEWTON_sphere_sample(...)
% We detect the available nargout and map outputs accordingly.
lambda0 = struct();
infoFull = struct();
nout = nargout('leigqNEWTON_sphere_sample');
if nout < 0
    % (should not happen here, but be defensive)
    nout = 6;
end
if nout >= 6
    [lamAllFull, lamSamplesFull, lambda0, clsFull, sphFull, infoFull] = ...
        leigqNEWTON_sphere_sample(A, ...
            'Collect',Collect, ...
            'RunsMax',RunsMax, ...
            'Restarts',Restarts, ...
            'Seed0',Seed0, ...
            'UniqueTol',UniqueTol, ...
            'MaxSpheres',MaxSpheres);
elseif nout == 5
    [lamAllFull, lamSamplesFull, lambda0, clsFull, sphFull] = ...
        leigqNEWTON_sphere_sample(A, ...
            'Collect',Collect, ...
            'RunsMax',RunsMax, ...
            'Restarts',Restarts, ...
            'Seed0',Seed0, ...
            'UniqueTol',UniqueTol, ...
            'MaxSpheres',MaxSpheres);
elseif nout == 4
    [lamAllFull, lamSamplesFull, clsFull, sphFull] = ...
        leigqNEWTON_sphere_sample(A, ...
            'Collect',Collect, ...
            'RunsMax',RunsMax, ...
            'Restarts',Restarts, ...
            'Seed0',Seed0, ...
            'UniqueTol',UniqueTol, ...
            'MaxSpheres',MaxSpheres);
else
    error('Unexpected nargout(leigqNEWTON_sphere_sample)=%d.', nout);
end

% residual certificates for printing + selection
resAllFull     = leigqNEWTON_cert_resMin(A, lamAllFull);
resSamplesFull = leigqNEWTON_cert_resMin(A, lamSamplesFull);

Ktot      = numel(lamAllFull);
Kdistinct = numel(lamSamplesFull);
fprintf('Stage 1 done.  Ktot=%d, Kdistinct=%d\n\n', Ktot, Kdistinct);

if PRINT_DEBUG
    fprintf('Stage 1 results (candidates):\n');
    print_lam_list(lamAllFull, resAllFull, 'lamAll (before refinement)');

    fprintf('\nStage 1 samples (distinct):\n');
    print_lam_list(lamSamplesFull, resSamplesFull, 'lamSamples (before refinement)');

    fprintf('\nStage 1 classification summary (cls):\n');
    print_cls_summary(clsFull);
else
    fprintf('Stage 1 residual summary (distinct samples):\n');
    fprintf('  med(resMin) ~ %.3e\n', median(double(resSamplesFull(:))));
    fprintf('  max(resMin) ~ %.3e\n', max(double(resSamplesFull(:))));
end

% Order detected spheres by radius (ascending) to match the Supplement narrative
sphFullOrd = reorder_spheres_by_radius(sphFull);

fprintf('\nSphere information returned by detection:\n');
    print_sphere_models(sphFullOrd, PRINT_DEBUG);

% Identify isolated and sphere indices in the DISTINCT sample list
clsNumFull  = local_cls_to_numeric(clsFull);
idxIsoFull  = find(clsNumFull==0);
idxSphereFull = find(clsNumFull>0);

% Reorder (and if needed, downselect) to match the Supplement narrative:
%   Sphere #1 (smaller radius) then Sphere #2, then isolated points.
% If more than Collect were returned, keep the best (lowest-residual) points.
targetSph = [17 13];  % desired inliers on sphere #1 and #2
targetIso = 2;        % desired isolated points

Ksel = min(Collect, Kdistinct);

% pick isolated (best residuals if more than needed)
idxIsoPick = idxIsoFull(:);
if ~isempty(idxIsoPick)
    [~,ordIso] = sort(resSamplesFull(idxIsoPick),'ascend');
    idxIsoPick = idxIsoPick(ordIso);
    idxIsoPick = idxIsoPick(1:min(targetIso, numel(idxIsoPick)));
end

% pick per-sphere inliers (best residuals)
idxSphPick = [];
for s = 1:min(2, numel(sphFullOrd))
    if ~isfield(sphFullOrd(s),'inliers') || isempty(sphFullOrd(s).inliers)
        continue
    end
    inl = sphFullOrd(s).inliers(:);
    [~,ordInl] = sort(resSamplesFull(inl),'ascend');
    inl = inl(ordInl);
    take = min(targetSph(s), numel(inl));
    idxSphPick = [idxSphPick; inl(1:take)]; %#ok<AGROW>
end

% desired order: sphere #1, sphere #2, isolated
idxKeep = unique([idxSphPick; idxIsoPick], 'stable');

% fill up to Ksel with best remaining (by residual)
if numel(idxKeep) < Ksel
    rest = setdiff((1:Kdistinct).', idxKeep, 'stable');
    [~,ord2] = sort(resSamplesFull(rest),'ascend');
    rest = rest(ord2);
    idxKeep = [idxKeep; rest(1:min(Ksel-numel(idxKeep), numel(rest)))];
elseif numel(idxKeep) > Ksel
    idxKeep = idxKeep(1:Ksel);
end

if Kdistinct ~= Ksel
    fprintf('\nNOTE: sampler returned K=%d distinct; using K=%d after selection/reordering.\n\n', Kdistinct, Ksel);
end

lamSamples0 = lamSamplesFull(idxKeep);
resSamples0 = resSamplesFull(idxKeep);
cls0        = clsNumFull(idxKeep);

% rebuild sphere models with indices re-mapped to the selected list
sph0 = remap_spheres_to_subset(sphFullOrd, idxKeep);

K0 = numel(lamSamples0);

fprintf('\nStage 1 (selected for Supplement): K=%d distinct candidates\n', K0);
if PRINT_DEBUG
    print_lam_list(lamSamples0, resSamples0, 'lamAll0 (selected distinct candidates)');
end

% Indices by label after selection
labels = unique(cls0(:)).';
idxIso = find(cls0==0);
idxSph1 = [];
idxSph2 = [];
if numel(labels(labels>0)) >= 1
    idxSph1 = find(cls0==labels(labels>0));
end

% Use sphere models (preferred) to define inliers for sphere-#1 and sphere-#2
idxSphereByModel = cell(1,numel(sph0));
for s=1:numel(sph0)
    if isfield(sph0(s),'inliers')
        idxSphereByModel{s} = sph0(s).inliers(:);
    else
        idxSphereByModel{s} = [];
    end
end

if numel(idxSphereByModel) >= 1, idxSph1 = idxSphereByModel{1}; end
if numel(idxSphereByModel) >= 2, idxSph2 = idxSphereByModel{2}; end

fprintf('\nDetected spherical clusters (selected set):\n');
print_sphere_models(sph0, PRINT_DEBUG);

fprintf('\nIsolated/unclassified DISTINCT samples (before refinement):\n');
print_lam_list(lamSamples0(idxIso), resSamples0(idxIso), 'isolated');

% ------------------- Supplement-style summary block (Stage 1) --------------
fprintf('\n--- Supplement-style summary (from Stage 1 candidates) ---\n');

normA = local_norm2_quat_left(A);
medRes = median(resSamples0);
maxRes = max(resSamples0);

fprintf('Ktot=%d distinct candidates.\n', K0);
if numel(sph0) >= 1
    fprintf('Sphere #1:  c1=%s,  r1=%.15g,  #inliers=%d\n', qfmt(sph0(1).center), sph0(1).radius, numel(sph0(1).inliers));
end
if numel(sph0) >= 2
    fprintf('Sphere #2:  c2=%s,  r2=%.15g,  #inliers=%d\n', qfmt(sph0(2).center), sph0(2).radius, numel(sph0(2).inliers));
end

if ~isempty(idxIso)
    fprintf('Isolated/outliers: %d\n', numel(idxIso));
    for t=1:numel(idxIso)
        fprintf('  iso #%d:  %s\n', t, qfmt(lamSamples0(idxIso(t))));
    end
end

fprintf('Residual certificates (Stage 1): medRes~%.3e, maxRes~%.3e\n', medRes, maxRes);
fprintf('||A||_2 (real left-embedding) ~ %.5g\n', normA);
fprintf('Normalized: medRes/||A||_2~%.3e, maxRes/||A||_2~%.3e\n', medRes/normA, maxRes/normA);

% Representative inliers (Stage 1)
if ~isempty(idxSph1)
    ia = idxSph1(1);
    ib = idxSph1(min(2,numel(idxSph1)));
    fprintf('Sphere #1 representative inliers: %s  and  %s\n', qfmt(lamSamples0(ia)), qfmt(lamSamples0(ib)));
end
if ~isempty(idxSph2)
    ia = idxSph2(1);
    ib = idxSph2(min(2,numel(idxSph2)));
    fprintf('Sphere #2 representative inliers: %s  and  %s\n', qfmt(lamSamples0(ia)), qfmt(lamSamples0(ib)));
end

% Parameterization hints (detected fixed coordinate)
if numel(sph0) >= 1 && ~isempty(idxSph1)
    [fixedDir, fixedVal] = local_detect_fixed_dir(lamSamples0(idxSph1));
    fprintf('Sphere #1 supporting affine 3-space: fixed %s-coefficient ~ %.15g\n', fixedDir, fixedVal);
end
if numel(sph0) >= 2 && ~isempty(idxSph2)
    [fixedDir, fixedVal] = local_detect_fixed_dir(lamSamples0(idxSph2));
    fprintf('Sphere #2 supporting affine 3-space: fixed %s-coefficient ~ %.15g\n', fixedDir, fixedVal);
end

% ===========================================================================
% Stage 2: refinement / polish
% ===========================================================================
fprintf('\nStage 2: refinement/polish of candidates ...\n');

verbRefine = 0;
if PRINT_DEBUG
    verbRefine = 1;
end

[lamRef, vRef, certRef] = leigqNEWTON_refine_batch(A, lamSamples0, ...
    'TargetResMin',TargetResMin, ...
    'Mode','auto', ...
    'DoPolish',true, ...
    'TolX',TolNewton, ...
    'TolFun',TolNewton, ...
    'Verbose',verbRefine);

% certRef should be a struct with vector field resMin; be defensive
resRef = NaN(size(lamRef));
if isstruct(certRef) && isfield(certRef,'resMin')
    resRef = double(certRef.resMin(:));
elseif iscell(certRef)
    % older variants: certRef{k}.resMin
    for k=1:numel(certRef)
        if isstruct(certRef{k}) && isfield(certRef{k},'resMin')
            resRef(k) = double(certRef{k}.resMin);
        end
    end
end

fprintf('Stage 2 done.\n\n');

fprintf('Stage 2 results (refined candidates):\n');
if PRINT_DEBUG
    print_lam_list(lamRef, resRef, 'lamRef (after Stage 2)');
end

fprintf('\nStage 2: residual summary (refined candidates)\n');
fprintf('  med(resMin) ~ %.3e\n', median(resRef));
fprintf('  max(resMin) ~ %.3e\n', max(resRef));
fprintf('  ||A||_2 (real left-embedding) ~ %.5g\n', normA);

% For convenience, keep the refined DISTINCT samples aligned with Stage 1
lamSamplesRef = lamRef;

fprintf('\nRefined DISTINCT samples (mapped to Stage 2 refinements):\n');
if PRINT_DEBUG
    print_lam_list(lamSamplesRef, resRef, 'lamSamplesRef (after Stage 2)');
else
    % Keep the output compact: show only the refined isolated/outlier candidates.
    if ~isempty(idxIso)
        fprintf('Isolated/outlier refined candidates:\n');
        for kk = 1:numel(idxIso)
            k = idxIso(kk);
            fprintf('  iso#%d: %s   |  res=%.3e\n', kk, qfmt(lamSamplesRef(k)), resRef(k));
        end
    end
end

% ===========================================================================
% Stage 3: sphere refit AFTER refinement (cluster-wise)
% ===========================================================================
fprintf('\nSphere fit AFTER refinement (cluster-wise) ...\n');

sphRef = struct([]);
for s=1:numel(sph0)
    idx = sph0(s).inliers(:);
    if isempty(idx)
        continue;
    end

    [sp, meta] = fitSPHEREfromLambdas(lamSamplesRef(idx), 'SnapInteger',true, 'Verbose',false);

    % map inliers back to global indices
    idxIn = idx(sp.inliers(:));

    sphRef(end+1).center = sp.center; %#ok<SAGROW>
    sphRef(end).radius   = sp.radius;
    sphRef(end).inliers  = idxIn(:);
    sphRef(end).meta     = meta;

    [medDev, maxDev] = local_sphere_consistency(lamSamplesRef(idxIn), sp.center, sp.radius);

    fprintf('Sphere #%d refined fit: c=%s, r=%.15g, #inliers=%d\n', ...
        numel(sphRef), qfmt(sp.center), sp.radius, numel(idxIn));
    fprintf('  consistency (on refined inliers): median | ||lam-c|| - r | ~ %.3e,  max ~ %.3e\n', medDev, maxDev);
end

% reorder refined spheres by radius to match Stage 1 ordering
sphRef = reorder_spheres_by_radius(sphRef);

fprintf('\nRefined sphere parametrizations (sphRef):\n');
print_sphere_models(sphRef, PRINT_DEBUG);

% Show representative refined inliers for each sphere
for s=1:numel(sphRef)
    idx = sphRef(s).inliers(:);
    if isempty(idx)
        continue;
    end
    a = idx(1);
    b = idx(min(2,numel(idx)));
    fprintf('\nTwo representative refined inliers on Sphere #%d:\n', s);
    fprintf('  lambda_a ~ %s\n', qfmt(lamSamplesRef(a)));
    fprintf('  lambda_b ~ %s\n', qfmt(lamSamplesRef(b)));
end

fprintf('\nDone.\n\n');

% ------------------- workspace summary ------------------------------------
fprintf('Variables created by this demo (inspect in the workspace):\n');
fprintf('  A                : the 6x6 quaternion matrix (AA6s2)\n');
fprintf('  lamSamples0,resSamples0 : selected DISTINCT candidates and residuals (Stage 1)\n');
fprintf('  cls0              : numeric classification labels (0=isolated, 1/2=two spheres)\n');
fprintf('  sph0              : sphere models detected from Stage 1 candidates (two spheres expected)\n');
fprintf('  lamRef,resRef     : refined candidates and refined residual certificates (Stage 2)\n');
fprintf('  sphRef            : spheres refit from refined candidates (Stage 3)\n');


%% ========================================================================
% Local functions
% ========================================================================

function A = build_AA6s2()
% Matrix from Supplement-B example “two spheres” (integer quaternion matrix).

q = @(a,b,c,d) quaternion(a,b,c,d);
z = q(0,0,0,0);

A = repmat(z,6,6);

A(1,:) = [ q(-14,-12, 12, 19),  q(  0, -2, -3, -5),  z, ...
           q(  6,  2,  6, -5),  q( -6, -4,  7,  4),  q( 12,  4, -4,-18) ];

A(2,:) = [ q( -6, -2, 18, 11),  q(  4, -4,  1, -5),  z, ...
           q( 18,  6,-22,-25),  q(  6,  4,  1, -4),  q( 12,  4,-20,-22) ];

A(3,:) = [ q( -7,-12, 15, 16),  q(  0, -2, -3, -5),  q( -6, 10,  0, -4), ...
           q(  7, 12, -7,-12),  q( -5,  6,  2, -3),  q( 12,  4,-12,-18) ];

A(4,:) = [ q(  0,  0,  0,  2),  q(  0,  0,  0, -2),  z, ...
           q(  4, -6, -2,-10),  z,                   q(  0,  0,  0, -4) ];

A(5,:) = [ q(  6,  2,-18,-11),  q(  0,  2,  3,  7),  z, ...
           q(-18, -6, 22, 25),  q( -2, -6,  3,  6),  q(-12, -4, 20, 22) ];

A(6,:) = [ q( -6, -2, -2,  7),  q(  0, -2, -3, -3),  z, ...
           q( -6, -2, 14, 15),  q( -6, -4,  3,  4),  q(  4, -6,  6, -4) ];

end

function print_lam_list(lam, res, titleStr)
if nargin < 3, titleStr = 'lambda list'; end
fprintf('--- %s ---\n', titleStr);
K = numel(lam);
for k=1:K
    if nargin >= 2 && ~isempty(res)
        fprintf('%4d)  %s   |  res=%.3e\n', k, qfmt(lam(k)), res(k));
    else
        fprintf('%4d)  %s\n', k, qfmt(lam(k)));
    end
end
end

function print_cls_summary(cls)
cls = local_cls_to_numeric(cls);
cls = cls(:);
labels = unique(cls).';
fprintf('  K=%d, unique labels = {', numel(cls));
for t=1:numel(labels)
    if t>1, fprintf(', '); end
    fprintf('%d', labels(t));
end
fprintf('}\n');
for t=1:numel(labels)
    fprintf('    label %d: %d\n', labels(t), nnz(cls==labels(t)));
end
end

function clsNum = local_cls_to_numeric(cls)
% Robust conversion of classification labels to numeric vector.
if isnumeric(cls)
    clsNum = cls;
    return;
end

% Some variants return a struct array or cell array.
if iscell(cls)
    try
        clsNum = cellfun(@double, cls);
    catch
        clsNum = zeros(numel(cls),1);
        for k=1:numel(cls)
            try, clsNum(k)=double(cls{k}); end
        end
    end
    return;
end

if isstruct(cls)
    % Try common fields
    f = fieldnames(cls);
    if ismember('label',f)
        clsNum = [cls.label].';
        return;
    end
    if ismember('cls',f)
        clsNum = [cls.cls].';
        return;
    end
end

% Fallback
try
    clsNum = double(cls);
catch
    clsNum = zeros(numel(cls),1);
end
end

function sphOrd = reorder_spheres_by_radius(sph)
if isempty(sph)
    sphOrd = sph;
    return;
end
r = zeros(numel(sph),1);
for s=1:numel(sph)
    if isfield(sph(s),'radius')
        r(s) = double(sph(s).radius);
    else
        r(s) = NaN;
    end
end
[~,ord] = sort(r,'ascend');
sphOrd = sph(ord);
end

function print_sphere_models(sph, showIdx)
if nargin < 2
    showIdx = false;
end
if isempty(sph)
    fprintf('  (none)\n');
    return;
end
for s=1:numel(sph)
    c = sph(s).center;
    r = sph(s).radius;
    nin = 0;
    if isfield(sph(s),'inliers') && ~isempty(sph(s).inliers)
        nin = numel(sph(s).inliers);
    end
    fprintf('  sph%d: center=%s   radius=%.15g\n', s, qfmt(c), r);
    if isfield(sph(s),'inliers') && ~isempty(sph(s).inliers)
        idx = sph(s).inliers(:).';
        fprintf('  sph%d: inliers count = %d\n', s, numel(idx));
        if showIdx
            fprintf('  sph%d: inliers idx(1:%d) = [', s, numel(idx));
            fprintf('%d ', idx);
            fprintf(']\n');
        end
    end
    fprintf('\n');
end
end

function sphSub = remap_spheres_to_subset(sph, idxKeep)
% idxKeep are indices into the original lamSamplesFull.
% We need to map inlier indices to positions in the reduced list.

sphSub = sph;
for s=1:numel(sph)
    if ~isfield(sph(s),'inliers') || isempty(sph(s).inliers)
        sphSub(s).inliers = [];
        continue;
    end
    inl = sph(s).inliers(:);
    [tf,loc] = ismember(inl, idxKeep);
    sphSub(s).inliers = loc(tf);
end
end

function s = qfmt(q)
% Compact display with MATLAB's quaternion formatting.
try
    s = char(q);
catch
    try
        % Aerospace quaternion: string(q) may work
        s = char(string(q));
    catch
        s = '<quaternion>';
    end
end
% Make spacing slightly tighter
s = strrep(s,'  ',' ');
end

function n2 = local_norm2_quat_left(A)
% Spectral norm of the real left-embedding of a quaternion matrix.
% This is a standard real embedding H^{n x n} -> R^{4n x 4n}.

n = size(A,1);
AR = zeros(4*n,4*n);
for r=1:n
    for c=1:n
        [a0,a1,a2,a3] = parts(A(r,c));
        % left multiplication matrix of a0 + a1 i + a2 j + a3 k
        L = [ a0, -a1, -a2, -a3;
              a1,  a0, -a3,  a2;
              a2,  a3,  a0, -a1;
              a3, -a2,  a1,  a0];
        rr = (4*(r-1)+1):(4*r);
        cc = (4*(c-1)+1):(4*c);
        AR(rr,cc) = L;
    end
end
n2 = norm(AR,2);
end

function [fixedDir, fixedVal] = local_detect_fixed_dir(lam)
% Heuristic: find which imaginary coordinate is (approximately) constant
% across a set of quaternions.

[a0,a1,a2,a3] = parts(lam);
S = [std(double(a1(:))), std(double(a2(:))), std(double(a3(:)))];
[~,k] = min(S);
dirs = {'i','j','k'};
fixedDir = dirs{k};
vals = [mean(double(a1(:))), mean(double(a2(:))), mean(double(a3(:)))];
fixedVal = vals(k);
end

function [medDev, maxDev] = local_sphere_consistency(lam, c, r)
% Compute | ||lam-c|| - r | over a set of quaternions.

[lc0,lc1,lc2,lc3] = parts(lam);
[cc0,cc1,cc2,cc3] = parts(c);

d0 = double(lc0) - double(cc0);
d1 = double(lc1) - double(cc1);
d2 = double(lc2) - double(cc2);
d3 = double(lc3) - double(cc3);

dist = sqrt(d0.^2 + d1.^2 + d2.^2 + d3.^2);
dev  = abs(dist - double(r));

medDev = median(dev(:));
maxDev = max(dev(:));
end

