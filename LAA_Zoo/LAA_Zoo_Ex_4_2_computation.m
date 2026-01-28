%LAA_ZOO_EX_4_2_COMPUTATION  Zoo demo — computation companion for Supplement B (Section 4 Example 2).
%
%   - LaTeX block:  Supplement B, Section 4 Example 2
%   - Script:       LAA_Zoo_Ex_4_2_computation.m
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
%   >> LAA_Zoo_Ex_4_2_computation
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

clc;
clearvars;

UniqueTol     = 1.0e-8;
TolNewton     = 1.0e-12;   % printed only (NOT passed to sphere_detect / sphere_sample)
TargetResMin  = 1.0e-16;

Collect       = 20;
RunsMax       = 120;
Restarts      = 1500;      % sphere_sample option name is 'Restarts'
Seed0         = 24680;

DoStage3      = true;

fprintf('=== LAA demo: Zoo Example 4--2 (4x4, spherical eigenvalues; AAsph2) ===\n');
fprintf('UniqueTol=%.1e | TolNewton=%.1e | TargetResMin=%.1e\n', UniqueTol, TolNewton, TargetResMin);
fprintf('Sampling: Collect=%d distinct | RunsMax=%d | Restarts=%d | Seed0=%d\n\n', ...
    Collect, RunsMax, Restarts, Seed0);

A = build_AAsph2();

% ---------------- Stage 1 ----------------
fprintf('Stage 1: sphere sampling/detection ...\n');
[lamAll, lamSamples, lambda0, cls, sph0, info0] = leigqNEWTON_sphere_detect( ...
    A, ...
    'Collect', Collect, ...
    'RunsMax', RunsMax, ...
    'Restarts', Restarts, ...
    'Seed0', Seed0, ...
    'UniqueTol', UniqueTol);

Ktot      = numel(lamAll);
Kdistinct = numel(lamSamples);
fprintf('Stage 1 done.  Ktot=%d, Kdistinct=%d\n\n', Ktot, Kdistinct);

resAll0     = leigqNEWTON_cert_resMin(A, lamAll);
resSamples0 = leigqNEWTON_cert_resMin(A, lamSamples);

fprintf('Stage 1 results (candidates):\n');
print_q_list('lamAll (before refinement)', lamAll, resAll0);

fprintf('\nStage 1 samples (distinct):\n');
print_q_list('lamSamples (before refinement)', lamSamples, resSamples0);

fprintf('\nStage 1 classification summary (cls):\n');
print_cls_summary(cls);

fprintf('\nSphere information returned by detection:\n');
print_sphere_struct('sph0', sph0);

[idxSphere, idxIso] = split_sphere_isolated(cls);

if ~isempty(idxIso)
    fprintf('\nUnclassified DISTINCT samples (before refinement):\n');
    print_q_list('unclassified', lamSamples(idxIso), resSamples0(idxIso));
end

% ---------------- Stage 2 ----------------
fprintf('\nStage 2: refinement/polish of candidates ...\n');

if exist('leigqNEWTON_refine_batch','file') ~= 2
    error('This demo assumes leigqNEWTON_refine_batch is on path.');
end

lamRef = lamAll;
resRef = nan(size(resAll0));

tStage2 = tic;
for k = 1:Ktot
    tk = tic;

    % refine_batch signature: [lambdaRef, vRef, cert, info]
    [lam_k, ~, cert_k] = leigqNEWTON_refine_batch( ...
        A, lamAll(k), ...
        'TargetResMin', TargetResMin, ...
        'Verbose', 0);

    lamRef(k) = lam_k(1);

    if ~isfield(cert_k,'resMin') || ~isnumeric(cert_k.resMin)
        error('refine_batch returned unexpected cert.resMin for k=%d', k);
    end
    resRef(k) = cert_k.resMin(1);

    fprintf('  Refinement %d/%d ... done (dt=%.2fs, r~%.3e)\n', k, Ktot, toc(tk), resRef(k));
end
fprintf('Stage 2 done. elapsed=%.1fs\n\n', toc(tStage2));

fprintf('Stage 2 results (refined candidates):\n');
print_q_list('lamRef (after Stage 2)', lamRef, resRef);

fprintf('\nStage 2: residual summary (refined candidates)\n');
fprintf('  med(resMin) ~ %.3e\n', median(resRef));
fprintf('  max(resMin) ~ %.3e\n', max(resRef));

% Map DISTINCT samples to their refined counterparts (using UniqueTol in R^4)
lamSamplesRef = map_samples_to_refined(lamSamples, lamAll, lamRef, UniqueTol);
resSamplesRef = leigqNEWTON_cert_resMin(A, lamSamplesRef);

fprintf('\nRefined DISTINCT samples (mapped to Stage 2 refinements):\n');
print_q_list('lamSamplesRef (after Stage 2)', lamSamplesRef, resSamplesRef);

if has_center_radius(sph0) && ~isempty(idxSphere)
    c = sph0.center;
    r = sph0.radius;
    qmag = qabs(lamSamplesRef(idxSphere) - c);
    dev  = builtin('abs', qmag - r);
fprintf('\nSphere consistency on refined inliers (using sph0 center/radius):\n');
    fprintf('  median | ||lam-c|| - r | ~ %.3e\n', median(dev));
    fprintf('  max    | ||lam-c|| - r | ~ %.3e\n', max(dev));
end

if ~isempty(idxIso)
    fprintf('\nUnclassified DISTINCT samples (refined):\n');
    print_q_list('unclassified', lamSamplesRef(idxIso), resSamplesRef(idxIso));
end

% ---------------- Stage 3 ----------------
if DoStage3
    fprintf('\nStage 3: sphere validation (optional) ...\n');
    fprintf('  Validating/refining detected sphere (may take long) ');

    beatPeriod = 15;
    beat = make_heartbeat_timer(beatPeriod);
    start(beat.timerObj);
    tStage3 = tic;

    [A2, lamAll2, resAll2, lamSamples2, resSamples2, sph2, conf2, out2] = ...
        leigqNEWTON_sphere_validate( ...
            A, lamRef, lamSamplesRef, ...
            'TargetRes', TargetResMin, ...
            'Verbose', 1);

    stop(beat.timerObj);
    delete(beat.timerObj);
    fprintf(' done (%.1fs)\n\n', toc(tStage3));

    fprintf('Sphere parametrization AFTER validation:\n');
    print_sphere_struct('sph2', sph2);

    fprintf('\nStage 3 results (validated/refined):\n');
    print_q_list('lamAll2 (after Stage 3)', lamAll2, resAll2);

    % push to workspace
    assignin('base','Asph',A);
    assignin('base','lamAll',lamAll);
    assignin('base','lamSamples',lamSamples);
    assignin('base','cls',cls);
    assignin('base','sph0',sph0);
    assignin('base','info0',info0);
    assignin('base','resAll0',resAll0);
    assignin('base','resSamples0',resSamples0);
    assignin('base','lamRef',lamRef);
    assignin('base','resRef',resRef);
    assignin('base','A2',A2);
    assignin('base','lamAll2',lamAll2);
    assignin('base','resAll2',resAll2);
    assignin('base','lamSamples2',lamSamples2);
    assignin('base','resSamples2',resSamples2);
    assignin('base','sph2',sph2);
    assignin('base','conf2',conf2);
    assignin('base','out2',out2);
else
    fprintf('\nStage 3 skipped (DoStage3=false).\n');
end

fprintf('\nDone.\n');

% =====================================================================
% local helpers
% =====================================================================

function A = build_AAsph2()
q = @(a,b,c,d) quaternion(a,b,c,d);
A = repmat(q(0,0,0,0), 4,4);

A(1,1) = q(  2,  5, -5,  6);
A(1,2) = q( 12, -3, -5,  4);
A(1,3) = q(  8, -1, -1, -2);
A(1,4) = q( 20, -4,  2,  2);

A(2,1) = q(  0,  0,  4,  0);
A(2,2) = q( 10,  4, -2,  4);
A(2,3) = q(  0,  0,  4,  0);
A(2,4) = q(  0,  0,  8,  0);

A(3,1) = q(  8, -1,  3, -2);
A(3,2) = q( 28, -5, -3,  0);
A(3,3) = q(  2,  5, -1,  6);
A(3,4) = q( 20, -4, 10,  2);

A(4,1) = q(  0,  0, -4,  0);
A(4,2) = q(-20,  4, -6, -2);
A(4,3) = q(  0,  0, -4,  0);
A(4,4) = q(-10,  8,-16,  2);
end

function print_q_list(titleStr, lam, res)
fprintf('--- %s ---\n', titleStr);
for k = 1:numel(lam)
    fprintf(' %3d)  %s   |  res=%.3e\n', k, fmt_q(lam(k)), res(k));
end
end

function s = fmt_q(q)
[w,x,y,z] = parts(q);
s = sprintf('%+.12g %+.12gi %+.12gj %+.12gk', w, x, y, z);
s = strrep(s,'+ -','- ');
end

function print_cls_summary(cls)
%PRINT_CLS_SUMMARY  Print a compact summary of cluster labels (no UNIQUE dependency).
%
% UNIQUE calls SORT internally, and some installations shadow MATLAB's SORT
% on the path. We therefore call BUILTIN('sort',...) explicitly and compute
% unique labels via DIFF.
%
% In leigqNEWTON_sphere_sample/_detect, cls is K-by-1 int32 labels for lamSamples:
%   0 = unclassified (not on any detected sphere)
%   s = sphere index (1,2,...) for inliers of detected sphere(s)

if isempty(cls)
    fprintf('  (cls empty)\n');
    return;
end

% If a struct is passed accidentally, try extracting a numeric label field.
if isstruct(cls)
    cand = {'cls','class','label','labels','id','idx'};
    found = false;
    for i = 1:numel(cand)
        if isfield(cls, cand{i}) && isnumeric(cls.(cand{i}))
            cls = cls.(cand{i});
            found = true;
            break;
        end
    end
    if ~found
        fprintf('  cls: struct with fields: %s\n', strjoin(fieldnames(cls).', ', '));
        return;
    end
end

try
    lbl = double(cls(:));
catch
    fprintf('  cls: unsupported type %s\n', class(cls));
    return;
end

lbls = builtin('sort', lbl);
isNew = [true; diff(lbls) ~= 0];
u = lbls(isNew);

fprintf('  K=%d, unique labels = {', numel(lbl));
for i = 1:numel(u)
    if i>1, fprintf(', '); end
    fprintf('%g', u(i));
end
fprintf('}\n');

for i = 1:numel(u)
    fprintf('    label %g: %d\n', u(i), sum(lbl == u(i)));
end
end

function tf = has_center_radius(sph)
tf = isstruct(sph) && isfield(sph,'center') && isfield(sph,'radius') && ~isempty(sph.center) && ~isempty(sph.radius);
end

function print_sphere_struct(name, sph)
%PRINT_SPHERE_STRUCT Print sphere struct (ASCII only, robust).
%
% Notes:
%   - Avoids fprintf repetition bugs when sph.inliers is a vector.
%   - Prints count + first indices (up to 30).

    if ~isstruct(sph) || isempty(fieldnames(sph))
        fprintf('  %s: (empty)\n', name);
        return;
    end

    if has_center_radius(sph)
        fprintf('  %s: center=%s   radius=%.12g\n', name, fmt_q(sph.center), sph.radius);
    else
        fprintf('  %s: (no center/radius fields found)\n', name);
    end

    if isfield(sph,'inliers')
        inl = sph.inliers;
        if islogical(inl)
            idx = find(inl);
        else
            idx = inl(:);
        end

        if isempty(idx)
            fprintf('  %s: inliers count = 0\n', name);
        else
            fprintf('  %s: inliers count = %d\n', name, numel(idx));
            m = min(numel(idx), 30);
            fprintf('  %s: inliers idx(1:%d) = %s\n', name, m, mat2str(idx(1:m).'));
            if numel(idx) > m
                fprintf('  %s: ... (%d more)\n', name, numel(idx)-m);
            end
        end
    end
end

function [idxSphere, idxIso] = split_sphere_isolated(cls)
%SPLIT_SPHERE_ISOLATED Split DISTINCT samples into sphere inliers vs unclassified.
%
% cls is K-by-1 int32 (labels for lamSamples):
%   idxSphere = indices where cls>0  (inliers of detected sphere(s))
%   idxIso    = indices where cls==0 (remaining/unclassified samples)

idxSphere = [];
idxIso    = [];
if isempty(cls), return; end

if isstruct(cls)
    % Try common numeric field names
    cand = {'cls','class','label','labels','id','idx'};
    for i = 1:numel(cand)
        if isfield(cls, cand{i}) && isnumeric(cls.(cand{i}))
            cls = cls.(cand{i});
            break;
        end
    end
end

try
    c = double(cls(:));
catch
    return;
end

idxSphere = find(c > 0);
idxIso    = find(c == 0);
end

function lamSamplesRef = map_samples_to_refined(lamSamples, lamAll, lamRef, UniqueTol)
if isempty(lamSamples)
    lamSamplesRef = lamSamples;
    return;
end
Ws = quat_parts_matrix(lamSamples);
Wa = quat_parts_matrix(lamAll);

idx = nan(numel(lamSamples),1);
for i = 1:numel(lamSamples)
    d2 = sum((Wa - Ws(i,:)).^2, 2);
    [m, j] = min(d2);
    if sqrt(m) <= UniqueTol
        idx(i) = j;
    end
end

if all(isfinite(idx))
    lamSamplesRef = lamRef(idx);
else
    lamSamplesRef = lamSamples;
end
end

function M = quat_parts_matrix(qv)
qv = qv(:);
M = zeros(numel(qv),4);
for k = 1:numel(qv)
    [w,x,y,z] = parts(qv(k));
    M(k,:) = [w x y z];
end

end

function mag = qabs(q)
%QABS Magnitude of quaternion array (robust; avoids calling abs()).
%   mag = sqrt(w^2 + x^2 + y^2 + z^2) computed via parts().
    q = q(:);
    mag = zeros(size(q));
    for ii = 1:numel(q)
        [w,x,y,z] = parts(q(ii));
        mag(ii) = sqrt(w*w + x*x + y*y + z*z);
    end
end
function beat = make_heartbeat_timer(periodSeconds)
beat.count = 0;
beat.timerObj = timer( ...
    'ExecutionMode','fixedSpacing', ...
    'Period', periodSeconds, ...
    'TimerFcn', @tick);
    function tick(~,~)
        beat.count = beat.count + 1;
        fprintf('.');
        if mod(beat.count,60)==0
            fprintf('\n  ');
        end
    end
end
