function [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = leigqNEWTON_sphere_refine(cases, k, varargin)
%LEIGQNEWTON_SPHERE_REFINE  Refine/validate sphere candidates and decide sphere-vs-artifact (advanced engine).
%
%   This routine is the computational engine behind LEIGQNEWTON_SPHERE_VALIDATE.
%   In the public release, LEIGQNEWTON_SPHERE_VALIDATE is a thin wrapper that
%   forwards all inputs/outputs to this function.
%
%   The function refines sampled left-eigenvalue candidates, computes residual
%   certificates, and (when a sphere model is provided) performs optional grid-based
%   checks and refitting to decide whether a detected sphere is likely genuine
%   or an artifact.
%
%   IMPORTANT (sphere model required for a sphere verdict)
%   -----------------------------------------------------
%   A sphere decision requires a sphere model in the input case, typically C.sph,
%   produced by LEIGQNEWTON_SPHERE_DETECT or LEIGQNEWTON_SPHERE_SAMPLE.
%   If no sphere model is available in the input, the routine can still refine
%   eigenvalues and certificates, but it will return sph=[] and the verdict will be
%   NO_SPHERE_MODEL (i.e., "nothing to validate").
%
% Syntax (recommended)
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_refine(C, Name,Value,...)
%
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_refine(cases, k, Name,Value,...)
%
% Convenience positional form (refinement-only; no sphere verdict unless sph provided via C)
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_refine(A, lamAll, lamSamples, Name,Value,...)
%
% Inputs
%   C / cases     Case container (struct, struct array, or cell array).
%                Each case should contain at least:
%                   C.A          (n-by-n quaternion)
%                   C.lamAll     (M-by-1 quaternion candidates)
%                   C.lamSamples (K-by-1 quaternion DISTINCT samples; may be [])
%                For sphere validation/decision, also provide:
%                   C.sph        sphere model struct (from sphere_detect/sample)
%                Additional fields are allowed (e.g., cls, info, seed, idx, origin).
%
%   k            Case index when passing a list (struct array / cell array).
%                Omit k when passing a single case struct C.
%
% Positional-form inputs
%   A            n-by-n quaternion matrix.
%   lamAll       M-by-1 quaternion, candidate list.
%   lamSamples   K-by-1 quaternion, DISTINCT samples (may be []).
%                NOTE: this form has no place to pass a sphere model; use a case struct
%                with field C.sph if you want a sphere verdict.
%
% Name,Value options (selected)
%   'RefineArgs'        Cell array of Name,Value passed to LEIGQNEWTON_REFINE_AUTO.
%                       Example: {'DoPolish',false,'TargetRes',1e-14}
%   'TargetRes'         Target residual used by the refinement/certification logic
%                       (default 1e-14).
%   'VerifySphere'      Enable grid-based sphere checks when a sphere model is present
%                       (default true).
%   'SphereIndex'       Which sphere model to validate when multiple are present
%                       (default 1).
%   'SphereTolRel'      Relative tolerance used in sphere validation (default 5e-3).
%   'GridTheta','GridPhi'  Grid resolution for sphere checks (defaults 8 and 16).
%   'GridMaxPoints'     Cap on grid points (default 256).
%   'MinGridForDecision' Minimum grid size for a conclusive decision (default 64).
%   'RefineGrid'        Refine grid points via the eigen-solver (default true).
%   'RefitSphere'       Refit sphere after refinement (default true).
%   'Verbose'           0/1/2 (default 1).
%
% Outputs
%   A                 The matrix of the selected case (same as input A).
%   lamAll            Refined candidate list.
%   lamSamples        Refined DISTINCT sample list (may be []).
%   resAll            Residual certificates for lamAll (via LEIGQNEWTON_CERT_RESMIN).
%   resSamples        Residual certificates for lamSamples.
%   sph               Updated/refitted sphere model (empty if none/invalid/missing).
%   conf              Confidence struct, including the decision/verdict fields.
%   out               Detailed diagnostics (grid statistics, inliers, etc.).
%
% Examples
%   % 1) Detect -> refine/validate (recommended workflow)
%   [lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_detect(A, ...
%       'Collect',8,'RunsMax',12,'Restarts',40,'Seed0',1);
%   C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'cls',cls,'sph',sph,'info',info);
%   [A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_validate(C, ...
%       'GridTheta',4,'GridPhi',4,'GridMaxPoints',32,'Verbose',1);
%
%   % 2) Sample -> refine/validate (sphere model provided by sphere_sample)
%   [lamAll, lamS, sph, infoS] = leigqNEWTON_sphere_sample(A, ...
%       'Collect',8,'RunsMax',12,'Restarts',40,'Seed0',1);
%   C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'sph',sph);
%   [A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_refine(C, ...
%       'GridTheta',4,'GridPhi',4,'GridMaxPoints',32,'Verbose',1);
%
%   % 3) Refinement-only (positional form; sphere verdict requires C.sph, so sph will be empty)
%   [A2, lamAll2, resAll, lamS2, resS, sph2, conf] = leigqNEWTON_sphere_refine(A, lamAll, lamSamples);
%
% See also leigqNEWTON_sphere_validate, leigqNEWTON_sphere_detect,
%          leigqNEWTON_sphere_sample, leigqNEWTON_refine_auto,
%          leigqNEWTON_cert_resMin


if nargin >= 2 && (ischar(k) || (isstring(k) && isscalar(k)))
    varargin = [{k}, varargin];
    k = [];
end

% If user passed (A, lamAll, lamSamples?, Name,Value,...):
if nargin >= 2 && ~(isstruct(cases) || iscell(cases))
    A_in = cases;
    lamAll_in = k;
    lamSamples_in = [];
    nv = varargin;
    if ~isempty(nv) && ~(ischar(nv{1}) || (isstring(nv{1}) && isscalar(nv{1})))
        lamSamples_in = nv{1};
        nv = nv(2:end);
    end
    cases = struct('A', A_in, 'lamAll', lamAll_in, 'lamSamples', lamSamples_in);
    k = [];
    varargin = nv;
end

p = inputParser;
p.addParameter('SphereIndex', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('RefineArgs', {'Verbose',0}, @(c)iscell(c));

p.addParameter('TargetRes', 1e-14, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SphereTolRel', 5e-3, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('VerifySphere', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('GridTheta', 8, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('GridPhi', 16, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
p.addParameter('GridMaxPoints', 400, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('RefineGrid', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('RefitSphere', true, @(x)islogical(x)&&isscalar(x));

p.addParameter('CollapseTol', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinGridForDecision', 64, @(x)isnumeric(x)&&isscalar(x)&&x>=8);

p.addParameter('Verbose', 1, @(x)isnumeric(x)&&isscalar(x));

p.parse(varargin{:});
opt = p.Results;

% ---------------- prerequisites ----------------
if exist('leigqNEWTON_refine_auto','file') ~= 2
    error('leigqNEWTON_sphere_refine:Missing','leigqNEWTON_refine_auto.m not found on the MATLAB path.');
end
if exist('leigqNEWTON_cert_resMin','file') ~= 2
    error('leigqNEWTON_sphere_refine:Missing','leigqNEWTON_cert_resMin.m not found on the MATLAB path.');
end

% ---------------- extract case ----------------
[C, meta] = local_pick_case(cases, k);
A = local_get_any(C, {'A'}, []);
lamAll0 = local_get_any(C, {'lamAll','lambdaAll','lambda'}, []);
lamSamples0 = local_get_any(C, {'lamSamples','lambdaSamples'}, []);
sph0 = local_select_sphere(C, opt.SphereIndex);

% base residuals (no refinement)
resAll0 = leigqNEWTON_cert_resMin(A, lamAll0);
resSamples0 = leigqNEWTON_cert_resMin(A, lamSamples0);

% augment sphere with stable sampler
sph = sph0;
if ~isempty(sph)
    if isfield(sph,'sampler'), sph.samplerOriginal = sph.sampler; else, sph.samplerOriginal = []; end
    sph.samplerStable = @(theta,phi)local_sphere_sampler(sph, theta, phi);
end

% init outputs/out struct
out = struct();
out.meta = meta;
out.lamAll0 = lamAll0;
out.resAll0 = resAll0;
out.lamSamples0 = lamSamples0;
out.resSamples0 = resSamples0;
out.sph0 = sph;

if opt.Verbose
    fprintf('=== leigqNEWTON_sphere_refine: idx=%g, seed=%g, origin=%s ===\n', meta.idx, meta.seedUsed, string(meta.origin));
    fprintf('Ktot=%d, Kdistinct=%d, spheres=%d\n', numel(lamAll0), numel(lamSamples0), ~isempty(sph));
    fprintf('Initial residuals: lamAll  med=%.2e  max=%.2e | lamSamples med=%.2e  max=%.2e\n', ...
        median(resAll0(:)), max(resAll0(:)), median(resSamples0(:)), max(resSamples0(:)));
end

% ---------------- refine all candidates ----------------
[~, lamAll_pol, ~, resAll_pol, ~, infoAll] = leigqNEWTON_refine_auto(A, lamAll0, opt.RefineArgs{:});
[~, lamSamples_pol, ~, resSamples_pol, ~, infoS] = leigqNEWTON_refine_auto(A, lamSamples0, opt.RefineArgs{:});

lamAll = reshape(lamAll_pol, size(lamAll0));
resAll = reshape(resAll_pol, size(lamAll0));
lamSamples = reshape(lamSamples_pol, size(lamSamples0));
resSamples = reshape(resSamples_pol, size(lamSamples0));

out.refineInfoAll = infoAll;
out.refineInfoSamples = infoS;

if opt.Verbose
    fprintf('Refined residuals: lamAll  med=%.2e  max=%.2e | lamSamples med=%.2e  max=%.2e\n', ...
        median(resAll(:)), max(resAll(:)), median(resSamples(:)), max(resSamples(:)));
end

% ---------------- eigenvalue certification summary ----------------
eigen = struct();
eigen.TargetRes = opt.TargetRes;
eigen.passAllFrac = mean(resAll(:) <= opt.TargetRes);
eigen.passSamplesFrac = mean(resSamples(:) <= opt.TargetRes);
eigen.medAll = median(resAll(:));
eigen.maxAll = max(resAll(:));
eigen.medSamples = median(resSamples(:));
eigen.maxSamples = max(resSamples(:));

% ---------------- inlier analysis + optional refit ----------------
out.inlier = struct();
out.sphFit = struct([]);
sphUse = sph; % may become refit
inlier = struct('N',0);

if ~isempty(sph) && isfield(sph,'inliers') && ~isempty(sph.inliers) && ~isempty(lamSamples0)
    inl = unique(round(sph.inliers(:)));
    inl = inl(inl>=1 & inl<=numel(lamSamples0));
    out.inlier.indices = inl;

    if ~isempty(inl)
        lamIn0 = lamSamples0(inl);
        lamIn1 = lamSamples(inl);

        inlier.N = numel(inl);
        inlier.res0 = resSamples0(inl);
        inlier.res1 = resSamples(inl);

        % devRel w.r.t original sphere model
        inlier.devRel0_sph0 = local_sphere_dev_rel(lamIn0, sph);
        inlier.devRel1_sph0 = local_sphere_dev_rel(lamIn1, sph);

        % optional refit to refined inliers
        if opt.RefitSphere
            [sphFit, fitStats] = local_refit_sphere_from_lam(lamIn1, sph);
            out.sphFit = sphFit;
            out.fitStats = fitStats;

            sphUse = sphFit;
            inlier.devRel1_sphFit = local_sphere_dev_rel(lamIn1, sphFit);
            inlier.devRel0_sphFit = local_sphere_dev_rel(lamIn0, sphFit);
        end

        if opt.Verbose
            fprintf('Inliers (n=%d): res med=%.2e max=%.2e\n', inlier.N, median(inlier.res1), max(inlier.res1));
            if opt.RefitSphere && ~isempty(out.sphFit)
                fprintf('Inliers on-sphere (refit): devRel med=%.2e max=%.2e\n', ...
                    median(inlier.devRel1_sphFit(:)), max(inlier.devRel1_sphFit(:)));
            else
                fprintf('Inliers on-sphere (orig):  devRel med=%.2e max=%.2e\n', ...
                    median(inlier.devRel1_sph0(:)), max(inlier.devRel1_sph0(:)));
            end
        end
    end
end
out.inlierStats = inlier;

% ---------------- grid verification ----------------
grid = struct('N0',0,'N1',0);
out.grid = grid;

if opt.VerifySphere && ~isempty(sphUse)
    [theta, phi, gridLam0] = local_make_grid(sphUse, opt.GridTheta, opt.GridPhi, opt.GridMaxPoints);

    out.grid.theta = theta;
    out.grid.phi = phi;
    out.grid.lam0 = gridLam0;

    % (A) residuals at sphere points without refinement (continuum evidence)
    r0 = leigqNEWTON_cert_resMin(A, gridLam0(:));
    out.grid.res0 = reshape(r0, size(gridLam0));
    out.grid.devRel0 = reshape(local_sphere_dev_rel(gridLam0(:), sphUse), size(gridLam0));
    grid.N0 = numel(r0);

    % (B) refinement from sphere points (collapse diagnostic)
    if opt.RefineGrid
        [~, lamPolG, ~, resPolG, ~, infoG] = leigqNEWTON_refine_auto(A, gridLam0(:), opt.RefineArgs{:});
        gridLam1 = reshape(lamPolG, size(gridLam0));
        r1 = reshape(resPolG, size(gridLam0));

        out.grid.lam1 = gridLam1;
        out.grid.res1 = r1;
        out.grid.devRel1 = reshape(local_sphere_dev_rel(gridLam1(:), sphUse), size(gridLam0));
        out.grid.refineInfo = infoG;
        grid.N1 = numel(r1);

        % collapse clustering (count distinct attractors)
        [K, sizes] = local_cluster_quat(gridLam1(:), opt.CollapseTol);
        out.collapse = struct('K',K,'sizes',sizes,'tol',opt.CollapseTol);
    else
        out.grid.lam1 = [];
        out.grid.res1 = [];
        out.grid.devRel1 = [];
        out.collapse = struct('K',NaN,'sizes',[],'tol',opt.CollapseTol);
    end
end
out.gridN0 = grid.N0;
out.gridN1 = grid.N1;

% ---------------- metrics + sphere decision ----------------
metrics = local_metrics(opt, out, eigen);
sphere = local_sphere_decision(opt, out, metrics);

% ---------------- pack conf ----------------
conf = struct();
conf.eigen = eigen;
conf.sphere = sphere;
conf.metrics = metrics;

% ---------------- reporting ----------------
if opt.Verbose
    fprintf('Eigenvalue certification (TargetRes=%.1e): lamAll pass=%.0f%% | lamSamples pass=%.0f%%\n', ...
        opt.TargetRes, 100*eigen.passAllFrac, 100*eigen.passSamplesFrac);

    fprintf('Sphere decision: %s  reliability=%s (%.2f)\n', sphere.verdict, sphere.reliabilityTag, sphere.reliabilityScore);
    fprintf('  Why: %s\n', sphere.why);

    if metrics.gridN0 > 0
        fprintf('Grid (no refine): N=%d  passRes=%.0f%%  passSphere=%.0f%%  passBoth=%.0f%%\n', ...
            metrics.gridN0, 100*metrics.gridPassRes0, 100*metrics.gridPassSphere0, 100*metrics.gridPassBoth0);
    end
    if metrics.gridN1 > 0
        fprintf('Grid (refined):  N=%d  passRes=%.0f%%  passSphere=%.0f%%  passBoth=%.0f%%\n', ...
            metrics.gridN1, 100*metrics.gridPassRes1, 100*metrics.gridPassSphere1, 100*metrics.gridPassBoth1);
        if isfield(out,'collapse') && isfinite(out.collapse.K)
            fprintf('Collapse clustering: K=%d attractor(s) (tol=%.1e). Largest clusters: %s\n', ...
                out.collapse.K, out.collapse.tol, local_top_sizes(out.collapse.sizes));
        end
    end
end

end

% =====================================================================
%                               helpers
% =====================================================================

function [C, meta] = local_pick_case(cases, k)
if nargin < 2 || isempty(k)
    if iscell(cases) && numel(cases)==1 && isstruct(cases{1})
        C = cases{1};
    elseif isstruct(cases) && isscalar(cases)
        C = cases;
    else
        error('leigqNEWTON_sphere_refine:Input','If cases is a list, you must provide k.');
    end
else
    if iscell(cases)
        C = cases{k};
    elseif isstruct(cases) && ~isscalar(cases)
        C = cases(k);
    elseif isstruct(cases) && isscalar(cases)
        C = cases;
    else
        error('leigqNEWTON_sphere_refine:Input','cases must be cell/struct array or a case struct.');
    end
end

meta = struct();
meta.idx      = local_get_any(C, {'idx','Idx','index'}, NaN);
meta.seedUsed = local_get_any(C, {'seedUsed','seed','SeedUsed'}, NaN);
meta.type     = local_get_any(C, {'type','Type'}, '');
meta.density  = local_get_any(C, {'density','Density'}, NaN);
meta.triPart  = local_get_any(C, {'triPart','TriPart'}, '');
meta.origin   = local_get_any(C, {'origin'}, '');
meta.userNote = local_get_any(C, {'userNote'}, '');
meta.info     = local_get_any(C, {'info'}, struct());
end

function sph = local_select_sphere(C, sidx)
sph = struct([]);
try
    sphAll = local_get_any(C, {'sph'}, struct([]));
    if isempty(sphAll), return; end
    s = min(max(1, round(sidx)), numel(sphAll));
    sph = sphAll(s);
catch
    sph = struct([]);
end
end

function v = local_get_any(S, names, default)
v = default;
try
    if ~isstruct(S) || ~isscalar(S), return; end
    fns = fieldnames(S);
    if ischar(names) || isstring(names), names = cellstr(names); end
    for i = 1:numel(names)
        nm = char(names{i});
        if isfield(S, nm), v = S.(nm); return; end
        hit = find(strcmpi(fns, nm), 1, 'first');
        if ~isempty(hit), v = S.(fns{hit}); return; end
    end
catch
    v = default;
end
end

function L = local_sphere_sampler(model, theta, phi)
if ~isfield(model,'basis4x3') || ~isfield(model,'radius')
    error('leigqNEWTON_sphere_refine:sampler','Sphere model missing basis4x3/radius.');
end

theta = double(theta);
phi   = double(phi);
[theta, phi] = local_compat(theta, phi);

u1 = cos(phi).*sin(theta);
u2 = sin(phi).*sin(theta);
u3 = cos(theta);
U = [u1(:)'; u2(:)'; u3(:)'];         % 3 x N

B  = double(model.basis4x3);
r  = double(model.radius);

useAff = isfield(model,'p0') && isfield(model,'center3') && ~isempty(model.p0) && ~isempty(model.center3);
if useAff
    p0 = double(model.p0(:));
    c3 = double(model.center3(:));
    P4 = p0 + B*(c3 + r*U);
else
    c4 = double(model.center4(:));
    P4 = c4 + B*(r*U);
end

w = reshape(P4(1,:), size(theta));
x = reshape(P4(2,:), size(theta));
y = reshape(P4(3,:), size(theta));
z = reshape(P4(4,:), size(theta));
L = quaternion(w,x,y,z);
end

function [a,b] = local_compat(a,b)
if isscalar(a) && ~isscalar(b)
    a = a + zeros(size(b));
elseif ~isscalar(a) && isscalar(b)
    b = b + zeros(size(a));
elseif ~isequal(size(a), size(b))
    error('leigqNEWTON_sphere_refine_v4:sampler','theta and phi must have same size or be scalar-expandable.');
end
end

function devRel = local_sphere_dev_rel(lam, sph)
if isempty(lam) || isempty(sph) || ~isfield(sph,'basis4x3') || ~isfield(sph,'p0') || ~isfield(sph,'center3') || ~isfield(sph,'radius')
    devRel = NaN(size(lam));
    return;
end

[w,x,y,z] = parts(lam);
P4 = [w(:)'; x(:)'; y(:)'; z(:)'];

p0 = double(sph.p0(:));
B  = double(sph.basis4x3);
c3 = double(sph.center3(:));
r  = double(sph.radius);
rEff = max(r, eps);

X3 = (B' * B) \ (B' * (P4 - p0));  % LS coords
D  = X3 - c3;
dist = sqrt(sum(D.^2,1));

dev = abs(dist - r);
devRel = reshape(dev./rEff, size(w));
end

function [sphFit, stats] = local_refit_sphere_from_lam(lam, sphBase)
stats = struct();
sphFit = sphBase;

if isempty(lam) || ~isfield(sphBase,'p0') || ~isfield(sphBase,'basis4x3')
    sphFit = struct([]); return;
end

[w,x,y,z] = parts(lam);
P4 = [w(:)'; x(:)'; y(:)'; z(:)'];

p0 = double(sphBase.p0(:));
B  = double(sphBase.basis4x3);

X3 = (B' * B) \ (B' * (P4 - p0));  % LS coords

X = X3';
N = size(X,1);

b = -sum(X.^2,2);
M = [X, ones(N,1)];
a = M \ b;

c = -0.5 * a(1:3);
rho2 = sum(c.^2) - a(4);
rho2 = max(rho2, 0);
r = sqrt(rho2);

sphFit.center3 = c(:);
sphFit.radius  = r;

center4 = p0 + B * sphFit.center3;
sphFit.center4 = center4(:).';
sphFit.center  = quaternion(center4(1), center4(2), center4(3), center4(4));

dist = sqrt(sum((X - c(:).').^2,2));
dev = abs(dist - r);
stats.N = N;
stats.devAbsMed = median(dev);
stats.devAbsMax = max(dev);
stats.devRelMed = stats.devAbsMed / max(r, eps);
stats.devRelMax = stats.devAbsMax / max(r, eps);
end

function [theta, phi, L] = local_make_grid(sph, nTheta, nPhi, maxPts)
theta = linspace(0, pi, round(nTheta));
phi   = linspace(0, 2*pi, round(nPhi));
[T,P] = ndgrid(theta, phi);
L = local_sphere_sampler(sph, T, P);

N = numel(L);
if N > maxPts
    stride = ceil(N/maxPts);
    L = L(1:stride:end);
    theta = [];
    phi = [];
end
end

function [K, sizes] = local_cluster_quat(lam, tol)
% Greedy clustering in R^4 by component distance <= tol (absolute).
if isempty(lam)
    K = 0; sizes = [];
    return;
end
[w,x,y,z] = parts(lam(:));
V = [double(w), double(x), double(y), double(z)];

N = size(V,1);
lab = zeros(N,1);
K = 0;
for i = 1:N
    if lab(i) == 0
        K = K + 1;
        d = sqrt(sum((V - V(i,:)).^2, 2));
        lab(d <= tol) = K;
    end
end
sizes = accumarray(lab, 1);
sizes = sort(sizes,'descend');
end

function s = local_top_sizes(sizes)
if isempty(sizes)
    s = '[]'; return;
end
m = min(numel(sizes), 5);
s = mat2str(sizes(1:m).');
end

function metrics = local_metrics(opt, out, eigen)
metrics = struct();

% eigen
metrics.eigenPassAll = eigen.passAllFrac;
metrics.eigenPassSamples = eigen.passSamplesFrac;

% inliers (refit if available)
metrics.inlN = 0; metrics.inlPassRes = NaN; metrics.inlPassSphere = NaN; metrics.inlPassBoth = NaN;
if isfield(out,'inlierStats') && isfield(out.inlierStats,'N') && out.inlierStats.N > 0
    metrics.inlN = out.inlierStats.N;
    passRes = out.inlierStats.res1(:) <= opt.TargetRes;

    if isfield(out.inlierStats,'devRel1_sphFit')
        dev = out.inlierStats.devRel1_sphFit(:);
    else
        dev = out.inlierStats.devRel1_sph0(:);
    end
    passSphere = dev <= opt.SphereTolRel;

    metrics.inlPassRes = mean(passRes);
    metrics.inlPassSphere = mean(passSphere);
    metrics.inlPassBoth = mean(passRes & passSphere);
end

% grid
metrics.gridN0 = 0; metrics.gridPassRes0 = NaN; metrics.gridPassSphere0 = NaN; metrics.gridPassBoth0 = NaN;
metrics.gridN1 = 0; metrics.gridPassRes1 = NaN; metrics.gridPassSphere1 = NaN; metrics.gridPassBoth1 = NaN;

if isfield(out,'grid') && isfield(out.grid,'res0') && ~isempty(out.grid.res0)
    r0 = out.grid.res0(:);
    d0 = out.grid.devRel0(:);
    metrics.gridN0 = numel(r0);
    metrics.gridPassRes0 = mean(r0 <= opt.TargetRes);
    metrics.gridPassSphere0 = mean(d0 <= opt.SphereTolRel);
    metrics.gridPassBoth0 = mean((r0 <= opt.TargetRes) & (d0 <= opt.SphereTolRel));
end

if isfield(out,'grid') && isfield(out.grid,'res1') && ~isempty(out.grid.res1)
    r1 = out.grid.res1(:);
    d1 = out.grid.devRel1(:);
    metrics.gridN1 = numel(r1);
    metrics.gridPassRes1 = mean(r1 <= opt.TargetRes);
    metrics.gridPassSphere1 = mean(d1 <= opt.SphereTolRel);
    metrics.gridPassBoth1 = mean((r1 <= opt.TargetRes) & (d1 <= opt.SphereTolRel));
end

% collapse
if isfield(out,'collapse')
    metrics.collapseK = out.collapse.K;
else
    metrics.collapseK = NaN;
end
end

function sphere = local_sphere_decision(opt, out, m)
sphere = struct();
sphere.verdict = 'INCONCLUSIVE';
sphere.why = '';
sphere.reliabilityTag = 'LOW';
sphere.reliabilityScore = 0.2;

if isempty(out.sph0)
    sphere.verdict = 'NO_SPHERE_MODEL';
    sphere.why = 'No sphere model stored in the case.';
    sphere.reliabilityTag = 'HIGH';
    sphere.reliabilityScore = 1.0;
    return;
end

% guard: enough grid?
hasGrid0 = (m.gridN0 >= opt.MinGridForDecision);
hasGrid1 = (m.gridN1 >= opt.MinGridForDecision);

% Collapse signature:
collapseStrong = hasGrid0 && hasGrid1 && ...
    (m.gridPassSphere0 > 0.95) && (m.gridPassRes0 < 0.05) && ...
    (m.gridPassRes1 > 0.95) && (m.gridPassSphere1 < 0.05);

% True sphere signature (primary: no-refine grid evidence):
trueStrong = hasGrid0 && (m.gridPassBoth0 > 0.50);

% Secondary (weaker): refined grid stays on sphere too
trueWeak = hasGrid1 && (m.gridPassBoth1 > 0.50) && (m.gridPassSphere1 > 0.50);

if trueStrong || trueWeak
    sphere.verdict = 'SPHERE_CONFIRMED';
    if trueStrong
        sphere.why = 'Many grid points on the fitted sphere already have small residuals without refinement (continuum evidence).';
        sphere.reliabilityScore = min(1, (m.gridPassBoth0 - 0.50) / 0.50);
    else
        sphere.why = 'After refinement, many points remain on one sphere and have small residuals (weaker than no-refine evidence).';
        sphere.reliabilityScore = min(0.8, (m.gridPassBoth1 - 0.50) / 0.50);
    end
    sphere.reliabilityTag = local_rel_tag(sphere.reliabilityScore);
    return;
end

if collapseStrong
    sphere.verdict = 'ARTIFACT_COLLAPSE';
    sphere.why = sprintf(['Grid points ON the fitted sphere do not certify as eigenvalues (no-refine passRes≈0%%), ' ...
                          'but refining those starts converges to certified eigenvalues that leave the sphere (passSphere≈0%%). ' ...
                          'This is a strong collapse/attractor artifact.']);
    % reliability from how extreme the signature is
    r1 = min(1, (0.05 - m.gridPassRes0)/0.05);
    r2 = min(1, (m.gridPassRes1 - 0.95)/0.05);
    r3 = min(1, (0.05 - m.gridPassSphere1)/0.05);
    sphere.reliabilityScore = max(0, min([r1,r2,r3]));
    sphere.reliabilityTag = local_rel_tag(sphere.reliabilityScore);
    return;
end

% Otherwise inconclusive; attempt to explain
if hasGrid0
    if m.gridPassRes0 < 0.05
        sphere.why = 'Sphere points fail residual test without refinement; this does not support a continuum at TargetRes.';
    elseif m.gridPassBoth0 < 0.20
        sphere.why = 'Only a small fraction of sphere points certify without refinement; evidence is weak at TargetRes/SphereTolRel.';
    else
        sphere.why = 'Mixed evidence; adjust tolerances or increase grid for a clearer decision.';
    end
else
    sphere.why = 'Insufficient grid size for a reliable decision; increase GridTheta/GridPhi or GridMaxPoints.';
end
sphere.reliabilityTag = 'LOW';
sphere.reliabilityScore = 0.2;
end

function tag = local_rel_tag(score)
if score >= 0.85
    tag = 'HIGH';
elseif score >= 0.55
    tag = 'MED';
else
    tag = 'LOW';
end
end
