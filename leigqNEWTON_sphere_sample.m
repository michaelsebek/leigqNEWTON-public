function [lambdaAll, lambdaSamples, lambda0, class, sph, info] = leigqNEWTON_sphere_sample(A, varargin)
%LEIGQNEWTON_SPHERE_SAMPLE  Sample left eigenvalues and detect candidate eigen-spheres (advanced engine).
%
%   [lambdaAll, lambdaSamples, lambda0, cls, sph, info] = ...
%       leigqNEWTON_sphere_sample(A, Name,Value,...)
%
%   This routine repeatedly calls LEIGQNEWTON from many initial guesses (controlled
%   by Seeds and RunsMax) to collect left-eigenvalue candidates. It then keeps a
%   DISTINCT subset for clustering and, when enabled, fits sphere models in R^4 to
%   detect potential 2-sphere families of left eigenvalues.
%
%   Most users should call LEIGQNEWTON_SPHERE_DETECT instead (same computation,
%   more user-facing defaults/docs). To strengthen or reject a detected sphere,
%   use LEIGQNEWTON_SPHERE_VALIDATE (which calls LEIGQNEWTON_SPHERE_REFINE).
%
% Inputs
%   A : n-by-n quaternion matrix.
%
% Name,Value options (selected)
%   'Collect'        Target number of DISTINCT samples to collect (default 5*n).
%   'RunsMax'        Maximum number of solver runs used to reach 'Collect' (default 50).
%   'Seed0'          Base seed; run r uses seed Seed0+(r-1) (default 1).
%   'Restarts'       Restart budget passed to LEIGQNEWTON (default 500).
%   'UniqueTol'      Tolerance used to thin samples to a DISTINCT subset (default 1e-8).
%   'ResMax'         Discard candidates with resMin(A,lambda) > ResMax (default inf).
%   'ZeroTol'        Threshold used to classify "certain" zeros (default 0).
%
%   'DetectSphere'   Enable sphere detection on the DISTINCT subset (default true).
%   'MinSphereCount' Minimum points required for a sphere attempt (default 10).
%   'RansacTrials'   Number of random trials in sphere fitting (default 200).
%   'TolSubspace'    Tolerance for the fitted affine 3-subspace (default 1e-6).
%   'TolRadius'      Tolerance for the fitted radius consistency (default 1e-6).
%   'MaxSpheres'     Maximum number of spheres to return (default 3).
%
%   'Report'         'off'|'summary'|'progress'|'full' (default 'off').
%                   NOTE: this controls INTERNAL sampler printing. In batch hunts
%                   use runExNEWTON / runExNEWTON_sparse which keep this silent by default.
%
% Outputs
%   lambdaAll        M-by-1 quaternion, counted candidate list (may include verified zeros).
%   lambdaSamples    K-by-1 quaternion, DISTINCT subset used for clustering/sphere detection.
%   lambda0          Struct with representative candidates and zero-related diagnostics.
%   cls              K-by-1 int32 cluster labels for lambdaSamples.
%   sph              Struct array with detected sphere model(s) (empty if none detected).
%   info             Struct with counts, timings, thresholds, and diagnostics.
%
% Examples
%   % Detect candidates (silent by default):
%   [lamAll, lamS, ~, cls, sph, info] = leigqNEWTON_sphere_sample(A, 'Collect', 30);
%
%   % Show internal sampler progress:
%   [lamAll, lamS] = leigqNEWTON_sphere_sample(A, 'Collect', 30, 'Report', 'progress');
%
% See also leigqNEWTON_sphere_detect, leigqNEWTON_sphere_validate,
%          leigqNEWTON_sphere_refine, leigqNEWTON
%
% Author: Michael Sebek (michael.sebek@fel.cvut.cz)
%   Version: v1.0
%
% Note:
%  This function is part of the public MATLAB toolbox leigqNEWTON accompanying the paper:
%     M. Sebek, "Computing Left Eigenvalues of Quaternion Matrices", submitted to
%     Linear Algebra and its Applications, 2026.
%  If you use this software in academic work, please cite the paper and/or the Zenodo archive:
%     https://doi.org/10.5281/zenodo.18410141
%

p = inputParser;

% collection
p.addParameter('Collect', [], @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('RunsMax', 50, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Seed0', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Restarts', 500, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('UniqueTol', 1e-8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('ResMax', inf, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('ZeroTol', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% sphere detection
p.addParameter('DetectSphere', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MinSphereCount', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=4);
p.addParameter('RansacTrials', 200, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('TolSubspace', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TolRadius', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MaxSpheres', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% reporting
p.addParameter('Report', 'off', @(x)ischar(x)||isstring(x)||islogical(x));

p.parse(varargin{:});
opt = p.Results;

reportMode = local_parse_report(opt.Report);

Aq = local_to_quat(A);
[n,m] = size(Aq);
if n ~= m
    error('leigqNEWTON_sphere_sample:Square', 'A must be square.');
end
if isempty(opt.Collect)
    opt.Collect = 5*n;
end

% --- sparse-friendly default tweak (heuristic), as in leigqNEWTON_sphere_sample ---
if ismember('UniqueTol', p.UsingDefaults)
    try
        [aw,ax,ay,az] = parts(Aq);
        nz = (aw~=0) | (ax~=0) | (ay~=0) | (az~=0);
        density = nnz(nz) / max(1, numel(nz));
        if density <= 0.20
            Af = sqrt( norm(aw,'fro')^2 + norm(ax,'fro')^2 + norm(ay,'fro')^2 + norm(az,'fro')^2 );
            scale = max(1, Af);
            opt.UniqueTol = max(opt.UniqueTol, 1e-4 * scale);
        end
    catch
        % keep default
    end
end

% ---------------- collect samples by repeated solver calls ----------------
% lambdaSamples: DISTINCT samples
lambdaSamples = local_qzeros(0,1);
resSamples    = zeros(0,1);
Psamples      = zeros(0,4);

% lambdaAll: counted list = (certain zeros with multiplicity once) + distinct nonzero samples
lambdaAll = local_qzeros(0,1);

lambda0 = struct();
lambda0.certain     = local_qzeros(0,1);
lambda0.certainMult = 0;
lambda0.newtonRaw   = local_qzeros(0,1);
lambda0.zeroTol     = opt.ZeroTol;

runCount = 0;
zeroMult_ref = [];
zeroMult_inconsistent = false;

if reportMode.progress
    fprintf('=== leigqNEWTON_sphere_sample: collecting samples (n=%d, targetDistinct=%d) ===\n', n, opt.Collect);
end

for r = 1:opt.RunsMax
    if numel(lambdaSamples) >= opt.Collect
        break;
    end
    runCount = runCount + 1;
    seed = opt.Seed0 + r - 1;

    % Call leigqNEWTON, requesting summary info if available.
    infoCell = [];
    try
        [lambda, ~, res, infoCell] = leigqNEWTON(Aq, ...
            'Num', n, ...
            'Restarts', opt.Restarts, ...
            'Seed', seed, ...
            'UniqueTol', opt.UniqueTol, ...
            'InfoLevel', 'summary');
    catch
        try
            [lambda, ~, res] = leigqNEWTON(Aq, ...
                'Num', n, ...
                'Restarts', opt.Restarts, ...
                'Seed', seed, ...
                'UniqueTol', opt.UniqueTol);
        catch
            [lambda, ~, res] = leigqNEWTON(Aq);
        end
    end
    lamRun = lambda(:);
    if isempty(res)
        resRun = zeros(size(lamRun));
    else
        resRun = res(:);
        if numel(resRun) ~= numel(lamRun)
            resRun = zeros(size(lamRun));
        end
    end

    % Detect how many VERIFIED zero-eigenvalues were accepted in the pre-pass.
    nZeroCertain = 0;
    try
        if iscell(infoCell) && ~isempty(infoCell) && isstruct(infoCell{1}) && isfield(infoCell{1},'nAcceptedZero')
            nZeroCertain = double(infoCell{1}.nAcceptedZero);
        end
    catch
        nZeroCertain = 0;
    end

    if isempty(zeroMult_ref)
        zeroMult_ref = nZeroCertain;
        lambda0.certainMult = nZeroCertain;
        if nZeroCertain > 0 && numel(lamRun) >= nZeroCertain
            lambda0.certain = lamRun(1:nZeroCertain);
            % Count these zeros ONCE in lambdaAll, WITH multiplicity.
            %
            % IMPORTANT (generator-level logic): verified zeros are NOT included
            % in lambdaSamples, because lambdaSamples is reserved for DISTINCT
            % NONZERO samples used for clustering / sphere detection.
            lambdaAll = [lambdaAll; lambda0.certain];
        end
    else
        if nZeroCertain ~= zeroMult_ref
            % Should not happen for a fixed A; track if it does.
            zeroMult_inconsistent = true;
        end
    end

    % Remove the certain-zero block from this run (we already accounted for it once).
    if nZeroCertain > 0 && numel(lamRun) >= nZeroCertain
        lamRun = lamRun(nZeroCertain+1:end);
        resRun = resRun(nZeroCertain+1:end);
    end

    % Residual filter
    keep = (resRun <= opt.ResMax);
    lamRun = lamRun(keep);
    resRun = resRun(keep);

    % Optional near-zero handling: record and exclude from DISTINCT sampling.
    if opt.ZeroTol > 0 && ~isempty(lamRun)
        qn = local_qnorm(lamRun);
        is0 = (qn <= opt.ZeroTol);
        if any(is0)
            lambda0.newtonRaw = [lambda0.newtonRaw; lamRun(is0)];
            lamRun = lamRun(~is0);
            resRun = resRun(~is0);
        end
    end

    % Append DISTINCT nonzero samples; also append newly added ones to lambdaAll.
    [lambdaSamples, resSamples, Psamples, added, lamAdded, ~] = local_append_unique(lambdaSamples, resSamples, Psamples, lamRun, resRun, opt.UniqueTol);
    if added > 0
        lambdaAll = [lambdaAll; lamAdded];
    end

    if reportMode.progress
        fprintf('  run %3d | Seed=%d | addedDistinct %2d | totalDistinct %3d | counted %3d\n', ...
            r, seed, added, numel(lambdaSamples), numel(lambdaAll));
    end
end

Kd = numel(lambdaSamples);
class = zeros(Kd,1,'int32');
sph = struct([]); %#ok<STRNU>

info = struct();
info.n = n;
info.runsUsed = runCount;
info.samplesDistinct = Kd;
info.samplesCounted  = numel(lambdaAll);
info.uniqueTol = opt.UniqueTol;
info.zeroTol = opt.ZeroTol;
info.zeroMultiplicity = lambda0.certainMult;
info.zeroNewtonHits  = numel(lambda0.newtonRaw);
info.zeroMultiplicityInconsistent = zeroMult_inconsistent;
info.resSamples = resSamples;

% Early exit if detection disabled or insufficient DISTINCT samples
if ~opt.DetectSphere || Kd < opt.MinSphereCount
    info.spheresFound = 0;
    info.remainingCount = Kd;
    info.note = 'Sphere detection disabled or not enough DISTINCT samples.';
    info.report = local_make_report(reportMode, info, sph);
    if reportMode.summary
        fprintf('%s', info.report);
    end
    return;
end

% ---------------- detect spheres by RANSAC fitting in R^4 ----------------
P = Psamples;  % Kd x 4
remaining = true(Kd,1);

sphereCount = 0;
while sphereCount < opt.MaxSpheres
    idxRem = find(remaining);
    if numel(idxRem) < opt.MinSphereCount
        break;
    end

    [model, inliersRem] = local_ransac_sphere(P(idxRem,:), opt);
    if isempty(model)
        break;
    end

    inliers = idxRem(inliersRem);
    if numel(inliers) < opt.MinSphereCount
        break;
    end

    sphereCount = sphereCount + 1;

    sph(sphereCount).center4  = model.center4;
    sph(sphereCount).center   = quaternion(model.center4(1), model.center4(2), model.center4(3), model.center4(4));
    sph(sphereCount).radius   = model.radius;
    sph(sphereCount).basis4x3 = model.basis4x3;
    sph(sphereCount).p0       = model.p0;
    sph(sphereCount).center3  = model.center3;
    sph(sphereCount).inliers  = inliers;
    sph(sphereCount).samples4 = lambdaSamples(inliers(1:min(4,end)));
    sph(sphereCount).sampler  = @(theta,phi) local_sphere_sampler(model, theta, phi);

    class(inliers) = int32(sphereCount);
    remaining(inliers) = false;
end

info.spheresFound = sphereCount;
info.remainingCount = nnz(remaining);
info.note = 'Sphere detection completed (DISTINCT samples).';
info.report = local_make_report(reportMode, info, sph);

if reportMode.summary
    fprintf('%s', info.report);
end

end

% =====================================================================
%                          helper functions
% =====================================================================

function reportMode = local_parse_report(r)
reportMode = struct('off',false,'summary',false,'progress',false,'full',false);
if islogical(r)
    if r
        r = 'summary';
    else
        r = 'off';
    end
end
r = lower(string(r));
switch r
    case "off"
        reportMode.off = true;
    case "summary"
        reportMode.summary = true;
    case "progress"
        reportMode.progress = true;
    case "full"
        reportMode.full = true;
        reportMode.progress = true;
        reportMode.summary = true;
    otherwise
        reportMode.off = true;
end
end

function Aq = local_to_quat(A)
if isa(A,'quaternion')
    Aq = A;
else
    Aq = quaternion(A);
end
end

function qz = local_qzeros(n,m)
qz = quaternion(zeros(n,m), zeros(n,m), zeros(n,m), zeros(n,m));
end

function P = local_q2r4(q)
% Return Nx4 real representation [w x y z]
[w,x,y,z] = parts(q);
P = [w(:), x(:), y(:), z(:)];
end

function qn = local_qnorm(q)
P = local_q2r4(q);
qn = sqrt(sum(P.^2,2));
end

function txt = local_make_report(reportMode, info, sph)
if reportMode.off
    txt = '';
    return;
end

txt = sprintf('=== leigqNEWTON_sphere_sample summary ===\n');
txt = [txt, sprintf(' n = %d\n', info.n)];
txt = [txt, sprintf(' runs used = %d\n', info.runsUsed)];
txt = [txt, sprintf(' distinct samples = %d\n', info.samplesDistinct)];
txt = [txt, sprintf(' counted K (incl. zero multiplicity) = %d\n', info.samplesCounted)];
if isfield(info,'zeroMultiplicity')
    txt = [txt, sprintf(' zero multiplicity (verified) = %d\n', info.zeroMultiplicity)];
end
if isfield(info,'zeroNewtonHits') && info.zeroNewtonHits > 0
    txt = [txt, sprintf(' near-zero Newton hits recorded = %d (ZeroTol=%g)\n', info.zeroNewtonHits, info.zeroTol)];
end
if isfield(info,'zeroMultiplicityInconsistent') && info.zeroMultiplicityInconsistent
    txt = [txt, sprintf(' WARNING: zero multiplicity inconsistent across runs.\n')];
end

txt = [txt, sprintf(' spheres found = %d\n', info.spheresFound)];
if info.spheresFound > 0 && ~isempty(sph)
    for s = 1:numel(sph)
        c4 = sph(s).center4;
        txt = [txt, sprintf('  sphere %d: center~(%.4g,%.4g,%.4g,%.4g), r~%.4g, inliers=%d\n', ...
            s, c4(1),c4(2),c4(3),c4(4), sph(s).radius, numel(sph(s).inliers))];
    end
end

txt = [txt, sprintf(' remaining (unclassified) distinct samples = %d\n', info.remainingCount)];
if isfield(info,'note')
    txt = [txt, sprintf(' note: %s\n', info.note)];
end
end
function q = local_sphere_sampler(model, theta, phi)
% model: fitted in R^4 with an affine 3-space basis; sample 2-sphere embedded in that space.
% theta,phi may be scalars or vectors of same size.

theta = theta(:);
phi   = phi(:);
if numel(theta) ~= numel(phi)
    error('local_sphere_sampler:AngleSize','theta and phi must have same number of elements.');
end

% Unit sphere in R^3
u = [cos(phi).*sin(theta), sin(phi).*sin(theta), cos(theta)];  % Mx3

% Map to affine 3-space in R^4
Y = model.center3(:).' + model.radius * u; %#ok<NASGU>
Y = model.center3(:).' + model.radius * u; % (kept for clarity)
Y = model.center3(:).' + model.radius * u; % (kept)

% Compute in R^4
P4 = model.p0 + ( (model.center3(:).' + model.radius*u) - model.center3(:).' ) * model.basis4x3.';

% Convert to quaternion
q = quaternion(P4(:,1), P4(:,2), P4(:,3), P4(:,4));
end

function [lamOut, resOut, Pout, added, lamAdded, resAdded] = local_append_unique(lamOld, resOld, Pold, lamNew, resNew, tol)
% Append elements of lamNew/resNew that are at distance > tol from ALL existing samples
% (in R^4). Returns also the list of newly added samples.

lamOut = lamOld;
resOut = resOld;
Pout   = Pold;
added  = 0;
lamAdded = local_qzeros(0,1);
resAdded = zeros(0,1);

if isempty(lamNew)
    return;
end

Pnew = local_q2r4(lamNew);

for k = 1:size(Pnew,1)
    pk = Pnew(k,:);

    if isempty(Pout)
        isUnique = true;
    else
        d2 = sum((Pout - pk).^2, 2);
        isUnique = all(d2 > tol^2);
    end

    if isUnique
        lamOut = [lamOut; lamNew(k)];
        resOut = [resOut; resNew(k)];
        Pout   = [Pout; pk];
        lamAdded = [lamAdded; lamNew(k)];
        resAdded = [resAdded; resNew(k)];
        added = added + 1;
    end
end
end

function [modelBest, inliersBest] = local_ransac_sphere(P, opt)
% Fit a 2-sphere embedded in an affine 3D subspace of R^4 using RANSAC.
% Returns model with fields: p0, basis4x3, center3, center4, radius.

N = size(P,1);
if N < 4
    modelBest = [];
    inliersBest = [];
    return;
end

bestInliers = 0;
modelBest = [];
inliersBest = [];

for t = 1:opt.RansacTrials
    % Choose 4 random points to define an affine 3-space (if possible)
    idx = randperm(N,4);
    Q = P(idx,:);

    % Affine subspace: p0 + span(U), with U 4x3 orthonormal
    p0 = Q(1,:);
    V = (Q(2:4,:) - p0).';  % 4x3

    % Orthonormal basis for span(V)
    [U,~,~] = svd(V,'econ');
    U = U(:,1:3);

    % Subspace distance for all points
    X = (P - p0)';        % 4xN
    proj = U*(U'*X);      % 4xN
    distSub = sqrt(sum((X - proj).^2,1));

    inSub = distSub <= opt.TolSubspace;
    idxSub = find(inSub);
    if numel(idxSub) < opt.MinSphereCount
        continue;
    end

    % Coordinates in the 3D subspace
    Y = (U'*(P(idxSub,:) - p0)')';  % (#sub)x3

    % Fit 2-sphere in R^3: ||y - c|| = r
    [c3, r0, ok] = local_fit_sphere_3d(Y);
    if ~ok
        continue;
    end

    % Inlier test by radius consistency
    rr = sqrt(sum((Y - c3').^2,2));
    inRad = abs(rr - r0) <= opt.TolRadius;

    inliers = idxSub(inRad);
    nIn = numel(inliers);

    if nIn > bestInliers
        bestInliers = nIn;

        center4 = p0 + (U*c3)';

        modelBest = struct();
        modelBest.p0 = p0;
        modelBest.basis4x3 = U;
        modelBest.center3 = c3;
        modelBest.center4 = center4;
        modelBest.radius  = r0;

        inliersBest = inliers;
    end
end

if bestInliers < opt.MinSphereCount
    modelBest = [];
    inliersBest = [];
end
end

function [c, r, ok] = local_fit_sphere_3d(Y)
% Algebraic sphere fit in R^3. Returns center c (3x1), radius r.
% Solve: ||y||^2 = 2*c'*y + d, then r^2 = ||c||^2 + d.

ok = false;
c = zeros(3,1);
r = NaN;

m = size(Y,1);
if m < 4
    return;
end

A = [2*Y, ones(m,1)];
b = sum(Y.^2,2);

% least squares
x = A\b;
c = x(1:3);
d = x(4);

r2 = sum(c.^2) + d;
if r2 <= 0
    return;
end

r = sqrt(r2);
ok = true;
end
