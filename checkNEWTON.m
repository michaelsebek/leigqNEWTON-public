function [out, cases, S] = checkNEWTON(A, lambda, v, varargin)
%CHECKNEWTON  Quick diagnostic report for quaternion left eigenpairs (leigqNEWTON toolbox).
%
%   out = CHECKNEWTON(A)
%   out = CHECKNEWTON(A, lambda)
%   out = CHECKNEWTON(A, lambda, v)
%   [out, cases, S] = CHECKNEWTON(...)
%   out = CHECKNEWTON(..., Name,Value,...)
%
% CHECKNEWTON prints a compact diagnostic report for quaternionic LEFT
% eigenpairs (lambda, v) of a square quaternion matrix A:
%       A*v = lambda*v,   lambda acts on the LEFT, v ~= 0.
%
% The function depends ONLY on:
%   - MATLAB + its toolboxes (notably the built-in quaternion class),
%   - the leigqNEWTON toolbox functions shipped in this folder.
%
% What it prints
%   1) A summary line (matrix size, how inputs were obtained, median/max
%      absolute and relative residuals).
%   2) Per-eigenvalue residuals with warnings ("recommend refine").
%   3) An "interestingness" assessment:
%        - standard: Kdistinct = n,
%        - nonstandard: Kdistinct > n,
%        - nonstandard: Kdistinct < n (suggest to search longer),
%        - optional sphere report (via leigqNEWTON_sphere_sample).
%   4) A short reminder where to find values in the outputs.
%
% Inputs
%   A      : n-by-n quaternion matrix.
%   lambda : (optional) K-by-1 quaternion array of candidate LEFT eigenvalues.
%            If omitted or empty, CHECKNEWTON calls leigqNEWTON first.
%   v      : (optional) n-by-K quaternion array, eigenvectors as columns.
%            If omitted/empty, CHECKNEWTON computes "best" vectors for the
%            provided lambda by minimizing resMin(A,lambda) using
%            leigqNEWTON_cert_resMin.
%
% Name,Value options (selected)
%   'SolveArgs'    : cell array of Name,Value forwarded to leigqNEWTON when
%                   lambda is not provided.
%                   Default {'SolveProfile','default','Seed',1}.
%   'CleanTol'     : tolerance used by qcleanNEWTON for presentation.
%                   Default 1e-12.
%   'UniqueTol'    : tolerance for counting distinct eigenvalues (4D distance
%                   on cleaned lambdas). Default 1e-6.
%   'WarnAbs'      : absolute warning threshold (applied to BOTH resMin and
%                   resPair). Default 1e-10.
%   'WarnRel'      : relative warning threshold (applied to BOTH resMin and
%                   resPair). Default 1e-10.
%   'PrintDigits'  : N used in qroundNEWTON(...,N,'significant') for printing.
%                   Default 4.
%   'SphereCheck'  : 'auto'|'on'|'off'.
%                   - 'auto' (default): run sphere sampling only when
%                     Kdistinct > n.
%                   - 'on' : always run sphere sampling.
%                   - 'off': never run sphere sampling.
%   'SphereCollect': number of solver batches used by leigqNEWTON_sphere_sample.
%                   Default 5*n.
%   'SphereSeed0'  : base seed for sphere sampling. Default 1.
%   'Verbose'      : 0/1 (default 1). If 0, suppresses most printing.
%
% Outputs
%   out   : struct with all key data (inputs, computed lambdas/vectors,
%           residuals, distinct counts, optional sphere report).
%   cases : cell array compatible with sphere tooling (last outputs are
%           provided in the same spirit as TESTSPHERE*). Typically {C}.
%   S     : summary struct (again in the spirit of TESTSPHERE*).
%
% Useful one-liners
%   out = checkNEWTON(A)
%   out = checkNEWTON(A, lam)
%   out = checkNEWTON(A, lam, V)
%   out = checkNEWTON(A, [], [], 'SolveArgs', {'SolveProfile','reliable','Seed',1})
%   out = checkNEWTON(A, lam, [], 'SphereCheck','on')
%
% See also: leigqNEWTON, leigqNEWTON_refine_batch, leigqNEWTON_cert_resMin,
%           leigqNEWTON_cert_resPair, leigqNEWTON_sphere_sample

% ---------------- normalize positional inputs ----------------
if nargin < 2
    lambda = [];
end
if nargin < 3
    v = [];
end

% Allow CHECKNEWTON(A, Name,Value,...) (lambda omitted)
if nargin >= 2 && (ischar(lambda) || (isstring(lambda) && isscalar(lambda)))
    varargin = [{lambda}, {v}, varargin];
    lambda = [];
    v = [];
end

% ---------------- parse options ----------------
p = inputParser;
p.addParameter('SolveArgs', {'SolveProfile','default','Seed',1}, @(c)iscell(c));
p.addParameter('CleanTol', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('UniqueTol', 1e-6, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('WarnAbs', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('WarnRel',  1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('PrintDigits', 4, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SphereCheck', 'auto', @(s)ischar(s)||isstring(s));
p.addParameter('SphereCollect', [], @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SphereSeed0', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Verbose', 1, @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

opt.SphereCheck = lower(string(opt.SphereCheck));
if isempty(opt.SphereCollect)
    try
        nTmp = size(A,1);
    catch
        nTmp = 1;
    end
    opt.SphereCollect = 5*double(nTmp);
end

if exist('quaternion','class') ~= 8
    error('checkNEWTON:NoQuaternionClass', 'MATLAB class "quaternion" not found. Requires Aerospace Toolbox.');
end

Aq = quaternion(A);
[n,m] = size(Aq);
if n ~= m
    error('checkNEWTON:Square', 'A must be square.');
end

% ---------------- compute missing lambdas/vectors ----------------
meta = struct();
meta.lambdaProvided = ~isempty(lambda);
meta.vProvided      = ~isempty(v);
meta.lambdaComputed = false;
meta.vComputed      = false;

lam = [];
V = [];
lamc = [];
resSolver = [];
infoSolver = [];

if isempty(lambda)
    meta.lambdaComputed = true;
    [lam, V, resSolver, infoSolver] = leigqNEWTON(Aq, opt.SolveArgs{:});
else
    lam = quaternion(lambda);
    lam = lam(:);
end

K = numel(lam);
if K == 0
    error('checkNEWTON:NoLambda', 'No eigenvalue candidates available (lambda list is empty).');
end

% Handle eigenvectors
if isempty(v)
    meta.vComputed = true;
    % Best vectors for each lambda via resMin minimizer
    [resMinAbs, xMin] = leigqNEWTON_cert_resMin(Aq, lam, 'ResidualNormalized', false);
    V = xMin;
else
    V = quaternion(v);
    if size(V,1) ~= n
        error('checkNEWTON:BadV', 'Eigenvector array V must have %d rows.', n);
    end
    % If V is a single vector but multiple lambdas, replicate as a fallback
    if size(V,2) == 1 && K > 1
        V = repmat(V, 1, K);
    end
    if size(V,2) ~= K
        error('checkNEWTON:BadVCols', 'Eigenvector array V must have K columns (K = numel(lambda)).');
    end
    % Still compute resMinAbs for reporting
    resMinAbs = leigqNEWTON_cert_resMin(Aq, lam, 'ResidualNormalized', false);
end

% If lambda came from solver, lamc/resSolver already computed there.
if isempty(lamc)
    lamc = qcleanNEWTON(lam, opt.CleanTol);
end

% ---------------- residuals (abs + relative) ----------------
A_fro = local_qfro(Aq);
% NOTE: MATLAB's built-in quaternion does not implement ABS.
% We therefore compute the magnitude elementwise from PARTS.
lam_abs = local_qabs(lam);

vnorm = local_qnorm2_cols(V);

resPairAbs = local_res_pair_abs(Aq, lam, V);

den = max(1, (A_fro + lam_abs(:)).*max(1, vnorm(:)));
resPairRel = resPairAbs(:) ./ den;

% resMinAbs already computed. Use the unit vector norm if available; otherwise vnorm.
try
    xnorm = local_qnorm2_cols(V);
catch
    xnorm = ones(K,1);
end

denMin = max(1, (A_fro + lam_abs(:)).*max(1, xnorm(:)));
resMinRel = resMinAbs(:) ./ denMin;

% ---------------- distinct count ----------------
[lamSamples, cls] = local_unique_quats(lamc, opt.UniqueTol);
Kdistinct = numel(lamSamples);

% ---------------- printing ----------------
if opt.Verbose
    fprintf('=== checkNEWTON ===\n');
    fprintf('Matrix: %dx%d quaternion\n', n, n);

    if meta.lambdaProvided
        srcLam = 'provided';
    else
        srcLam = 'computed';
    end
    if meta.vProvided
        srcV = 'provided';
    else
        srcV = 'computed';
    end

    fprintf('Lambdas: %d (%s).  Vectors: %d columns (%s).\n', K, srcLam, size(V,2), srcV);

    fprintf('Certificate/residual summary (median | max):\n');
    fprintf('  resMin abs: %8.2e | %8.2e    resMin rel: %8.2e | %8.2e\n', ...
        median(resMinAbs(:)), max(resMinAbs(:)), median(resMinRel(:)), max(resMinRel(:)));
    fprintf('  resPair abs:%8.2e | %8.2e    resPair rel:%8.2e | %8.2e\n', ...
        median(resPairAbs(:)), max(resPairAbs(:)), median(resPairRel(:)), max(resPairRel(:)));

    fprintf('\nPer-eigenvalue residuals:\n');
    fprintf('  #  lambda (cleaned/rounded)                          resMin(abs/rel)        resPair(abs/rel)\n');
    for k = 1:K
        lamPrint = qroundNEWTON(lamc(k), opt.PrintDigits, "significant");
        tag = '';
        if (resMinAbs(k) > opt.WarnAbs) || (resPairAbs(k) > opt.WarnAbs) || (resMinRel(k) > opt.WarnRel) || (resPairRel(k) > opt.WarnRel)
            tag = '  [LOW ACCURACY: try leigqNEWTON_refine_batch]';
        end
        fprintf('  %2d  %-48s  %9.2e/%9.2e   %9.2e/%9.2e%s\n', ...
            k, local_qstr(lamPrint), resMinAbs(k), resMinRel(k), resPairAbs(k), resPairRel(k), tag);
    end

    fprintf('\nInterestingness:\n');
    if Kdistinct == n
        fprintf('  Standard: Kdistinct = n = %d.\n', n);
    elseif Kdistinct > n
        fprintf('  Nonstandard: Kdistinct = %d > n = %d (more left eigenvalues than dimension).\n', Kdistinct, n);
    else
        fprintf('  Nonstandard: Kdistinct = %d < n = %d. Consider running longer: leigqNEWTON(A,''SolveProfile'',''reliable'').\n', Kdistinct, n);
    end
end

% ---------------- optional sphere sampling ----------------
sphere = struct();
sphere.ran = false;
sphere.lambdaAll = quaternion.empty(0,1);
sphere.lambdaSamples = quaternion.empty(0,1);
sphere.cls = int32([]);
sphere.sph = struct([]);
sphere.info = struct();
sphere.conf = struct();

runSphere = false;
if opt.SphereCheck == "on"
    runSphere = true;
elseif opt.SphereCheck == "auto"
    runSphere = (Kdistinct > n);
end

if runSphere
    sphere.ran = true;
    try
        [lamAllS, lamS, lambda0S, clsS, sphS, infoS] = leigqNEWTON_sphere_sample(Aq, ...
            'Collect', opt.SphereCollect, 'Seed0', opt.SphereSeed0, 'Report', 'off');
        sphere.lambdaAll = lamAllS;
        sphere.lambdaSamples = lamS;
        sphere.lambda0 = lambda0S;
        sphere.cls = clsS;
        sphere.sph = sphS;
        sphere.info = infoS;

        % Validate each detected sphere model (quick, no verbose)
        if ~isempty(sphS)
            Ccase = struct('A', Aq, 'lamAll', lamAllS, 'lamSamples', lamS, 'sph', sphS, 'cls', clsS, 'lambda0', lambda0S, 'info', infoS);
            confAll = cell(1, numel(sphS));
            for si = 1:numel(sphS)
                [~,~,~,~,~,~, confSi, outSi] = leigqNEWTON_sphere_refine(Ccase, [], ...
                    'SphereIndex', si, 'Verbose', 0, 'RefineArgs', {'Mode','none','DoPolish',false,'Verbose',0});
                confAll{si} = struct('sphereIndex', si, 'conf', confSi, 'out', outSi);
            end
            sphere.conf = confAll;
        end

        if opt.Verbose
            if isempty(sphS)
                fprintf('  Sphere check: ran leigqNEWTON_sphere_sample, no sphere model detected.\n');
            else
                fprintf('  Sphere check: detected %d sphere model(s).\n', numel(sphS));
                for si = 1:numel(sphS)
                    sph = sphS(si);
                    fprintf('    Sphere #%d: radius=%.6g  (basis4x3 present=%d)\n', si, local_getfield(sph,'radius',NaN), isfield(sph,'basis4x3'));
                    if isfield(sph,'p0')
                        fprintf('      p0   = %s\n', local_qstr(qroundNEWTON(quaternion(sph.p0(1),sph.p0(2),sph.p0(3),sph.p0(4)), opt.PrintDigits, "significant")));
                    end
                    if isfield(sph,'center3')
                        c3 = sph.center3(:); c3 = double(c3(:));
                        fprintf('      center3 = [%g %g %g]^T\n', c3(1), c3(2), c3(3));
                    end
                    if isfield(sphere,'conf') && numel(sphere.conf) >= si
                        c = sphere.conf{si}.conf;
                        if isfield(c,'sphere')
                            fprintf('      verdict = %s  reliability=%s (%.2f)\n', string(c.sphere.verdict), string(c.sphere.reliabilityTag), double(c.sphere.reliabilityScore));
                        end
                    end
                end

                % Sphere membership (lambdaSamples indices; best-effort)
                try
                    used = false(numel(lamS),1);
                    for si = 1:numel(sphS)
                        if isfield(sphS(si),'inliers') && ~isempty(sphS(si).inliers)
                            ii = sphS(si).inliers(:);
                            ii = ii(ii>=1 & ii<=numel(lamS));
                            used(ii) = true;
                            fprintf('      Sphere #%d: inliers = %d distinct samples\n', si, numel(ii));
                            jj = ii(1:min(numel(ii),6));
                            for t = 1:numel(jj)
                                ltmp = qroundNEWTON(qcleanNEWTON(lamS(jj(t)), opt.CleanTol), opt.PrintDigits, "significant");
                                fprintf('        lamSamples(%d) = %s\n', jj(t), local_qstr(ltmp));
                            end
                        end
                    end
                    iso = find(~used);
                    if ~isempty(iso)
                        fprintf('      Isolated (not assigned to any sphere): %d distinct samples\n', numel(iso));
                    end
                catch
                    % best-effort only
                end

                % Assign provided lambdas to detected spheres (best-effort)
                if ~isempty(lam)
                    assign = local_assign_to_spheres(lam, sphS, 5e-3);
                    nOn = sum(assign > 0);
                    fprintf('    Your %d provided/computed lambdas: %d are within SphereTolRel=5e-3 of a detected sphere model.\n', K, nOn);
                end
            end
        end
    catch ME
        sphere.error = ME;
        if opt.Verbose
            fprintf('  Sphere check: FAILED (%s).\n', ME.message);
        end
    end
else
    if opt.Verbose
        fprintf('  Sphere check: not run. Use checkNEWTON(...,''SphereCheck'',''on'') to force.\n');
    end
end

if opt.Verbose
    fprintf('\nOutputs:\n');
    fprintf('  out           : main struct (out.lambda, out.v, out.resMinAbs, out.resPairAbs, ...).\n');
    fprintf('  cases{1}      : simple case struct with fields A, lamAll, lamSamples, sph, info (sphere tooling style).\n');
    fprintf('  S             : summary struct (n, Ktot, Kdistinct, spheresFound, ...).\n');
end

% ---------------- build outputs ----------------
out = struct();
out.A = Aq;
out.n = n;
out.lambda = lam;
out.v = V;
out.lambdac = lamc;
out.resMinAbs = resMinAbs(:);
out.resMinRel = resMinRel(:);
out.resPairAbs = resPairAbs(:);
out.resPairRel = resPairRel(:);
out.Ktot = K;
out.Kdistinct = Kdistinct;
out.lambdaSamples = lamSamples;
out.cls = cls;
out.solver = struct('resSolver', resSolver, 'info', infoSolver, 'meta', meta);
out.options = opt;
out.sphere = sphere;

C = struct();
C.idx = NaN;
C.seedUsed = NaN;
C.A = Aq;
C.type = 'user';
C.density = local_density(Aq);
C.triPart = local_tri_part(Aq);
C.lamAll = lam;
C.lamSamples = lamSamples;
C.lambda0 = struct();
C.cls = cls;
C.sph = local_getfield(sphere,'sph', struct([]));
C.info = struct('Ktot',K,'Kdistinct',Kdistinct,'n',n,'sphereRan',sphere.ran);

cases = {C};
S = struct();
S.n = n;
S.nMat = 1;
S.Type = "user";
S.Ktot = K;
S.Kdistinct = Kdistinct;
S.spheresFound = numel(C.sph);
S.elapsed = NaN;
S.info = C.info;

end

% =====================================================================
% helpers (stand-alone)
% =====================================================================

function s = local_qstr(q)
% String for a quaternion scalar (robust, readable)
try
    if ~isa(q,'quaternion'), q = quaternion(q); end
    [w,x,y,z] = parts(q);
    s = sprintf('% .6g %+.6gi %+.6gj %+.6gk', double(w), double(x), double(y), double(z));
catch
    s = '<q?>';
end
end

function f = local_getfield(S, name, default)
try
    if isstruct(S) && isfield(S,name)
        f = S.(name);
    else
        f = default;
    end
catch
    f = default;
end
end

function Af = local_qfro(A)
[w,x,y,z] = parts(A);
Af = sqrt( norm(w,'fro')^2 + norm(x,'fro')^2 + norm(y,'fro')^2 + norm(z,'fro')^2 );
Af = double(Af);
end

function qa = local_qabs(q)
% Elementwise quaternion magnitude (works for scalars, vectors, and arrays).
% For q = w + x*i + y*j + z*k, returns sqrt(w^2+x^2+y^2+z^2) componentwise.
try
    q = quaternion(q);
    [w,x,y,z] = parts(q);
    qa = sqrt( double(w).^2 + double(x).^2 + double(y).^2 + double(z).^2 );
catch
    % Fallback for numeric inputs
    qa = abs(double(q));
end
qa = double(qa);
end

function vn = local_qnorm2_cols(V)
% 2-norm per column for quaternion vectors (or matrices treated column-wise)
[w,x,y,z] = parts(V);
K = size(V,2);
vn = zeros(K,1);
for k = 1:K
    vn(k) = sqrt( norm(w(:,k),2)^2 + norm(x(:,k),2)^2 + norm(y(:,k),2)^2 + norm(z(:,k),2)^2 );
end
vn = double(vn);
end

function rp = local_res_pair_abs(A, lam, V)
K = numel(lam);
rp = zeros(K,1);
for k = 1:K
    % MATLAB's built-in quaternion does not support matrix mtimes.
    % Use the stand-alone real-embedding implementation shipped with the toolbox.
    r = qmtimesNEWTON(A, V(:,k)) - qmtimesNEWTON(lam(k), V(:,k));
    [w,x,y,z] = parts(r);
    rp(k) = sqrt( norm(w,2)^2 + norm(x,2)^2 + norm(y,2)^2 + norm(z,2)^2 );
end
rp = double(rp);
end

function [lamU, cls] = local_unique_quats(lam, tol)
% Distinct quaternions by 4D Euclidean distance on components.
lam = quaternion(lam);
lam = lam(:);
K = numel(lam);
if K == 0
    lamU = lam;
    cls = int32([]);
    return;
end
[w,x,y,z] = parts(lam);
P = [double(w(:)), double(x(:)), double(y(:)), double(z(:))];
keep = false(K,1);
cls = zeros(K,1,'int32');
rep = zeros(0,4);
for i = 1:K
    pi = P(i,:);
    if isempty(rep)
        rep = pi;
        keep(i) = true;
        cls(i) = 1;
    else
        d2 = sum((rep - pi).^2,2);
        [dmin, j] = min(d2);
        if sqrt(dmin) <= tol
            cls(i) = int32(j);
        else
            rep = [rep; pi];
            keep(i) = true;
            cls(i) = int32(size(rep,1));
        end
    end
end
lamU = lam(keep);
end

function dens = local_density(A)
try
    [aw,ax,ay,az] = parts(A);
    nz = (aw~=0) | (ax~=0) | (ay~=0) | (az~=0);
    dens = nnz(nz) / max(1, numel(nz));
catch
    dens = NaN;
end
end

function tp = local_tri_part(A)
% best-effort: 'upper'/'lower'/'none'
try
    [aw,ax,ay,az] = parts(A);
    nz = (aw~=0) | (ax~=0) | (ay~=0) | (az~=0);
    if istriu(nz)
        tp = 'upper';
    elseif istril(nz)
        tp = 'lower';
    else
        tp = 'none';
    end
catch
    tp = 'none';
end
end

function asg = local_assign_to_spheres(lam, sphArr, tolRel)
% Assign each lambda to the closest sphere model if devRel <= tolRel.
lam = quaternion(lam);
lam = lam(:);
K = numel(lam);
asg = zeros(K,1,'int32');
if isempty(sphArr)
    return;
end

bestDev = inf(K,1);
for si = 1:numel(sphArr)
    sph = sphArr(si);
    devRel = local_sphere_dev_rel(lam, sph);
    hit = devRel <= tolRel & devRel < bestDev;
    asg(hit) = int32(si);
    bestDev(hit) = devRel(hit);
end
end

function devRel = local_sphere_dev_rel(lam, sph)
% Mirror of refineSPHERE local_sphere_dev_rel (best-effort).
try
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

    X3 = (B' * B) \ (B' * (P4 - p0));
    D  = X3 - c3;
    dist = sqrt(sum(D.^2,1));
    dev = abs(dist - r);
    devRel = reshape(dev./rEff, size(w));
catch
    devRel = NaN(size(lam));
end
end
