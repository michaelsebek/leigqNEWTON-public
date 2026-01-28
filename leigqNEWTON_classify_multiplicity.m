function R = leigqNEWTON_classify_multiplicity(A, lambdaList, varargin)
%LEIGQNEWTON_CLASSIFY_MULTIPLICITY  Deduplicate solver hits and estimate geometric multiplicity (quaternion).
%
%   R = leigqNEWTON_classify_multiplicity(A, lambdaList)
%   R = leigqNEWTON_classify_multiplicity(A, lambdaList, V)                  % positional V
%   R = leigqNEWTON_classify_multiplicity(A, lambdaList, 'Vectors', V)
%   R = leigqNEWTON_classify_multiplicity(..., Name,Value,...)
%
% Clusters duplicate eigenvalue hits (by distance in R^4 of quaternion parts),
% estimates geometric multiplicity from nullity of chi(A-lambda I),
% and optionally clusters eigenvectors by direction up to right scaling.

% --- accept positional V
posVectors = [];
if ~isempty(varargin)
    first = varargin{1};
    isTextKey = ischar(first) || (isstring(first) && isscalar(first));
    if ~isTextKey
        posVectors = first;
        varargin(1) = [];
    end
end

% --- parse Name-Value
p = inputParser;
p.addParameter('Vectors', [], @(x) true);
p.addParameter('LambdaTol', 1e-10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('NullTol', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>=0));
p.addParameter('NullTolFactor', 1e3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('VectorTol', 1e-8, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Verbose', false, @(x)islogical(x)||isnumeric(x));
p.parse(varargin{:});
opt = p.Results;
opt.Verbose = logical(opt.Verbose);

if ~isempty(posVectors) && isempty(opt.Vectors)
    opt.Vectors = posVectors;
end

% --- checks
[n,mA] = size(A);
if n ~= mA, error('A must be square.'); end
if ~isa(lambdaList,'quaternion'), error('lambdaList must be MATLAB quaternion.'); end
lambdaList = lambdaList(:);
m = numel(lambdaList);

V = opt.Vectors;
haveVec = ~isempty(V);
if haveVec
    if ~isa(V,'quaternion') || ~isequal(size(V),[n,m])
        error('Vectors must be an n-by-m quaternion matrix (columns correspond to lambdaList).');
    end
end

% =========================================================================
% (1) CLUSTER lambdas robustly using real parts (NO quaternion array arithmetic)
% =========================================================================
lambdaTol = opt.LambdaTol;

[aL,bL,cL,dL] = parts(lambdaList);
aL = aL(:); bL = bL(:); cL = cL(:); dL = dL(:);

lambdaUnique = quaternion.empty(0,1);
groups = {};
aU = []; bU = []; cU = []; dU = [];

for k = 1:m
    if isempty(aU)
        lambdaUnique(1,1) = lambdaList(k);
        groups{1,1} = k;
        aU = aL(k); bU = bL(k); cU = cL(k); dU = dL(k);
        continue;
    end

    da = aL(k) - aU;
    db = bL(k) - bU;
    dc = cL(k) - cU;
    dd = dL(k) - dU;
    dist = sqrt(da.^2 + db.^2 + dc.^2 + dd.^2);

    [dmin, j] = min(dist);

    if dmin <= lambdaTol
        groups{j,1}(end+1) = k;
    else
        lambdaUnique(end+1,1) = lambdaList(k); %#ok<AGROW>
        groups{end+1,1} = k; %#ok<AGROW>
        aU(end+1,1) = aL(k); %#ok<AGROW>
        bU(end+1,1) = bL(k); %#ok<AGROW>
        cU(end+1,1) = cL(k); %#ok<AGROW>
        dU(end+1,1) = dL(k); %#ok<AGROW>
    end
end

u = numel(lambdaUnique);
hitCount = cellfun(@numel, groups);

% =========================================================================
% (2) Multiplicity via nullity of chi(A - lambda I)
% =========================================================================
geomMultC = zeros(u,1);
geomMultH = zeros(u,1);
vecDirs   = nan(u,1);
vecGroups = cell(u,1);

for j = 1:u
    lam = lambdaUnique(j);

    B = A;
    B(1:n+1:end) = B(1:n+1:end) - lam;
    ChiB = local_qmat2chi(B);

    s = svd(ChiB);

    if isempty(opt.NullTol)
        tol = max(size(ChiB)) * eps(max(s)) * opt.NullTolFactor;
    else
        tol = opt.NullTol;
    end

    nullC = sum(s <= tol);
    geomMultC(j) = nullC;
    geomMultH(j) = round(nullC/2);

    % =========================================================================
    % (3) Optional: group eigenvectors by direction up to right scaling
    % =========================================================================
    if haveVec
        idxHits = groups{j};

        reps = {};
        dirGroups = {};

        for t = 1:numel(idxHits)
            ii = idxHits(t);
            x = V(:,ii);

            if isempty(reps)
                reps{1} = x;
                dirGroups{1} = ii;
                continue;
            end

            [isSame, jBest] = local_same_direction(reps, x, opt.VectorTol);
            if isSame
                dirGroups{jBest}(end+1) = ii; %#ok<AGROW>
            else
                reps{end+1} = x; %#ok<AGROW>
                dirGroups{end+1} = ii; %#ok<AGROW>
            end
        end

        vecDirs(j) = numel(reps);
        vecGroups{j} = dirGroups;
    end

    if opt.Verbose
        fprintf('lambda #%d: hits=%d, nullC=%d => geomMultH~%d', ...
            j, hitCount(j), geomMultC(j), geomMultH(j));
        if haveVec, fprintf(', vecDirs=%d', vecDirs(j)); end
        fprintf('\n');
    end
end

% --- output
R = struct();
R.lambdaUnique = lambdaUnique;
R.groups       = groups;
R.hitCount     = hitCount(:);
R.geomMultC    = geomMultC;
R.geomMultH    = geomMultH;
R.vecDirs      = vecDirs;
R.vecGroups    = vecGroups;
R.tols = struct('LambdaTol', opt.LambdaTol, 'NullTol', opt.NullTol, ...
                'NullTolFactor', opt.NullTolFactor, 'VectorTol', opt.VectorTol);

end

% ========================================================================
% Local helpers
% ========================================================================

function Iq = local_qeye(n)
Iq = quaternion(eye(n), zeros(n), zeros(n), zeros(n));
end

function Chi = local_qmat2chi(Q)
[a,b,c,d] = parts(Q);
alpha = complex(a,b);
beta  = complex(c,d);
Chi = [alpha, beta; -conj(beta), conj(alpha)];
end

function nrm = local_qvecnorm(x)
[a,b,c,d] = parts(x(:));
ax = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
nrm = sqrt(sum(ax.^2));
end

function [tf, jBest] = local_same_direction(reps, x, tol)
tf = false; jBest = 0;

nx = local_qvecnorm(x);
if nx == 0, return; end

bestErr = inf; bestJ = 0;

for j = 1:numel(reps)
    r = reps{j};
    nr = local_qvecnorm(r);
    if nr == 0, continue; end

    num = quaternion(0,0,0,0);
    den = 0;

    for k=1:numel(r)
        rk = r(k);
        xk = x(k);
        num = num + local_qconj(rk) * xk;

        [ar,br,cr,dr] = parts(rk);
        den = den + (ar^2 + br^2 + cr^2 + dr^2);
    end

    if den <= 0, continue; end

    q = num / den;
    err = local_qvecnorm(x - r*q) / (nx + eps);

    if err < bestErr
        bestErr = err;
        bestJ = j;
    end
end

if bestErr <= tol
    tf = true;
    jBest = bestJ;
end
end

function qc = local_qconj(q)
[a,b,c,d] = parts(q);
qc = quaternion(a, -b, -c, -d);
end
