function [lambdaRef, vRef, cert, info] = leigqNEWTON_refine_batch(A, lambda0, varargin)
%LEIGQNEWTON_REFINE_BATCH  Refine + certify a list of eigenvalue candidates.
%
%   [lambdaRef, vRef] = leigqNEWTON_refine_batch(A, lambda0)
%   [lambdaRef, vRef, cert] = leigqNEWTON_refine_batch(A, lambda0, Name,Value,...)
%
% Batch wrapper that turns coarse candidates lambda0 into refined/certified
% eigenvalues (and optionally polished eigenvectors). This is the recommended
% “production” post-processing step after leigqNEWTON.
%
% Inputs
%   A       : n-by-n quaternion
%   lambda0 : K-by-1 quaternion candidates (raw or cleaned)
%
% Name,Value options (selected)
%   'Mode'         : 'auto'|'rand+fminsearch'|'fminsearch'|'none' (default 'auto')
%   'TargetResMin' : target certificate value (default 1e-13)
%   'DoPolish'     : true/false (default true)
%
% Outputs
%   lambdaRef : K-by-1 quaternion (refined)
%   vRef      : n-by-K quaternion (associated vectors; columns)
%   cert      : struct with fields such as:
%                 resMin, resPair, lambdaClean, acceptedMask, timings, ...
%
% Useful one-liners
%   [lamR,vR,cert] = leigqNEWTON_refine_batch(A, lam);
%   [lamR,vR,cert] = leigqNEWTON_refine_batch(A, lam, 'DoPolish',true,'TargetResMin',1e-14);
%
% See also: leigqNEWTON_refine_auto, leigqNEWTON_refine_lambda, leigqNEWTON_cert_resPair

p = inputParser;
p.FunctionName = mfilename;

addParameter(p,'Mode','auto',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetResMin',1e-13,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Method','svd',@(s)ischar(s)||isstring(s));

addParameter(p,'MaxIter',600,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolX',1e-16,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolFun',1e-16,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Tol',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0)); % alias: sets TolX/TolFun/TolResPolish

addParameter(p,'Seed',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
addParameter(p,'Verbose',1,@(x)isnumeric(x)&&isscalar(x));

addParameter(p,'Radii',[1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4],@(x)isnumeric(x)&&isvector(x)&&all(x>0));
addParameter(p,'NRand',2000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'PrintEvery',200,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'NStarts',20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Sigma0',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'DoPolish',[],@(x)islogical(x)&&isscalar(x) || isempty(x));
addParameter(p,'TolResPolish',1e-15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MaxIterPolish',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);

parse(p,varargin{:});
opt = p.Results;

% Backward/shortcut alias: 'Tol' sets TolX/TolFun/TolResPolish (unless explicitly provided).
if ~isempty(opt.Tol)
    using = p.UsingDefaults;
    if ismember('TolX', using), opt.TolX = opt.Tol; end
    if ismember('TolFun', using), opt.TolFun = opt.Tol; end
    if ismember('TolResPolish', using), opt.TolResPolish = opt.Tol; end
end

mode = char(opt.Mode);
target = opt.TargetResMin;

if ~isempty(opt.Seed)
    rng(opt.Seed,'twister');
end

% ---------- shape ----------
sz0 = size(lambda0);
lambda0v = lambda0(:);
K = numel(lambda0v);

n = size(A,1);
if size(A,2) ~= n
    error('A must be square.');
end

% Decide polishing default
if isempty(opt.DoPolish)
    opt.DoPolish = (exist('leigqNEWTON_refine_polish','file')==2);
end

% fminsearch options
fsopt = optimset('Display','off', 'TolX',opt.TolX, 'TolFun',opt.TolFun, 'MaxIter',opt.MaxIter);

% ---------- helpers ----------
resmin_scalar = @(lam) leigqNEWTON_cert_resMin(A, lam, 'Method', opt.Method); % 1 output => scalar

q2x = @(q) local_q2x(q);
x2q = @(x) quaternion(x(1),x(2),x(3),x(4));

% ---------- main loop refine each candidate ----------
lambdaBest = lambda0v;
info = struct;
info.perK = repmat(struct(), K, 1);

if opt.Verbose
    fprintf('=== %s: start (n=%d, K=%d, mode=%s, target=%.1e) ===\n', ...
        mfilename, n, K, mode, target);
end

for k = 1:K
    t0 = tic;

    lam0 = lambda0v(k);
    f0 = resmin_scalar(lam0);

    if opt.Verbose
        fprintf('\n[k=%d/%d] start r_min=%.3e\n', k, K, f0);
    end

    % pick mode automatically if requested
    modeK = mode;
    if strcmpi(modeK,'auto')
        if f0 <= target
            modeK = 'none';
        elseif f0 > 1e-10
            modeK = 'multistart';
        else
            modeK = 'rand+fminsearch';
        end
    end

    lamB = lam0;
    fB = f0;

    nEvalRand = 0;
    nEvalFmin = 0;
    fminOutput = [];
    fminExitflag = NaN;

    % ---- mode execution ----
    if strcmpi(modeK,'none')
        % keep as is

    elseif strcmpi(modeK,'fminsearch')
        x0 = q2x(lam0);
        obj = @(x) resmin_scalar(x2q(x));
        [xB, fB, fminExitflag, fminOutput] = fminsearch(obj, x0, fsopt);
        nEvalFmin = fminOutput.funcCount;
        lamB = x2q(xB);

        if opt.Verbose
            fprintf('[k=%d] fminsearch done: r_min=%.3e (exitflag=%g)\n', k, fB, fminExitflag);
        end

    elseif strcmpi(modeK,'rand+fminsearch')
        % random probing around current best
        [lamB, fB, nEvalRand] = local_rand_probe(resmin_scalar, lam0, f0, opt, k, K);

        % local optimization from best random
        x0 = q2x(lamB);
        obj = @(x) resmin_scalar(x2q(x));
        if opt.Verbose
            fprintf('[k=%d] fminsearch starting from r_min=%.3e ...\n', k, fB);
            fsopt2 = optimset(fsopt,'Display','iter');
        else
            fsopt2 = fsopt;
        end

        [xB, fB, fminExitflag, fminOutput] = fminsearch(obj, x0, fsopt2);
        nEvalFmin = fminOutput.funcCount;
        lamB = x2q(xB);

        if opt.Verbose
            fprintf('[k=%d] after fminsearch: r_min=%.3e (exitflag=%g)\n', k, fB, fminExitflag);
        end

    elseif strcmpi(modeK,'multistart')
        % multistart fminsearch in R^4
        obj = @(x) resmin_scalar(x2q(x));
        xBest = q2x(lam0);
        fBest = f0;
        nEvalFmin = 0;

        if opt.Verbose
            fprintf('[k=%d] multistart: NStarts=%d, Sigma0=%.3g\n', k, opt.NStarts, opt.Sigma0);
        end

        for s = 1:opt.NStarts
            if s == 1
                xs = q2x(lam0);
            else
                xs = q2x(lam0) + opt.Sigma0*randn(4,1);
            end

            if opt.Verbose
                fprintf('  start %2d/%2d ... (best=%.3e)\n', s, opt.NStarts, fBest);
            end

            [xS, fS, efS, outS] = fminsearch(obj, xs, fsopt);
            nEvalFmin = nEvalFmin + outS.funcCount;

            if fS < fBest
                fBest = fS;
                xBest = xS;
                if opt.Verbose
                    fprintf('    improved: r_min=%.3e (exitflag=%g)\n', fBest, efS);
                end
                if fBest <= target
                    if opt.Verbose
                        fprintf('    target reached -> early stop.\n');
                    end
                    break;
                end
            end
        end

        lamB = x2q(xBest);
        fB = fBest;

        if opt.Verbose
            fprintf('[k=%d] multistart done: r_min=%.3e\n', k, fB);
        end

    else
        error('Unknown Mode: %s', modeK);
    end

    lambdaBest(k) = lamB;

    % store info per k
    info.perK(k).mode = modeK;
    info.perK(k).resMin0 = f0;
    info.perK(k).resMinBest = fB;
    info.perK(k).nEvalRand = nEvalRand;
    info.perK(k).nEvalFmin = nEvalFmin;
    info.perK(k).fminExitflag = fminExitflag;
    info.perK(k).fminOutput = fminOutput;
    info.perK(k).dt = toc(t0);

    if opt.Verbose
        fprintf('[k=%d] DONE: r_min %.3e -> %.3e   (dt=%.2fs)\n', ...
            k, f0, fB, info.perK(k).dt);
    end
end

% reshape back
lambdaRef = reshape(lambdaBest, sz0);

% ---------- final certification for all lambdas at once ----------
[resMin, xMin_all] = leigqNEWTON_cert_resMin(A, lambdaRef, 'Method', opt.Method);

% respair using returned minimizers (works well as a certificate for the pair)
resPair = zeros(K,1);
for k = 1:K
    resPair(k) = leigqNEWTON_cert_resPair(A, lambdaBest(k), xMin_all(:,k));
end

% optional polish
resInf = NaN(K,1);
if opt.DoPolish
    for k = 1:K
        v0 = xMin_all(:,k);
        [lamP, vP, rk] = leigqNEWTON_refine_polish(A, lambdaBest(k), v0, ...
            'TolRes', opt.TolResPolish, 'MaxIter', opt.MaxIterPolish, 'Verbose', 0);
        lambdaBest(k) = lamP; %#ok<AGROW>
        xMin_all(:,k) = vP;   %#ok<AGROW>
        resInf(k) = rk.resInf;
    end
    lambdaRef = reshape(lambdaBest, sz0);
    % refresh cert after polish (resMin won't necessarily decrease further; respair/resInf is the focus)
    [resMin, xMin_all] = leigqNEWTON_cert_resMin(A, lambdaRef, 'Method', opt.Method);
    for k = 1:K
        resPair(k) = leigqNEWTON_cert_resPair(A, lambdaBest(k), xMin_all(:,k));
    end
end

vRef = xMin_all;

% distances in R^4 (useful to show "distinct")
D = zeros(K);
for i = 1:K
    for j = 1:K
        D(i,j) = local_qdist(lambdaBest(i), lambdaBest(j));
    end
end
D2 = D + diag(inf(K,1));
[minOff, idx] = min(D2(:));
[iMin, jMin] = ind2sub([K K], idx);

% output cert struct
cert = struct();
cert.resMin = resMin(:);
cert.resPair = resPair(:);
cert.resInf = resInf(:);
cert.dist = D;
cert.minDist = minOff;
cert.minDistPair = [iMin, jMin];

if opt.Verbose
    fprintf('\n=== %s: summary ===\n', mfilename);
    for k = 1:K
        fprintf('k=%d  r_min=%.3e  respair=%.3e', k, cert.resMin(k), cert.resPair(k));
        if ~isnan(cert.resInf(k))
            fprintf('  resInf(polish)=%.3e', cert.resInf(k));
        end
        fprintf('\n');
    end
    fprintf('min pairwise distance = %.3e between (%d,%d)\n', cert.minDist, cert.minDistPair(1), cert.minDistPair(2));
end

end

% ===== local helpers =====

function x = local_q2x(q)
[a,b,c,d] = parts(q);
x = [a;b;c;d];
end

function d = local_qdist(p,q)
r = p - q;
[da,db,dc,dd] = parts(r);
d = sqrt(da.^2 + db.^2 + dc.^2 + dd.^2);
d = double(d);
end

function [lamB, fB, nEval] = local_rand_probe(resmin_scalar, lam0, f0, opt, k, K)
lamB = lam0;
fB = f0;
nEval = 0;

for ir = 1:numel(opt.Radii)
    r = opt.Radii(ir);
    if opt.Verbose
        fprintf('[k=%d/%d] rand probe: radius=%g, NRand=%d\n', k, K, r, opt.NRand);
    end
    tR = tic;

    for t = 1:opt.NRand
        % random direction in R^4
        g = randn(4,1);
        ng = norm(g);
        if ng == 0, continue; end
        x = local_q2x(lamB) + r*(g/ng);

        lam = quaternion(x(1),x(2),x(3),x(4));
        f = resmin_scalar(lam);
        nEval = nEval + 1;

        if f < fB
            fB = f;
            lamB = lam;
        end

        if opt.Verbose && (mod(t,opt.PrintEvery)==0)
            fprintf('  t=%5d/%5d  best=%.3e  dt=%.2fs\n', t, opt.NRand, fB, toc(tR));
        end
    end

    if opt.Verbose
        fprintf('  done radius=%g -> best=%.3e (dt=%.2fs)\n', r, fB, toc(tR));
    end
end
end
