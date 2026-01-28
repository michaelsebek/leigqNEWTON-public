function [lambda_best, lambda_polish, resMin_best, resMin_polished, resInf_polish, info] = leigqNEWTON_refine_auto(A, lambda0, varargin)
%LEIGQNEWTON_REFINE_AUTO  Policy-driven refinement for one or many eigenvalue candidates.
%
%   [lambdaBest, lambdaPolish, rBest, rPolish, resInfPolish, info] = ...
%       leigqNEWTON_refine_auto(A, lambda0, Name,Value,...)
%
% Given candidate(s) lambda0, chooses a refinement strategy (certificate minimization,
% multistart, and optional Newton polishing) to obtain improved eigenvalues with
% certificate values near machine precision when possible.
%
% Inputs
%   A       : n-by-n quaternion
%   lambda0 : scalar or K-by-1 quaternion candidates
%
% Name,Value options (selected)
%   'Mode'         : 'auto'|'none'|'rand+fminsearch'|'fminsearch'|'multistart'|'none' (default 'auto')
%   'TargetResMin' : target for resMin(A,lambda) (default 1e-13)
%   'DoPolish'     : true/false (default true)
%   'Verbose'      : 0/1/2 (default 0)
%
% Outputs
%   lambdaBest     : refined lambdas (K-by-1 quaternion)
%   lambdaPolish   : polished lambdas (K-by-1 quaternion; may equal lambdaBest)
%   rBest          : resMin values for lambdaBest
%   rPolish        : resPair values after polishing (if DoPolish)
%   resInfPolish   : infinity-norm residual after polishing (diagnostic)
%   info           : struct (per-candidate details; chosen modes; timings)
%
% Useful one-liners
%   [lamB,lamP] = leigqNEWTON_refine_auto(A, lam0, 'TargetResMin',1e-14);
%   [lamB,lamP] = leigqNEWTON_refine_auto(A, lam0, 'Mode','none','DoPolish',false);
%
% See also: leigqNEWTON_refine_lambda, leigqNEWTON_refine_polish, leigqNEWTON_cert_resMin

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'Mode','auto',@(s)ischar(s)||isstring(s));

addParameter(p,'TargetResMin',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0) || local_is_auto(x));
addParameter(p,'TargetResMinFactor',5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TargetResMinMin',1e-16,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TargetResMinMax',1e-8,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'HardThresh',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0) || local_is_auto(x));
addParameter(p,'HardThreshFactor',1e3,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'HardThreshMin',1e-8,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'ForceRefine',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Method','svd',@(s)ischar(s)||isstring(s));

addParameter(p,'Radii',[1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4],@(x)isnumeric(x)&&isvector(x)&&all(x>0));
addParameter(p,'NRand',2000,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'PrintEvery',200,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'NStarts',20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Sigma0',0.3,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'MaxIter',600,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolX',1e-16,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TolFun',1e-16,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'Tol',[],@(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0)); % alias: sets TolX/TolFun/TolResPolish

addParameter(p,'DoPolish',[],@(x)islogical(x)&&isscalar(x) || isempty(x));
addParameter(p,'TolResPolish',1e-15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'MaxIterPolish',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(p,'Seed',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
addParameter(p,'Verbose',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'UseExternalRefiner',false,@(x)islogical(x)&&isscalar(x));

parse(p,varargin{:});
opt = p.Results;

% Backward/shortcut alias: 'Tol' sets TolX/TolFun/TolResPolish (unless explicitly provided).
if ~isempty(opt.Tol)
    using = p.UsingDefaults;
    if ismember('TolX', using), opt.TolX = opt.Tol; end
    if ismember('TolFun', using), opt.TolFun = opt.Tol; end
    if ismember('TolResPolish', using), opt.TolResPolish = opt.Tol; end
end

if ~isempty(opt.Seed)
    rng(opt.Seed,'twister');
end

if isempty(opt.DoPolish)
    opt.DoPolish = (exist('leigqNEWTON_refine_polish','file')==2);
end

useExt = opt.UseExternalRefiner && (exist('leigqNEWTON_refine_lambda','file')==2);

% ---------------- basic checks ----------------
n = size(A,1);
if size(A,2) ~= n
    error('A must be square.');
end

lambda0v = lambda0(:);
K = numel(lambda0v);

% ---------------- AUTO thresholds (scale-aware defaults) ----------------
normScale = local_norm_scale(A);  % >= 1, componentwise Frobenius norm in R^4
targetResMinEff = local_effective_target(opt.TargetResMin, opt.TargetResMinFactor, opt.TargetResMinMin, opt.TargetResMinMax, normScale);
hardThreshEff   = local_effective_hard(opt.HardThresh, opt.HardThreshFactor, opt.HardThreshMin, targetResMinEff);

% ---------------- prealloc outputs ----------------
q0 = quaternion(0,0,0,0);
lambda_best   = repmat(q0, K, 1);
lambda_polish = repmat(q0, K, 1);

resMin_best      = zeros(K,1);
resMin_polished  = zeros(K,1);
resInf_polish    = NaN(K,1);

% info struct
info = struct();
info.resMin_orig = zeros(K,1);
info.mode        = strings(K,1);
info.dt_refine   = zeros(K,1);
info.dt_polish   = zeros(K,1);
info.nEvalRand   = zeros(K,1);
info.nEvalFmin   = zeros(K,1);
info.fminExitflag= NaN(K,1);

% scalar diagnostics for reproducibility
info.normScale = normScale;
info.TargetResMinEff = targetResMinEff;
info.HardThreshEff = hardThreshEff;
info.TargetResMinUser = opt.TargetResMin;
info.HardThreshUser = opt.HardThresh;
info.ForceRefine = opt.ForceRefine;

% ---------------- Step 1: compute original certificates (vectorized) ----------------
[resMin_orig, xMin_orig] = leigqNEWTON_cert_resMin(A, lambda0v, 'Method', opt.Method);
info.resMin_orig(:) = resMin_orig(:);

if opt.Verbose
    fprintf('=== %s: n=%d, K=%d ===\n', mfilename, n, K);
    fprintf('TargetResMinEff=%.1e, HardThreshEff=%.1e, normScale=%.3g, Method=%s', targetResMinEff, hardThreshEff, normScale, string(opt.Method));
end

% fminsearch options
fsopt = optimset('Display','off', 'MaxIter',opt.MaxIter, 'TolX',opt.TolX, 'TolFun',opt.TolFun);
if opt.Verbose >= 2
    fsopt = optimset(fsopt,'Display','iter');
end

% ---------------- Step 2: refine each candidate ----------------
for k = 1:K
    t0 = tic;

    lam0 = lambda0v(k);
    r0 = resMin_orig(k);

    if opt.Verbose
        fprintf('\n[k=%d/%d] orig r_min = %.3e\n', k, K, r0);
    end

    % choose mode
modeUser = local_norm_mode(opt.Mode);
if modeUser == "auto"
    if ~opt.ForceRefine && r0 <= targetResMinEff
        mode = "none";
    elseif r0 > hardThreshEff
        mode = "multistart";
    else
        mode = "rand+fmin";
    end
else
    mode = modeUser;
end
info.mode(k) = mode;

    lamB = lam0;
    vB   = xMin_orig(:,k);
    rB   = r0;

    nEvalRand = 0;
    nEvalFmin = 0;
    ef = NaN;

    switch mode
        case "none"
            % keep lam0, but recompute (r_min, v_min) just to be safe
            [rB, vB] = leigqNEWTON_cert_resMin(A, lam0, 'Method', opt.Method);

        case "rand+fmin"
            if useExt
                % Use your existing external routine (if you prefer it)
                [lamB, vB, out] = leigqNEWTON_refine_lambda(A, lam0, ...
                    'NRand', opt.NRand, ...
                    'Radii', opt.Radii, ...
                    'UseFmin', true, ...
                    'MaxIter', opt.MaxIter, ...
                    'Verbose', max(0,opt.Verbose-1));
                rB = out.resMinBest;
                % eval counts unknown here
            else
                % Internal robust variant: random probing then fminsearch
                [lamB, rB, nEvalRand] = local_rand_probe(@(lam) local_resmin_scalar(A, lam, opt.Method), lam0, r0, opt, k, K);

                % local minimization from best random point
                obj = @(x) local_resmin_scalar(A, local_x2q(x), opt.Method);
                x0 = local_q2x(lamB);
                [xB, rB, ef, outFS] = fminsearch(obj, x0, fsopt);
                nEvalFmin = outFS.funcCount;
                lamB = local_x2q(xB);

                % get minimizing vector for lamB
                [rB, vB] = leigqNEWTON_cert_resMin(A, lamB, 'Method', opt.Method);
            end

                case "fmin"
            % Single local optimization (one fminsearch from the current lambda)
            obj = @(x) local_resmin_scalar(A, local_x2q(x), opt.Method);

            xs = local_q2x(lam0);
            fBest = r0;
            xBest = xs;
            ef = NaN;

            [xS, fS, efS, outS] = fminsearch(obj, xs, fsopt);
            nEvalFmin = nEvalFmin + outS.funcCount;

            if fS < fBest
                fBest = fS;
                xBest = xS;
                ef = efS;
            end

            lamB = local_x2q(xBest);
            [rB, vB] = leigqNEWTON_cert_resMin(A, lamB, 'Method', opt.Method);

case "multistart"
            obj = @(x) local_resmin_scalar(A, local_x2q(x), opt.Method);

            xBest = local_q2x(lam0);
            fBest = r0;

            for s = 1:opt.NStarts
                if s == 1
                    xs = local_q2x(lam0);
                else
                    xs = local_q2x(lam0) + opt.Sigma0*randn(4,1);
                end

                if opt.Verbose
                    fprintf('  multistart %2d/%2d (best=%.3e)\n', s, opt.NStarts, fBest);
                end

                [xS, fS, efS, outS] = fminsearch(obj, xs, fsopt);
                nEvalFmin = nEvalFmin + outS.funcCount;

                if fS < fBest
                    fBest = fS;
                    xBest = xS;
                    ef = efS;

                    if opt.Verbose
                        fprintf('    improved -> %.3e\n', fBest);
                    end

                    if fBest <= targetResMinEff
                        if opt.Verbose
                            fprintf('    target reached -> early stop.\n');
                        end
                        break;
                    end
                end
            end

            lamB = local_x2q(xBest);
            % get certified r_min + minimizing vector
            [rB, vB] = leigqNEWTON_cert_resMin(A, lamB, 'Method', opt.Method);

        otherwise
            error('Internal error: unknown mode "%s".', mode);
    end

    info.dt_refine(k) = toc(t0);
    info.nEvalRand(k) = nEvalRand;
    info.nEvalFmin(k) = nEvalFmin;
    info.fminExitflag(k) = ef;

    lambda_best(k)  = lamB;
    resMin_best(k)  = rB;

    if opt.Verbose
        fprintf('[k=%d] refined: r_min=%.3e  (mode=%s, dt=%.2fs)\n', k, rB, mode, info.dt_refine(k));
    end

    % ---------------- Step 3: polish (optional) ----------------
    if opt.DoPolish
        tp = tic;
        [lamP, vP, rk] = leigqNEWTON_refine_polish(A, lamB, vB, ...
            'TolRes', opt.TolResPolish, 'MaxIter', opt.MaxIterPolish, 'Verbose', max(0,opt.Verbose-1));
        info.dt_polish(k) = toc(tp);

        lambda_polish(k) = lamP;
        resInf_polish(k) = rk.resInf;

        if opt.Verbose
            fprintf('[k=%d] polish: resInf=%.3e (dt=%.2fs)\n', k, rk.resInf, info.dt_polish(k));
        end
    else
        lambda_polish(k) = lamB;
        resInf_polish(k) = NaN;
    end
end

% ---------------- Step 4: final certificate for polished lambdas (vectorized) ----------------
resMin_polished = leigqNEWTON_cert_resMin(A, lambda_polish, 'Method', opt.Method);

% return shapes consistent with lambda0 input
lambda_best   = reshape(lambda_best,   size(lambda0));
lambda_polish = reshape(lambda_polish, size(lambda0));
resMin_best      = reshape(resMin_best,     size(lambda0));
resMin_polished  = reshape(resMin_polished, size(lambda0));
resInf_polish    = reshape(resInf_polish,   size(lambda0));

end

% ===== local helpers =====

function f = local_resmin_scalar(A, lam, method)
% 1-output call => scalar
f = leigqNEWTON_cert_resMin(A, lam, 'Method', method);
end

function x = local_q2x(q)
[a,b,c,d] = parts(q);
x = [a;b;c;d];
end

function q = local_x2q(x)
q = quaternion(x(1),x(2),x(3),x(4));
end

function [lamB, fB, nEval] = local_rand_probe(resmin_scalar, lam0, f0, opt, k, K)
lamB = lam0;
fB = f0;
nEval = 0;

xB = local_q2x(lamB);

for ir = 1:numel(opt.Radii)
    r = opt.Radii(ir);

    if opt.Verbose
        fprintf('  [k=%d/%d] rand: radius=%g, NRand=%d\n', k, K, r, opt.NRand);
    end

    tR = tic;

    for t = 1:opt.NRand
        g = randn(4,1);
        ng = norm(g);
        if ng == 0, continue; end
        x = xB + r*(g/ng);

        lam = local_x2q(x);
        f = resmin_scalar(lam);
        nEval = nEval + 1;

        if f < fB
            fB = f;
            lamB = lam;
            xB = x;
        end

        if opt.Verbose && mod(t,opt.PrintEvery)==0
            fprintf('    t=%5d/%5d  best=%.3e  dt=%.2fs\n', t, opt.NRand, fB, toc(tR));
        end
    end

    if opt.Verbose
        fprintf('    done radius=%g -> best=%.3e (dt=%.2fs)\n', r, fB, toc(tR));
    end
end
end


function s = local_norm_scale(A)
% Scale proxy for setting "AUTO" thresholds.
% We intentionally use a cheap, robust proxy: ||A||_F computed componentwise in R^4.
try
    [w,x,y,z] = parts(A);
    s = sqrt(sum(w(:).^2 + x(:).^2 + y(:).^2 + z(:).^2));
    if ~isfinite(s) || s <= 0
        s = 1;
    end
catch
    s = 1;
end
s = max(1, s);
end

function m = local_norm_mode(mode)
%LOCAL_NORM_MODE  Normalize the user-facing Mode option to internal labels.
%
% Accepted values (case-insensitive):
%   'auto' (default)   : choose based on r0 and thresholds
%   'none'             : no refinement
%   'rand+fmin'        : random probe + fminsearch (internal default when r0 is moderate)
%   'multistart'       : multistart fminsearch around lambda0 (internal default when r0 is large)
%   'fmin'             : one fminsearch from lambda0
%
% Synonyms are accepted, e.g. 'rand+fminsearch', 'fminsearch'.

mode = string(mode);
mode = lower(strtrim(mode));

switch mode
    case {""}
        m = "auto";

    case {"auto"}
        m = "auto";

    case {"none","off","no","skip"}
        m = "none";

    case {"rand+fmin","rand+fminsearch","randfmin","rand+fminsearch"}
        m = "rand+fmin";

    case {"multistart","multi","ms"}
        m = "multistart";

    case {"fmin","fminsearch","local"}
        m = "fmin";

    otherwise
        error('Unknown Mode "%s". Use ''auto'', ''none'', ''rand+fminsearch'', ''fminsearch'', or ''multistart''.', mode);
end
end


function tf = local_is_auto(x)
tf = isempty(x);
if tf, return; end
if isstring(x) || ischar(x)
    xs = lower(string(x));
    tf = (xs == "auto") || (xs == "default");
else
    tf = false;
end
end

function tEff = local_effective_target(tUser, factor, tMin, tMax, normScale)
% AUTO: tEff = clip( factor * eps(normScale), tMin, tMax )
if local_is_auto(tUser)
    tEff = factor * eps(normScale);
else
    tEff = double(tUser);
end
tEff = max(tMin, min(tMax, tEff));
end

function hEff = local_effective_hard(hUser, hFactor, hMin, tEff)
% AUTO: hEff = max(hMin, hFactor * tEff)
if local_is_auto(hUser)
    hEff = max(hMin, hFactor * tEff);
else
    hEff = double(hUser);
end
end
