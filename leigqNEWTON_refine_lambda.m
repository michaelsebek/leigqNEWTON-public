function [lambdaBest, vBest, out] = leigqNEWTON_refine_lambda(A, lambda0, varargin)
%LEIGQNEWTON_REFINE_LAMBDA  Refine a single eigenvalue candidate (lambda-only workflow).
%
%   [lambdaBest, vBest] = leigqNEWTON_refine_lambda(A, lambda0)
%   [lambdaBest, vBest, cert] = leigqNEWTON_refine_lambda(A, lambda0, Name,Value,...)
%
% Refines one candidate left eigenvalue lambda0 by minimizing the certificate-like
% residual resMin(A,lambda). Optionally, it can also polish the resulting eigenpair.
%
% Inputs
%   A       : n-by-n quaternion
%   lambda0 : scalar quaternion (initial candidate)
%
% Name,Value options (selected)
%   'Mode'        : 'auto'|'rand+fminsearch'|'fminsearch'|'none' (default 'auto')
%   'TargetResMin': stopping target for resMin (default 1e-13)
%   'DoPolish'    : true/false (default true)
%   'Seed'        : RNG seed for random starts (default [])
%
% Outputs
%   lambdaBest : scalar quaternion
%   vBest      : n-by-1 quaternion (best associated vector; may be from polishing)
%   cert       : struct with fields such as resMin, resPair, v0Used, etc.
%
% Useful one-liners
%   [lam1,v1] = leigqNEWTON_refine_lambda(A, lam(1));
%   [lam1,v1] = leigqNEWTON_refine_lambda(A, lam(1), 'DoPolish',false);
%
% See also: leigqNEWTON_cert_resMin, leigqNEWTON_refine_polish, leigqNEWTON_refine_batch

opt.Radii   = [];
opt.NRand   = 2000;
opt.UseFmin = true;
opt.TolFun  = 1e-16;
opt.TolX    = 1e-16;
opt.MaxIter = 600;
opt.Verbose = 1;
opt.Seed    = [];
opt.Method  = 'svd';

opt = parseOpts(opt, varargin{:});

if ~isempty(opt.Seed)
    rng(opt.Seed);
end

% Ensure quaternion scalar
if ~isa(lambda0,'quaternion')
    lambda0 = quaternion(lambda0,0,0,0);
end
if ~isscalar(lambda0)
    error('lambda0 must be a scalar quaternion.');
end

[a0,b0,c0,d0] = parts(lambda0);
x0 = [a0;b0;c0;d0];

% auto radii: scale with |lambda0| and ||A|| (very mild)
lamMag = sqrt(a0^2+b0^2+c0^2+d0^2);
if isempty(opt.Radii)
    base = max(1, lamMag);
    opt.Radii = base * [1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5];
end

obj = @(x) resmin_scalar(A, x, opt.Method);

f0 = obj(x0);
xBest = x0;
fBest = f0;

if opt.Verbose
    fprintf('refine: start resMin=%.3e\n', f0);
end

hist = struct('stage',{},'radius',{},'k',{},'resMin',{},'a',{},'b',{},'c',{},'d',{});
pushHist('init', 0, 0, f0, x0);

% ---- (1) randomized local search ----
for ir = 1:numel(opt.Radii)
    r = opt.Radii(ir);
    for k = 1:opt.NRand
        dx = r * randn(4,1);
        x  = xBest + dx;
        fx = obj(x);
        pushHist('rand', r, k, fx, x);
        if fx < fBest
            fBest = fx;
            xBest = x;
        end
    end
    if opt.Verbose
        fprintf('refine: after rand r=%.3g -> best resMin=%.3e\n', r, fBest);
    end
end

% ---- (2) Nelder-Mead refinement ----
if opt.UseFmin
    if opt.Verbose
        fprintf('refine: fminsearch...\n');
    end
    fopts = optimset('Display','off','TolFun',opt.TolFun,'TolX',opt.TolX,'MaxIter',opt.MaxIter);
    [xF, fF] = fminsearch(obj, xBest, fopts);
    pushHist('fmin', 0, 0, fF, xF);
    if fF < fBest
        fBest = fF;
        xBest = xF(:);
    end
    if opt.Verbose
        fprintf('refine: after fminsearch -> best resMin=%.3e\n', fBest);
    end
end

lambdaBest = quaternion(xBest(1), xBest(2), xBest(3), xBest(4));

% compute best v (xMin) for the final lambda
[resMinBest, vBest] = leigqNEWTON_cert_resMin(A, lambdaBest, 'Method', opt.Method);

out = struct();
out.lambda0 = lambda0;
out.resMin0 = f0;
out.lambdaBest = lambdaBest;
out.resMinBest = resMinBest;
out.hist = hist;
out.method = opt.Method;

    function pushHist(stage, radius, kk, fval, xvec)
        S.stage = stage;
        S.radius = radius;
        S.k = kk;
        S.resMin = fval;
        S.a = xvec(1); S.b = xvec(2); S.c = xvec(3); S.d = xvec(4);
        hist(end+1) = S; %#ok<AGROW>
    end
end

% ===== objective: scalar resMin for given x=[a;b;c;d] =====
function f = resmin_scalar(A, x, method)
lam = quaternion(x(1),x(2),x(3),x(4));
f = leigqNEWTON_cert_resMin(A, lam, 'Method', method); % 1 output => scalar
end

function opt = parseOpts(opt, varargin)
if mod(numel(varargin),2)~=0, error('Name-value pairs expected.'); end
for k=1:2:numel(varargin)
    name = varargin{k};
    val  = varargin{k+1};
    if ~isfield(opt,name), error('Unknown option: %s',name); end
    opt.(name) = val;
end
end
