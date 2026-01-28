function [lambda, v, res, info] = leigqNEWTON_refine_polish(A, lambda0, v0, varargin)
%LEIGQNEWTON_REFINE_POLISH  Newton polish a quaternion LEFT eigenpair (lambda,v).
%
%   [lambda, v, res] = leigqNEWTON_refine_polish(A, lambda0, v0)
%   [lambda, v, res, info] = leigqNEWTON_refine_polish(A, lambda0, v0, Name,Value,...)
%   [...] = leigqNEWTON_refine_polish(A, lambda0, [])     % if v0 is empty, an initial v0 is constructed
%
% Performs a Newton-type local refinement of an eigenpair candidate for
%     A*v = lambda*v  (LEFT eigenvalue).
% This is the “make it machine-precision” step used after a coarse solver.
%
% Inputs
%   A       : n-by-n quaternion
%   lambda0 : scalar quaternion (initial guess)
%   v0      : n-by-1 quaternion (initial eigenvector guess; may be [])
%
% Name,Value options (selected)
%   'MaxIter'            : Newton iterations (default 20–50, depending on internal defaults)
%   'Tol'                : stopping tolerance on residual
%   'Backtrack'          : true/false line search (default true)
%   'ResidualNormalized' : true/false (default false). Controls reported res.
%
% Outputs
%   lambda : scalar quaternion (polished)
%   v      : n-by-1 quaternion (polished, normalized/gauged)
%   res    : double (final residual norm)
%   info   : struct (iteration history, alphas, residuals; if requested)
%
% Useful one-liners
%   [lamP,vP,resP] = leigqNEWTON_refine_polish(A, lam(1), V(:,1));
%   [lamP,vP,resP] = leigqNEWTON_refine_polish(A, lam(1), []);
%
% See also: leigqNEWTON_init_vec, leigqNEWTON_refine_lambda, leigqNEWTON_cert_resPair

opt.Side = 'left';
opt.TolRes = 1e-14;
opt.TolStep = 1e-14;
opt.MaxIter = 20;
opt.Damping = true;
opt.Verbose = 0;

opt = parseOpts(opt, varargin{:});

n = size(A,1);
assert(size(A,2)==n,'A must be square.');

AR = qA_real_left(A);

% If v0 not provided, compute it from SVD nullvector
if isempty(v0)
    [v0, ninfo] = leigqNEWTON_init_vec(A, lambda0, 'Side', opt.Side);
    p = ninfo.pivot;
else
    [va,vb,vc,vd] = parts(v0);
    mag = sqrt(va.^2 + vb.^2 + vc.^2 + vd.^2);
    [~,p] = max(mag);
end

lambda = lambda0;
v = v0;

% initial gauge+normalize
v = gaugeFix(v,p,opt.Side);
v = normalizeV(v);

[rInf, rvec] = residualInf(AR, lambda, v, opt.Side);
hist.resInf = zeros(opt.MaxIter+1,1);
hist.stepInf = zeros(opt.MaxIter,1);
hist.alpha = zeros(opt.MaxIter,1);
hist.resInf(1) = rInf;

if opt.Verbose
    fprintf('polish: iter %2d  resInf=%.3e\n',0,rInf);
end

for it = 1:opt.MaxIter
    % Build Jacobian blocks
    [M, G] = buildMG(AR, lambda, v, opt.Side); % residual = M*vvec, d/dlam term = -G*dlam

    % Constraints: imag parts of pivot + norm(v)^2-1
    [va,vb,vc,vd] = parts(v);
    vvec = [va;vb;vc;vd];

    c = [ vb(p);
          vc(p);
          vd(p);
          (va.'*va + vb.'*vb + vc.'*vc + vd.'*vd) - 1 ];

    F = [rvec; c];

    % Constraint Jacobian wrt dv
    Cdv = zeros(4, 4*n);
    Cdv(1, n + p)     = 1;         % d(vb(p))
    Cdv(2, 2*n + p)   = 1;         % d(vc(p))
    Cdv(3, 3*n + p)   = 1;         % d(vd(p))
    Cdv(4, :)         = 2*vvec.';  % d(||v||^2-1)

    % Full Jacobian for unknowns [dv(4n); dlam(4)]
    J = [M, -G;
         Cdv, zeros(4,4)];

    dx = -J \ F;
    dv = dx(1:4*n);
    dl = dx(4*n+1:end);

    stepInf = norm(dx, inf);
    hist.stepInf(it) = stepInf;

    % Candidate update
    alpha = 1.0;
    lambda_new = lambda + quaternion(dl(1),dl(2),dl(3),dl(4));
    v_new = addToV(v, dv, n);

    % Gauge+normalize after update
    v_new = gaugeFix(v_new,p,opt.Side);
    v_new = normalizeV(v_new);

    % Damping / backtracking
    if opt.Damping
        [rInf_new, rvec_new] = residualInf(AR, lambda_new, v_new, opt.Side);
        while rInf_new > rInf && alpha > 1/64
            alpha = alpha/2;
            lambda_new = lambda + quaternion(alpha*dl(1),alpha*dl(2),alpha*dl(3),alpha*dl(4));
            v_new = addToV(v, alpha*dv, n);
            v_new = gaugeFix(v_new,p,opt.Side);
            v_new = normalizeV(v_new);
            [rInf_new, rvec_new] = residualInf(AR, lambda_new, v_new, opt.Side);
        end
    else
        [rInf_new, rvec_new] = residualInf(AR, lambda_new, v_new, opt.Side);
    end

    hist.alpha(it) = alpha;

    lambda = lambda_new;
    v = v_new;
    rInf = rInf_new;
    rvec = rvec_new;

    hist.resInf(it+1) = rInf;

    if opt.Verbose
        fprintf('polish: iter %2d  resInf=%.3e  stepInf=%.3e  alpha=%.3g\n',...
            it,rInf,stepInf,alpha);
    end

    if rInf <= opt.TolRes || stepInf <= opt.TolStep
        break;
    end
end

res.resInf = rInf;
res.res2   = norm(rvec,2);

info.iter = it;
info.pivot = p;
info.hist = hist;
info.side = opt.Side;

end

% ---------------- helpers ----------------
function [rInf, rvec] = residualInf(AR, lambda, v, side)
n = size(AR,1)/4;
[va,vb,vc,vd] = parts(v);
vvec = [va;vb;vc;vd];

[a,b,c,d] = parts(lambda);
q = [a;b;c;d];
if strcmpi(side,'left')
    L = qLmat(q);
else
    L = qRmat(q);
end
M = AR - kron(L, eye(n));
rvec = M * vvec;

ra = rvec(1:n); rb = rvec(n+1:2*n); rc = rvec(2*n+1:3*n); rd = rvec(3*n+1:4*n);
rInf = max(sqrt(ra.^2 + rb.^2 + rc.^2 + rd.^2));
end

function [M,G] = buildMG(AR, lambda, v, side)
n = size(AR,1)/4;
[va,vb,vc,vd] = parts(v);

[a,b,c,d] = parts(lambda);
q = [a;b;c;d];

if strcmpi(side,'left')
    L = qLmat(q);
    % G maps dlam -> (dlam)*v  (component-stacked)
    G = [ va, -vb, -vc, -vd;
          vb,  va, -vd,  vc;
          vc,  vd,  va, -vb;
          vd, -vc,  vb,  va ];
else
    L = qRmat(q);
    % G maps dlam -> v*(dlam)
    G = [ va, -vb, -vc, -vd;
          vb,  va,  vd, -vc;
          vc, -vd,  va,  vb;
          vd,  vc, -vb,  va ];
end

M = AR - kron(L, eye(n));
end

function v = addToV(v, dv, n)
[va,vb,vc,vd] = parts(v);
va = va + dv(1:n);
vb = vb + dv(n+1:2*n);
vc = vc + dv(2*n+1:3*n);
vd = vd + dv(3*n+1:4*n);
v = quaternion(va,vb,vc,vd);
end

function v = normalizeV(v)
[va,vb,vc,vd] = parts(v);
nv = sqrt(sum(va.^2 + vb.^2 + vc.^2 + vd.^2));
if nv==0, return; end
v = v ./ nv;
end

function v = gaugeFix(v,p,side)
vp = v(p);

[a,b,c,d] = parts(vp);
m = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
if m==0, return; end

q = quaternion(a,-b,-c,-d) ./ m;   % conj(vp)/|vp|

if strcmpi(side,'left')
    v = v .* q;    % right-gauge keeps LEFT eigenproblem invariant
else
    v = q .* v;    % left-gauge keeps RIGHT eigenproblem invariant
end

[ar,~,~,~] = parts(v(p));
if ar < 0
    v = -v;
end
end

function AR = qA_real_left(A)
[A0,A1,A2,A3] = parts(A);
AR = [ A0, -A1, -A2, -A3;
       A1,  A0, -A3,  A2;
       A2,  A3,  A0, -A1;
       A3, -A2,  A1,  A0 ];
end

function L = qLmat(q)
a=q(1); b=q(2); c=q(3); d=q(4);
L = [ a, -b, -c, -d;
      b,  a, -d,  c;
      c,  d,  a, -b;
      d, -c,  b,  a ];
end

function R = qRmat(q)
a=q(1); b=q(2); c=q(3); d=q(4);
R = [ a, -b, -c, -d;
      b,  a,  d, -c;
      c, -d,  a,  b;
      d,  c, -b,  a ];
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

