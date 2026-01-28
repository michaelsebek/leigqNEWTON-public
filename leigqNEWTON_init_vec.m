function [v0, info] = leigqNEWTON_init_vec(A, lambda, varargin)
%LEIGQNEWTON_INIT_VEC  Construct an initial eigenvector guess for a quaternion LEFT eigenvalue candidate.
%
%   v0 = leigqNEWTON_init_vec(A, lambda)
%   [v0, info] = leigqNEWTON_init_vec(A, lambda, Name,Value,...)
%
% Given A and a candidate left eigenvalue lambda, constructs a reasonable initial
% right eigenvector guess v0 for use by polishing/refinement routines.
%
% Inputs
%   A      : n-by-n quaternion (or numeric; interpreted as real quaternion)
%   lambda : scalar quaternion
%
% Name,Value options (selected)
%   'Method' : 'svd'|'null'|'eigs' (default 'svd')   % how to obtain a null-like vector
%   'Tol'    : tolerance for rank/null decisions
%
% Outputs
%   v0   : n-by-1 quaternion (normalized)
%   info : struct (method used; diagnostics)
%
% Useful one-liners
%   v0 = leigqNEWTON_init_vec(A, lam(1));
%   v0 = leigqNEWTON_init_vec(A, lam(1), 'Method','svd');
%
% See also: leigqNEWTON_refine_polish, leigqNEWTON_cert_resPair

opt.Side = 'left';
opt.Pivot = [];
opt.Gauge = true;
opt.Normalize = true;

opt = parseOpts(opt, varargin{:});

n = size(A,1);
assert(size(A,2)==n,'A must be square.');

AR = qA_real_left(A); % 4n x 4n mapping for v -> A*v in component-stacked form

% Lambda to 4-vector
[a,b,c,d] = parts(lambda);
if strcmpi(opt.Side,'left')
    L = qLmat([a;b;c;d]);
else
    L = qRmat([a;b;c;d]);
end

% Real system matrix for residual: (A*v - lambda*v) or (A*v - v*lambda)
M = AR - kron(L, eye(n));  % 4n x 4n

% Smallest right singular vector of M
[U,S,V] = svd(M,'econ'); %#ok<ASGLU>
x = V(:,end);

[va,vb,vc,vd] = split4(x,n);
v0 = quaternion(va, vb, vc, vd);

% choose pivot
if isempty(opt.Pivot)
    mag = sqrt(va.^2 + vb.^2 + vc.^2 + vd.^2);
    [~,p] = max(mag);
else
    p = opt.Pivot;
end

% Gauge-fix (optional): make pivot component real-positive
if opt.Gauge
    v0 = gaugeFix(v0,p,opt.Side);
end

% Normalize (optional)
if opt.Normalize
    v0 = normalizeV(v0);
end

% Info
info.pivot = p;
info.sigmaMin = S(end,end); % smallest singular value from econ SVD

end

% ---------------- helpers ----------------
function AR = qA_real_left(A)
% component-stacked real embedding for quaternion matrix-vector product A*v.
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

function [v0,v1,v2,v3] = split4(x,n)
v0 = x(1:n);
v1 = x(n+1:2*n);
v2 = x(2*n+1:3*n);
v3 = x(3*n+1:4*n);
end

function v = gaugeFix(v,p,side)
% Make pivot component v(p) real-positive by multiplying by a unit quaternion.
vp = v(p);

[a,b,c,d] = parts(vp);
m = sqrt(a.^2 + b.^2 + c.^2 + d.^2);
if m==0, return; end

q = quaternion(a,-b,-c,-d) ./ m; % conj(vp)/|vp|, avoids abs/conj overloads

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

function v = normalizeV(v)
[va,vb,vc,vd] = parts(v);
nv = sqrt(sum(va.^2 + vb.^2 + vc.^2 + vd.^2));
if nv==0, return; end
v = v ./ nv;
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
