function [resPair, resPairRaw, vUnit, rPair, report] = leigqNEWTON_cert_resPair(A, lambda, v, varargin)
%LEIGQNEWTON_CERT_RESPAIR  Residual certificate for a quaternion LEFT eigenpair (lambda,v).
%
%   resPair = leigqNEWTON_cert_resPair(A, lambda, v)
%   [resPair, rVec] = leigqNEWTON_cert_resPair(A, lambda, v)
%   [resPair, rVec, info] = leigqNEWTON_cert_resPair(A, lambda, v, Name,Value,...)
%
% Computes ||A*v - lambda*v||_2 for one or many eigenpairs. If v is provided as
% n-by-K (columns), lambda may be scalar or K-by-1.
%
% Inputs
%   A      : n-by-n quaternion (or numeric; interpreted as real quaternion)
%   lambda : scalar quaternion or K-by-1 quaternion array
%   v      : n-by-1 or n-by-K quaternion array (eigenvectors in columns)
%
% Name,Value options (selected)
%   'NormalizeV'         : true/false (default true). If true, normalize each column of v.
%   'ResidualNormalized' : true/false (default false). If true, return a scale-invariant residual.
%
% Outputs
%   resPair : K-by-1 double     (residual norm(s))
%   rVec    : cell array        (optional residual vectors; for diagnostics)
%   info    : struct            (diagnostics; normalization; denominators if used)
%
% Useful one-liners
%   r = leigqNEWTON_cert_resPair(A, lam(1), V(:,1));
%   r = leigqNEWTON_cert_resPair(A, lam, V, 'NormalizeV',true);
%
% See also: leigqNEWTON_cert_resMin, leigqNEWTON_refine_polish, leigqNEWTON

if ~exist('quaternion','class')
    error('leigqNEWTON_cert_resPair:NoQuaternionClass', ...
        'MATLAB class "quaternion" not found. Requires Aerospace Toolbox.');
end

% ---------------- parse options ----------------
opt.NormalizeV         = true;
opt.ResidualNormalized = false;
opt.LambdaC            = [];
opt.UseLambda          = 'lambda';
opt.ResExpected        = [];

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('leigqNEWTON_cert_resPair:BadNV', 'Name/value inputs must come in pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('leigqNEWTON_cert_resPair:BadNV', 'Option names must be char or string.');
        end
        key = lower(char(name));
        switch key
            case {'normalizev','normalize','unitizev'}
                opt.NormalizeV = logical(val);
            case {'residualnormalized','normalized','scalefree','scaleinvariant'}
                opt.ResidualNormalized = logical(val);
            case {'lambdac','lamc','lambdaclean','cleanlambda','cleanedlambda'}
                opt.LambdaC = val;
            case {'uselambda','use','whichlambda','lamtype'}
                opt.UseLambda = lower(char(val));
            case {'resexpected','res0','resin','resinput','resprovided','resfromsolver','ressolver','resvec','resvector','res'}
                opt.ResExpected = double(val);
            otherwise
                % ignore unknown options (consistent with leigqNEWTON)
        end
    end
end

% ---------------- type normalization ----------------
Aq = local_to_quat(A);
[n,m] = size(Aq);
if n ~= m
    error('leigqNEWTON_cert_resPair:BadInput', 'A must be square.');
end

lamIn = local_to_quat(lambda);
lamIn = lamIn(:);

% choose which eigenvalue vector to use
useLam = opt.UseLambda;
if isempty(lamIn) && ~isempty(opt.LambdaC)
    useLam = 'lambdac';
end
switch useLam
    case {'lambda','lam','raw'}
        lamUse = lamIn;
    case {'lambdac','lamc','clean','cleaned'}
        if isempty(opt.LambdaC)
            error('leigqNEWTON_cert_resPair:MissingLambdaC', ...
                'UseLambda="lambdac" requested but LambdaC was not provided.');
        end
        lamUse = local_to_quat(opt.LambdaC);
        lamUse = lamUse(:);
    otherwise
        % fall back to input lambda (do not error; match vanilla style)
        lamUse = lamIn;
end

K = numel(lamUse);

Vq = local_to_quat(v);
if isempty(Vq)
    resPair = zeros(0,1);
    resPairRaw = zeros(0,1);
    vUnit = quaternion.empty(n,0);
    if nargout >= 4
        rPair = quaternion.empty(n,0);
    else
        rPair = [];
    end
    report = local_report_stub(n, K, useLam);
    return;
end
if size(Vq,1) ~= n
    error('leigqNEWTON_cert_resPair:BadInput', 'v must have n rows (n=size(A,1)).');
end
if size(Vq,2) == 1 && K > 1
    error('leigqNEWTON_cert_resPair:BadInput', ...
        'For vector lambda (K>1), v must be n-by-K with matching columns.');
end
if size(Vq,2) ~= K
    if K == 1
        % allow n-by-any with K=1? No; ambiguous.
        error('leigqNEWTON_cert_resPair:BadInput', 'v must be n-by-1 when lambda is scalar.');
    else
        error('leigqNEWTON_cert_resPair:BadInput', 'Number of columns in v must match numel(lambda).');
    end
end

% ---------------- real embedding precompute ----------------
LA   = local_left_block(Aq);  % 4n x 4n
aF   = norm(LA,'fro');

% pack v (all columns)
vR = local_pack_qmat(Vq);     % 4n x K

% normalize each column (enforce ||v||=1)
colNorm = sqrt(sum(vR.^2, 1));
colNorm(colNorm == 0) = 1;
if opt.NormalizeV
    vR = vR ./ colNorm;
    colNorm = ones(1,K);
end
vUnit = local_unpack_qmat(vR);

% A*v for all columns (dominant shared cost)
AvR = LA * vR;                % 4n x K

resPairRaw  = zeros(K,1);
resPairNorm = zeros(K,1);

if nargout >= 4
    rPairR = zeros(4*n, K);
else
    rPairR = [];
end

for kk = 1:K
    lam4 = local_pack_qscalar(local_to_quat_scalar(lamUse(kk)));
    Llam = local_Lmat(lam4(1), lam4(2), lam4(3), lam4(4));

    lamvR = local_apply_blockdiag(Llam, vR(:,kk));
    rR    = AvR(:,kk) - lamvR;
    resPairRaw(kk) = norm(rR);

    den = (aF + norm(Llam,'fro')) * colNorm(kk) + 1;
    resPairNorm(kk) = resPairRaw(kk) / den;

    if nargout >= 4
        rPairR(:,kk) = rR;
    end
end

if opt.ResidualNormalized
    resPair = resPairNorm;
else
    resPair = resPairRaw;
end

if nargout >= 4
    rPair = local_unpack_qmat(rPairR);
else
    rPair = [];
end

% ---------------- report / comparison ----------------
report = local_report_stub(n, K, useLam);
report.normalizeV = opt.NormalizeV;
report.resPairRaw = resPairRaw;
report.resPairNorm = resPairNorm;
report.resPairReturned = resPair;

if ~isempty(opt.ResExpected)
    re = opt.ResExpected(:);
    if numel(re) ~= K
        report.resExpected = re;
        report.warning = "ResExpected length does not match numel(lambda).";
    else
        report.resExpected = re;
        report.diffRaw  = resPairRaw  - re;
        report.diffNorm = resPairNorm - re;
        report.maxAbsDiffRaw  = max(abs(report.diffRaw));
        report.maxAbsDiffNorm = max(abs(report.diffNorm));

        % which flavor matches better (useful when solver ran with ResidualNormalized=true)
        if report.maxAbsDiffNorm < report.maxAbsDiffRaw
            report.bestMatch = 'normalized';
            report.bestMaxAbsDiff = report.maxAbsDiffNorm;
        else
            report.bestMatch = 'raw';
            report.bestMaxAbsDiff = report.maxAbsDiffRaw;
        end
    end
end

end

% =====================================================================
%                         LOCAL HELPERS
% =====================================================================

function Aq = local_to_quat(x)
if isa(x,'quaternion')
    Aq = x;
elseif isnumeric(x)
    Aq = quaternion(x, zeros(size(x)), zeros(size(x)), zeros(size(x)));
else
    error('leigqNEWTON_cert_resPair:Type', ...
        'Unsupported type "%s". Input must be quaternion or numeric.', class(x));
end
end

function q = local_to_quat_scalar(x)
q = local_to_quat(x);
if ~isscalar(q)
    error('leigqNEWTON_cert_resPair:BadInput', 'lambda entries must be scalar quaternions.');
end
end

function y = local_pack_qmat(Vq)
% Pack n-by-K quaternion matrix into 4n-by-K real matrix (stack [w;x;y;z] per row).
[n,K] = size(Vq);
y = zeros(4*n, K);
[W,X,Y,Zz] = parts(Vq);
for k = 1:K
    for i = 1:n
        ii = 4*(i-1);
        y(ii+1,k) = W(i,k);
        y(ii+2,k) = X(i,k);
        y(ii+3,k) = Y(i,k);
        y(ii+4,k) = Zz(i,k);
    end
end
end

function Vq = local_unpack_qmat(Y)
% Unpack 4n-by-K real matrix into n-by-K quaternion matrix.
[m,K] = size(Y);
n = m/4;
W = zeros(n,K); X = W; Yy = W; Zz = W;
for k = 1:K
    for i = 1:n
        ii = 4*(i-1);
        W(i,k)  = Y(ii+1,k);
        X(i,k)  = Y(ii+2,k);
        Yy(i,k) = Y(ii+3,k);
        Zz(i,k) = Y(ii+4,k);
    end
end
Vq = quaternion(W,X,Yy,Zz);
end

function y4 = local_pack_qscalar(q)
[w,x,y,z] = parts(q);
y4 = [w; x; y; z];
end

function Z = local_left_block(Q)
% Real embedding rho(Q) for quaternion matrix Q (n-by-n): 4n-by4n block matrix.
n = size(Q,1);
Z = zeros(4*n, 4*n);
[W,X,Y,Zz] = parts(Q);
for i=1:n
    ii = 4*(i-1)+1:4*i;
    for j=1:n
        jj = 4*(j-1)+1:4*j;
        Z(ii,jj) = local_Lmat(W(i,j), X(i,j), Y(i,j), Zz(i,j));
    end
end
end

function L = local_Lmat(a,b,c,d)
% 4-by-4 real matrix representing LEFT multiplication by q=a+bi+cj+dk.
L = [ a, -b, -c, -d;
      b,  a, -d,  c;
      c,  d,  a, -b;
      d, -c,  b,  a ];
end

function y = local_apply_blockdiag(L4, x)
% Apply kron(I_n,L4) to a packed 4n-by-1 vector x.
n = numel(x)/4;
y = zeros(4*n,1);
for i = 1:n
    ii = 4*(i-1)+1:4*i;
    y(ii) = L4 * x(ii);
end
end

function report = local_report_stub(n, K, useLam)
report = struct();
report.n = n;
report.K = K;
report.lambdaTypeUsed = useLam;
report.warning = "";
end
