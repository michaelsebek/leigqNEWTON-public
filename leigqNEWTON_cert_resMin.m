function [resMin, xMin, rMin, info] = leigqNEWTON_cert_resMin(A, lambda, varargin)
%LEIGQNEWTON_CERT_RESMIN  Minimum-residual certificate for a quaternion LEFT eigenvalue candidate.
%
%   resMin = leigqNEWTON_cert_resMin(A, lambda)
%   [resMin, xMin] = leigqNEWTON_cert_resMin(A, lambda)
%   [resMin, xMin, rMin, info] = leigqNEWTON_cert_resMin(A, lambda, Name,Value,...)
%
% Defines the certificate-like quantity
%   resMin(A,lambda) := min_{||x||_2=1} ||A*x - lambda*x||_2,
% computed via a real/complex embedding and a smallest-singular-value routine.
%
% Inputs
%   A      : n-by-n quaternion (or numeric; complex values use the (1,i) slice)
%   lambda : scalar quaternion or K-by-1 quaternion array
%
% Name,Value options (selected)
%   'Method'             : 'svd'|'svds'|'auto' (default 'svd')
%   'MaxIter'            : iterations for iterative methods (if used)
%   'Tol'                : tolerance for iterative methods (if used)
%   'ResidualNormalized' : if true, return scale-invariant resMin (default false)
%       denominator = norm(A,2)+abs(lambda), no +1 floor; etaMin(0,0)=0.
%       This is the same normalization as leigqNEWTON for a unit vector.
%       Core revision: LAA-R1-relative-2026-10-03.
%   'UseLambda'          : 'lambda'|'lambdac' (default 'lambda')  % choose which lambda vector to use
%   'LambdaC'            : cleaned lambdas (required if UseLambda='lambdac')
%
% Outputs
%   resMin : K-by-1 double   (certificate value for each lambda)
%   xMin   : n-by-K quaternion (minimizers; columns; returned if requested)
%   rMin   : n-by-K quaternion (defects of the returned unit minimizers)
%   info   : struct          (method statistics; timings; etc.)
%
% Useful one-liners
%   r = leigqNEWTON_cert_resMin(A, lam(1));
%   [r,x] = leigqNEWTON_cert_resMin(A, lam, 'ResidualNormalized',true);
%
% See also: leigqNEWTON_cert_resPair, leigqNEWTON_refine_lambda, leigqNEWTON

if ~exist('quaternion','class')
    error('leigqNEWTON_cert_resMin:NoQuaternionClass', ...
        'MATLAB class "quaternion" not found. Requires Aerospace Toolbox.');
end

% ---------------- parse options ----------------
opt.Method             = 'svd';
opt.Tol                = 1e-12;
opt.MaxIter            = 1000;
opt.ResidualNormalized = false;
opt.LambdaC            = [];
opt.UseLambda          = 'lambda';
opt.ResPair            = [];

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('leigqNEWTON_cert_resMin:BadNV', 'Name/value inputs must come in pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        if ~ischar(name) && ~isstring(name)
            error('leigqNEWTON_cert_resMin:BadNV', 'Option names must be char or string.');
        end
        key = lower(char(name));
        switch key
            case {'method'}
                opt.Method = lower(char(val));
            case {'tol','tolerance'}
                opt.Tol = double(val);
            case {'maxiter','maxit','iterations'}
                opt.MaxIter = double(val);
            case {'residualnormalized','normalized','scalefree','scaleinvariant'}
                opt.ResidualNormalized = logical(val);
            case {'lambdac','lamc','lambdaclean','cleanlambda','cleanedlambda'}
                opt.LambdaC = val;
            case {'uselambda','use','whichlambda','lamtype'}
                opt.UseLambda = lower(char(val));
            case {'respair','res','resvec','resvector','resexpected','resfromsolver','ressolver'}
                opt.ResPair = double(val);
            otherwise
                % ignore unknown options (consistent with leigqNEWTON)
        end
    end
end

% ---------------- type normalization ----------------
Aq = local_to_quat(A);
[n,m] = size(Aq);
if ndims(Aq) ~= 2 || n ~= m || n == 0
    error('leigqNEWTON_cert_resMin:BadInput', 'A must be square.');
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
            error('leigqNEWTON_cert_resMin:MissingLambdaC', ...
                'UseLambda="lambdac" requested but LambdaC was not provided.');
        end
        lamUse = local_to_quat(opt.LambdaC);
        lamUse = lamUse(:);
    otherwise
        lamUse = lamIn;
end

K = numel(lamUse);

% ---- precompute rho(A) once ----
LA = local_left_block(Aq);
if any(~isfinite(LA(:)))
    error('leigqNEWTON_cert_resMin:NonfiniteInput','A must be finite.');
end
a2 = norm(LA,2);
if ~isfinite(a2)
    error('leigq:ScaleOverflow','The matrix norm is not finite; rescale A.');
end
In = eye(n);

wantX = (nargout >= 2);
wantR = (nargout >= 3);

resMinRaw  = zeros(K,1);
resMinNorm = zeros(K,1);
denominators = zeros(K,1);

% ---- SAFE preallocation for quaternion arrays (avoids xMin(n,K)=quaternion) ----
if wantX
    Z = zeros(n,K);
    xMin = quaternion(Z,Z,Z,Z);
else
    xMin = [];
end
if wantR
    Z = zeros(n,K);
    rMin = quaternion(Z,Z,Z,Z);
else
    rMin = [];
end

methodUsed = strings(K,1);

for kk = 1:K
    lam4 = local_pack_qscalar(local_to_quat_scalar(lamUse(kk)));
    if any(~isfinite(lam4))
        error('leigqNEWTON_cert_resMin:NonfiniteInput','lambda must be finite.');
    end
    Llam = local_Lmat(lam4(1), lam4(2), lam4(3), lam4(4));
    M    = LA - kron(In, Llam);  % rho(A - lambda I)

    % decide method
    meth = opt.Method;
    if strcmp(meth,'auto')
        if size(M,1) >= 512
            meth = 'svds';
        else
            meth = 'svd';
        end
    end

    if strcmp(meth,'svds')
        try
            optsS.tol    = opt.Tol;
            optsS.maxit  = opt.MaxIter;
            [~,s,V,flag] = svds(M, 1, 'smallest', optsS);
            if flag ~= 0 || isempty(s) || ~isfinite(s(1,1)) || ...
                    any(~isfinite(V(:))) || norm(V(:,1)) == 0
                error('leigq:SVDSNotConverged','SVDS did not return a finite converged triplet.');
            end
            if wantX, vR = V(:,1); else, vR = []; end
            sigmaMin = s(1,1);
        catch
            meth = 'svd';
            [~,S,V] = svd(M,'econ');
            sigmaMin = S(end,end);
            if wantX
                vR = V(:,end);
            else
                vR = [];
            end
        end
    elseif strcmp(meth,'svd')
        [~,S,V] = svd(M,'econ');
        sigmaMin = S(end,end);
        if wantX
            vR = V(:,end);
        else
            vR = [];
        end
    else
        error('leigqNEWTON_cert_resMin:BadMethod', ...
            'Unknown Method="%s". Use svd|svds|auto.', opt.Method);
    end

    resMinRaw(kk) = sigmaMin;
    methodUsed(kk) = string(meth);

    % matches leigqNEWTON's denominator with ||x||=1
    [resMinNorm(kk), denominators(kk)] = ...
        leigqNEWTON_relres(resMinRaw(kk), a2, norm(lam4), 1);

    if wantX
        if isempty(vR)
            [~,~,V] = svd(M,'econ');
            vR = V(:,end);
        end
        nv = norm(vR);
        if nv == 0
            vR = randn(4*n,1);
            vR = vR / norm(vR);
        else
            vR = vR / nv;
        end
        xMin(:,kk) = local_unpack_qvec(vR);
        if wantR
            rMin(:,kk) = local_unpack_qvec(M*vR);
        end
    end
end

if opt.ResidualNormalized
    resMin = resMinNorm;
else
    resMin = resMinRaw;
end

if nargout >= 4
    info = struct();
    info.n = n;
    info.K = K;
    info.lambdaTypeUsed = useLam;
    info.lambdaUsed = lamUse;
    info.embeddingSize = 4*n;
    info.methodUsed = methodUsed;
    info.resMinRaw = resMinRaw;
    info.resMinNorm = resMinNorm;
    info.resMinReturned = resMin;
    info.normA2 = a2;
    info.denominators = denominators;
    info.normalization = 'norm(A,2)+abs(lambda)';

    if ~isempty(opt.ResPair)
        rp = opt.ResPair(:);
        info.resPairProvided = rp;
        if numel(rp) ~= K
            info.boundWarning = 'ResPair length does not match numel(lambda).';
        else
            info.boundCheckRaw.maxViolation  = max(resMinRaw  - rp);
            info.boundCheckRaw.nViolations  = sum((resMinRaw  - rp) > 0);
            info.boundCheckRaw.ok = (info.boundCheckRaw.maxViolation <= 0);

            info.boundCheckNorm.maxViolation = max(resMinNorm - rp);
            info.boundCheckNorm.nViolations = sum((resMinNorm - rp) > 0);
            info.boundCheckNorm.ok = (info.boundCheckNorm.maxViolation <= 0);
        end
    end
else
    info = [];
end

end

% =====================================================================
%                         LOCAL HELPERS
% =====================================================================

function Aq = local_to_quat(x)
if isa(x,'quaternion')
    Aq = x;
elseif isnumeric(x)
    Aq = quaternion(double(real(x)), double(imag(x)), zeros(size(x)), zeros(size(x)));
else
    error('leigqNEWTON_cert_resMin:Type', ...
        'Unsupported type "%s". Input must be quaternion or numeric.', class(x));
end
end

function q = local_to_quat_scalar(x)
q = local_to_quat(x);
if ~isscalar(q)
    error('leigqNEWTON_cert_resMin:BadInput', 'lambda entries must be scalar quaternions.');
end
end

function y4 = local_pack_qscalar(q)
[w,x,y,z] = parts(q);
y4 = [w; x; y; z];
end

function xq = local_unpack_qvec(y)
% y is 4n-by-1 real; unpack to n-by-1 quaternion vector
Y = reshape(y, 4, []);
xq = quaternion(Y(1,:).', Y(2,:).', Y(3,:).', Y(4,:).');
end

function Z = local_left_block(Q)
% Real embedding rho(Q) for quaternion matrix Q (n-by-n): 4n-by-4n block matrix.
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
