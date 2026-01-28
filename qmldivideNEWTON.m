function X = qmldivideNEWTON(A,B,tol)
%QMLDIVIDENEWTON  Solve quaternion linear systems A\B via real embedding (stand-alone).
%
%   X = qmldivideNEWTON(A,B)
%   X = qmldivideNEWTON(A,B,TOL)
%
% Solves the quaternion linear system
%       A * X = B
% in the least-squares / minimum-norm sense, analogously to MATLAB's backslash
% operator (\). This is useful because MATLAB's built-in quaternion class does
% not implement matrix mldivide for quaternion arrays.
%
% The solver converts the quaternion system to a real system using the standard
% 4-by-4 block embedding R(·) for LEFT multiplication:
%       R(A) * [X0;X1;X2;X3] = [B0;B1;B2;B3],
% where A = A0 + A1*i + A2*j + A3*k (componentwise), and similarly for X,B.
%
% Inputs
%   A, B : Quaternion or numeric arrays.
%          Supported types:
%            * quaternion arrays (MATLAB class quaternion),
%            * real numeric arrays,
%            * complex numeric arrays (embedded into the (1,i)-slice).
%          Complex embedding is: z = a + b*i  ->  a + b*i + 0*j + 0*k.
%
%   TOL  : (optional) nonnegative scalar.
%          If provided, the solution is computed via a pseudoinverse:
%               X = pinv(R(A),TOL) * R(B).
%          If omitted or empty, MATLAB's backslash is used on the real system:
%               R(A) \ R(B).
%
% Output
%   X    : Quaternion solution array.
%
% Examples (one-liners)
%   % Square system:
%   A = quaternion(randn(4),randn(4),randn(4),randn(4));
%   B = quaternion(randn(4,2),randn(4,2),randn(4,2),randn(4,2));
%   X = qmldivideNEWTON(A,B);
%
%   % Least-squares with tolerance:
%   X = qmldivideNEWTON(A,B,1e-12);
%
% See also: qmrdivideNEWTON, qmtimesNEWTON, pinv, mldivide, quaternion, parts

% Original author: M. Sebek (2024/2025)
% NEWTON toolbox refactor / packaging: 2026

if nargin < 2
    error('qmldivideNEWTON:NotEnoughInputs','Not enough input arguments.');
end

usePinv = (nargin >= 3) && ~isempty(tol);
if usePinv
    validateattributes(tol, {'double','single'}, {'real','scalar','nonnegative'}, mfilename, 'tol', 3);
end

A = local_any2quat(A);
B = local_any2quat(B);

[m,n]   = size(A);
[mB,~]  = size(B);
if mB ~= m
    error('qmldivideNEWTON:DimMismatch','A and B must have the same number of rows.');
end

RA = local_real_embed_left(A);
BR = local_stack_parts(B);

if usePinv
    XR = pinv(RA, double(tol)) * BR;
else
    XR = RA \ BR;
end

X = local_unstack_to_quat(XR, n);

end

% -------------------------------------------------------------------------
% Local helpers (kept inside this file to stay stand-alone)
% -------------------------------------------------------------------------

function Q = local_any2quat(X)
    if isa(X,'quaternion')
        Q = X;
        return
    end
    if ~isnumeric(X)
        error('qmldivideNEWTON:BadType','Input must be quaternion or numeric.');
    end
    if isreal(X)
        Z = zeros(size(X));
        Q = quaternion(double(X), double(Z), double(Z), double(Z));
    else
        Z = zeros(size(X));
        Q = quaternion(double(real(X)), double(imag(X)), double(Z), double(Z));
    end
end

function RA = local_real_embed_left(A)
    [a0,a1,a2,a3] = parts(A);
    RA = [ a0, -a1, -a2, -a3;
           a1,  a0, -a3,  a2;
           a2,  a3,  a0, -a1;
           a3, -a2,  a1,  a0];
end

function BR = local_stack_parts(B)
    [b0,b1,b2,b3] = parts(B);
    BR = [b0; b1; b2; b3];
end

function X = local_unstack_to_quat(XR, n)
    x0 = XR(1:n, :);
    x1 = XR(n+1:2*n, :);
    x2 = XR(2*n+1:3*n, :);
    x3 = XR(3*n+1:4*n, :);
    X  = quaternion(x0,x1,x2,x3);
end
