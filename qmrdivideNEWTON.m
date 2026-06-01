function X = qmrdivideNEWTON(A,B,tol)
%QMRDIVIDENEWTON  Solve quaternion right-division A/B via real embedding (stand-alone).
%
%   X = qmrdivideNEWTON(A,B)
%   X = qmrdivideNEWTON(A,B,TOL)
%
% Solves for X in the quaternion matrix equation
%       X * B = A
% in the least-squares / minimum-norm sense, analogously to MATLAB''s mrdivide
% operator (/). This is useful because MATLAB''s built-in quaternion class does
% not implement matrix mrdivide for quaternion arrays.
%
% Method
%   We use the identity for quaternion matrices (conjugate transpose):
%       (X*B)'' = B'' * X''
% which holds because quaternion conjugation reverses the multiplication order.
% Therefore we solve
%       B'' * Y = A'',   where Y = X'',
% using QMLDIVIDENEWTON, and return X = Y''.
%
%   IMPORTANT for the stand-alone bundle:
%   this implementation does NOT rely on MATLAB''s quaternion CTRANSPOSE (''),
%   because that overload may be unavailable or incomplete. The conjugate
%   transpose is formed explicitly from quaternion components.
%
% Inputs
%   A, B : Quaternion or numeric arrays.
%          Supported types:
%            * quaternion arrays (MATLAB class quaternion),
%            * real numeric arrays,
%            * complex numeric arrays (embedded into the (1,i)-slice).
%
%   TOL  : (optional) nonnegative scalar.
%          If provided and non-empty, it is forwarded to QMLDIVIDENEWTON.
%
% Output
%   X    : Quaternion solution array.
%
% Examples (one-liners)
%   A = quaternion(randn(3),randn(3),randn(3),randn(3));
%   B = quaternion(randn(3),randn(3),randn(3),randn(3));
%   X = qmrdivideNEWTON(A,B);
%
%   % With tolerance (pseudoinverse-based):
%   X = qmrdivideNEWTON(A,B,1e-12);
%
% See also: qmldivideNEWTON, qmtimesNEWTON, mrdivide, quaternion

% Original author: M. Sebek (2024/2025)
% NEWTON toolbox refactor / packaging: 2026

if nargin < 2
    error('qmrdivideNEWTON:NotEnoughInputs','Not enough input arguments.');
end

if nargin < 3
    tol = [];
end

A = local_any2quat(A);
B = local_any2quat(B);

At = local_ctranspose_quat(A);
Bt = local_ctranspose_quat(B);

if isempty(tol)
    Y = qmldivideNEWTON(Bt, At);
else
    Y = qmldivideNEWTON(Bt, At, tol);
end

X = local_ctranspose_quat(Y);

end

function Q = local_any2quat(X)
    if isa(X,'quaternion')
        Q = X;
        return
    end
    if ~isnumeric(X)
        error('qmrdivideNEWTON:BadType','Input must be quaternion or numeric.');
    end
    if isreal(X)
        Z = zeros(size(X));
        Q = quaternion(double(X), double(Z), double(Z), double(Z));
    else
        Z = zeros(size(X));
        Q = quaternion(double(real(X)), double(imag(X)), double(Z), double(Z));
    end
end

function Qt = local_ctranspose_quat(Q)
% Quaternion conjugate transpose formed explicitly from components.
    [w,x,y,z] = parts(Q);
    Qt = quaternion(w.', -x.', -y.', -z.');
end
