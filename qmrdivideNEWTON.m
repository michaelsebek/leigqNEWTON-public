function X = qmrdivideNEWTON(A,B,tol)
%QMRDIVIDENEWTON  Solve quaternion right-division A/B via real embedding (stand-alone).
%
%   X = qmrdivideNEWTON(A,B)
%   X = qmrdivideNEWTON(A,B,TOL)
%
% Solves for X in the quaternion matrix equation
%       X * B = A
% in the least-squares / minimum-norm sense, analogously to MATLAB's mrdivide
% operator (/). This is useful because MATLAB's built-in quaternion class does
% not implement matrix mrdivide for quaternion arrays.
%
% Method
%   We use the identity for quaternion matrices (conjugate transpose):
%       (X*B)' = B' * X'
% which holds because quaternion conjugation reverses the multiplication order.
% Therefore we solve
%       B' * Y = A',   where Y = X',
% using QMLDIVIDENEWTON, and return X = Y'.
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

if isempty(tol)
    X = (qmldivideNEWTON(B', A'))';
else
    X = (qmldivideNEWTON(B', A', tol))';
end

end
