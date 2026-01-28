function C = qmtimesNEWTON(A,B,tol)
%QMTIMESNEWTON  Quaternion matrix product A*B via real embedding (stand-alone).
%
%   C = qmtimesNEWTON(A,B)
%   C = qmtimesNEWTON(A,B,TOL)
%
% Computes the quaternion matrix product C = A*B even though MATLAB's built-in
% quaternion class does not implement matrix multiplication (MTIMES) for
% quaternion arrays.
%
% Inputs
%   A, B : Matrices to multiply.
%          Supported types:
%            * quaternion arrays (MATLAB class quaternion),
%            * real numeric arrays,
%            * complex numeric arrays (embedded into the (1,i)-slice).
%          Complex embedding is: z = a + b*i  ->  a + b*i + 0*j + 0*k.
%
%   TOL  : (optional) nonnegative scalar.
%          If TOL>0, a magnitude-based "zeroing" post-pass is applied to the
%          output: entries with |C_ij| < TOL*scale are set to exact zero, where
%             scale = minNonzeroNorm(A) * minNonzeroNorm(B).
%          If omitted or empty, no zeroing is applied.
%
% Output
%   C    : Quaternion matrix product (quaternion array).
%
% Notes
%   * The core computation uses a standard 4x4 real block embedding for left
%     multiplication, typically realized as a single BLAS-backed real matrix
%     multiplication. A low-memory fallback computes the component formula.
%   * This function is stand-alone (no dependencies on your other toolboxes).
%   * This generates the algebraic product in H, not an elementwise product.
%
% Examples (one-liners)
%   A = quaternion(randn(3),randn(3),randn(3),randn(3));
%   B = quaternion(randn(3),randn(3),randn(3),randn(3));
%   C = qmtimesNEWTON(A,B);
%
%   % Complex embedding (i-slice):
%   C = qmtimesNEWTON(eye(3), 1+2i);
%
%   % Optional zeroing:
%   C = qmtimesNEWTON(A,B,1e-12);
%
% See also: qmldivideNEWTON, qmrdivideNEWTON, quaternion, parts

% Original author: M. Sebek (2024/2025)
% NEWTON toolbox refactor / packaging: 2026

if nargin < 2
    error('qmtimesNEWTON:NotEnoughInputs','Not enough input arguments.');
end
if nargin < 3 || isempty(tol)
    tol = 0;
else
    validateattributes(tol, {'double','single'}, {'real','scalar','nonnegative'}, mfilename, 'tol', 3);
end

A = local_any2quat(A);
B = local_any2quat(B);

[m,n] = size(A);
[p,q] = size(B);

% Scalar cases (preserve order!)
if isscalar(A)
    C = local_times_scalar_left(A,B);
    C = local_zero_if_requested(C, tol, local_min_nonzero_norm(A)*local_min_nonzero_norm(B));
    return
elseif isscalar(B)
    C = local_times_scalar_right(A,B);
    C = local_zero_if_requested(C, tol, local_min_nonzero_norm(A)*local_min_nonzero_norm(B));
    return
end

if n ~= p
    error('qmtimesNEWTON:InnerDimensions','Inner matrix dimensions must agree.');
end

% Split parts
[a0,a1,a2,a3] = parts(A);
[b0,b1,b2,b3] = parts(B);

% Heuristic: build the 4x4 real block embedding unless it would be too large.
% (keeps memory use bounded for large matrices)
bytesRA = 8 * 16 * double(m) * double(n);  % doubles in RA times 8 bytes
useBlock = (bytesRA <= 1e8);               % ~100 MB threshold

if useBlock
    % Real block embedding for left multiplication:
    % vec([C0;C1;C2;C3]) = R(A) * vec([B0;B1;B2;B3])
    RA = [ a0, -a1, -a2, -a3;
           a1,  a0, -a3,  a2;
           a2,  a3,  a0, -a1;
           a3, -a2,  a1,  a0];

    Bst = [b0; b1; b2; b3];
    Cst = RA * Bst;

    c0 = Cst(1:m, :);
    c1 = Cst(m+1:2*m, :);
    c2 = Cst(2*m+1:3*m, :);
    c3 = Cst(3*m+1:4*m, :);
else
    % Component formula (less memory, more BLAS calls)
    c0 =  a0*b0 - a1*b1 - a2*b2 - a3*b3;
    c1 =  a0*b1 + a1*b0 + a2*b3 - a3*b2;
    c2 =  a0*b2 - a1*b3 + a2*b0 + a3*b1;
    c3 =  a0*b3 + a1*b2 - a2*b1 + a3*b0;
end

C = quaternion(c0,c1,c2,c3);

% Optional zeroing
scale = local_min_nonzero_norm(A) * local_min_nonzero_norm(B);
C = local_zero_if_requested(C, tol, scale);

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
        error('qmtimesNEWTON:BadType','Input must be quaternion or numeric.');
    end
    if isreal(X)
        Z = zeros(size(X));
        Q = quaternion(double(X), double(Z), double(Z), double(Z));
    else
        Z = zeros(size(X));
        Q = quaternion(double(real(X)), double(imag(X)), double(Z), double(Z));
    end
end

function C = local_times_scalar_left(s,B)
    % C = s * B   (s is 1x1 quaternion)
    [s0,s1,s2,s3] = parts(s);
    [b0,b1,b2,b3] = parts(B);
    c0 = s0.*b0 - s1.*b1 - s2.*b2 - s3.*b3;
    c1 = s0.*b1 + s1.*b0 + s2.*b3 - s3.*b2;
    c2 = s0.*b2 - s1.*b3 + s2.*b0 + s3.*b1;
    c3 = s0.*b3 + s1.*b2 - s2.*b1 + s3.*b0;
    C = quaternion(c0,c1,c2,c3);
end

function C = local_times_scalar_right(A,s)
    % C = A * s   (s is 1x1 quaternion)
    [a0,a1,a2,a3] = parts(A);
    [s0,s1,s2,s3] = parts(s);
    c0 = a0.*s0 - a1.*s1 - a2.*s2 - a3.*s3;
    c1 = a0.*s1 + a1.*s0 + a2.*s3 - a3.*s2;
    c2 = a0.*s2 - a1.*s3 + a2.*s0 + a3.*s1;
    c3 = a0.*s3 + a1.*s2 - a2.*s1 + a3.*s0;
    C = quaternion(c0,c1,c2,c3);
end

function mn = local_min_nonzero_norm(Q)
    [q0,q1,q2,q3] = parts(Q);
    mag = sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    mag = mag(:);
    mag = mag(mag>0);
    if isempty(mag)
        mn = 0;
    else
        mn = min(mag);
    end
end

function Q = local_zero_if_requested(Q,tol,scale)
    if tol <= 0 || scale <= 0
        return
    end
    thr = double(tol) * double(scale);
    [q0,q1,q2,q3] = parts(Q);
    mag = sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
    mask = (mag < thr);
    if any(mask(:))
        q0(mask) = 0; q1(mask) = 0; q2(mask) = 0; q3(mask) = 0;
        Q = quaternion(q0,q1,q2,q3);
    end
end
