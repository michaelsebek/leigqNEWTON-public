function Y = qroundNEWTON(X, N, mode)
%QROUNDNEWTON  Round quaternion elements to a specified number of digits.
%
%   Y = qroundNEWTON(X)
%   Y = qroundNEWTON(X, N)
%   Y = qroundNEWTON(X, N, "significant")
%
% Rounds quaternion values **component-wise**, mirroring MATLAB's ROUND for
% real/complex arrays.
%
% For a quaternion
%     X = w + x*i + y*j + z*k,
% qroundNEWTON returns (component-wise)
%     round(w, ...) + round(x, ...)*i + round(y, ...)*j + round(z, ...)*k.
%
% Inputs
%   X    : quaternion array (any size) OR numeric array.
%          - If X is quaternion, each component is rounded separately.
%          - If X is numeric, qroundNEWTON calls MATLAB ROUND(X, ...).
%
%   N    : (optional) integer scalar specifying the rounding precision.
%          - If omitted or empty, N = 0 (round to nearest integer).
%          - If N > 0, round to N decimal digits.
%          - If N < 0, round to powers of ten (e.g., N = -1 rounds to tens).
%
%   mode : (optional) "significant" (string/char). If provided, rounds to
%          N significant digits (component-wise for quaternions).
%
% Output
%   Y : same type/shape as X (quaternion if X is quaternion; numeric otherwise).
%
% Notes
%   * This is a component-wise operation; it does NOT preserve unit quaternions.
%     If X represents rotations, re-normalize as needed.
%
% Examples (copy/paste)
%   q  = quaternion( 1.2345, -3.7654, 0.0123, 4.9876 );
%   qroundNEWTON(q)                    % integer rounding
%   qroundNEWTON(q, 2)                 % 2 decimal digits
%   qroundNEWTON(q, -1)                % to tens
%   qroundNEWTON(q, 3, "significant")  % 3 significant digits
%
%   % Numeric inputs fall back to MATLAB's round:
%   qroundNEWTON(pi, 4)
%
% Legacy mapping
%   qround  ->  qroundNEWTON
%
% See also: round, quaternion, parts, qcleanNEWTON.

% ---- defaults & argument normalization ----
if nargin < 2 || isempty(N)
    N = 0;
end

if nargin < 3
    mode = [];
end

useSignificant = false;
if ~isempty(mode)
    if isstring(mode) || ischar(mode)
        modeStr = string(mode);
        if isscalar(modeStr) && lower(modeStr) == "significant"
            useSignificant = true;
        else
            error('qroundNEWTON:Mode', 'Unsupported mode. Only "significant" is supported.');
        end
    else
        error('qroundNEWTON:ModeType', 'Mode must be a string or character vector.');
    end
end

% ---- quaternion vs numeric ----
if isa(X,'quaternion')
    [w,x,y,z] = parts(X);

    if useSignificant
        Y = quaternion( round(w, N, "significant"), ...
                        round(x, N, "significant"), ...
                        round(y, N, "significant"), ...
                        round(z, N, "significant") );
    else
        if N == 0
            Y = quaternion( round(w), round(x), round(y), round(z) );
        else
            Y = quaternion( round(w, N), round(x, N), round(y, N), round(z, N) );
        end
    end
else
    % numeric fallback behaves like MATLAB ROUND
    if useSignificant
        Y = round(X, N, "significant");
    else
        if N == 0
            Y = round(X);
        else
            Y = round(X, N);
        end
    end
end

end
