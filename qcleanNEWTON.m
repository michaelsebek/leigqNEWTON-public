function A = qcleanNEWTON(A, tol)
%QCLEANNEWTON  Component-wise cleanup (zeroing) of quaternion arrays.
%
%   Y = qcleanNEWTON(X)
%   Y = qcleanNEWTON(X, TOL)
%
% Cleans ("denoises") a quaternion array by setting small real and imaginary
% components to zero **component-wise**.
%
% For each quaternion entry
%     x = w + x*i + y*j + z*k,
% the function zeroes any component whose magnitude is smaller than TOL:
%     w := 0 if |w| < TOL,  x := 0 if |x| < TOL,  etc.
%
% This operation improves readability (and may promote sparsity), but reduces
% numerical precision by at most the chosen tolerance.
%
% Inputs
%   X   : quaternion array (any size). Numeric inputs are converted via
%         quaternion(X) (purely real quaternions).
%   TOL : scalar double in [0,1]. If omitted or empty, qcleanNEWTON uses
%         the global tolerance PGLOBAL.ZEROING when available; otherwise a
%         built-in default of 1e-12 is used.
%
% Output
%   Y : quaternion array, same size as X.
%
% Notes
%   * Cleaning is performed component-wise. An entry may retain some non-zero
%     components while others are cleared, so the entry itself is not
%     necessarily set to 0.
%   * If X represents rotations (unit quaternions), component-wise cleaning
%     will generally destroy unit length; re-normalize if needed.
%
% Examples (copy/paste)
%   q  = quaternion(1, 1e-14, -2e-13, 3);
%   qcleanNEWTON(q)                 % uses default tolerance
%
%   A  = quaternion(randn(3), 1e-14*randn(3), randn(3), 0*randn(3));
%   A2 = qcleanNEWTON(A, 1e-12);    % explicit tolerance
%
%   % Numeric inputs are converted to purely real quaternions:
%   qcleanNEWTON(eye(3), 1e-12)
%
% Legacy mapping
%   qclean  ->  qcleanNEWTON
%
% See also: quaternion, parts, qroundNEWTON.

% ----------------------- tolerance -----------------------
if nargin == 0
    error('Not enough input arguments.');
end

if nargin < 2 || isempty(tol)
    % Mirror the legacy behaviour when Polynomial Toolbox globals exist,
    % but remain fully stand-alone.
    tolDefault = 1e-12;
    try
        % If present, use the same global as the legacy qclean.
        global PGLOBAL;
        if ~isempty(PGLOBAL) && isstruct(PGLOBAL) && isfield(PGLOBAL,'ZEROING')
            tol = PGLOBAL.ZEROING;
        else
            tol = tolDefault;
        end
    catch
        tol = tolDefault;
    end
else
    if ~isa(tol,'double') || numel(tol) ~= 1 || ~isreal(tol) || tol < 0 || tol > 1
        error('Invalid tolerance.');
    end
end

% ----------------------- input cast -----------------------
try
    A = quaternion(A);
catch
    error('Invalid 1st argument.');
end

% ----------------------- cleanup -----------------------
[wA,xA,yA,zA] = parts(A);

wA(abs(wA) < tol) = 0;
xA(abs(xA) < tol) = 0;
yA(abs(yA) < tol) = 0;
zA(abs(zA) < tol) = 0;

A = quaternion(wA,xA,yA,zA);

end
