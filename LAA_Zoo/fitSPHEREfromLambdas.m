function [sp, meta] = fitSPHEREfromLambdas(lambda, varargin)
%fitSPHEREfromLambdas  Fit a 2-sphere hypothesis to quaternion samples (diagnostic utility).
%
%   [sp, meta] = fitSPHEREfromLambdas(lambda, 'Name',Value,...)
%
% Fits a sphere to a set of quaternion samples lambda (Nx1 quaternion), intended for
% "spherical components" diagnostics in the left spectrum.
%
% The fit is performed in a best-fit 3D affine subspace of R^4 obtained by PCA/SVD.
% A least-squares sphere fit is then solved in R^3 and lifted back to R^4.
%
% Inputs
%   lambda : vector of quaternions (candidates, typically sampled/refined eigenvalues)
%
% Name-value options
%   'TolInlier'     : inlier threshold on | ||lam-c|| - r |  (default: 1e-8)
%   'SnapInteger'   : snap fitted center/radius close to integers if very near (default: false)
%   'SnapTol'       : snapping tolerance (default: 1e-10)
%   'Verbose'       : true/false (default: true)
%
% Outputs
%   sp.center   : fitted center (quaternion)
%   sp.radius   : fitted radius (double)
%   sp.inliers  : logical inliers mask (Nx1 logical)
%   sp.dev      : deviations di = ||lam-c|| - r  (Nx1 double)
%   sp.basis    : 4x3 basis of fitted subspace (columns)
%   sp.mean4    : mean (1x4) used for centering
%
% Notes
%   - This is a diagnostic helper intended for the Zoo examples. It is not part of the
%     leigqNEWTON core solver API.
%   - The private "experiment driver" toolbox may contain a richer implementation; this
%     public version is self-contained and uses only MATLAB quaternion.
%
% See also: leigqNEWTON_sphere_detect, leigqNEWTON_sphere_refine, quaternion

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'TolInlier',1e-8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'SnapInteger',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'SnapTol',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
opt = p.Results;

meta = struct();
meta.opt = opt;

% Ensure quaternion is available
try
    quaternion(0,0,0,0);
catch ME
    error('MATLAB quaternion class not available (requires Aerospace Toolbox). (%s)', ME.message);
end

lambda = lambda(:);
n = numel(lambda);
if n < 4
    error('Need at least 4 samples to fit a 2-sphere robustly (got %d).', n);
end

% Convert quaternion samples to points in R^4
X = local_q2r4(lambda); % n x 4
mu = mean(X,1);
X0 = X - mu;

% Best-fit 3D subspace in R^4 (PCA via SVD)
[~,S,V] = svd(X0, 'econ'); %#ok<ASGLU>
basis = V(:,1:3);          % 4 x 3
Y = X0 * basis;            % n x 3 (coordinates)

% Least-squares sphere fit in R^3:
% Solve [-2Y, 1] * [c; alpha] = -||Y||^2, with alpha = c^T c - r^2
b = -sum(Y.^2,2);
A = [-2*Y, ones(n,1)];
u = A \ b;                  % 4x1
c3 = u(1:3);
alpha = u(4);
r2 = max(dot(c3,c3) - alpha, 0);
r = sqrt(r2);

% Lift center back to R^4
c4 = mu(:) + basis*c3(:);   % 4x1
center = quaternion(c4(1), c4(2), c4(3), c4(4));

% Deviations and inliers (in original R^4 norm)
d = local_dev_r4(X, c4(:), r);
inliers = abs(d) <= opt.TolInlier;

% Optional snapping
if opt.SnapInteger
    [center, r] = local_snap_integer(center, r, opt.SnapTol);
    % recompute deviations with snapped parameters
    c4s = local_q2r4(center);
    d = local_dev_r4(X, c4s(:), r);
    inliers = abs(d) <= opt.TolInlier;
end

sp = struct();
sp.center  = center;
sp.radius  = r;
sp.inliers = inliers;
sp.dev     = d;
sp.basis   = basis;
sp.mean4   = mu;

meta.singular_values = diag(S);
meta.n = n;
meta.inlier_count = nnz(inliers);
meta.dev_median = median(abs(d(inliers)));
meta.dev_max    = max(abs(d(inliers)));

if opt.Verbose
    fprintf('fitSPHEREfromLambdas: n=%d, inliers=%d, r=%.15g\n', n, nnz(inliers), r);
    fprintf('  median |dev| (inliers) = %.3e, max |dev| (inliers) = %.3e\n', meta.dev_median, meta.dev_max);
end
end

% -------------------------------------------------------------------------
function X = local_q2r4(q)
% Return n x 4 matrix [a b c d]
try
    [a,b,c,d] = parts(q);
    X = [a(:), b(:), c(:), d(:)];
catch
    % fallback for older releases: compact(q) -> n x 4
    C = compact(q);
    X = C(:,1:4);
end
end

function dev = local_dev_r4(X, c4, r)
% dev_i = ||X_i - c|| - r
D = X - reshape(c4(:).',1,4);
dist = sqrt(sum(D.^2,2));
dev = dist - r;
end

function [c, r] = local_snap_integer(c, r, tol)
% Snap quaternion center and radius to nearest integer if within tol
X = local_q2r4(c); % 1x4
Xi = round(X);
if max(abs(X - Xi)) <= tol
    c = quaternion(Xi(1), Xi(2), Xi(3), Xi(4));
end
ri = round(r);
if abs(r - ri) <= tol
    r = ri;
end
end
