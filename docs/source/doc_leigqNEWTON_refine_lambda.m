%% leigqNEWTON_refine_lambda — refine a single candidate eigenvalue
% Improves a candidate lambda by minimizing the eigenvalue-only certificate resMin(A,lambda).
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

%% Setup
hasQuat = true;
try
    quaternion(0,0,0,0);
catch
    hasQuat = false;
end
if ~hasQuat
    disp('This toolbox requires MATLAB''s built-in quaternion class (quaternion(w,x,y,z)).');
    disp('Examples in this page are skipped.');
    return;
end

if exist('leigqNEWTON_refine_lambda','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'leigqNEWTON_refine_lambda.m'),'file')
            addpath(rootGuess);
        end
    end
end

if exist('leigqNEWTON_refine_lambda','file') ~= 2
    error('leigqNEWTON_refine_lambda not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |[lam1, v1, out] = leigqNEWTON_refine_lambda(A, lambda0)|

%% Example
% Start from a solver hit and refine lambda via resMin.
q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);
qj = quaternion(0,0,1,0);
A  = [ q0, qi;
       qj, q1 ];

lambda = leigqNEWTON(A,'SolveProfile','default','Seed',1);
lam0   = lambda(1);

r0 = leigqNEWTON_cert_resMin(A, lam0);
[lam1, v1, out] = leigqNEWTON_refine_lambda(A, lam0);
r1 = leigqNEWTON_cert_resMin(A, lam1);
[r0, r1]

%% See also
% leigqNEWTON_refine_batch, leigqNEWTON_cert_resMin
