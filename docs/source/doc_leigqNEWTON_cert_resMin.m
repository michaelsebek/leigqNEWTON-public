%% leigqNEWTON_cert_resMin — eigenvalue-only residual certificate
% Compute
%   resMin(A,lambda) = min_{||x||=1} ||A*x - lambda*x||
% and (optionally) the minimizer vector.
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

%% Setup (requirements + path)
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

if exist('leigqNEWTON_cert_resMin','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'leigqNEWTON_cert_resMin.m'),'file')
            addpath(rootGuess);
        end
    end
end

if exist('leigqNEWTON_cert_resMin','file') ~= 2
    error('leigqNEWTON_cert_resMin not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |r = leigqNEWTON_cert_resMin(A, lambda)|
% * |[r,x] = leigqNEWTON_cert_resMin(A, lambda)|
% * |... = leigqNEWTON_cert_resMin(..., 'ResidualNormalized',true/false)|

%% Example 1: Huang–So Example 2.5 (2×2, known eigenvalues)
q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);

A = [ q0, q1 + qi;
      q1 - qi, q0 ];

lamTrue = [ quaternion(sqrt(2),0,0,0);
            quaternion(-sqrt(2),0,0,0) ];

[r,x] = leigqNEWTON_cert_resMin(A, lamTrue, 'ResidualNormalized',false);
disp(r);

%% Example 2: compute best vectors for a candidate set (Pan–Ng 4×4)
q = @(w,x,y,z) quaternion(w,x,y,z);
a = q(-2, 1, 1, 4);
b = q( 2, 4, 1, 1);
c = q( 1, 3, 2, 2);
d = q(-1, 2, 2, 3);
A = [a b c d;
     d a b c;
     c d a b;
     b c d a];

lambda = leigqNEWTON(A,'SolveProfile','reliable','Seed',1);
[resMinAbs, X] = leigqNEWTON_cert_resMin(A, lambda, 'ResidualNormalized',false);
size(X)
median(resMinAbs), max(resMinAbs)

%% See also
% leigqNEWTON_cert_resPair, leigqNEWTON_refine_batch, checkNEWTON
