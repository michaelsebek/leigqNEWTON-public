%% leigqNEWTON_refine_batch — refine + certify a list of candidates
% Refine a list of candidate left eigenvalues and compute residual certificates.
% This is the recommended post-processing after calling |leigqNEWTON|.

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

if exist('leigqNEWTON_refine_batch','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'leigqNEWTON_refine_batch.m'),'file')
            addpath(rootGuess);
        end
    end
end

if exist('leigqNEWTON_refine_batch','file') ~= 2
    error('leigqNEWTON_refine_batch not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |[lamR,VR,cert] = leigqNEWTON_refine_batch(A, lam)|
% * |... = leigqNEWTON_refine_batch(..., 'DoPolish',true,'TargetResMin',1e-14)|

%% Notes
% * |cert| is a struct; common fields include |resMin| (lambda-only) and |resPair| (pair residual).
% * If |DoPolish| is true, each accepted candidate is polished to near machine precision.

%% Example 1: refine solver output (2×2)
q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);
qj = quaternion(0,0,1,0);
A  = [ q0, qi;
       qj, q1 ];

lambda1 = leigqNEWTON(A,'Seed',1,'SolveProfile','default');
[lamR1,VR1,cert1] = leigqNEWTON_refine_batch(A, lambda1, 'DoPolish',true);
median(cert1.resMin), max(cert1.resMin)

%% Example 2: Pan–Ng 4×4 (start from solver hits)
q = @(w,x,y,z) quaternion(w,x,y,z);
a = q(-2, 1, 1, 4);
b = q( 2, 4, 1, 1);
c = q( 1, 3, 2, 2);
d = q(-1, 2, 2, 3);
A = [a b c d;
     d a b c;
     c d a b;
     b c d a];

lambda2 = leigqNEWTON(A,'Seed',1,'SolveProfile','reliable','K',2*size(A,1));
[lamR2,VR2,cert2] = leigqNEWTON_refine_batch(A, lambda2, 'DoPolish',true,'TargetResMin',1e-14);
median(cert2.resMin), max(cert2.resMin)

%% See also
% leigqNEWTON, leigqNEWTON_refine_polish, leigqNEWTON_cert_resMin, checkNEWTON
