%% leigqNEWTON_quickstart
% A short page of minimal example calls (Live Script friendly).
% You can copy/paste one‑liners.

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

if exist('leigqNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'leigqNEWTON.m'),'file')
            addpath(rootGuess);
        end
    end
end

if exist('leigqNEWTON','file') ~= 2
    error('leigqNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Create a small demo matrix (Pan–Ng 4×4 circulant)
q = @(w,x,y,z) quaternion(w,x,y,z);
a = q(-2, 1, 1, 4);
b = q( 2, 4, 1, 1);
c = q( 1, 3, 2, 2);
d = q(-1, 2, 2, 3);
A = [a b c d;
     d a b c;
     c d a b;
     b c d a];

%% Solve (candidates)
% Request more hits than n to reduce the chance of missing distinct eigenvalues.
[lambda, V, res, info, lambdaU, VU, resU] = leigqNEWTON(A, 'K', 2*size(A,1), 'Seed',1, 'SolveProfile','reliable');
median(resU), max(resU)

%% Refine + certify (recommended)
[lambdaR, VR, cert] = leigqNEWTON_refine_batch(A, lambdaU, 'DoPolish',true, 'TargetResMin',1e-14);
median(cert.resMin), max(cert.resMin)

%% Certificates only (eigenvalue-only vs eigenpair)
rMin  = leigqNEWTON_cert_resMin(A, lambdaR(1));
rPair = leigqNEWTON_cert_resPair(A, lambdaR(1), VR(:,1));
[rMin rPair]

%% One-call report
out = checkNEWTON(A, lambdaR, VR, 'SphereCheck','off');

%% See also
% doc_GettingStarted, doc_leigqNEWTON, doc_RefinementAndCertificates
