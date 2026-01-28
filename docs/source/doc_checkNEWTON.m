%% checkNEWTON — one-call diagnostic workflow
% Print a compact report, compute missing lambdas/vectors, and optionally trigger sphere diagnostics.
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

%% Requirements
% This toolbox relies on MATLAB's built-in quaternion class.
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

% Ensure toolbox is on the path (add root only, if needed).
if exist('checkNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'checkNEWTON.m'),'file')
            addpath(rootGuess);
        end
    end
end

if exist('checkNEWTON','file') ~= 2
    error('checkNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |out = checkNEWTON(A)|
% * |out = checkNEWTON(A, lambda)|
% * |out = checkNEWTON(A, lambda, V)|
% * |out = checkNEWTON(A, [], [], 'SolveArgs', {'SolveProfile','reliable','Seed',1}, 'SphereCheck','auto')|

%% Notes
% * If |lambda| is omitted, CHECKNEWTON calls |leigqNEWTON|.
% * If |V| is omitted, CHECKNEWTON computes best vectors by minimizing |resMin|.
% * For presentation, it cleans and rounds lambdas (use 'CleanTol'/'PrintDigits').

%% Examples
% These use fixed small matrices from Huang–So (LAA 323, 2001) for reproducibility.
q0 = quaternion(0,0,0,0);
q1 = quaternion(1,0,0,0);
qi = quaternion(0,1,0,0);
qj = quaternion(0,0,1,0);
qk = quaternion(0,0,0,1);

%% Example 1: simplest call (solver + report)
% Huang–So Example 2.6
A = [ q0, qi;
      qj, q1 ];
out = checkNEWTON(A);

%% Example 2: lambda-only workflow
% (you already have candidates; CHECKNEWTON will compute vectors and report residuals)
[lambda] = leigqNEWTON(A,'SolveProfile','fast','Seed',1);
out2 = checkNEWTON(A, lambda, []);

%% Example 3: sphere diagnostics (Huang–So Example 2.7)
A27 = [ quaternion(2,0,0,0),  qi;
       -qi,                 quaternion(2,0,0,0) ];

% With SphereCheck='on', checkNEWTON runs leigqNEWTON_sphere_sample internally.
out3 = checkNEWTON(A27, [], [], 'SphereCheck','on', 'SphereCollect',20, 'SphereSeed0',1);

%% See also
% |leigqNEWTON|, |leigqNEWTON_refine_batch|, |leigqNEWTON_sphere_sample|
