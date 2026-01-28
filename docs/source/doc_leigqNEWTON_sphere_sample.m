%% leigqNEWTON_sphere_sample
% Sampling engine for sphere-hunting: repeatedly runs |leigqNEWTON| and optionally
% applies simple sphere detection on the DISTINCT samples.
%
% This page runs a deliberately small "doc profile" (seconds) and prints a short
% summary. For sphere validation, see |doc_leigqNEWTON_sphere_validate|.

%% NOTE (documentation profile)
% This live-script source is designed to be *fast and quiet* by default.
% It demonstrates the API of the documented function without running long
% sphere-hunting loops. For more reliable (but slower) settings, see the
% "RELIABLE profile" section near the end of this page.
%
% The examples require MATLAB's built-in quaternion class.


%% Requirements and path
hasQuat = true;
try, quaternion(0,0,0,0); catch, hasQuat = false; end
if ~hasQuat
    disp("This toolbox requires MATLAB's built-in quaternion class (quaternion(w,x,y,z)).");
    return;
end
if exist('leigqNEWTON_sphere_sample','file') ~= 2
    error('leigqNEWTON_sphere_sample not found on the MATLAB path. Add the toolbox root folder.');
end

%% Test matrix with a known eigen-sphere (Huang–So style 2x2 example)
% A = [2, i; -i, 2] has a spherical component with center ~ 2 and radius ~ 1.
qi = quaternion(0,1,0,0);
A  = [ quaternion(2,0,0,0),  qi;
      -qi,                 quaternion(2,0,0,0) ];

%% Doc profile (fast)
Collect   = 8;     % DISTINCT samples requested
RunsMax   = 12;
Restarts  = 40;
Seed0     = 1;

[lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_sample( ...
    A, 'Collect',Collect, 'RunsMax',RunsMax, 'Restarts',Restarts, 'Seed0',Seed0, ...
    'DetectSphere',true, 'Report','summary');

% Show the detected sphere model (if any)
if isempty(sph)
    disp('No sphere detected in this short run (this can happen). Increase Collect/RunsMax/Restarts.');
else
    disp('Detected sphere candidate(s):');
    disp(sph);
end

%% RELIABLE profile (optional, slower)
% Set RUN_SPHERE_SAMPLE_RELIABLE = true in the base workspace and re-run.
RUN_SPHERE_SAMPLE_RELIABLE = false;
try
    RUN_SPHERE_SAMPLE_RELIABLE = evalin('base','RUN_SPHERE_SAMPLE_RELIABLE');
catch
end
if RUN_SPHERE_SAMPLE_RELIABLE
    disp('RELIABLE profile: this may take minutes (many solver restarts).');
    [lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_sample( ...
        A, 'Collect', max(30, 5*size(A,1)), 'RunsMax', 80, 'Restarts', 400, 'Seed0', 1, ...
        'DetectSphere', true, 'Report','progress'); %#ok<ASGLU>
    disp(sph);
end

%% See also
% leigqNEWTON_sphere_detect, leigqNEWTON_sphere_validate, leigqNEWTON_sphere_refine, leigqNEWTON
