%% leigqNEWTON_sphere_detect
% User-facing entry point for detecting candidate eigen-spheres.
%
% Internally this is a thin wrapper around |leigqNEWTON_sphere_sample| with a
% default progress report (unless you set |'Report','off'|).
%
% This page uses a short "doc profile" and prints what was found.

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
if exist('leigqNEWTON_sphere_detect','file') ~= 2
    error('leigqNEWTON_sphere_detect not found on the MATLAB path. Add the toolbox root folder.');
end

%% Test matrix with a known eigen-sphere (2x2)
qi = quaternion(0,1,0,0);
A  = [ quaternion(2,0,0,0),  qi;
      -qi,                 quaternion(2,0,0,0) ];

%% Doc profile (fast)
Collect  = 10;
RunsMax  = 20;
Restarts = 80;

[lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_detect( ...
    A, 'Collect',Collect, 'RunsMax',RunsMax, 'Restarts',Restarts, ...
    'Seed0',1, 'Report','summary');

fprintf('Distinct samples: %d, counted list: %d, spheres found: %d\n', numel(lamS), numel(lamAll), ~isempty(sph));
if ~isempty(sph)
    disp(sph);
end

%% RELIABLE profile (optional, slower)
RUN_SPHERE_DETECT_RELIABLE = false;
try
    RUN_SPHERE_DETECT_RELIABLE = evalin('base','RUN_SPHERE_DETECT_RELIABLE');
catch
end
if RUN_SPHERE_DETECT_RELIABLE
    disp('RELIABLE profile: this may take minutes.');
    [lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_detect( ...
        A, 'Collect', 40, 'RunsMax', 120, 'Restarts', 600, 'Seed0', 1, 'Report','progress'); %#ok<ASGLU>
    disp(sph);
end

%% See also
% leigqNEWTON_sphere_sample, leigqNEWTON_sphere_validate, leigqNEWTON_sphere_refine
