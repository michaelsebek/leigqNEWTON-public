%% Examples — Huang–So (LAA 323, 2001) explicit 2×2 cases
% Reproduce the explicit 2×2 examples (Examples 2.5–2.7) from:
%   L. Huang and W. So, "On left eigenvalues of a quaternionic matrix,"
%   Linear Algebra and its Applications 323 (2001) 105–116.
%   DOI: 10.1016/S0024-3795(00)00246-9
%
% The shipped script |examples/ExNEWTON_1_HuangSo.m| builds the matrices using
% MATLAB built-in |quaternion| class and runs |checkNEWTON| plus the sphere
% utilities for Example 2.7.
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

%% Requirements
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

%% Run the shipped script
root = fileparts(which('leigqNEWTON'));
ex   = fullfile(root,'examples','ExNEWTON_1_HuangSo.m');
if exist(ex,'file')
    run(ex);
else
    error('Example not found: %s', ex);
end

%% What to look for
% * Example 2.5 and 2.6: two isolated left eigenvalues, agreement with the paper.
% * Example 2.7: an infinite left spectrum (sphere family). The script runs
%   |leigqNEWTON_sphere_sample| and then |leigqNEWTON_sphere_validate| (which calls
%   |leigqNEWTON_sphere_refine| in the public release) and prints the fitted sphere model
%   sphere model and a confidence report.

%% See also
% doc_SphereHunting, leigqNEWTON_sphere_sample, leigqNEWTON_sphere_refine, checkNEWTON