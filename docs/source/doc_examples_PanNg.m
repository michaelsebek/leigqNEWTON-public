%% Examples — Pan–Ng (2024) quaternion circulant 4×4
% Reproduce the dense 4×4 quaternion circulant matrix example and run solve+refine.
% The shipped script is |examples/ExNEWTON_3.m|.
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
ex   = fullfile(root,'examples','ExNEWTON_3.m');
if exist(ex,'file')
    run(ex);
else
    error('Example not found: %s', ex);
end

%% Notes
% * The script builds the Pan–Ng 4×4 matrix explicitly, runs |checkNEWTON|,
%   then (optionally) refines/polishes the returned candidates.
% * If you want a faster run, switch SolveProfile 'reliable' → 'default' in the script.

%% See also
% checkNEWTON, leigqNEWTON_refine_batch, doc_RefinementAndCertificates
