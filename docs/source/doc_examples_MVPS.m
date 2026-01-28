%% Examples — MVPS topological 3×3 cases
% Reproduce explicit 3×3 examples from Macías‑Virgos & Pereira‑Sáez (MVPS)
% using the shipped script |examples/ExNEWTON_2_MVPS.m|.
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
ex   = fullfile(root,'examples','ExNEWTON_2_MVPS.m');
if exist(ex,'file')
    run(ex);
else
    error('Example not found: %s', ex);
end

%% What to look for
% * The script prints the matrices (cleaned) and runs |checkNEWTON| on each case.
% * Some cases are non-generic (e.g., more than n distinct left eigenvalues).
% * If you want to increase solver budgets, edit the script and change SolveProfile
%   (e.g., 'default' → 'reliable') or pass larger trial counts.

%% See also
% checkNEWTON, doc_RefinementAndCertificates, leigqNEWTON_refine_batch
