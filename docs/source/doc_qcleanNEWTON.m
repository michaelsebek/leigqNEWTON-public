%% qcleanNEWTON — componentwise cleanup (zeroing) of quaternion arrays
% Zero small components of a quaternion array to improve readability and promote sparsity.

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

if exist('qcleanNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'qcleanNEWTON.m'),'file')
            addpath(rootGuess);
            rehash toolboxcache
        end
    end
end

if exist('qcleanNEWTON','file') ~= 2
    error('qcleanNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |Y = qcleanNEWTON(X)|
% * |Y = qcleanNEWTON(X, tol)|

%% Notes
% * Cleaning is **componentwise**: it does not preserve unit quaternions.
% * If X represents rotations, re-normalize after cleaning.

%% Example 1: a single quaternion (before/after)
q  = quaternion(1, 1e-14, -2e-13, 3);
q2 = qcleanNEWTON(q, 1e-12);

disp('Before:'); disp(q);
disp('After (tol=1e-12):'); disp(q2);

%% Example 2: a matrix (before/after)
rng(1);
A  = quaternion(randn(3), 1e-14*randn(3), randn(3), zeros(3));
A2 = qcleanNEWTON(A, 1e-12);

disp('A (original):');  disp(A);
disp('A2 (cleaned):');  disp(A2);

% Quantify how much was removed in the tiny components
[~,bx,~,~] = parts(A);
[~,bx2,~,~] = parts(A2);
fprintf('Max |x|-component before: %.3e   after: %.3e\n', max(abs(bx(:))), max(abs(bx2(:))));

%% See also
% qroundNEWTON, parts
