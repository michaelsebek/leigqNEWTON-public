%% qroundNEWTON — componentwise rounding of quaternion arrays
% Round quaternion components to improve readability (for logs, tables, diagnostics).

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

if exist('qroundNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'qroundNEWTON.m'),'file')
            addpath(rootGuess);
            rehash toolboxcache
        end
    end
end

if exist('qroundNEWTON','file') ~= 2
    error('qroundNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |Y = qroundNEWTON(X)|
% * |Y = qroundNEWTON(X, ndigits)|

%% Example 1: a single quaternion
q  = quaternion(1.234567, -0.0001234, 2, -3.1415926);
q2 = qroundNEWTON(q, 3);

disp('Before:'); disp(q);
disp('After (ndigits=3):'); disp(q2);

%% Example 2: a matrix
rng(2);
A  = quaternion(randn(3),randn(3),randn(3),randn(3));
A2 = qroundNEWTON(A, 2);

disp('A (original):'); disp(A);
disp('A2 (rounded to 2 digits):'); disp(A2);

%% See also
% qcleanNEWTON, parts
