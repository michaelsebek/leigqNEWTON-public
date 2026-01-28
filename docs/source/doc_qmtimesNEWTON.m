%% qmtimesNEWTON — quaternion matrix multiply helper
% Compute quaternion matrix products using a real embedding (toolbox-controlled helper).
%
% Some MATLAB installations do not support the built-in operator |A*B| for quaternion arrays.
% This helper provides a stand-alone, reproducible alternative used throughout leigqNEWTON.

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

if exist('qmtimesNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'qmtimesNEWTON.m'),'file')
            addpath(rootGuess);
            rehash toolboxcache
        end
    end
end
if exist('qmtimesNEWTON','file') ~= 2
    error('qmtimesNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |C = qmtimesNEWTON(A,B)|

%% Example
rng(1);
A = quaternion(randn(2),randn(2),randn(2),randn(2));
B = quaternion(randn(2),randn(2),randn(2),randn(2));

C = qmtimesNEWTON(A,B);

disp('A:'); disp(A);
disp('B:'); disp(B);
disp('C = qmtimesNEWTON(A,B):'); disp(C);

% If your MATLAB supports |A*B| for quaternion arrays, compare results.
try
    Cbuiltin = A*B;
    D = C - Cbuiltin;
    [w,x,y,z] = parts(D);
    maxAbsDiff = max(abs([w(:);x(:);y(:);z(:)]));
    fprintf('Max abs difference vs built-in A*B: %.3e\n', maxAbsDiff);
catch
    disp('Note: built-in quaternion A*B is not available in this MATLAB configuration.');
    disp('      qmtimesNEWTON provides this functionality using a real embedding.');
end

%% See also
% qmldivideNEWTON, qmrdivideNEWTON
