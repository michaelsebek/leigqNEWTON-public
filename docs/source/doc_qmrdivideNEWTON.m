%% qmrdivideNEWTON — right matrix division B/A for quaternion matrices
% Solve X*A = B for X using a real embedding.
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

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

if exist('qmrdivideNEWTON','file') ~= 2
    thisFile = mfilename('fullpath');
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,'qmrdivideNEWTON.m'),'file')
            addpath(rootGuess);
            rehash toolboxcache
        end
    end
end
if exist('qmrdivideNEWTON','file') ~= 2
    error('qmrdivideNEWTON not found on the MATLAB path. Add the toolbox root folder.');
end

%% Syntax
% * |X = qmrdivideNEWTON(B,A)|

%% Example
rng(1);
n = 3;
A = quaternion(randn(n),randn(n),randn(n),randn(n));
% Make A better conditioned for the demo:
A = A + quaternion(5,0,0,0)*eye(n);

B = quaternion(randn(n),randn(n),randn(n),randn(n));
X = qmrdivideNEWTON(B,A);

disp('A:'); disp(A);
disp('B:'); disp(B);
disp('X = B/A (computed):'); disp(X);

% Check residual: ||X*A - B|| (componentwise max-norm)
R = qmtimesNEWTON(X,A) - B;
[w,x,y,z] = parts(R);
maxAbsRes = max(abs([w(:);x(:);y(:);z(:)]));
fprintf('Residual max-abs component: %.3e\n', maxAbsRes);

%% See also
% qmldivideNEWTON, qmtimesNEWTON
