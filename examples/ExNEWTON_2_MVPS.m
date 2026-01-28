%% ExNEWTON_2_MVPS.m (UPDATED for leigqNEWTON toolbox)
% Quaternion matrices from explicit examples in:
%   Macias-Virgos and Pereira-Saez,
%   "A topological approach to left eigenvalues of quaternionic matrices".
%
% This script builds the matrices exactly as in the paper examples (no changes)
% and runs the leigqNEWTON diagnostics.
%
% Requirements
%   - MATLAB class quaternion (Aerospace Toolbox).
%   - leigqNEWTON toolbox on the MATLAB path.
%
% -------------------------------------------------------------------------
clc;

% --- Sanity checks --------------------------------------------------------
try
    quaternion(0,0,0,0);
catch
    error(['This script requires MATLAB''s built-in quaternion class. ' ...
           'Ensure Aerospace Toolbox is installed and quaternion(w,x,y,z) works.']);
end

if exist('leigqNEWTON','file') ~= 2
    error(['leigqNEWTON not found on the MATLAB path. ' ...
           'Add the leigqNEWTON toolbox folder, e.g. addpath(genpath(<toolboxRoot>)).']);
end

% --- basis elements and handy constants ----------------------------------
q0 = quaternion(0,0,0,0);   % 0
q1 = quaternion(1,0,0,0);   % 1
qi = quaternion(0,1,0,0);   % i
qj = quaternion(0,0,1,0);   % j
qk = quaternion(0,0,0,1);   % k

% Pretty-print tolerance.
tolShow = 1e-12;

%% Example 19
% A = [ i  0  0
%       k  j  0
%      -3i 2k k ]
A19 = [  qi,    q0,   q0;
         qk,    qj,   q0;
      -3*qi,  2*qk,   qk ];

%% Example 38
% A = [ 0      i      1
%       3i-k   0      1
%       k   -1+j+k    0 ]
A38 = [  q0,           qi,          q1;
       3*qi - qk,      q0,          q1;
         qk,     (-q1 + qj + qk),   q0 ];

%% Example 51 (uses the matrix from Example 38)
% In the paper: pi_A = -i, and B = A - pi_A*Id = A + i*Id
B51 = A38 + qi*eye(3);

% B^{-1} = (1/10) * [ 4i-2k   -4i+2k    0
%                    -1-3i+8j-6k   1+3i-3j+k   -5j-5k
%                    11+i-8j-8k   -1-i+3j+3k   -5j+5k ]
scale = 0.1;
Binv51 = scale .* [  (4*qi - 2*qk),              (-4*qi + 2*qk),                 q0;
                     (-q1 - 3*qi + 8*qj - 6*qk), (q1 + 3*qi - 3*qj + qk),     (-5*qj - 5*qk);
                     (11*q1 + qi - 8*qj - 8*qk), (-q1 - qi + 3*qj + 3*qk),    (-5*qj + 5*qk) ];

%% Example 52
% A = [  j    1    0
%       2i   -k    1
%     2-i-2j  -1-j+k  -i-k ]
A52 = [  qj,                 q1,              q0;
       2*qi,               (-qk),             q1;
     (2*q1 - qi - 2*qj), (-q1 - qj + qk),   (-qi - qk) ];

%% Example 54 (p,q,r arbitrary; p,q nonzero in the paper)
% A = [ 0    -j    i
%      -1+j   j    k
%        p    q    r ]
% NOTE: The paper allows arbitrary p,q,r (with p,q nonzero).
% Keep the placeholders as NaN until you fill them in.
p54 = quaternion(NaN,NaN,NaN,NaN);   % <-- fill in
q54 = quaternion(NaN,NaN,NaN,NaN);   % <-- fill in (nonzero)
r54 = quaternion(NaN,NaN,NaN,NaN);   % <-- fill in
A54 = [  q0,      (-qj),     qi;
       (-q1+qj),   qj,       qk;
         p54,      q54,      r54 ];

%% Example 55
% A = [  k     0    0
%      3i-j   -i    i
%      1-2k    j   -j ]
A55 = [  qk,         q0,    q0;
       3*qi - qj,   (-qi),  qi;
       (q1 - 2*qk),  qj,   (-qj) ];

%% Example 56
% A = [ -i-j   0    0
%        k    -i    i
%       1-i    j   -j ]
A56 = [ (-qi - qj),  q0,    q0;
          qk,       (-qi),  qi;
        (q1 - qi),   qj,   (-qj) ];

%% Overview / diagnostics --------------------------------------------------

fprintf('\nMVPS topological examples (leigqNEWTON)\n');
fprintf('=============================================================\n');

% Show the matrices (cleaned) for quick inspection.
disp('A19 =');  disp(qcleanNEWTON(A19, tolShow));
disp('A38 =');  disp(qcleanNEWTON(A38, tolShow));
disp('B51 =');  disp(qcleanNEWTON(B51, tolShow));
disp('Binv51 =');disp(qcleanNEWTON(Binv51, tolShow));
disp('A52 =');  disp(qcleanNEWTON(A52, tolShow));
disp('A54 =');  disp(qcleanNEWTON(A54, tolShow));
disp('A55 =');  disp(qcleanNEWTON(A55, tolShow));
disp('A56 =');  disp(qcleanNEWTON(A56, tolShow));

% Quick sanity: verify that Binv51 is an inverse of B51 (matrix multiply via NEWTON helper).
try
    I3 = q1*eye(3);
    E = qmtimesNEWTON(B51, Binv51) - I3;
    fprintf('\nInverse check for Example 51: max|B51*Binv51 - I| = %.3e\n', local_qabsmax(E));
catch ME
    warning('Inverse check for Example 51 failed: %s', ME.message);
end

% Run checkNEWTON on each fully-specified matrix.
items = {
    struct('name','Example 19', 'A', A19)
    struct('name','Example 38', 'A', A38)
    struct('name','Example 51 (B)', 'A', B51)
    struct('name','Example 52', 'A', A52)
    struct('name','Example 54', 'A', A54)  % may contain NaNs
    struct('name','Example 55', 'A', A55)
    struct('name','Example 56', 'A', A56)
};

for t = 1:numel(items)
    name = items{t}.name;
    A = items{t}.A;

    fprintf('\n-------------------------------------------------------------\n');
    fprintf('%s\n', name);

    if local_hasNaN(A)
        fprintf('Skipped: matrix contains NaN placeholders (fill p,q,r first).\n');
        continue;
    end

    % Use a deterministic seed; increase SolveProfile if you want more robustness.
    out = checkNEWTON(A, [], [], 'SolveArgs', {'SolveProfile','reliable','Seed',1}); %#ok<NASGU>
end

fprintf('\nDone.\n');

%% ------------------------------------------------------------------------
%% Local helpers

function tf = local_hasNaN(Q)
%LOCAL_HASNAN True if quaternion array has any NaN in its components.
[w,x,y,z] = parts(Q);
tf = any(isnan(w(:)) | isnan(x(:)) | isnan(y(:)) | isnan(z(:)));
end

function m = local_qabsmax(Q)
%LOCAL_QABSMAX Max elementwise quaternion magnitude in an array.
[w,x,y,z] = parts(Q);
mags = sqrt(w.^2 + x.^2 + y.^2 + z.^2);
m = max(mags(:));
end
