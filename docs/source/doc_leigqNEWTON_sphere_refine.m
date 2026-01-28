%% leigqNEWTON_sphere_refine
% Refine/validate a *given* sphere candidate for left eigenvalues.
% Convert to Live Script: File → Save As… → Live Script (*.mlx)
%
% Notes on naming:
%   * |leigqNEWTON_sphere_refine| is the main implementation (engine).
%   * |leigqNEWTON_sphere_validate| is the recommended user-facing entry point.
%     In the public toolbox, VALIDATE is a thin wrapper that forwards all inputs
%     to REFINE:
%         leigqNEWTON_sphere_validate(varargin{:}) == leigqNEWTON_sphere_refine(varargin{:})
%
% IMPORTANT (inputs):
%   Both functions require a *sphere model* in the inputs, typically as a field
%   |cases.sph| inside the case struct. If no sphere model is supplied, the
%   verdict is "NO_SPHERE_MODEL" (this is not an error; it just means: nothing
%   to validate).
%
% Typical pipeline (recommended for real problems):
%   [lamAll,lamS,lam0,cls,sph,info] = leigqNEWTON_sphere_detect(A, ...);
%   C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'cls',cls,'sph',sph,'info',info);
%   [A2,lamAll2,resAll,lamS2,resS,sph2,conf] = leigqNEWTON_sphere_validate(C, ...);

%% Setup (requirements + path)
hasQuat = true;
try
    quaternion(0,0,0,0);
catch
    hasQuat = false;
end
if ~hasQuat
    disp('This toolbox requires the built-in MATLAB quaternion class.');
    return
end
if exist('leigqNEWTON_sphere_refine','file') ~= 2
    error('Sphere functions not found on the MATLAB path. Add the toolbox root folder.');
end

%% Example (documentation-sized, fast)
% For documentation, we use a tiny 2-by-2 example with a known spherical component.
% We provide a *pre-built* candidate sphere model (as if it came from ..._sphere_detect).
%
% The run is intentionally fast (coarse grid, no polishing). For a more reliable
% decision, increase GridTheta/GridPhi and allow polishing (see the options below).

qi = quaternion(0,1,0,0);
A = [ quaternion(2,0,0,0),  qi;
     -qi,                  quaternion(2,0,0,0) ];

% Candidate samples near the sphere (small perturbation for illustration)
s2 = 1/sqrt(2);
lamS = [ quaternion(2, 1, 0, 0);
         quaternion(2,-1, 0, 0);
         quaternion(2, 0, 1, 0);
         quaternion(2, 0,-1, 0);
         quaternion(2, 0, 0, 1);
         quaternion(2, 0, 0,-1);
         quaternion(2, s2, s2, 0);
         quaternion(2,-s2,-s2, 0) ];
% Light perturbation (so refinement has something to do)
epsPert = 1e-3;
lamAll = lamS + quaternion(epsPert, -epsPert, 0, epsPert);

% Pre-built sphere model in affine form: p0 + B*(c3 + r*u),  u on unit S^2 in R^3
B  = [0 0 0; 1 0 0; 0 1 0; 0 0 1];     % maps R^3 -> imag(i,j,k) part
p0 = [2;0;0;0];                        % 4x1
c3 = [0;0;0];                          % 3x1
r  = 1;

sph = struct();
sph.p0       = p0;
sph.basis4x3  = B;
sph.center3   = c3;
sph.radius    = r;
sph.center4   = p0 + B*c3;
sph.center    = quaternion(sph.center4(1), sph.center4(2), sph.center4(3), sph.center4(4));
sph.inliers   = (1:numel(lamS)).';

fprintf('Note:');
fprintf('\n  Refinement is a lengthy process and may take several minutes.');
fprintf('\n  Therefore, a shortened version with inconclusive results');
fprintf('\n  is called for documentation purposes.\n');


% Assemble minimal case struct
C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'sph',sph,'idx',NaN,'seedUsed',NaN,'origin',"doc_sphere_refine");

% Documentation profile: fast + quiet (no polishing), coarse validation grid
RefineArgs = {'DoPolish', false, 'Verbose', 0, 'MaxIter', 120, 'TolX', 1e-12, 'TolFun', 1e-12};
[A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_refine( ...
    C, 'Verbose', 1, 'RefineArgs', RefineArgs, 'GridTheta', 4, 'GridPhi', 8, 'MinGridForDecision', 32); %#ok<ASGLU>

%% Nicely printed summaries
fprintf('resAll: median = %.3e, max = %.3e\n', median(resAll), max(resAll));
fprintf('resS  : median = %.3e, max = %.3e\n', median(resS),   max(resS));

disp('sph2 (refined model):');
disp(sph2)

disp('conf.sphere (decision + reliability):');
disp(conf.sphere)

%% See also
% leigqNEWTON_sphere_validate, leigqNEWTON_sphere_detect, leigqNEWTON_sphere_sample
