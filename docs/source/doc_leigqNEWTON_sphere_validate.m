%% leigqNEWTON_sphere_validate
% Validate/refine a detected eigen-sphere candidate (recommended entry point).
% Convert to Live Script: File → Save As… → Live Script (*.mlx)
%
% IMPORTANT:
%   In the public toolbox, this function is a thin wrapper that forwards all
%   inputs to |leigqNEWTON_sphere_refine|. Therefore, the calling convention and
%   outputs are the same.
%
%   Both VALIDATE and REFINE expect a *case struct* that contains a sphere model
%   in |cases.sph|. If no sphere model is supplied, the verdict will be
%   "NO_SPHERE_MODEL".
%
% Typical pipeline (recommended):
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
if exist('leigqNEWTON_sphere_validate','file') ~= 2
    error('Sphere functions not found on the MATLAB path. Add the toolbox root folder.');
end

%% Example (documentation-sized, fast)
% We reuse the same small 2-by-2 example as in doc_leigqNEWTON_sphere_refine.
% Here we call the user-facing entry point VALIDATE, which forwards to REFINE.

qi = quaternion(0,1,0,0);
A = [ quaternion(2,0,0,0),  qi;
     -qi,                  quaternion(2,0,0,0) ];

s2 = 1/sqrt(2);
lamS = [ quaternion(2, 1, 0, 0);
         quaternion(2,-1, 0, 0);
         quaternion(2, 0, 1, 0);
         quaternion(2, 0,-1, 0);
         quaternion(2, 0, 0, 1);
         quaternion(2, 0, 0,-1);
         quaternion(2, s2, s2, 0);
         quaternion(2,-s2,-s2, 0) ];
epsPert = 1e-3;
lamAll = lamS + quaternion(epsPert, -epsPert, 0, epsPert);

B  = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
p0 = [2;0;0;0];
c3 = [0;0;0];
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
fprintf('\n  Validation is a lengthy process and may take several minutes.');
fprintf('\n  Therefore, a shortened version with inconclusive results');
fprintf('\n  is called for documentation purposes.\n');

C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'sph',sph,'idx',NaN,'seedUsed',NaN,'origin',"doc_sphere_validate");

RefineArgs = {'DoPolish', false, 'Verbose', 0, 'MaxIter', 120, 'TolX', 1e-12, 'TolFun', 1e-12};
[A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_validate( ...
    C, 'Verbose', 1, 'RefineArgs', RefineArgs, 'GridTheta', 4, 'GridPhi', 8, 'MinGridForDecision', 32); %#ok<ASGLU>

%% Nicely printed summaries
fprintf('resAll: median = %.3e, max = %.3e\n', median(resAll), max(resAll));
fprintf('resS  : median = %.3e, max = %.3e\n', median(resS),   max(resS));

disp('sph2 (refined model):');
disp(sph2)

disp('conf.sphere (decision + reliability):');
disp(conf.sphere)

%% See also
% leigqNEWTON_sphere_refine, leigqNEWTON_sphere_detect, leigqNEWTON_sphere_sample
