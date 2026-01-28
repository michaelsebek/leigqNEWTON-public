function leigqNEWTON_quickstart
%LEIGQNEWTON_QUICKSTART  Minimal workflow examples for the LEIGQNEWTON toolbox.
%
% This file is intentionally “copy/paste friendly”: each example is a
% one-liner (or two) that you can recall quickly.
%
% -------------------------------------------------------------------------
% Example 0: Create a random quaternion matrix (demo only)
%   n = 5; A = quaternion(randn(n),randn(n),randn(n),randn(n));
%
% -------------------------------------------------------------------------
% Example 1: Solve for candidate left eigenvalues
%   [lam,V,res,info,lamU,VU,resU] = leigqNEWTON(A,'Num',2*size(A,1),'Seed',1,'Verbose',false);
%
% Example 2: Refine + certify a list of candidates (recommended)
%   [lamR,vR,cert] = leigqNEWTON_refine_batch(A,lam);
%
% Example 3: Refine a single lambda (lambda-only workflow)
%   [lam1,v1] = leigqNEWTON_refine_lambda(A,lam(1));
%
% Example 4: Polish an eigenpair (when you already have v0)
%   [lamP,vP,resP] = leigqNEWTON_refine_polish(A,lam(1),v(:,1));
%
% Example 5: Certificates only (useful for acceptance thresholds)
%   rMin  = leigqNEWTON_cert_resMin(A,lam(1));
%   rPair = leigqNEWTON_cert_resPair(A,lam(1),v(:,1));
%
% Example 6: Sphere hunting (experimental; many samples + detection)
%   [lamAll,lamS,lam0,cls,sph,info] = leigqNEWTON_sphere_sample(A,'Collect',250,'Seed0',1,'Report','progress');
%   [A2,lamAll2,resAll,lamS2,resS,sph2,conf,out] = leigqNEWTON_sphere_refine(A,lamAll,lamS);  % positional overload
% caseStruct = struct('A',A,'lamAll',lamAll,'lamSamples',lamS); [A2,lamAll2,resAll,lamS2,resS,sph2,conf,out] = leigqNEWTON_sphere_refine(caseStruct,[]);
%
% -------------------------------------------------------------------------
% Legacy name mapping (original names you used earlier)
%   leigqVANILA                   -> leigqNEWTON
%   leigqVANILA_polish            -> leigqNEWTON_refine_polish
%   leigqVANILAresmin             -> leigqNEWTON_cert_resMin
%   leigqVANILArespair            -> leigqNEWTON_cert_resPair
%   leigqVANILA_refineAuto(_v2)   -> leigqNEWTON_refine_auto
%   leigqVANILA_refineLambdas     -> leigqNEWTON_refine_batch
%   leigqVANILAsphEXP2            -> leigqNEWTON_sphere_sample
%   refineSPHERE                  -> leigqNEWTON_sphere_refine
%
% Tip: Type "help leigqNEWTON" or "help leigqNEWTON_refine_batch" for full API docs.
%
% No runtime code needed; this is a documentation script.
end
