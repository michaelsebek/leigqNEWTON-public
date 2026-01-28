function varargout = leigqNEWTON_sphere_validate(varargin)
%LEIGQNEWTON_SPHERE_VALIDATE  Validate sphere candidates and decide 
%                               sphere-vs-artifact (user-facing entry point).
%                             Alias of leigqNEWTON_sphere_refine 
%                               (same signature and outputs).
%   
%   In the public release, LEIGQNEWTON_SPHERE_VALIDATE is a thin wrapper that
%   forwards all inputs/outputs to LEIGQNEWTON_SPHERE_REFINE.
%
%   Use this function as the recommended entry point for sphere-hunting: it refines
%   sampled left-eigenvalue candidates, computes residual certificates, and (when a
%   sphere model is provided) performs optional grid-based checks and refitting to
%   decide whether a detected sphere is likely genuine or an artifact.
%
%   IMPORTANT (sphere model required for a sphere verdict)
%   -----------------------------------------------------
%   A sphere decision requires a sphere model in the input case, typically C.sph,
%   produced by LEIGQNEWTON_SPHERE_DETECT or LEIGQNEWTON_SPHERE_SAMPLE.
%   If no sphere model is available in the input, the routine can still refine
%   eigenvalues and certificates, but it will return sph=[] and the verdict will be
%   NO_SPHERE_MODEL (i.e., "nothing to validate").
%
% Syntax (recommended)
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_validate(C, Name,Value,...)
%
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_validate(cases, k, Name,Value,...)
%
% Convenience positional form (refinement-only; no sphere verdict unless sph provided via C)
%   [A, lamAll, resAll, lamSamples, resSamples, sph, conf, out] = ...
%       leigqNEWTON_sphere_validate(A, lamAll, lamSamples, Name,Value,...)
%
% Inputs
%   C / cases     Case container (struct, struct array, or cell array).
%                Each case should contain at least:
%                   C.A          (n-by-n quaternion)
%                   C.lamAll     (M-by-1 quaternion candidates)
%                   C.lamSamples (K-by-1 quaternion DISTINCT samples; may be [])
%                For sphere validation/decision, also provide:
%                   C.sph        sphere model struct (from sphere_detect/sample)
%                Additional fields are allowed (e.g., cls, info, seed, idx, origin).
%
%   k            Case index when passing a list (struct array / cell array).
%                Omit k when passing a single case struct C.
%
% Outputs
%   A                 The matrix of the selected case (same as input A).
%   lamAll            Refined candidate list.
%   lamSamples        Refined DISTINCT sample list (may be []).
%   resAll            Residual certificates for lamAll (via LEIGQNEWTON_CERT_RESMIN).
%   resSamples        Residual certificates for lamSamples.
%   sph               Updated/refitted sphere model (empty if none/invalid/missing).
%   conf              Confidence struct, including decision/verdict fields.
%   out               Detailed diagnostics (grid statistics, inliers, etc.).
%
% Examples
%   % Detect -> validate (recommended workflow)
%   [lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_detect(A, ...
%       'Collect',8,'RunsMax',12,'Restarts',40,'Seed0',1);
%   C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'cls',cls,'sph',sph,'info',info);
%   [A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_validate(C, ...
%       'GridTheta',4,'GridPhi',4,'GridMaxPoints',32,'Verbose',1);
%
%   % Sample -> validate (sphere model provided by sphere_sample)
%   [lamAll, lamS, sph, infoS] = leigqNEWTON_sphere_sample(A, ...
%       'Collect',8,'RunsMax',12,'Restarts',40,'Seed0',1);
%   C = struct('A',A,'lamAll',lamAll,'lamSamples',lamS,'sph',sph);
%   [A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_validate(C);
%
% Notes
%   For the full list of Name,Value options and detailed behavior, see
%   LEIGQNEWTON_SPHERE_REFINE (this function forwards to it).
%
% See also leigqNEWTON_sphere_refine, leigqNEWTON_sphere_detect,
%          leigqNEWTON_sphere_sample, leigqNEWTON_refine_auto,
%          leigqNEWTON_cert_resMin

[varargout{1:nargout}] = leigqNEWTON_sphere_refine(varargin{:});

end
