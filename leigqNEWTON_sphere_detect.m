function [lambdaAll, lambdaSamples, lambda0, class, sph, info] = leigqNEWTON_sphere_detect(A, varargin)
%LEIGQNEWTON_SPHERE_DETECT  Detect candidate eigen-spheres by sampling left eigenvalues.
%
%   [lambdaAll, lambdaSamples, lambda0, cls, sph, info] = ...
%       leigqNEWTON_sphere_detect(A, Name,Value,...)
%
%   This is the recommended user-facing entry point for "sphere detection".
%   It repeatedly calls LEIGQNEWTON to obtain many left-eigenvalue candidates and
%   looks for clusters that lie on a 2-sphere in R^4.
%
%   Internally, this is a thin wrapper around LEIGQNEWTON_SPHERE_SAMPLE.
%
% Inputs
%   A : n-by-n quaternion matrix.
%
% Name,Value options
%   This wrapper accepts the same options as LEIGQNEWTON_SPHERE_SAMPLE.
%   The most important ones are:
%     'Collect'      target number of DISTINCT samples (default 5*n)
%     'RunsMax'      maximum solver runs used to reach 'Collect' (default 50)
%     'Seed0'        base seed (default 1)
%     'UniqueTol'    distinctness tolerance (default 1e-8)
%     'DetectSphere' true/false (default true)
%     'Report'       'off'|'summary'|'progress'|'full'
%
%   Default reporting
%     If you do not specify 'Report', this wrapper enables 'Report'='progress'
%     so you can see that sampling proceeds.
%
% Outputs
%   lambdaAll        all collected candidates (counted list)
%   lambdaSamples    DISTINCT subset used for clustering/sphere detection
%   lambda0          diagnostics (e.g., zero handling, representatives)
%   cls              cluster labels for lambdaSamples
%   sph              detected sphere model(s) (empty if none detected)
%   info             counts, timings, thresholds
%
% Examples
%   % Quick detection with progress printing:
%   [lamAll,lamS,~,cls,sph,info] = leigqNEWTON_sphere_detect(A,'Collect',30);
%
%   % Silent mode:
%   [lamAll,lamS] = leigqNEWTON_sphere_detect(A,'Collect',30,'Report','off');
%
% Next step
%   Use LEIGQNEWTON_SPHERE_VALIDATE to confirm/refine a detected sphere.
%
% See also leigqNEWTON_sphere_sample, leigqNEWTON_sphere_validate, leigqNEWTON

% Inject default Report='progress' unless user provided Report explicitly.
hasReport = false;
for i = 1:2:numel(varargin)
    if ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i}))
        if strcmpi(string(varargin{i}), "Report")
            hasReport = true;
            break;
        end
    end
end
if ~hasReport
    varargin = [varargin, {'Report','progress'}];
end

[lambdaAll, lambdaSamples, lambda0, class, sph, info] = leigqNEWTON_sphere_sample(A, varargin{:});
end
