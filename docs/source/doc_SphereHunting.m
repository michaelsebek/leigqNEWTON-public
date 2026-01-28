%% leigqNEWTON — Sphere Hunting (experimental)
% Convert to Live Script: File → Save As… → Live Script (*.mlx)
%
% Sphere hunting tries to detect and validate *spherical components* in the left spectrum.
% Typical workflow:
%   1) sample many left-eigenvalue candidates (Newton runs)
%   2) detect/fit candidate spheres from DISTINCT samples
%   3) refine + validate (grid test) to decide sphere-vs-artifact
%
% This page is written to be *documentation-friendly*:
%   - Default profile aims to finish quickly (seconds), yet produce a reliable decision.
%   - A heavier profile is provided for research runs.

%% 0) Setup (requirements + path)
rng(0,'twister'); 
hasQuat = true;
try
    quaternion(0,0,0,0);
catch
    hasQuat = false;
end
if ~hasQuat
    disp("This toolbox requires MATLAB's built-in quaternion class (Aerospace Toolbox).");
    disp("Examples in this page are skipped.");
    return;
end

% Ensure toolbox root is on the path (add only the root, no genpath).
if exist("leigqNEWTON_sphere_detect","file") ~= 2
    thisFile = mfilename("fullpath");
    if ~isempty(thisFile)
        rootGuess = fileparts(fileparts(fileparts(thisFile))); % .../docs/source -> toolbox root
        if exist(fullfile(rootGuess,"leigqNEWTON_sphere_detect.m"),"file")
            addpath(rootGuess);
            rehash toolboxcache
        end
    end
end
if exist("leigqNEWTON_sphere_detect","file") ~= 2
    error("Sphere functions not found on the MATLAB path. Add the toolbox root folder.");
end

%% 1) Choose a profile (recommended: "DOC")
% Profiles control the *work budget* (RunsMax×Restarts) and the *validation grid* size.
% If you previously observed very long runtimes, reduce RunsMax and/or Restarts.

Profile = "DOC";   % "DOC" (default) | "FAST" (smoke) | "FULL" (research)

switch upper(Profile)
    case "FAST"
        % Very quick smoke test. May produce "INCONCLUSIVE" if the validation grid is too small.
        Collect   = 8;
        RunsMax   = 12;
        Restarts  = 40;
        GridTheta = 6;
        GridPhi   = 12;
        GridMaxPoints = 60;   % <-- may fall below MinGridForDecision (see below)
        MinGridForDecision = 64;

        TargetRes = 1e-12;    % relaxed for smoke tests
        RefineArgs = {"Mode","auto","TargetResMin",1e-12,"Verbose",0};

    case "DOC"
        % Documentation-friendly: usually seconds, and the grid is large enough
        % to avoid the "INSUFFICIENT GRID SIZE" / "INCONCLUSIVE" outcome.
        Collect   = 16;
        RunsMax   = 30;
        Restarts  = 80;
        GridTheta = 10;
        GridPhi   = 20;
        GridMaxPoints = 200;  % >= MinGridForDecision
        MinGridForDecision = 64;

        TargetRes = 1e-13;    % certificate threshold
        RefineArgs = {"Mode","auto","TargetResMin",1e-13,"Verbose",0};

    case "FULL"
        % Research run (can be slow). Use when you need more distinct samples.
        Collect   = 30;
        RunsMax   = 80;
        Restarts  = 200;
        GridTheta = 12;
        GridPhi   = 24;
        GridMaxPoints = 400;
        MinGridForDecision = 128;

        TargetRes = 1e-14;
        RefineArgs = {"Mode","auto","TargetResMin",1e-14,"Verbose",0};

    otherwise
        error("Unknown Profile: %s", Profile);
end

approxRestarts = RunsMax*Restarts;
approxGrid = min(GridTheta*GridPhi, GridMaxPoints);
fprintf("Profile=%s | Collect=%d | RunsMax=%d | Restarts=%d | approxRestarts=%d | grid~%d pts\n", ...
    Profile, Collect, RunsMax, Restarts, approxRestarts, approxGrid);

%% 2) A small test matrix with a known spherical component (Huang–So style)
% This 2-by-2 example is intentionally small so that documentation runs are fast.
qi = quaternion(0,1,0,0);
A  = [ quaternion(2,0,0,0),  qi;
      -qi,                  quaternion(2,0,0,0) ];

%% 3) Detect: collect samples and fit sphere candidates (from DISTINCT samples)
tDetect = tic;
fprintf("\nDetection is fast:");
fprintf("\n=== leigqNEWTON_sphere_detect ===\n");
[lamAll, lamS, lam0, cls, sph, info] = leigqNEWTON_sphere_detect( ...
    A, ...
    "Collect",  Collect, ...
    "RunsMax",  RunsMax, ...
    "Restarts", Restarts, ...
    "Seed0",    1, ...
    "Report",   "off"); 

fprintf("Detect done in %.2fs: Ktot=%d  Kdistinct=%d  spheres=%d\n", ...
    toc(tDetect), numel(lamAll), numel(lamS), numel(sph));

%% 4) Validate / refine (recommended entry point)
% If you ever see:
%   "Sphere decision: INCONCLUSIVE ... Insufficient grid size ..."
% increase GridTheta/GridPhi and/or GridMaxPoints (and keep GridMaxPoints>=MinGridForDecision).
fprintf("\nValidation/refining is slow:");

C = struct();
C.A = A;
C.lamAll = lamAll;
C.lamSamples = lamS;
C.cls = cls;
C.sph = sph;
C.info = info;
C.idx = 1;
C.seedUsed = 1;
C.origin = "doc_SphereHunting";

tVal = tic;
fprintf("\n=== leigqNEWTON_sphere_validate ===\n");
[A2, lamAll2, resAll, lamS2, resS, sph2, conf, out] = leigqNEWTON_sphere_validate( ...
    C, ...
    "Verbose", 1, ...
    "TargetRes", TargetRes, ...
    "RefineArgs", RefineArgs, ...
    "GridTheta", GridTheta, ...
    "GridPhi", GridPhi, ...
    "GridMaxPoints", GridMaxPoints, ...
    "MinGridForDecision", MinGridForDecision, ...
    "RefineGrid", true); 

fprintf("Validate done in %.2fs.\n", toc(tVal));

%% 5) Quick interpretation hooks
% - resAll / resS: certificate-like values (smaller is better)
% - sph2: sphere model (empty if no sphere validated)
% - conf.sphere: decision + reliability score

% Residual summary
medRes = median(resAll);
maxRes = max(resAll);
fprintf("resAll: median = %.3e, max = %.3e\n", medRes, maxRes);

% Sphere model (struct or empty)
fprintf("sph2:\n");
disp(sph2)

% Sphere confidence (sub-struct)
fprintf("conf.sphere:\n");
disp(conf.sphere)

% Tip: if Max residual is above TargetRes, either:
%   - increase refinement effort (RefineArgs / TargetResMin), or
%   - relax TargetRes for documentation runs.
