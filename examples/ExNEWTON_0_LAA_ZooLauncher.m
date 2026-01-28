function ExNEWTON_0_LAA_ZooLauncher()
%ExNEWTON_0_LAA_ZooLauncher  Run the optional LAA_Zoo companion bundle (if present).
%
% This helper is provided for the PUBLIC distribution where we avoid addpath(genpath(...)).
% It adds the LAA_Zoo folder to the path (without genpath), runs its overview (if any),
% and prints guidance if the bundle is not found.
%
% Typical use:
%   root = fileparts(which("leigqNEWTON"));
%   addpath(root); addpath(fullfile(root,"examples")); rehash toolboxcache
%   ExNEWTON_0_LAA_ZooLauncher
%
% The launcher looks for:
%   (A) <root>/LAA_Zoo        (recommended when bundled inside leigqNEWTON_public)
%   (B) sibling folder        <parent-of-root>/LAA_Zoo
%
% See also: PACKAGEVerifierNEWTON_fromContents

root = fileparts(which("leigqNEWTON"));
if isempty(root)
    error('leigqNEWTON not found on path. Add the toolbox root folder first.');
end

candA = fullfile(root,'LAA_Zoo');
candB = fullfile(fileparts(root),'LAA_Zoo');

if exist(candA,'dir')==7
    zooDir = candA;
elseif exist(candB,'dir')==7
    zooDir = candB;
else
    fprintf('LAA_Zoo bundle not found.\n');
    fprintf('Expected either:\n  %s\n  %s\n', candA, candB);
    fprintf('If you downloaded LAA_Zoo separately, place it next to the toolbox root or inside it.\n');
    return;
end

addpath(zooDir);
rehash toolboxcache

% Try to run an overview script if it exists
ov1 = fullfile(zooDir,'LAA_Zoo_Bundle_Overview.m');
ov2 = fullfile(zooDir,'ReadMe.md');

if exist(ov1,'file')==2
    fprintf('Running Zoo overview: %s\n', ov1);
    run(ov1);
elseif exist(ov2,'file')==2
    fprintf('Zoo is on path. Please open/read:\n  %s\n', ov2);
    fprintf('Then run one of the scripts LAA_Zoo_Ex_*_computation.\n');
else
    fprintf('Zoo is on path. Run one of the scripts LAA_Zoo_Ex_*_computation.\n');
end
end
