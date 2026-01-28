%% LEIGQNEWTON — Getting Started
% This Live Script is a copy/paste-friendly introduction to the leigqNEWTON toolbox.

%% 0) Setup (add toolbox to path)
% If you are in the toolbox root folder, this is enough:
addpath(genpath(pwd));

%% 1) Create a small demo quaternion matrix (reproducible)
rng(1);
n = 5;
A = quaternion(randn(n),randn(n),randn(n),randn(n));
display(A);

%% 2) Solve: get candidate left eigenvalues
% One-liner (recommended):
[lam,V,lamc,res,info] = leigqNEWTON(A,'SolveProfile','default','Seed',1);
[lambda, V, res, info, lambdaU, VU, resU] = leigqNEWTON(A, 'SolveProfile','default', 'Seed',1);

fprintf("\nAll lambda (lambda):\n");
display(lambda);

fprintf("\nDistinct lambda (lambdaU):\n");
display(lambdaU);

% Quality summary (lower is better):
%median(resU), max(resU)
medRes = median(resU);
maxRes = max(resU);
fprintf("resMin: median = %.3e, max = %.3e\n", medRes, maxRes);

%% 3) Refine + certify a batch of candidates (recommended workflow)
[lamR,VR,cert] = leigqNEWTON_refine_batch(A, lam,'Verbose',0);

fprintf("\nRefined lambda (lamR):\n");
display(lamR);

% Quality summary (lower is better):
%median(resU), max(resU)
medRes = median(cert.resMin);
maxRes = max(cert.resMin);
fprintf("resMin: median = %.3e, max = %.3e\n", medRes, maxRes);

%% 4) Certificates post hoc computation (smaller is better):
fprintf('cert.resMin: median = %.3e, max = %.3e', median(cert.resMin), max(cert.resMin));%% 4) Certificates only (eigenvalue-only vs eigenpair residual)
rMin1  = leigqNEWTON_cert_resMin(A, lamR(1));
rPair1 = leigqNEWTON_cert_resPair(A, lamR(1), VR(:,1));
disp([rMin1, rPair1]);


%% Spherical eigenvalues (advanced)
% Some quaternion matrices have infinitely many left eigenvalues forming a sphere.
% Detecting/validating such spherical families is an advanced (and typically slower)
% workflow. For a dedicated, reproducible tutorial, see:
%
%   doc_SphereHunting
%
% That page demonstrates the sphere-sampling/validation pipeline based on
% LEIGQNEWTON_SPHERE_SAMPLE / _DETECT / _VALIDATE / _REFINE and discusses the
% speed–reliability trade-offs (the run time can range from tens of seconds to
% minutes depending on settings and hardware).
%
% Tip: from the Command Window you can run:
%   doc_SphereHunting
