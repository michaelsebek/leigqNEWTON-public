%% LEIGQNEWTON — API Reference (short map)
% Convert to Live Script: File → Save As… → Live Script (*.mlx)

%% Core solver
%   leigqNEWTON
help leigqNEWTON

%% Certification
%   leigqNEWTON_cert_resMin
%   leigqNEWTON_cert_resPair
help leigqNEWTON_cert_resMin
help leigqNEWTON_cert_resPair

%% Refinement
%   leigqNEWTON_init_vec
%   leigqNEWTON_refine_polish
%   leigqNEWTON_refine_lambda
%   leigqNEWTON_refine_auto
%   leigqNEWTON_refine_batch
help leigqNEWTON_refine_batch

%% Sphere module
%   leigqNEWTON_sphere_sample
%   leigqNEWTON_sphere_refine
help leigqNEWTON_sphere_refine

%% Legacy mapping (original names)
% leigqVANILA                   -> leigqNEWTON
% leigqVANILAresmin             -> leigqNEWTON_cert_resMin
% leigqVANILArespair            -> leigqNEWTON_cert_resPair
% leigqVANILA_nullvec           -> leigqNEWTON_init_vec
% leigqVANILA_polish            -> leigqNEWTON_refine_polish
% leigqVANILA_refineLambdaMinRes-> leigqNEWTON_refine_lambda
% leigqVANILA_refineAuto(_v2)   -> leigqNEWTON_refine_auto
% leigqVANILA_refineLambdas     -> leigqNEWTON_refine_batch
% leigqVANILAsphEXP2            -> leigqNEWTON_sphere_sample
% refineSPHERE                  -> leigqNEWTON_sphere_refine
