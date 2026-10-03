function G = leigqNEWTON_polish_coupling(va,vb,vc,vd)
% Coupling d(lambda*v)/d(lambda), component-stacked real coordinates.
% For fixed v this is RIGHT multiplication of deltaLambda by each v_i.
G = [va, -vb, -vc, -vd;
     vb,  va,  vd, -vc;
     vc, -vd,  va,  vb;
     vd,  vc, -vb,  va];
end
