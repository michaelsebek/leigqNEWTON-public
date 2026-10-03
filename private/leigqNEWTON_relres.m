function [eta, denominator] = leigqNEWTON_relres(raw, normA, absLambda, normV)
%LEIGQNEWTON_RELRES Internal scalar scale-invariant eigenpair residual.
% eta = raw / ((normA+absLambda)*normV), without any absolute floor.
% Invalid/nonfinite inputs and a zero vector produce Inf, never success.
% For A=lambda=0 and a nonzero vector, the exact zero defect has eta=0.
% Exponent splitting avoids intermediate overflow/underflow in the ratio.
eta = Inf;
denominator = NaN;
if ~all(cellfun(@(x)isnumeric(x)&&isreal(x)&&isscalar(x), ...
        {raw,normA,absLambda,normV}))
    error('leigq:BadResidualInput','Residual metric arguments must be real numeric scalars.');
end
if any(~isfinite([raw,normA,absLambda,normV])) || ...
        raw < 0 || normA < 0 || absLambda < 0 || normV <= 0
    return;
end
scale = max(normA,absLambda);
if scale == 0
    denominator = 0;
    if raw == 0, eta = 0; end
    return;
end
t = 1 + min(normA,absLambda)/scale;
[fs,es] = log2(scale);
[fv,ev] = log2(normV);
denominator = pow2(fs*fv*t,es+ev);
if raw == 0
    eta = 0;
else
    [fr,er] = log2(raw);
    eta = pow2((fr/fs)/fv/t,er-es-ev);
end
end
