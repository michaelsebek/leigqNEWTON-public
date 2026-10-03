function [lambda, v, res, info] = leigqNEWTON_refine_polish(A, lambda0, v0, varargin)
%LEIGQNEWTON_REFINE_POLISH Local Newton polish of A*v=lambda*v (LEFT).
% [lambda,v,res,info] = leigqNEWTON_refine_polish(A,lambda0,v0,Name,Value,...)
% Empty v0 selects a smallest right singular vector via leigqNEWTON_init_vec.
% TolRes (alias Tol), default 1e-14, applies to the relative TWO-norm residual
%   eta = norm(A*v-lambda*v)/((norm(A,2)+abs(lambda))*norm(v)).
% ToleranceMode='absolute' explicitly selects a raw TWO-norm stopping test.
% TolStep (default 1e-14) detects stagnation, never certifies convergence.
% MaxIter=20, Damping=true (alias Backtrack), Verbose=0; Side='left' only.
% ResidualNormalized=false only controls the added res.returned field.
% res remains a STRUCT with the existing raw .resInf and .res2 fields;
% .relative, .converged, .returned and .toleranceMode are added.
% info.converged is true only if the final residual passes the requested test.
% Core revision: LAA-R1-relative-2026-10-03. Not a historical benchmark rerun.

opt = struct('Side','left','TolRes',1e-14,'TolStep',1e-14,'MaxIter',20, ...
    'Damping',true,'Verbose',0,'ToleranceMode','relative','ResidualNormalized',false);
opt = parseOpts(opt,varargin{:});
if ~strcmpi(opt.Side,'left')
    error('leigq:PolishSide','This routine supports left eigenpairs only.');
end
opt.ToleranceMode = lower(char(opt.ToleranceMode));
if ~any(strcmp(opt.ToleranceMode,{'relative','absolute'}))
    error('leigq:BadToleranceMode','ToleranceMode must be relative or absolute.');
end
for name={'TolRes','TolStep','MaxIter'}
    value=opt.(name{1});
    if ~(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=0)
        error('leigq:BadOption','%s must be finite and nonnegative.',name{1});
    end
end
opt.MaxIter = floor(opt.MaxIter);
if isnumeric(A), A=quaternion(real(A),imag(A),zeros(size(A)),zeros(size(A))); end
if isnumeric(lambda0), lambda0=quaternion(real(lambda0),imag(lambda0),zeros(size(lambda0)),zeros(size(lambda0))); end
n = size(A,1);
if ndims(A)~=2 || n==0 || size(A,2)~=n || ~isscalar(lambda0)
    error('leigq:BadInput','A must be nonempty square and lambda0 scalar.');
end
AR = qA_real_left(A);
[a,b,c,d]=parts(lambda0);
if any(~isfinite(AR(:))) || any(~isfinite([a;b;c;d]))
    error('leigq:NonfiniteInput','A and lambda0 must be finite.');
end
normA = norm(AR,2);
if ~isfinite(normA), error('leigq:ScaleOverflow','Rescale A.'); end
if isempty(v0)
    [v0,ninfo] = leigqNEWTON_init_vec(A,lambda0,'Side','left');
    p=ninfo.pivot;
else
    if isnumeric(v0), v0=quaternion(real(v0),imag(v0),zeros(size(v0)),zeros(size(v0))); end
    v0=v0(:);
    if numel(v0)~=n, error('leigq:BadV0','v0 must have n entries.'); end
    [va,vb,vc,vd]=parts(v0);
    vv=[va;vb;vc;vd];
    if any(~isfinite(vv)) || norm(vv)==0 || ~isfinite(norm(vv))
        error('leigq:BadV0','v0 must be finite and nonzero.');
    end
    [~,p]=max(hypot(hypot(va,vb),hypot(vc,vd)));
end
lambda=lambda0;
v=normalizeV(gaugeFix(v0,p,'left'));
[rInf,rvec]=residualInf(AR,lambda,v,'left');
eta=relativeMetric(rvec,normA,lambda,v);
hist=struct('resInf',rInf,'resRelative',eta,'stepInf',zeros(0,1),'alpha',zeros(0,1));
iter=0; reason='iteration limit';
for it=1:opt.MaxIter
    if accepted(rvec,eta,opt), reason='residual'; break; end
    [M,G]=buildMG(AR,lambda,v,'left');
    [va,vb,vc,vd]=parts(v); vvec=[va;vb;vc;vd];
    [a,b,c,d]=parts(lambda); lambdaScale=max(normA,norm([a;b;c;d]));
    if lambdaScale==0, lambdaScale=1; end
    cvec=[vb(p);vc(p);vd(p);vvec.'*vvec-1];
    C=zeros(4,4*n); C(1,n+p)=1; C(2,2*n+p)=1; C(3,3*n+p)=1; C(4,:)=2*vvec.';
    J=[M/lambdaScale,-G; C,zeros(4,4)];
    dx= -J\[rvec/lambdaScale;cvec];
    if any(~isfinite(dx)), reason='nonfinite correction'; break; end
    dv=dx(1:4*n); dl=lambdaScale*dx(4*n+1:end);
    stepInf=norm(dx,inf); alpha=1; found=false;
    while alpha>=1/64
        ln=lambda+quaternion(alpha*dl(1),alpha*dl(2),alpha*dl(3),alpha*dl(4));
        vn=normalizeV(gaugeFix(addToV(v,alpha*dv,n),p,'left'));
        [rn,rvn]=residualInf(AR,ln,vn,'left');
        en=relativeMetric(rvn,normA,ln,vn);
        if isfinite(en) && isfinite(rn) && (~opt.Damping || norm(rvn)<=norm(rvec))
            found=true; break;
        end
        alpha=alpha/2;
    end
    if ~found, reason='line search'; break; end
    lambda=ln; v=vn; rInf=rn; rvec=rvn; eta=en; iter=it;
    hist.resInf(end+1,1)=rInf; hist.resRelative(end+1,1)=eta;
    hist.stepInf(end+1,1)=stepInf; hist.alpha(end+1,1)=alpha;
    if opt.Verbose, fprintf('polish: iter %d  eta=%.3e  raw2=%.3e\n',it,eta,norm(rvec)); end
    if accepted(rvec,eta,opt), reason='residual'; break; end
    if alpha*stepInf<=opt.TolStep, reason='step stagnation'; break; end
end
converged=accepted(rvec,eta,opt);
if converged, reason='residual'; end
res=struct('resInf',rInf,'res2',norm(rvec),'relative',eta, ...
    'converged',converged,'toleranceMode',opt.ToleranceMode);
if opt.ResidualNormalized, res.returned=eta; else, res.returned=res.res2; end
info=struct('iter',iter,'pivot',p,'hist',hist,'side','left','converged',converged, ...
    'reason',reason,'normA2',normA,'coreRevision','LAA-R1-relative-2026-10-03');
end

function eta=relativeMetric(rvec,normA,lambda,v)
[a,b,c,d]=parts(lambda); [va,vb,vc,vd]=parts(v);
eta=leigqNEWTON_relres(norm(rvec),normA,norm([a;b;c;d]),norm([va;vb;vc;vd]));
end

function ok=accepted(rvec,eta,opt)
ok=isfinite(eta)&&all(isfinite(rvec));
if strcmp(opt.ToleranceMode,'relative'), ok=ok&&(eta<=opt.TolRes);
else, ok=ok&&(norm(rvec)<=opt.TolRes); end
end

% ---------------- helpers ----------------
function [rInf, rvec] = residualInf(AR, lambda, v, side)
n = size(AR,1)/4;
[va,vb,vc,vd] = parts(v);
vvec = [va;vb;vc;vd];

[a,b,c,d] = parts(lambda);
q = [a;b;c;d];
if strcmpi(side,'left')
    L = qLmat(q);
else
    L = qRmat(q);
end
M = AR - kron(L, eye(n));
rvec = M * vvec;

ra = rvec(1:n); rb = rvec(n+1:2*n); rc = rvec(2*n+1:3*n); rd = rvec(3*n+1:4*n);
rInf = max(hypot(hypot(ra,rb),hypot(rc,rd)));
end

function [M,G] = buildMG(AR, lambda, v, side) %#ok<INUSD>
n = size(AR,1)/4;
[va,vb,vc,vd] = parts(v);
[a,b,c,d] = parts(lambda);
M = AR - kron(qLmat([a;b;c;d]), eye(n));
% d(lambda*v)/d(lambda) is RIGHT multiplication by each entry of v.
G = leigqNEWTON_polish_coupling(va,vb,vc,vd);
end

function v = addToV(v, dv, n)
[va,vb,vc,vd] = parts(v);
va = va + dv(1:n);
vb = vb + dv(n+1:2*n);
vc = vc + dv(2*n+1:3*n);
vd = vd + dv(3*n+1:4*n);
v = quaternion(va,vb,vc,vd);
end

function v = normalizeV(v)
[va,vb,vc,vd] = parts(v);
nv = norm([va;vb;vc;vd]);
if nv==0, return; end
v = local_qscale(v, 1/nv);
end

function v = gaugeFix(v,p,side)
vp = v(p);

[a,b,c,d] = parts(vp);
m = norm([a;b;c;d]);
if m==0, return; end

q = quaternion(a/m,-b/m,-c/m,-d/m);   % conj(vp)/|vp|, avoids conj/rdivide overloads

v = qmtimesNEWTON(v, q); % right gauge preserves A*v=lambda*v

[ar,~,~,~] = parts(v(p));
if ar < 0
    v = -v;
end
end

function AR = qA_real_left(A)
[A0,A1,A2,A3] = parts(A);
AR = [ A0, -A1, -A2, -A3;
       A1,  A0, -A3,  A2;
       A2,  A3,  A0, -A1;
       A3, -A2,  A1,  A0 ];
end

function L = qLmat(q)
a=q(1); b=q(2); c=q(3); d=q(4);
L = [ a, -b, -c, -d;
      b,  a, -d,  c;
      c,  d,  a, -b;
      d, -c,  b,  a ];
end

function R = qRmat(q)
a=q(1); b=q(2); c=q(3); d=q(4);
R = [ a, -b, -c, -d;
      b,  a,  d, -c;
      c, -d,  a,  b;
      d,  c, -b,  a ];
end

function opt = parseOpts(opt, varargin)
if mod(numel(varargin),2)~=0, error('leigq:BadArgs','Name-value pairs expected.'); end
fields = fieldnames(opt);
for k=1:2:numel(varargin)
    key = lower(char(varargin{k}));
    if strcmp(key,'tol'), key='tolres'; end
    if strcmp(key,'backtrack'), key='damping'; end
    j = find(strcmpi(fields,key),1);
    if isempty(j), error('leigq:BadOption','Unknown polish option: %s',key); end
    opt.(fields{j}) = varargin{k+1};
end
end

function v = local_qscale(v, s)
% Scale quaternion array by a real scalar without relying on quaternion RDIVIDE/TIMES.
[w,x,y,z] = parts(v);
v = quaternion(s*w, s*x, s*y, s*z);
end

