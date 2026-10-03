function report = test_leigqNEWTON_revision()
%TEST_LEIGQNEWTON_REVISION Targeted regression tests for the LAA relative-residual patch.
% report = test_leigqNEWTON_revision;
% No historical benchmark, no large matrices, no third-party quaternion toolbox.
% Run in a fresh MATLAB session after replacing the core. Original wrappers stay untouched.
% The test requires the main solver and certificates from THIS directory on the path.
root = fileparts(mfilename('fullpath'));
expected = {'leigqNEWTON','leigqNEWTON_cert_resPair','leigqNEWTON_cert_resMin', ...
    'leigqNEWTON_refine_polish','leigqNEWTON_init_vec','checkNEWTON'};
for k=1:numel(expected)
    actual = which(expected{k});
    if isempty(actual) || ~strcmpi(fileparts(actual),root)
        error('leigq:test:Path','%s resolves outside the test directory: %s',expected{k},actual);
    end
end
if exist('quaternion','class') ~= 8
    error('leigq:test:Quaternion','The MATLAB quaternion class is required.');
end
rng0=rng; cleanup=onCleanup(@()rng(rng0)); %#ok<NASGU>
cases = { ...
 'relative scalar ratio and extreme finite denominators', @test_ratio; ...
 'scalar pair residual and quaternion order', @test_scalar_pair; ...
 'pair scale invariance (1e-100 to 1e100)', @test_pair_scaling; ...
 'minimal residual normalization (SVD)', @test_minimal; ...
 'minimal residual iterative SVD/fallback', @test_minimal_iterative; ...
 'zero matrix certificates', @test_zero_certificates; ...
 'zero vectors and mismatched lists rejected', @test_bad_vectors; ...
 'nonfinite input and invalid tolerances rejected', @test_invalid; ...
 'reviewer tiny-matrix false acceptance prevented', @test_tiny; ...
 'explicit absolute compatibility criterion', @test_absolute; ...
 'reporting option does not alter convergence', @test_output_mode; ...
 'scaled Newton correction with prescribed start', @test_newton_scaling; ...
 'RefineV preserves accepted certificates', @test_refine_vector; ...
 'complex nullspace conversion; early return', @test_null_one; ...
 'two right-H independent null vectors', @test_null_two; ...
 'zero matrix with diagonal shortcut disabled', @test_zero_core; ...
 'small nonzero diagonal entries are not zeros', @test_diagonal_rank; ...
 'approximate diagonals cannot use exact shortcut', @test_near_diagonal; ...
 'zero branch cannot bypass common acceptance', @test_null_safety; ...
 'explicit trial cap is respected', @test_budget; ...
 'seeded call restores RNG state', @test_rng; ...
 'numeric complex input maps to (1,i) slice', @test_complex; ...
 'polish Jacobian uses right multiplication', @test_coupling; ...
 'polish scaled noncommuting example', @test_polish; ...
 'polish exact zero stops without a singular solve', @test_polish_zero; ...
 'small polish step is not convergence', @test_polish_stagnation; ...
 'unsupported right gauge is rejected', @test_right_guard; ...
 'checkNEWTON normalization for nonunit vectors', @test_check; ...
 'info histories distinguish raw/relative residuals', @test_info ...
};
report=struct();
report.origin='Native MATLAB targeted regression, not a historical benchmark';
report.coreRevision='LAA-R1-relative-2026-10-03';
report.matlabVersion=version;
report.computer=computer;
report.date=datestr(now,30);
report.root=root;
report.sources=sourceHashes(root,expected);
report.checks=repmat(struct('name','','ok',false,'seconds',0,'message',''),size(cases,1),1);
t0=tic;
fprintf('leigqNEWTON relative-residual regression\nSource: %s\n',root);
for k=1:size(cases,1)
    t=tic; ok=false; msg='';
    try
        cases{k,2}(); ok=true;
    catch ME
        msg=sprintf('%s: %s',ME.identifier,ME.message);
        if ~isempty(ME.stack), msg=sprintf('%s (line %d in %s)',msg,ME.stack(1).line,ME.stack(1).name); end
    end
    report.checks(k)=struct('name',cases{k,1},'ok',ok,'seconds',toc(t),'message',msg);
    if ok, tag='OK'; else, tag='FAIL'; end
    fprintf('%2d/%2d %-4s %s\n',k,size(cases,1),tag,cases{k,1});
    if ~ok, fprintf('       %s\n',msg); end
end
report.elapsedSeconds=toc(t0);
report.total=numel(report.checks);
report.passed=sum([report.checks.ok]);
report.failed=report.total-report.passed;
if report.failed==0, report.status='OK'; else, report.status='FAILED'; end
fprintf('OVERALL: %s (%d/%d), %.2f s\n',report.status,report.passed,report.total,report.elapsedSeconds);
if nargout==0 && report.failed>0
    error('leigq:test:Failed','%d regression groups failed.',report.failed);
end
end

function test_ratio()
assert(leigqNEWTON_relres(1,1,0,1)==1);
assert(leigqNEWTON_relres(0,0,0,1)==0);
assert(isinf(leigqNEWTON_relres(0,0,0,0)));
assert(isinf(leigqNEWTON_relres(NaN,1,1,1)));
near(leigqNEWTON_relres(1e308,1e308,1e308,1),0.5,1e-14);
near(leigqNEWTON_relres(1e308,1e308,0,1e308),1e-308,1e-14);
near(leigqNEWTON_relres(1e-308,1e-308,0,1e-308),1e308,1e-14);
end

function test_scalar_pair()
A=q(1,2,3,4); l=q(-2,1,2,-1); v=q(2,-3,1,4);
[eta,raw,~,def,rep]=leigqNEWTON_cert_resPair(A,l,v,'NormalizeV',false,'ResidualNormalized',true);
w=prodq(A,v)-prodq(l,v);
near(raw,qn(w),1e-14); near(qn(def-w),0,1e-12);
near(eta,qn(A-l)/(qn(A)+qn(l)),1e-14);
near(rep.normA2,qn(A),1e-14);
end

function test_pair_scaling()
A=q([1,2;-1,3],[0,1;2,0],[1,-1;0,2],[0,2;1,-1]);
l=q(.7,.2,-.4,.3); v=q([1;2],[.2;-.5],[.3;.7],[-.1;.6]);
e=leigqNEWTON_cert_resPair(A,l,v,'NormalizeV',false,'ResidualNormalized',true);
for s=[1e-100,1e-16,1,1e16,1e100]
    ee=leigqNEWTON_cert_resPair(scale(A,s),scale(l,s),v,'NormalizeV',false,'ResidualNormalized',true);
    near(ee,e,2e-13);
end
for s=[1e-100,1e100]
    ee=leigqNEWTON_cert_resPair(A,l,scale(v,s),'NormalizeV',false,'ResidualNormalized',true);
    near(ee,e,2e-13);
end
end

function test_minimal()
for s=[1e-100,1e-16,1,1e16,1e100]
    A=q(s*diag([2,5])); l=q(s);
    [e,v,~,info]=leigqNEWTON_cert_resMin(A,l,'ResidualNormalized',true);
    near(e,1/6,2e-13); near(info.resMinRaw,s,2e-13);
    ep=leigqNEWTON_cert_resPair(A,l,v,'ResidualNormalized',true);
    near(e,ep,2e-13);
end
end

function test_minimal_iterative()
[e,v]=leigqNEWTON_cert_resMin(q(diag([2,5])),q(1),'Method','svds','ResidualNormalized',true);
near(e,1/6,1e-10);
near(leigqNEWTON_cert_resPair(q(diag([2,5])),q(1),v,'ResidualNormalized',true),e,1e-10);
end

function test_zero_certificates()
A=q(zeros(2)); v=q([1;2]);
assert(leigqNEWTON_cert_resPair(A,q(0),v,'ResidualNormalized',true)==0);
assert(leigqNEWTON_cert_resMin(A,q(0),'ResidualNormalized',true)==0);
end

function test_bad_vectors()
throws(@()leigqNEWTON_cert_resPair(q(1),q(1),q(0)),'leigqNEWTON_cert_resPair:BadVector');
throws(@()leigqNEWTON_cert_resPair(q(1),q(1),[]),'leigqNEWTON_cert_resPair:BadVector');
throws(@()leigqNEWTON_cert_resPair(q(1),q([1;2]),q(1)),'leigqNEWTON_cert_resPair:BadInput');
[r,raw,v,~,info]=leigqNEWTON_cert_resPair(q(eye(2)),[],[]);
assert(isempty(r)&&isempty(raw)&&isempty(v)&&info.K==0);
end

function test_invalid()
throws(@()leigqNEWTON(q(NaN)),'leigq:NonfiniteInput');
throws(@()leigqNEWTON(q(Inf)),'leigq:NonfiniteInput');
throws(@()leigqNEWTON(zeros(0)),'leigq:BadInput');
throws(@()leigqNEWTON(q(1),'Tol',NaN),'leigq:BadTol');
throws(@()leigqNEWTON(q(1),'Tol',Inf),'leigq:BadTol');
throws(@()leigqNEWTON(q(1),'Tol',-1),'leigq:BadTol');
throws(@()leigqNEWTON(q(1),'ToleranceMode','other'),'leigq:BadToleranceMode');
throws(@()leigqNEWTON(q(1),'TriangularShortcut','off','V0',q(0),'Trials',1),'leigq:BadV0');
throws(@()leigqNEWTON_cert_resMin(q(1),q(Inf)),'leigqNEWTON_cert_resMin:NonfiniteInput');
end

function args=tiny_args()
args={'Num',1,'Trials',1,'MaxIter',1,'Tol',1e-10,'Lambda0',q(0),'V0',q(1), ...
 'Damping',.05,'RefineV',false,'TriangularShortcut','off','TriangularInit',false,'AutoExtend',false};
end
function test_tiny()
args=tiny_args();
for normalized=[false,true]
    [l,~,~,i]=leigqNEWTON(q(1e-16),args{:},'ResidualNormalized',normalized);
    assert(isempty(l)&&i{1}.nAccepted==0);
end
assert(leigqNEWTON_cert_resPair(q(1e-16),q(0),q(1),'ResidualNormalized',true)==1);
end
function test_absolute()
args=tiny_args();
[l,~,~,i]=leigqNEWTON(q(1e-16),args{:},'ToleranceMode','absolute');
assert(numel(l)==1&&qn(l)==0&&strcmp(i{1}.toleranceMode,'absolute'));
assert(i{2}.rFinalRelative==1);
end

function [A,l,v]=hs_example()
A=q([0,1;1,0],[0,1;-1,0]);
l=q(sqrt(2));
v=q([.5;1/sqrt(2)],[.5;0]);
end
function args=newton_args(l,v)
args={'Num',1,'Trials',1,'MaxIter',30,'Tol',1e-12,'Lambda0',l,'V0',v, ...
 'RefineV',false,'TriangularShortcut','off','TriangularInit',false,'InfoLevel','full','AutoExtend',false};
end
function test_output_mode()
[A,l,v]=hs_example(); s=1e5; A=scale(A,s); l=scale(l,.9*s);
v=v+q([.01;-.02],[.02;.01],[.01;-.03],[.02;.01]);
a=newton_args(l,v);
[l1,V1,r1]=leigqNEWTON(A,a{:},'ResidualNormalized',true);
[l2,V2,r2]=leigqNEWTON(A,a{:},'ResidualNormalized',false);
assert(numel(l1)==1&&numel(l2)==1); assert(qn(l1-l2)==0&&qn(V1-V2)==0);
near(r1,leigqNEWTON_cert_resPair(A,l1,V1,'NormalizeV',false,'ResidualNormalized',true),1e-10);
near(r2,leigqNEWTON_cert_resPair(A,l2,V2,'NormalizeV',false,'ResidualNormalized',false),1e-10);
end
function test_newton_scaling()
[A,l,v]=hs_example();
v=v+q([.01;-.02],[.02;.01],[.01;-.03],[.02;.01]);
for s=[1e-100,1e-16,1,1e16,1e100]
    a=newton_args(scale(l,.9*s),v);
    [lr,vr,~,info]=leigqNEWTON(scale(A,s),a{:});
    assert(numel(lr)==1,'Scaled prescribed-start trial did not converge.');
    eta=leigqNEWTON_cert_resPair(scale(A,s),lr,vr,'ResidualNormalized',true);
    assert(isfinite(eta)&&eta<=1.01e-12);
    near(qn(scale(lr,1/s)-l),0,1e-9);
    assert(strcmp(info{1}.coreRevision,'LAA-R1-relative-2026-10-03'));
end
end
function test_refine_vector()
[A,l,v]=hs_example(); a=newton_args(scale(l,.9),v+q([.01;-.02],[.02;.01],[.01;-.03],[.02;.01]));
[lr,vr,r]=leigqNEWTON(A,a{:},'RefineV',true);
assert(numel(lr)==1&&r<=1e-12);
e=leigqNEWTON_cert_resPair(A,lr,vr,'NormalizeV',false,'ResidualNormalized',true);
near(r,e,1e-9);
end

function A=null_example(n)
W=zeros(n);W(1,1)=1;Y=zeros(n);Y(1,2)=1;Z=zeros(n);
if n>2, Z(1,3)=1; end
A=q(W,zeros(n),Y,Z);
end
function test_null_one()
A=null_example(2);
[l,V,r,i]=leigqNEWTON(A,'Num',1,'Trials',1,'UseNullFallbackLA',false,'Tol',1e-12);
assert(numel(l)==1&&qn(l)==0&&i{1}.nAcceptedZero==1&&i{1}.nRuns==0);
assert(r<=1e-12&&qn(prodq(A,V))<1e-12);
near(i{2}.rFinalRaw,qn(prodq(A,V)),1e-10);
end
function test_null_two()
A=null_example(3);
[l,V,r,i,lU]=leigqNEWTON(A,'Num',2,'Trials',1,'UseNullFallbackLA',false, ...
 'VerifyZeroNull',false,'Tol',1e-12);
assert(numel(l)==2&&numel(lU)==1&&i{1}.nAcceptedZero==2&&i{1}.nRuns==0);
assert(all(r<=1e-12)&&qn(prodq(A,V))<1e-12);
assert(rank(cembed(V),1e-10)==4,'Returned null vectors are not right-H independent.');
end
function test_zero_core()
[l,V,r,i,lU]=leigqNEWTON(q(zeros(3)),'TriangularShortcut','off','Num',3);
assert(numel(l)==3&&numel(lU)==1&&all(r==0)&&i{1}.nAcceptedZero==3&&i{1}.nRuns==0);
assert(rank(cembed(V),1e-10)==6);
end
function test_diagonal_rank()
[l,~,r,i]=leigqNEWTON(q(diag([1e-16,2e-16])));
assert(numel(l)==2&&all(r==0)&&i{1}.rankA==2&&i{1}.df0==0&&i{1}.nAcceptedZero==0);
end
function test_near_diagonal()
[l,~,~,i]=leigqNEWTON(q([1,1e-6;1e-6,2]),'Num',1,'Trials',1,'MaxIter',1, ...
 'Damping',.05,'Tol',1e-14,'TriTol',1,'Lambda0',q(1),'V0',q([1;0]),'RefineV',false);
assert(isempty(l)&&i{1}.nRuns==1);
end
function test_null_safety()
[l,~,~,i]=leigqNEWTON(q([1,1e-14;0,1e-14]),'Num',1,'Trials',1,'MaxIter',1, ...
 'Tol',1e-18,'RankTolFactor',100,'VerifyZeroNull',false,'UseNullFallbackLA',false, ...
 'ZeroNullTol',1e99,'Damping',.05,'Lambda0',q(0),'V0',q([1;1]),'RefineV',false);
assert(i{1}.nAcceptedZero==0);
if ~isempty(l), assert(qn(l)>0); end
end
function test_budget()
[~,~,~,i]=leigqNEWTON(q([1,.2;.1,2]),'Num',2,'Trials',1,'MaxTrials',1, ...
 'MaxIter',1,'Tol',1e-14,'Damping',.05,'Lambda0',q(0),'V0',q([1;1]),'RefineV',false);
assert(i{1}.nRuns<=1&&i{1}.maxTrials==1);
end
function test_rng()
rng(726,'twister'); before=rng;
leigqNEWTON(q([1,.2;.1,2]),'Seed',27,'Num',1,'Trials',1,'MaxIter',2,'RefineV',false);
after=rng; assert(isequal(before,after));
end
function test_complex()
[l,V,r]=leigqNEWTON([1+2*1i,0;0,-2+3*1i]);
assert(all(r==0)&&qn(l-q([1;-2],[2;3]))==0);
e=leigqNEWTON_cert_resPair([1+2*1i,0;0,-2+3*1i],l,V,'ResidualNormalized',true);
assert(all(e==0));
end

function test_coupling()
v=q([1;2],[-3;4],[5;6],[7;-8]); dl=q(2,-3,5,7);
[a,b,c,d]=parts(v); G=leigqNEWTON_polish_coupling(a,b,c,d);
[e,f,g,h]=parts(dl); expected=prodq(dl,v); [a,b,c,d]=parts(expected);
assert(norm(G*[e;f;g;h]-[a;b;c;d])==0,'Wrong side in polish derivative.');
end
function test_polish()
[A,l,v]=hs_example();
v=v+q([.01;-.02],[.02;.01],[.01;-.03],[.02;.01]);
for s=[1e-16,1,1e16]
    l0=scale(l+q(-.1,0,.02,-.03),s);
    [lp,vp,r,info]=leigqNEWTON_refine_polish(scale(A,s),l0,v, ...
        'TolRes',1e-12,'TolStep',1e-15,'MaxIter',20);
    assert(isstruct(r)&&info.converged&&r.converged&&r.relative<=1e-12);
    e=leigqNEWTON_cert_resPair(scale(A,s),lp,vp,'ResidualNormalized',true);
    assert(e<=1.01e-12); near(qn(scale(lp,1/s)-l),0,1e-8);
end
end
function test_polish_zero()
[l,v,r,i]=leigqNEWTON_refine_polish(q(zeros(2)),q(0),q([1;0]));
assert(i.iter==0&&i.converged&&r.relative==0&&qn(l)==0&&qn(v)>0);
end
function test_polish_stagnation()
[A,l,v]=hs_example();
[~,~,r,i]=leigqNEWTON_refine_polish(A,scale(l,.5),v+q([.1;-.2]), ...
 'TolRes',0,'TolStep',1e6,'MaxIter',2);
assert(i.converged==(r.relative==0));
assert(~i.converged,'The test must remain inexact after one damped Newton step.');
assert(strcmp(i.reason,'step stagnation'));
end
function test_right_guard()
throws(@()leigqNEWTON_refine_polish(q(1),q(1),q(1),'Side','right'),'leigq:PolishSide');
throws(@()leigqNEWTON_init_vec(q(1),q(1),'Side','right'),'leigq:RightGauge');
end
function test_check()
out=checkNEWTON(q(2),q(1),q(3),'Verbose',0,'SphereCheck','off');
near(out.resMinAbs,1,1e-14); near(out.resPairAbs,3,1e-14);
near(out.resMinRel,1/3,1e-14); near(out.resPairRel,1/3,1e-14);
end
function test_info()
[A,l,v]=hs_example();a=newton_args(scale(l,.9),v);
[lr,vr,r,i]=leigqNEWTON(A,a{:}); assert(numel(lr)==1);
run=i{2}; assert(isfield(run,'rFinalRaw')&&isfield(run,'rFinalRelative'));
assert(numel(run.resHistRaw)==numel(run.resHistRelative));
near(run.rFinalRelative,r,1e-10);
raw=leigqNEWTON_cert_resPair(A,lr,vr,'NormalizeV',false);
near(run.rFinalRaw,raw,1e-10);
end

function Q=q(a,b,c,d)
if nargin<2, b=zeros(size(a)); end
if nargin<3, c=zeros(size(a)); end
if nargin<4, d=zeros(size(a)); end
Q=quaternion(a,b,c,d);
end
function B=scale(A,s)
[a,b,c,d]=parts(A); B=q(s*a,s*b,s*c,s*d);
end
function n=qn(A)
[a,b,c,d]=parts(A); n=norm([a(:);b(:);c(:);d(:)]);
end
function C=prodq(A,B)
% Independent component multiplication, not qmtimesNEWTON.
[a,b,c,d]=parts(A);[e,f,g,h]=parts(B);
if isscalar(A) || isscalar(B)
    C=q(a.*e-b.*f-c.*g-d.*h,a.*f+b.*e+c.*h-d.*g, ...
        a.*g-b.*h+c.*e+d.*f,a.*h+b.*g-c.*f+d.*e);
else
    C=q(a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e);
end
end
function X=cembed(A)
[a,b,c,d]=parts(A); u=a+1i*b; v=c+1i*d; X=[u,v;-conj(v),conj(u)];
end
function near(actual,expected,tol)
% Use a relative check except for an exact-zero target.
if expected==0, err=abs(actual); else, err=abs((actual-expected)/expected); end
assert(isfinite(err)&&err<=tol,sprintf('actual %.17g, expected %.17g, error %.3g',actual,expected,err));
end
function throws(fun,id)
try
    fun();
catch ME
    assert(strcmp(ME.identifier,id),sprintf('Expected %s, got %s: %s',id,ME.identifier,ME.message));
    return;
end
error('leigq:test:ExpectedError','Expected exception %s was not raised.',id);
end
function out=sourceHashes(root,names)
out=struct();
for k=1:numel(names)
    f=fullfile(root,[names{k},'.m']);
    fid=fopen(f,'rb'); assert(fid>=0); bytes=fread(fid,Inf,'*uint8'); fclose(fid);
    try
        md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
        md.update(typecast(bytes,'int8')); h=typecast(md.digest(),'uint8');
        out.(names{k})=lower(reshape(dec2hex(h,2).',1,[]));
    catch
        out.(names{k})='SHA-256 unavailable (no JVM); source path recorded';
    end
end
end
