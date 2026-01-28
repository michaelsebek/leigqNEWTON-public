function [lambda, V, res, info, lambdaU, VU, resU] = leigqNEWTON(A, varargin)
%LEIGQNEWTON  Left eigenpairs of a quaternionic matrix via a Newton-type solver (stand-alone).
%
%   [lambda, V, res] = leigqNEWTON(A)
%   [lambda, V, res, info] = leigqNEWTON(A, Name,Value,...)
%   [lambda, V, res, info, lambdaU, VU, resU] = leigqNEWTON(A, Name,Value,...)
%   [...] = leigqNEWTON(A, Num, Name,Value,...)          %% convenience positional Num
%
% Computes quaternionic LEFT eigenpairs (lambda, V) of a square quaternion matrix A:
%       A*v = lambda*v,   v ~= 0,
% where lambda is a quaternion acting on the LEFT.
%
% This implementation is self-contained and uses only:
%   - MATLAB built-in quaternion class (Aerospace Toolbox),
%   - standard MATLAB linear algebra on REAL/COMPLEX embeddings built locally.
%
% -------------------------------------------------------------------------
% Outputs (types and meaning)
%
%   lambda   : K-by-1 quaternion
%       The K eigenvalue "hits" accepted by the solver (K is controlled by Num).
%       These may contain duplicates (the solver does not enforce uniqueness).
%
%   V        : n-by-K quaternion
%       Corresponding right eigenvectors as columns (A*V(:,k)  approx  lambda(k)*V(:,k)).
%
%   res      : K-by-1 double
%       Final residual norms reported for each accepted eigenpair.
%       If ResidualNormalized=true (default), res is scale-invariant:
%           res = ||A*v - lambda*v|| / den(A,lambda,v),
%       where den(...) is a mild scale factor (see local_residual_den).
%       If ResidualNormalized=false, res is the raw Euclidean norm.
%
%   info     : cell array (only if requested)
%       info{1} is a summary struct; info{2:end} are per-trial structs when InfoLevel='full'.
%       In addition, info{1}.distinct describes how lambdaU was formed (if requested).
%
%   lambdaU  : Ku-by-1 quaternion  (only if requested)
%       A DISTINCT representative set extracted from lambda using a tolerance.
%       Each group of near-identical hits contributes one representative, chosen by:
%         (i) significantly better residual (if present), otherwise
%         (ii) "prettiness" (prefers exact integers/zeros over small numerical tails).
%       By design, lambdaU is a subset of the returned lambda (representatives are picked
%       from existing hits; no additional cleaning is applied).
%
%   VU       : n-by-Ku quaternion  (only if requested)
%   resU     : Ku-by-1 double      (only if requested)
%       The vectors/residuals corresponding to lambdaU.
%
% -------------------------------------------------------------------------
% Profiles (recommended)
% LEIGQNEWTON is restart-based; success depends strongly on the TRIAL budget and
% Newton iterations per trial. SolveProfile presets apply only to options you did
% NOT specify explicitly:
%   'fast'     : small budgets; quick scans (may return fewer eigenpairs).
%   'default'  : balanced; typically succeeds for generic dense matrices.
%   'reliable' : larger budgets + auto-extension; aims to return (almost) all eigenpairs.
%
% -------------------------------------------------------------------------
% Options (Name,Value)  (aliases accepted; unknown options are ignored)
%
%   Desired count:
%     'Num'/'K'/'NEigs'/'NEV'         : number of eigenpairs requested (default n)
%     (positional) Num               : leigqNEWTON(A, 9, ...) sets Num=9
%
%   Trial budgets / restarts:
%     'Trials'/'Restarts'/'NTrials'  : initial trial budget (profile-dependent default)
%     'TrialsFactor'                : convenience Trials := TrialsFactor*n (if Trials not set)
%     'MaxTrials'/'MaxRestarts'      : hard cap on total trials (also disables AutoExtend)
%     'AutoExtend'                  : true/false (profile-dependent default)
%     'MaxTrialsCap'/'TrialsCap'    : cap used by AutoExtend
%     'ExtendBy'/'ExtendFactor'     : how much to increase budget when auto-extending
%
%   Newton / acceptance:
%     'MaxIter'/'MaxIt'              : Newton iterations per trial (profile-dependent default)
%     'Tol'/'ResTol'                 : convergence tolerance (raw residual is used internally)
%     'Damping'/'Alpha'/'StepSize'   : initial step size alpha in (0,1] (default 1)
%     'Backtrack'/'LineSearch'       : backtracking line search (default true)
%     'MinAlpha'/'AlphaMin'          : minimum alpha in backtracking (default 1/64)
%
%   Distinct representative extraction (for lambdaU):
%     'DistinctTol'/'DistinctTolAbs'/'UniqueTol' : absolute grouping tolerance (default max(1e-6, sqrt(Tol)))
%     'DistinctTolRel'                          : relative grouping tolerance (default 0)
%     'DistinctResFactor'                       : residual "order-of-magnitude" factor (default 10)
%       Note: The distinct-grouping tolerance is intentionally more tolerant than Tol by default.
%       You still see all raw hits in lambda; lambdaU only groups near-duplicates into representatives.
%
%   Reproducibility / logging:
%     'Seed'/'RNG'                   : RNG seed (default [])
%     'Verbose'/'Disp'               : 0/1 (default 0)
%     'InfoLevel'/'InfoMode'         : 'summary'|'full' OR numeric 0/1/2 (2='full')
%
%   Robustness / post-processing:
%     'VerifyZeroNull'/'ZeroNullVerify' : verify Aq*x approx 0 in the zero-eigenvalue pre-pass (default true)
%     'ZeroNullTol'                     : tolerance for the zero-eig pre-pass (raw units). [] -> auto
%     'UseNullFallbackLA'/'FallbackLA'  : if verification fails, recompute nullspace using real 4n embedding (default true)
%     'ResidualNormalized'              : if true, res output is scale-invariant (default true)
%     'RefineV'                         : recompute eigenvectors for final lambdas via real least-squares/SVD (default true)
%
%   Initialization / triangular handling:
%     'Lambda0'/'Lam0'/'InitLambda'  : initial lambda guess (scalar quaternion)
%     'V0'/'X0'/'InitVec'            : initial eigenvector guess (n-by-1 quaternion)
%     'TriangularInit'/'TriInit'     : if true, deterministic diagonal/basis seeding for (detected) triangular A (default true)
%     'TriTol'/'TriangularTol'       : triangular detection tolerance (default 0)
%     'IsTriangular'/'TriangularFlag'/'ForceTriangular' : user override for triangular flag (true/false)
%     'TriangularShortcut'/'TriShortcut' : 'diag' (default) shortcut for diagonal A; 'off' forces Newton
%
% Notes on triangular matrices:
%   If TriangularInit is true and A is treated as triangular, the solver uses
%   deterministic diagonal/basis seeding for the first n trials. This often leads
%   to very small iteration counts (sometimes 0), which is expected.
%   To benchmark "generic" random initialization on triangular A, set:
%       'TriangularInit', false.
%
% -------------------------------------------------------------------------
% Useful one-liners (copy/paste)
%   % Default (balanced):
%   [lam,V,res] = leigqNEWTON(A);
%
%   % Also return distinct representatives (lambdaU):
%   [lam,V,res,info,lamU,VU,resU] = leigqNEWTON(A,'SolveProfile','reliable');
%   % Add a larger distinct-grouping tolerance explicitly: 
%   [lam,V,res,info,lamU,VU,resU] = leigqNEWTON(A,'SolveProfile','reliable','DistinctTolAbs', 1e-5);  % or 3e-5, 1e-4 if you want stronger merging
%   % If you want it scale-aware (relative as well)
%   [lam,V,res,info,lamU,VU,resU] = leigqNEWTON(A, 'SolveProfile','reliable','DistinctTolAbs', 1e-5,'DistinctTolRel', 1e-8);
% 
%   % Rule of thumb: Start with DistinctTolAbs = 10*sqrt(Tol) (or 30*sqrt(Tol) for very defective/triangular-ish cases).
%     This affects only deduplication into lamU. You still see all raw hits in lam, so being more tolerant here is usually safe.
%
%   % Request exactly 9 eigenpairs (two equivalent forms):
%   [lam,V,res] = leigqNEWTON(A, 9, 'SolveProfile','reliable');
%   [lam,V,res] = leigqNEWTON(A,'Num',9,'SolveProfile','reliable');
%
%   % Reproducible run:
%   [lam,V,res] = leigqNEWTON(A,'Seed',1,'SolveProfile','default');
%
% -------------------------------------------------------------------------
% Legacy name mapping
%   leigqVANILA  ->  leigqNEWTON
%
% See also: quaternion, parts, null, rank, leigqNEWTON_refine_polish, leigqNEWTON_cert_resMin

if ~exist('quaternion','class')
    error('leigq:NoQuaternionClass', ...
        'MATLAB class "quaternion" not found. Requires Aerospace Toolbox.');
end

% ---------------- type normalization ----------------
Aq = local_to_quat(A);
n = size(Aq,1);
if n ~= size(Aq,2)
    error('leigq:BadInput','A must be square.');
end

% ---------------- parse options (profile-aware) ----------------
opt = local_parse_opts(n, varargin{:});
%if ~isempty(opt.Seed), rng(opt.Seed); end
% To remove side-effects, replaced by
rngCleanup = [];
if ~isempty(opt.Seed)
    s0 = rng;
    rngCleanup = onCleanup(@() rng(s0));
    rng(double(opt.Seed),'twister');     % also pins the generator
end


% ---------------- triangular detection (for seeding) ----------------
isTriDetected = local_is_triangular(Aq, opt.TriTol);
if ~isempty(opt.IsTriangularUser)
    opt.IsTriangular = logical(opt.IsTriangularUser);
else
    opt.IsTriangular = isTriDetected;
end

% ---------------- output intent ----------------
wantInfo     = (nargout >= 4);
wantDistinct = (nargout >= 5);

% ---------------- optional diagonal shortcut ----------------
% By default, diagonal matrices are handled without running Newton (can be disabled).
if strcmpi(opt.TriangularShortcut,'diag') && local_is_diagonal(Aq, opt.TriTol)
    [lambda, V, res, info, lambdaU, VU, resU] = local_diagonal_shortcut(Aq, opt, wantInfo, wantDistinct);
    return;
end

% ---------------- allocate outputs ----------------
lambda  = quaternion.empty(0,1);
V       = quaternion.empty(n,0);
res     = zeros(0,1);


lambdaU  = quaternion.empty(0,1);
VU       = quaternion.empty(n,0);
resU     = zeros(0,1);
info     = [];
% ---------------- info allocation ----------------
if wantInfo
    runlog = cell(0,1);
else
    runlog = [];
end

targetK    = opt.Num;
trial      = 0;
nConverged = 0;

% ---------------- precompute real left-multiplication block for A ----------------
LA = local_left_block(Aq);  % 4n x 4n (double)
In = eye(n);

% Precompute a stable scaling for normalized residuals (real embedding)
opt.Afro = norm(LA,'fro');

% ============================================================
%  ZERO-EIGENVALUE PRE-PASS (rank/null based)
% ============================================================
rankA        = NaN;
df0          = 0;
nAcceptedZero = 0;

if targetK > 0
    try
        rankA = local_rank_quat(Aq, opt);  % quaternion rank via complex embedding
    catch ME
        if opt.Verbose
            fprintf('leigq: rank(A) failed (%s); skipping zero-eig pre-pass.\n', ME.message);
        end
        rankA = NaN;
    end

    if isfinite(rankA)
        df0 = n - rankA;
        df0 = max(0, min(n, round(df0)));

        if df0 > 0
            % Compute a quaternion basis for the right nullspace of A (i.e., Aq*x = 0).
            % Primary method: complex embedding (fast). We VERIFY candidates and can FALLBACK
            % to the real 4n embedding if needed.
            Xq = quaternion.empty(n,0);

            % --- 1) Fast method (complex embedding) ---
            try
                Xq = local_null_quat(Aq, df0, opt); % n-by-<=df0 quaternion (best effort)
            catch ME
                if opt.Verbose
                    fprintf('leigq: null(A) (complex) failed (%s); will try LA fallback.\n', ME.message);
                end
                Xq = quaternion.empty(n,0);
            end

            % --- 2) Verify and filter the candidates ---
            if opt.VerifyZeroNull && ~isempty(Xq)
                [Xq, okMask] = local_verify_zero_nullspace(LA, In, Xq, opt);
                if opt.Verbose && any(~okMask)
                    fprintf('leigq: rejected %d/%d zero-null candidates (residual too large).\n', ...
                        sum(~okMask), numel(okMask));
                end
            end

            % --- 3) Fallback: real 4n embedding nullspace, then select df0 independent quaternion columns ---
            if (size(Xq,2) < df0) && opt.UseNullFallbackLA
                try
                    XqLA = local_null_quat_LA(Aq, df0, LA);
                    if opt.VerifyZeroNull && ~isempty(XqLA)
                        [XqLA, okMaskLA] = local_verify_zero_nullspace(LA, In, XqLA, opt);
                        if opt.Verbose && any(~okMaskLA)
                            fprintf('leigq: LA-null rejected %d/%d candidates (residual too large).\n', ...
                                sum(~okMaskLA), numel(okMaskLA));
                        end
                    end
                    if size(XqLA,2) > size(Xq,2)
                        Xq = XqLA;
                    end
                catch ME
                    if opt.Verbose
                        fprintf('leigq: null(A) (LA) failed (%s); skipping zero-eig pre-pass.\n', ME.message);
                    end
                end
            end

            if ~isempty(Xq)
                dfFound = size(Xq,2);
                dfTake  = min([df0, dfFound, targetK]);

                lam0R = zeros(4,1);
                lam0  = quaternion(0,0,0,0);
                for k0 = 1:dfTake
                    xq = Xq(:,k0);

                    % gauge + normalize (packed-real form, consistent with Newton phase)
                    xR = local_pack_qvec(xq);
                    [xR, ~] = local_gauge(xR);

                    % residual (raw + optional normalized output)
                    [~, rRaw, rOut] = local_residual_metrics(LA, In, lam0R, xR, opt);

                    % Safety: do NOT accept a bogus nullspace vector
                    if rRaw > local_zero_accept_tol(LA, xR, opt)
                        if opt.Verbose
                            fprintf('leigq: skipping a zero-eig candidate, rawRes=%.3e too large.\n', rRaw);
                        end
                        continue;
                    end

                    xq = local_unpack_qvec(xR);

                    lambda(end+1,1)  = lam0;   %#ok<AGROW>
                    V(:,end+1)       = xq;     %#ok<AGROW>
                    res(end+1,1)     = rOut;   %#ok<AGROW>

                    nAcceptedZero = nAcceptedZero + 1;

                    if wantInfo
                        runlog{end+1,1} = local_make_run_stub_zero(k0, lam0, xq, rOut); %#ok<AGROW>
                    end

                    if opt.Verbose
                        fprintf('leigq: accepted zero eigenvalue #%d/%d (df0=%d), res=%.3e\n', ...
                            k0, dfTake, df0, rOut);
                    end

                    if numel(lambda) >= targetK
                        if wantInfo
                            info = local_finalize_info(opt, trial, nConverged, targetK, maxTrials, ...
                                rankA, df0, nAcceptedZero, runlog, lambda);
                        else
                            info = [];
                        end
                        if wantDistinct
                            [lambdaU, VU, resU, metaD] = local_distinct_representatives(lambda, V, res, opt);
                            if wantInfo && ~isempty(info)
                                info{1}.summary.nDistinct = numel(lambdaU);
                                info{1}.distinct = metaD;
                            end
                        else
                            lambdaU = quaternion.empty(0,1);
                            VU      = quaternion.empty(n,0);
                            resU    = zeros(0,1);
                        end
                        return;
                    end
                end
            end
        end
    end
end

% ============================================================
%  NEWTON TRIALS FOR REMAINING EIGENPAIRS
% ============================================================

% Trial budget logic:
%  - opt.Trials is the initial budget.
%  - If AutoExtend is true and Trials was NOT user-specified, we can increase
%    the budget up to opt.MaxTrialsCap.
maxTrials    = max(opt.Trials, targetK);
maxTrialsCap = max(maxTrials, opt.MaxTrialsCap);

while numel(lambda) < targetK
    found = false;

    while numel(lambda) < targetK && trial < maxTrials
        trial = trial + 1;

        % ---- initialize guess (triangular seeding if enabled) ----
        [xq, lamq, seedInfo] = local_init_guess(Aq, trial, opt);

        xR   = local_pack_qvec(xq);
        lamR = local_pack_qscalar(lamq);

        if wantInfo
            run = struct();
            run.trial      = trial;
            run.accepted   = false;
            run.success    = false;
            run.reason     = "";
            run.iters      = 0;
            run.rFinal     = NaN;

            run.seedX      = seedInfo.seedX;
            run.seedLambda = seedInfo.seedLambda;

            run.lambda0     = lamq;
            run.x0          = local_unpack_qvec(xR);
            run.x0norm      = norm(xR);

            run.lambdaFinal = quaternion.empty(0,1);
            run.xFinal      = quaternion.empty(n,0);

            if strcmpi(opt.InfoLevel,'full')
                run.lambdaHist      = quaternion.empty(0,1);
                run.resHist         = [];
                run.alphaHist       = [];
                run.gaugeIndexHist  = [];
                run.gaugeAbsHist    = [];
                run.xNormHist       = [];
                run.pivotHist       = quaternion.empty(0,1);
            end
        end

        % ---- refine one eigenpair via Newton ----
        [lamR, xR, rRaw, rOut, hist] = local_newton_one(LA, In, lamR, xR, opt);

        converged = (isfinite(rRaw) && rRaw <= opt.Tol);
        if converged
            nConverged = nConverged + 1;
        end

        if wantInfo
            run.success        = converged;
            run.iters          = hist.iters;
            run.rFinal         = rOut;
            run.rFinalRaw      = rRaw;

            if strcmpi(opt.InfoLevel,'full') && isfield(hist,'lambdaHist')
                run.lambdaHist     = hist.lambdaHist;
                run.resHist        = hist.resHist;
                run.alphaHist      = hist.alphaHist;

                run.gaugeIndexHist = hist.gaugeIndexHist;
                run.gaugeAbsHist   = hist.gaugeAbsHist;
                run.xNormHist      = hist.xNormHist;
                run.pivotHist      = hist.pivotHist;
            end

            run.lambdaFinal    = local_unpack_qscalar(lamR);
            run.xFinal         = local_unpack_qvec(xR);
        end

        if ~converged
            if wantInfo
                run.reason = "not converged";
                runlog{end+1,1} = run; %#ok<AGROW>
            end
            continue;
        end

        % ---- (optional) recompute v for the final lambda by solving a real least-squares / SVD problem ----
        if opt.RefineV
            try
                [xRref, rRawRef, rOutRef] = local_refine_vec(LA, In, lamR, xR, opt);
                if isfinite(rRawRef) && (rRawRef <= rRaw)
                    xR   = xRref;
                    rRaw = rRawRef;
                    rOut = rOutRef;
                end
            catch
                % If refinement fails, keep the Newton vector.
            end
        end

        if wantInfo
            run.rFinal     = rOut;
            run.rFinalRaw  = rRaw;
            run.xFinal     = local_unpack_qvec(xR);
            run.lambdaFinal= local_unpack_qscalar(lamR);
        end

        lamq = local_unpack_qscalar(lamR);
        found = true;
        break;
    end

    if ~found
        % --- AutoExtend logic (only if enabled AND Trials not explicitly set) ---
        if opt.AutoExtend && ~opt.TrialsWasSet && maxTrials < maxTrialsCap
            newMax = min(maxTrialsCap, maxTrials + max(ceil(opt.ExtendBy*n), 10));
            if newMax <= maxTrials
                break;
            end
            if opt.Verbose
                fprintf('leigq: AutoExtend increasing Trials %d -> %d (found %d/%d)\n', ...
                    maxTrials, newMax, numel(lambda), targetK);
            end
            maxTrials = newMax;
            continue;
        end
        break;
    end

    % ---- accept ----
    lambda(end+1,1)  = lamq;                    %#ok<AGROW>
    V(:,end+1)       = local_unpack_qvec(xR);   %#ok<AGROW>
    res(end+1,1)     = rOut;                    %#ok<AGROW>

    if wantInfo
        run.accepted = true;
        run.reason   = "accepted";
        runlog{end+1,1} = run; %#ok<AGROW>
    end

    if opt.Verbose
        fprintf('leigq: accepted #%d at trial %d/%d, res=%.3e', ...
            numel(lambda), trial, maxTrials, rOut);
    end
end

% ---- finalize info ----
if wantInfo
    info = local_finalize_info(opt, trial, nConverged, targetK, maxTrials, ...
        rankA, df0, nAcceptedZero, runlog, lambda);
else
    info = [];
end

% ---- build distinct representatives (optional) ----
if wantDistinct
    [lambdaU, VU, resU, distinctMeta] = local_distinct_representatives(lambda, V, res, opt);
    if wantInfo && ~isempty(info)
        info{1}.summary.nDistinct = numel(lambdaU);
        info{1}.distinct = distinctMeta;
    end
else
    lambdaU = quaternion.empty(0,1);
    VU      = quaternion.empty(n,0);
    resU    = zeros(0,1);
end

end

% =====================================================================
%                 INFO helpers
% =====================================================================
function info = local_finalize_info(opt, trial, nConverged, targetK, maxTrials, rankA, df0, nAcceptedZero, runlog, lambda)
summary = struct();
summary.nRuns         = trial;
summary.nAccepted     = numel(lambda);
summary.nConverged    = nConverged;
summary.targetK       = targetK;
summary.maxTrials     = maxTrials;
summary.opt           = opt;
summary.rankA         = rankA;
summary.df0           = df0;
summary.nAcceptedZero = nAcceptedZero;

% --- restart / iteration statistics (Newton trials only) ---
trials = trial;
summary.trialsTotal = trials;

trialIters = zeros(trials,1);
accTrialsPerEig = [];
accIters = [];
accTrialIdx = [];
countSince = 0;
for k=1:numel(runlog)
    r = runlog{k};
    if ~isstruct(r) || ~isfield(r,'trial') || r.trial<=0
        continue;
    end
    t = r.trial;
    if t>=1 && t<=trials
        if isfield(r,'iters') && ~isempty(r.iters)
            trialIters(t) = double(r.iters);
        end
    end
    countSince = countSince + 1;
    if isfield(r,'accepted') && r.accepted
        accTrialsPerEig(end+1,1) = countSince; %#ok<AGROW>
        accIters(end+1,1) = double(r.iters); %#ok<AGROW>
        accTrialIdx(end+1,1) = t; %#ok<AGROW>
        countSince = 0;
    end
end

summary.trialIters = trialIters;
summary.itersTotal = sum(trialIters);
summary.acceptedNewton = numel(accTrialsPerEig);
summary.acceptedTrialsPerEig = accTrialsPerEig;
summary.acceptedTrialIdx = accTrialIdx;
summary.acceptedIters = accIters;
summary.itersAcceptedTotal = sum(accIters);
summary.tailTrialsAfterLastAccept = countSince;
summary.restartsTotal = max(0, trials - summary.acceptedNewton);
summary.infoLevel = opt.InfoLevel;

info = [{summary}; runlog(:)];
end

function run = local_make_run_stub_zero(k0, lam0, xq, rRaw)
run = struct();
run.trial      = 0;
run.accepted   = true;
run.success    = true;
run.reason     = "zero-eig from null(A)";
run.iters      = 0;
run.rFinal     = rRaw;

run.seedX      = "null(A)";
run.seedLambda = "0";

run.lambda0     = lam0;
run.x0          = xq;
run.x0norm      = local_qvecnorm_q(xq);

run.lambdaFinal = lam0;
run.xFinal      = xq;

run.lambdaHist      = lam0;
run.resHist         = rRaw;
run.alphaHist       = [];
run.gaugeIndexHist  = NaN;
run.gaugeAbsHist    = NaN;
run.xNormHist       = NaN;
run.pivotHist       = quaternion.empty(0,1);

run.zeroIndex = k0;
end

% =====================================================================
%                 OPTION PARSING (PROFILE-AWARE)
% =====================================================================
function opt = local_parse_opts(n, varargin)

% We use NaN sentinels for profile-controlled defaults (easy to tune later).
opt.Profile     = 'default';
opt.Num         = n;

opt.Trials      = NaN;
opt.MaxIter     = NaN;
opt.Tol         = NaN;

opt.Damping     = 1.0;
opt.Backtrack   = true;
opt.MinAlpha    = 1/64;

opt.UniqueTol   = 1e-6;
opt.CleanTol    = 1e-12;

% Distinct representative selection (post-processing of raw hits)
opt.DistinctTolAbs    = opt.UniqueTol;
opt.DistinctTolRel    = 0;
opt.DistinctResFactor = 10;

opt.DistinctTolAbsWasSet = false;  % track explicit user setting

opt.RankTolFactor = 1e2;   % multiplier for rank/null tolerance in the ZERO-eigenvalue pre-pass

opt.Seed        = [];
opt.Verbose     = 0;

opt.InfoLevel   = 'full';   % 'full' | 'summary'
opt.RecordHist  = true;     % internal, derived from InfoLevel


%
% Robustness / post-processing (defaults chosen for public, stable behavior)
opt.VerifyZeroNull     = true;
opt.ZeroNullTol        = [];     % [] -> auto (based on opt.Tol and scaling)
opt.UseNullFallbackLA  = true;
opt.ResidualNormalized = true;   % affects output res only (convergence uses raw residual)
opt.RefineV            = true;   % refine eigenvectors after Newton convergence

opt.Lambda0     = [];
opt.V0          = [];

opt.TriangularInit = true;
opt.TriTol         = 0;

opt.TriangularShortcut = 'diag';
opt.TriangularShortcutWasSet = false;

opt.IsTriangularUser = [];   % user override for triangular flag (empty -> auto-detect)

opt.TrialsFactor = [];       % optional convenience: Trials = TrialsFactor*n
opt.MaxTrials    = [];       % optional hard cap on total trials (disables AutoExtend)

opt.AutoExtend   = NaN;     % profile-dependent
opt.MaxTrialsCap = NaN;     % profile-dependent
opt.ExtendBy     = NaN;     % profile-dependent (in multiples of n)

% Track which knobs the user explicitly set.
opt.TrialsWasSet  = false;
opt.MaxIterWasSet = false;
opt.TolWasSet     = false;

% Leading numeric K?
if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
    opt.Num = varargin{1};
    varargin = varargin(2:end);
end

% Leading struct?
if ~isempty(varargin) && isstruct(varargin{1})
    s = varargin{1};
    varargin = varargin(2:end);
    f = fieldnames(s);
    for k=1:numel(f)
        opt.(f{k}) = s.(f{k});
    end
end

if ~isempty(varargin)
    if mod(numel(varargin),2) ~= 0
        error('leigq:BadArgs','Options must be name/value pairs (or a struct).');
    end

    for k=1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        if ~(ischar(name) || (isstring(name) && isscalar(name)))
            continue;
        end
        key = lower(char(name));

        switch key
            case {'solveprofile','profile'}
                opt.Profile = char(val);

            case {'num','k','neigs','nev'}
                opt.Num = double(val);

            case {'trials','restarts','ntrials'}
                opt.Trials = double(val);
                opt.TrialsWasSet = true;

            case {'maxiter','maxit','iter'}
                opt.MaxIter = double(val);
                opt.MaxIterWasSet = true;

            case {'tol','restol','tolerance'}
                opt.Tol = double(val);
                opt.TolWasSet = true;

            case {'damping','alpha','stepsize'}
                opt.Damping = double(val);

            case {'backtrack','linesearch'}
                opt.Backtrack = logical(val);

            case {'minalpha','alphamin'}
                opt.MinAlpha = double(val);

            case {'uniquetol','duptol'}
                % Backward-compatibility alias: now used for distinct grouping.
                opt.UniqueTol = double(val);
                opt.DistinctTolAbs = double(val);
                opt.DistinctTolAbsWasSet = true;
                opt.DistinctTolAbsWasSet = true;

            case {'distincttol','distincttolabs','distincttol_abs','distincttolabsolute'}
                opt.DistinctTolAbs = double(val);
                opt.DistinctTolAbsWasSet = true;

            case {'distincttolrel','distincttol_rel','distincttolrelative'}
                opt.DistinctTolRel = double(val);

            case {'distinctresfactor','distinctres','distinctresdom'}
                opt.DistinctResFactor = double(val);

            case {'cleantol','clean','cleantolerance'}
                opt.CleanTol = double(val);

            case {'ranktolfactor','ranktol','nulltolfactor'}
                opt.RankTolFactor = double(val);

            case {'seed','rng'}
                opt.Seed = double(val);

            case {'verbose','disp'}
                opt.Verbose = double(val) ~= 0;

            case {'infolevel','infomode'}
                if isnumeric(val) && isscalar(val)
                    opt.InfoLevel = double(val);
                else
                    opt.InfoLevel = char(val);
                end

            case {'lambda0','lam0','initlambda'}
                opt.Lambda0 = val;

            case {'v0','x0','initvec'}
                opt.V0 = val;

            case {'triangularinit','triinit'}
                opt.TriangularInit = logical(val);

            case {'tritol','triangulartol'}
                opt.TriTol = double(val);

            case {'istriangular','triangularflag','forcetriangular'}
                % User override: treat A as triangular (or not), regardless of detection.
                opt.IsTriangularUser = logical(val);

            case {'triangularshortcut','trishortcut','diagshortcut'}
                opt.TriangularShortcut = val;
                opt.TriangularShortcutWasSet = true;

            case {'trialsfactor'}
                % Convenience option: set Trials = TrialsFactor*n (unless Trials explicitly set).
                opt.TrialsFactor = double(val);

            case {'maxtrials','maxrestarts','trialscap_hard'}
                % Hard cap on total trials (also disables AutoExtend).
                opt.MaxTrials = double(val);

            case {'autoextend','extendtrials'}
                opt.AutoExtend = logical(val);

            case {'maxtrialscap','trialscap'}
                opt.MaxTrialsCap = double(val);

            case {'extendby','extendfactor'}
                opt.ExtendBy = double(val);


            case {'verifyzeronull','zeronullverify'}
                opt.VerifyZeroNull = logical(val);

            case {'zeronulltol','zeronulltolerance'}
                opt.ZeroNullTol = double(val);

            case {'usenullfallbackla','fallbackla','nullfallbackla'}
                opt.UseNullFallbackLA = logical(val);

            case {'residualnormalized','normalizedresidual','normres'}
                opt.ResidualNormalized = logical(val);

            case {'refinev','refinevec','refinevector'}
                opt.RefineV = logical(val);

            otherwise
                % ignore unknown options intentionally
        end
    end
end

% ---- Apply profile defaults to any NaN knobs (dense matrices base) ----
prof = local_profile_defaults(n, opt.Profile);

if isnan(opt.Trials),      opt.Trials      = prof.Trials;      end
if isnan(opt.MaxIter),     opt.MaxIter     = prof.MaxIter;     end
if isnan(opt.Tol),         opt.Tol         = prof.Tol;         end
if isnan(opt.AutoExtend),  opt.AutoExtend  = prof.AutoExtend;  end
if isnan(opt.MaxTrialsCap),opt.MaxTrialsCap= prof.MaxTrialsCap;end
if isnan(opt.ExtendBy),    opt.ExtendBy    = prof.ExtendBy;    end

% ---- Convenience overrides: TrialsFactor / MaxTrials ----
if ~isempty(opt.TrialsFactor) && ~opt.TrialsWasSet
    opt.Trials = opt.TrialsFactor * n;
end
if ~isempty(opt.MaxTrials)
    % Hard cap disables auto-extension and clamps Trials.
    opt.AutoExtend = false;
    opt.Trials = min(opt.Trials, opt.MaxTrials);
    opt.MaxTrialsCap = min(opt.MaxTrialsCap, opt.MaxTrials);
end

% ---- Bounds / normalization ----
opt.Num      = max(0, min(n, round(opt.Num)));
opt.Trials   = max(1, round(opt.Trials));
opt.MaxIter  = max(1, round(opt.MaxIter));
opt.Tol      = max(0, opt.Tol);

opt.Damping  = max(0.05, min(1.0, opt.Damping));
opt.MinAlpha = max(1e-6, min(1.0, opt.MinAlpha));

opt.UniqueTol = max(0, opt.UniqueTol);
opt.CleanTol  = max(0, opt.CleanTol);
opt.DistinctTolAbs = max(0, opt.DistinctTolAbs);

% ---- Default for distinct grouping (lambdaU): more tolerant than Tol unless user overrides ----
% Rationale: on defective/triangular cases, eigenvalue error can scale like sqrt(residual).
if ~opt.DistinctTolAbsWasSet
    opt.DistinctTolAbs = max(opt.DistinctTolAbs, max(1e-6, sqrt(opt.Tol)));
end

% ---- Normalize TriangularShortcut ----
opt.TriangularShortcut = local_norm_triangular_shortcut(opt.TriangularShortcut);
opt.DistinctTolRel = max(0, opt.DistinctTolRel);
opt.DistinctResFactor = max(1, opt.DistinctResFactor);

opt.MaxTrialsCap = max(opt.Trials, round(opt.MaxTrialsCap));
opt.ExtendBy = max(0.1, opt.ExtendBy);

% ---- Info logging level ----
if isnumeric(opt.InfoLevel) && isscalar(opt.InfoLevel)
    lvl = double(opt.InfoLevel);
    if lvl >= 2
        opt.InfoLevel = 'full';
    else
        opt.InfoLevel = 'summary';
    end
else
    opt.InfoLevel = lower(strtrim(char(opt.InfoLevel)));
    if isempty(opt.InfoLevel), opt.InfoLevel = 'full'; end
end
if ~ismember(opt.InfoLevel, {'full','summary'})
    error('leigq:BadInfoLevel','InfoLevel must be ''full'' or ''summary'' (or numeric 0/1/2).');
end
opt.RecordHist = strcmpi(opt.InfoLevel,'full');


function mode = local_norm_triangular_shortcut(val)
% Normalize TriangularShortcut option to a small set of strings.
% Supported modes: 'diag' (diagonal shortcut), 'off'.

if islogical(val) && isscalar(val)
    if val, mode = 'diag'; else, mode = 'off'; end
    return
end
if isnumeric(val) && isscalar(val)
    if val~=0, mode = 'diag'; else, mode = 'off'; end
    return
end

s = lower(strtrim(char(string(val))));
if isempty(s)
    mode = 'diag';
    return
end

if any(strcmp(s, {'off','none','false','0'}))
    mode = 'off';
else
    % default on
    mode = 'diag';
end
end

end

function prof = local_profile_defaults(n, profile)
% ======================= DEFAULTS (TUNABLE) =======================
% Edit THIS function if you want to tune the built-in behavior.
%
% Dense-matrix oriented heuristics:
%   - Trials scales ~ O(n) to O(n log n) depending on profile.
%   - Larger n typically needs many more restarts to cover basins.
%   - Reliable profile enables AutoExtend.

p = lower(strtrim(char(profile)));
if isempty(p), p = 'default'; end

switch p
    case {'fast','quick'}
        prof.Trials      = max(2*n, 20);
        prof.MaxIter     = 30;
        prof.Tol         = 1e-9;
        prof.AutoExtend  = false;
        prof.MaxTrialsCap= prof.Trials;
        prof.ExtendBy    = 0.0;

    case {'reliable','robust','safe'}
        prof.Trials      = max(18*n, 200);   % high initial budget
        prof.MaxIter     = 80;
        prof.Tol         = 1e-10;
        prof.AutoExtend  = true;
        prof.MaxTrialsCap= max(60*n, 1500);  % cap for extension
        % prof.MaxTrialsCap= max(100*n, 1500);  % increased to help n = 8 for 'hermitian'
        %prof.ExtendBy    = 6.0;              % add ~6n trials each extension
        prof.ExtendBy    = 10.0;              % add ~6n trials each extension

    otherwise % {'default','balanced',...}
        prof.Trials      = max(8*n, 80);
        prof.MaxIter     = 50;
        prof.Tol         = 1e-10;
        prof.AutoExtend  = true;
        prof.MaxTrialsCap= max(24*n, 600);
        prof.ExtendBy    = 4.0;
end

% Slight bump for larger n (dense):
if n >= 64
    prof.Trials = ceil(1.25*prof.Trials);
end

end

% =====================================================================
%                 INITIALIZATION
% =====================================================================
function [xq, lamq, seedInfo] = local_init_guess(Aq, trial, opt)
n = size(Aq,1);
seedInfo = struct('seedX','', 'seedLambda','');

% User-provided initial guesses apply ONLY to trial 1.
if trial == 1 && ~isempty(opt.V0)
    xq = local_to_quat(opt.V0);
    xq = xq(:);
    if numel(xq) ~= n
        error('leigq:BadV0','V0/X0 must have length n.');
    end
    seedInfo.seedX = 'user';
else
    % Triangular seeding: first n trials use standard basis vectors.
    if opt.TriangularInit && opt.IsTriangular && trial <= n
        xq = quaternion(zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1));
        xq(trial) = quaternion(1,0,0,0);
        seedInfo.seedX = 'triangular-basis';
    else
        xq = local_rand_vec_quat(n);
        seedInfo.seedX = 'random';
    end
end

if trial == 1 && ~isempty(opt.Lambda0)
    lamq = local_to_quat(opt.Lambda0);
    lamq = lamq(1);
    seedInfo.seedLambda = 'user';
else
    if opt.TriangularInit && opt.IsTriangular && trial <= n
        lamq = Aq(trial,trial); % diagonal seed
        seedInfo.seedLambda = 'triangular-diagonal';
    else
        % Use Rayleigh quotient in packed-real form (cheap)
        lamR = local_rayleigh_left(local_left_block(Aq), local_pack_qvec(xq));
        lamq = local_unpack_qscalar(lamR);
        seedInfo.seedLambda = 'rayleigh';
    end
end
end

% =====================================================================
%                 CLEAN / UNIQUE
% =====================================================================
function lamc = local_clean_lambda(lam, tol)
% Component-wise "zeroing" of small components.
[w,x,y,z] = parts(lam);
w(abs(w)<tol) = 0;
x(abs(x)<tol) = 0;
y(abs(y)<tol) = 0;
z(abs(z)<tol) = 0;
lamc = quaternion(w,x,y,z);
end

function ok = local_is_new_lambda(lamc, list, tol)
if isempty(list)
    ok = true; return;
end
[w1,x1,y1,z1] = parts(lamc);
d = zeros(numel(list),1);
for k=1:numel(list)
    [w2,x2,y2,z2] = parts(list(k));
    dw = w1-w2; dx = x1-x2; dy = y1-y2; dz = z1-z2;
    d(k) = sqrt(dw*dw + dx*dx + dy*dy + dz*dz);
end
ok = all(d > tol);
end

% =====================================================================
%                 NEWTON CORE (one eigenpair)
% =====================================================================
function [lamR, xR, rRaw, rOut, hist] = local_newton_one(LA, In, lamR, xR, opt)
n = size(In,1);
record = isfield(opt,'RecordHist') && logical(opt.RecordHist);

if record
    hist.lambdaHist     = quaternion.empty(0,1);
    hist.resHist        = zeros(0,1);
    hist.resHistRaw     = zeros(0,1);
    hist.alphaHist      = zeros(0,1);
    hist.gaugeIndexHist = zeros(0,1);
    hist.gaugeAbsHist   = zeros(0,1);
    hist.xNormHist      = zeros(0,1);
    hist.pivotHist      = quaternion.empty(0,1);
else
    hist = struct();
end
hist.iters = 0;

% Gauge once at start (and log if requested)
if record
    [xR, j, pivAbs, xPreNorm, pivVal] = local_gauge_log(xR);
else
    [xR, j] = local_gauge(xR);
    pivAbs = NaN; xPreNorm = norm(xR); pivVal = quaternion(0,0,0,0);
end

[rR, rRaw, rOut] = local_residual_metrics(LA, In, lamR, xR, opt);

if record
    hist.lambdaHist(1,1)     = local_unpack_qscalar(lamR);
    hist.resHist(1,1)        = rOut;
    hist.resHistRaw(1,1)     = rRaw;
    hist.gaugeIndexHist(1,1) = j;
    hist.gaugeAbsHist(1,1)   = pivAbs;
    hist.xNormHist(1,1)      = xPreNorm;
    hist.pivotHist(1,1)      = pivVal;
end

for it = 1:opt.MaxIter
    if rRaw <= opt.Tol
        hist.iters = it-1;
        return;
    end

    [Mreal, b] = local_build_newton_system(LA, In, lamR, xR, rR, j);
    delta = Mreal \ b;

    dxR   = delta(1:4*n);
    dlamR = delta(4*n+1:end);

    alpha = opt.Damping;
    if opt.Backtrack
        alpha = local_backtrack(LA, In, lamR, xR, dxR, dlamR, alpha, opt.MinAlpha, rRaw);
    end

    xR   = xR   + alpha*dxR;
    lamR = lamR + alpha*dlamR;

    if record
        [xR, j, pivAbs, xPreNorm, pivVal] = local_gauge_log(xR);
    else
        [xR, j] = local_gauge(xR);
        pivAbs = NaN; xPreNorm = norm(xR); pivVal = quaternion(0,0,0,0);
    end

    [rR, rRaw, rOut] = local_residual_metrics(LA, In, lamR, xR, opt);

    if record
        hist.alphaHist(it,1)           = alpha;
        hist.lambdaHist(it+1,1)        = local_unpack_qscalar(lamR);
        hist.resHist(it+1,1)           = rOut;
        hist.resHistRaw(it+1,1)        = rRaw;
        hist.gaugeIndexHist(it+1,1)    = j;
        hist.gaugeAbsHist(it+1,1)      = pivAbs;
        hist.xNormHist(it+1,1)         = xPreNorm;
        hist.pivotHist(it+1,1)         = pivVal;
    end

    if opt.Verbose && (mod(it,10)==0)
        fprintf('  it=%d, res=%.3e\n', it, rOut);
    end
end

hist.iters = opt.MaxIter;
end

function alpha = local_backtrack(LA, In, lamR, xR, dxR, dlamR, alpha0, alphaMin, r0)
alpha = alpha0;
while alpha >= alphaMin
    xtR   = xR   + alpha*dxR;
    lamtR = lamR + alpha*dlamR;
    xtR   = local_gauge_only(xtR);
    rtR   = local_residual(LA, In, lamtR, xtR);
    r1    = norm(rtR);
    if r1 < r0
        return;
    end
    alpha = alpha/2;
end
end

function [Mreal, b] = local_build_newton_system(LA, In, lamR, xR, rR, j)
n = size(In,1);

Llam  = local_Lmat(lamR(1), lamR(2), lamR(3), lamR(4));
LAlam = LA - kron(In, Llam);
Rx    = local_right_stack_from_pack(xR);

Mcore = [LAlam, -Rx];
bcore = -rR;

% constraint: dot(xR, dxR) = 0
c1 = [xR.', zeros(1,4)];
b1 = 0;

% constraints: imag parts of dx_j = 0
idx0 = 4*(j-1);
c2 = zeros(1,4*n+4); c2(idx0+2) = 1;
c3 = zeros(1,4*n+4); c3(idx0+3) = 1;
c4 = zeros(1,4*n+4); c4(idx0+4) = 1;

Mreal = [Mcore; c1; c2; c3; c4];
b     = [bcore; b1; 0; 0; 0];
end

% =====================================================================
%                 GAUGE / NORMALIZATION
% =====================================================================
function [xR, j] = local_gauge(xR)
absx = local_block_abs(xR);
[~, j] = max(absx(:));

if absx(j) == 0
    j = 1;
else
    q = local_get_block(xR, j) / absx(j);
    s = local_qconj_pack(q);
    xR = local_right_scale_pack(xR, s);
end

nx = norm(xR);
if nx > 0
    xR = xR / nx;
end

w = xR(4*(j-1)+1);
if w < 0
    xR = -xR;
end
end

function [xR, j, pivAbs, xPreNorm, pivValQ] = local_gauge_log(xR)
absx = local_block_abs(xR);
[pivAbs, j] = max(absx(:));

if pivAbs == 0
    j = 1;
    pivValQ = quaternion(0,0,0,0);
    s = [1;0;0;0];
else
    pivPack = local_get_block(xR, j);
    pivValQ = quaternion(pivPack(1),pivPack(2),pivPack(3),pivPack(4));
    q = pivPack / pivAbs;
    s = local_qconj_pack(q);
end

xR = local_right_scale_pack(xR, s);
xPreNorm = norm(xR);

if xPreNorm > 0
    xR = xR / xPreNorm;
end

w = xR(4*(j-1)+1);
if w < 0
    xR = -xR;
end
end

function xR = local_gauge_only(xR)
[xR, ~] = local_gauge(xR);
end

% =====================================================================
%                 RESIDUALS / RAYLEIGH
% =====================================================================
function lamR = local_rayleigh_left(LA, xR)
AxR = LA * xR;
den = norm(xR)^2;
if den == 0
    lamR = zeros(4,1);
    return;
end
n = numel(xR)/4;
num = zeros(4,1);
for k = 1:n
    ak = local_get_block(AxR, k);
    xk = local_get_block(xR,  k);
    xkc = local_qconj_pack(xk);
    num = num + local_qmul_pack(ak, xkc);
end
lamR = num / den;
end

function rR = local_residual(LA, In, lamR, xR)
Llam  = local_Lmat(lamR(1), lamR(2), lamR(3), lamR(4));
rR = (LA - kron(In, Llam)) * xR;
end

function den = local_residual_den(LA, In, lamR, xR, opt)
% Scale used for normalized residual reporting (NOT used for convergence).
xnorm = norm(xR);
if xnorm == 0
    den = 1;
    return;
end

if isfield(opt,'Afro') && ~isempty(opt.Afro)
    aF = opt.Afro;
else
    aF = norm(LA,'fro');
end

Llam = local_Lmat(lamR(1), lamR(2), lamR(3), lamR(4));
den  = (aF + norm(Llam,'fro')) * xnorm + 1;
end

function [rR, rRaw, rOut] = local_residual_metrics(LA, In, lamR, xR, opt)
% rRaw is the raw Euclidean norm in the real embedding.
% rOut is what we REPORT in the output 'res' (optionally normalized).
rR   = local_residual(LA, In, lamR, xR);
rRaw = norm(rR);
if isfield(opt,'ResidualNormalized') && logical(opt.ResidualNormalized)
    den  = local_residual_den(LA, In, lamR, xR, opt);
    rOut = rRaw / den;
else
    rOut = rRaw;
end
end

function tolRaw = local_zero_accept_tol(LA, xR, opt)
% Acceptance tolerance for zero-eigenvalue pre-pass vectors (raw residual).
% If opt.ZeroNullTol is set, it is used directly (raw units). Otherwise, use
% opt.Tol with a mild scale-aware eps floor.
if isfield(opt,'ZeroNullTol') && ~isempty(opt.ZeroNullTol)
    tolRaw = double(opt.ZeroNullTol);
    return;
end
base = opt.Tol;
den  = (norm(LA,'fro')*norm(xR) + 1);
tolRaw = max(base, 1e3*eps*den);
end

function [XqOK, okMask] = local_verify_zero_nullspace(LA, In, Xq, opt)
% Verify Aq*x≈0 (in real embedding) for each candidate column of Xq.
m = size(Xq,2);
okMask = true(1,m);
lam0R = zeros(4,1);

for k = 1:m
    xR = local_pack_qvec(Xq(:,k));
    [xR, ~] = local_gauge(xR);
    [~, rRaw, ~] = local_residual_metrics(LA, In, lam0R, xR, opt);
    if rRaw > local_zero_accept_tol(LA, xR, opt)
        okMask(k) = false;
    else
        % store the gauged/normalized version back (for consistent output)
        Xq(:,k) = local_unpack_qvec(xR);
    end
end

XqOK = Xq(:, okMask);
end

function Xq = local_null_quat_LA(Aq, df0, LA)
% Nullspace via REAL 4n embedding: null(LA) gives a real basis of dimension 4*df0.
% We then select df0 right-H independent quaternion columns.
n = size(Aq,1);
if df0 <= 0
    Xq = quaternion.empty(n,0);
    return;
end

ZR = null(LA,'r');  % 4n-by-mR (orthonormal real basis)
mR = size(ZR,2);
if mR == 0
    Xq = quaternion.empty(n,0);
    return;
end

% Convert each real basis vector to a quaternion vector candidate
Xcand = quaternion.empty(n, mR);
for k = 1:mR
    Xcand(:,k) = local_unpack_qvec(ZR(:,k));
end

% Select df0 right-H independent columns (using complex embedding rank test)
keep = local_select_independent_quat_cols(Xcand, df0);
Xq   = Xcand(:, keep);
end

function keep = local_select_independent_quat_cols(Xq, kNeed)
% Select kNeed columns that are right-H independent (best effort).
% Uses complex embedding rank on columns.
n = size(Xq,1);
m = size(Xq,2);
keep = [];

if kNeed <= 0 || m == 0
    return;
end

% Complex embedding: q = a + b*j, with a,b in C (i is MATLAB's complex unit).
[W,X,Y,Z] = parts(Xq);
A = W + 1i*X;
B = Y + 1i*Z;

% Column embedding: [A; B] is 2n-by-m complex
C = [A; B];

% Greedy selection
for j = 1:m
    cand = [keep, j];
    if rank(C(:,cand)) == numel(cand)
        keep = cand;
    end
    if numel(keep) >= kNeed
        break;
    end
end

% If still short, take what we have (caller will min with dfFound)
end

function [xRref, rRawRef, rOutRef] = local_refine_vec(LA, In, lamR, xR, opt)
% Recompute v for the FINAL lambda by solving a real least-squares / SVD problem:
%   minimize || (LA - kron(I,L(lam))) x ||  subject to ||x||=1 and gauge constraints handled by local_gauge.
%
% We obtain x as the eigenvector corresponding to the smallest eigenvalue of M'*M.
n4 = size(LA,1);
Llam  = local_Lmat(lamR(1), lamR(2), lamR(3), lamR(4));
M     = LA - kron(In, Llam);

x0 = xR;
if norm(x0) == 0
    x0 = randn(n4,1);
end

xRref = x0 / norm(x0);

try
    G = M.'*M;              % symmetric PSD
    optsE.isreal = true;
    % NOTE (MATLAB compatibility): opts.issym is only honored when the first
    % input to EIGS is a function handle. For numeric matrices, setting opts.issym
    % triggers repetitive warnings ("Ignoring issym field...") in some MATLAB versions.
    optsE.tol    = 1e-12;
    optsE.maxit  = 1000;
    optsE.v0 = xRref;       % <- ADD THIS (already normalized)
    [vmin,~] = eigs(G, 1, 'smallestreal', optsE);
    xRref = vmin;
catch
    % fallback: full SVD (small problems) or last resort
    [~,~,V] = svd(M,'econ');
    xRref = V(:,end);
end

% Gauge + normalize to match the solver's convention
[xRref, ~] = local_gauge(xRref);

% Residual metrics for the refined vector
[~, rRawRef, rOutRef] = local_residual_metrics(LA, In, lamR, xRref, opt);
end


% =====================================================================
%                 REAL BLOCK EMBEDDING
% =====================================================================
function Z = local_left_block(Q)
n = size(Q,1);
Z = zeros(4*n, 4*n);
[W,X,Y,Zz] = parts(Q);
for i=1:n
    ii = 4*(i-1)+1:4*i;
    for j=1:n
        jj = 4*(j-1)+1:4*j;
        Z(ii,jj) = local_Lmat(W(i,j), X(i,j), Y(i,j), Zz(i,j));
    end
end
end

function Z = local_right_stack_from_pack(xR)
n = numel(xR)/4;
Z = zeros(4*n, 4);
for i=1:n
    ii = 4*(i-1)+1:4*i;
    blk = xR(ii);
    Z(ii,:) = local_Rmat(blk(1), blk(2), blk(3), blk(4));
end
end

function L = local_Lmat(a,b,c,d)
L = [ a, -b, -c, -d;
    b,  a, -d,  c;
    c,  d,  a, -b;
    d, -c,  b,  a ];
end

function R = local_Rmat(a,b,c,d)
R = [ a, -b, -c, -d;
    b,  a,  d, -c;
    c, -d,  a,  b;
    d,  c, -b,  a ];
end

% =====================================================================
%                 PACK / UNPACK
% =====================================================================
function y = local_pack_qvec(xq)
n = size(xq,1);
y = zeros(4*n,1);
[W,X,Y,Zz] = parts(xq);
for k=1:n
    idx = 4*(k-1);
    y(idx+1) = W(k);
    y(idx+2) = X(k);
    y(idx+3) = Y(k);
    y(idx+4) = Zz(k);
end
end

function xq = local_unpack_qvec(y)
n = numel(y)/4;
W = zeros(n,1); X = W; Y = W; Zz = W;
for k=1:n
    idx = 4*(k-1);
    W(k)  = y(idx+1);
    X(k)  = y(idx+2);
    Y(k)  = y(idx+3);
    Zz(k) = y(idx+4);
end
xq = quaternion(W,X,Y,Zz);
end

function y4 = local_pack_qscalar(q)
[w,x,y,z] = parts(q);
y4 = [w; x; y; z];
end

function q = local_unpack_qscalar(y4)
q = quaternion(y4(1), y4(2), y4(3), y4(4));
end

% =====================================================================
%                 RANDOM / TYPE
% =====================================================================
function xq = local_rand_vec_quat(n)
xq = quaternion(randn(n,1), randn(n,1), randn(n,1), randn(n,1));
end

function Aq = local_to_quat(x)
if isa(x,'quaternion')
    Aq = x;
elseif isnumeric(x)
    Aq = quaternion(x, zeros(size(x)), zeros(size(x)), zeros(size(x)));
else
    error('leigq:Type','Unsupported type "%s". Input must be quaternion or numeric.', class(x));
end
end

% =====================================================================
%                 TRIANGULAR DETECTION
% =====================================================================
function tf = local_is_triangular(Aq, tol)
% Detect strict upper or lower triangular (within tol) across all components.
[W,X,Y,Z] = parts(Aq);

if tol <= 0
    tfU = isequal(W, triu(W)) && isequal(X, triu(X)) && isequal(Y, triu(Y)) && isequal(Z, triu(Z));
    tfL = isequal(W, tril(W)) && isequal(X, tril(X)) && isequal(Y, tril(Y)) && isequal(Z, tril(Z));
else
    tfU = (norm(W - triu(W), 'fro') <= tol) && (norm(X - triu(X), 'fro') <= tol) && ...
        (norm(Y - triu(Y), 'fro') <= tol) && (norm(Z - triu(Z), 'fro') <= tol);
    tfL = (norm(W - tril(W), 'fro') <= tol) && (norm(X - tril(X), 'fro') <= tol) && ...
        (norm(Y - tril(Y), 'fro') <= tol) && (norm(Z - tril(Z), 'fro') <= tol);
end

tf = tfU || tfL;
end


function tf = local_is_diagonal(Aq, tol)
% Detect diagonal (within tol) across all components.
[W,X,Y,Z] = parts(Aq);

if tol <= 0
    tf = isequal(W, diag(diag(W))) && isequal(X, diag(diag(X))) && ...
         isequal(Y, diag(diag(Y))) && isequal(Z, diag(diag(Z)));
else
    tf = (norm(W - diag(diag(W)), 'fro') <= tol) && (norm(X - diag(diag(X)), 'fro') <= tol) && ...
         (norm(Y - diag(diag(Y)), 'fro') <= tol) && (norm(Z - diag(diag(Z)), 'fro') <= tol);
end
end


function [lambda, V, res, info, lambdaU, VU, resU] = local_diagonal_shortcut(Aq, opt, wantInfo, wantDistinct)
% Fast path for diagonal matrices (default behavior unless disabled).

n = size(Aq,1);
K = opt.Num;

lamAll = diag(Aq);
lambda = lamAll(1:K);

% Standard basis eigenvectors
I = eye(n);
Vfull = quaternion(I, zeros(n), zeros(n), zeros(n));
V = Vfull(:,1:K);

% Exact residuals for diagonal basis pairs
res = zeros(K,1);

% Rank / df0 are cheap for diagonal matrices.
absLam = local_quat_abs(lamAll);
zeroMask = absLam <= max(opt.CleanTol, 0);
df0 = sum(zeroMask);
rankA = n - df0;

% Count how many exact zeros are present in the RETURNED raw list.
nAcceptedZero = sum(zeroMask(1:K));

runlog = cell(0,1);
maxTrials = 0;
trial = 0;
nConverged = 0;

if wantInfo
    % Minimal per-eigenpair stubs (trial=0 indicates shortcut, not Newton)
    runlog = cell(K,1);
    for kk = 1:K
        r = struct();
        r.trial      = 0;
        r.accepted   = true;
        r.success    = true;
        r.reason     = "shortcut(diagonal)";
        r.iters      = 0;
        r.rFinal     = res(kk);
        r.rFinalRaw  = 0;
        r.seedX      = 'basis';
        r.seedLambda = 'diag';
        r.lambda0    = lambda(kk);
        r.x0         = V(:,kk);
        r.x0norm     = 1;
        r.lambdaFinal= lambda(kk);
        r.xFinal     = V(:,kk);
        runlog{kk,1} = r;
    end

    info = local_finalize_info(opt, trial, nConverged, K, maxTrials, rankA, df0, nAcceptedZero, runlog, lambda);
    info{1}.summary.shortcut = struct('used',true,'mode','diagonal','TriangularShortcut',opt.TriangularShortcut,'TriTol',opt.TriTol);
else
    info = [];
end

if wantDistinct
    [lambdaU, VU, resU, metaD] = local_distinct_representatives(lambda, V, res, opt);
    if wantInfo && ~isempty(info)
        info{1}.summary.nDistinct = numel(lambdaU);
        info{1}.distinct = metaD;
    end
else
    lambdaU = quaternion.empty(0,1);
    VU      = quaternion.empty(n,0);
    resU    = zeros(0,1);
end

end

function a = local_quat_abs(q)
% Euclidean magnitude of quaternion scalars (vectorized).
[w,x,y,z] = parts(q);
a = sqrt(w.^2 + x.^2 + y.^2 + z.^2);
end

% =====================================================================
%                 NORMS / PACK HELPERS
% =====================================================================
function r = local_qvecnorm_q(xq)
ax = local_block_abs(local_pack_qvec(xq));
r = sqrt(sum(ax(:).^2));
end

function absx = local_block_abs(xR)
n = numel(xR)/4;
absx = zeros(n,1);
for k=1:n
    blk = local_get_block(xR,k);
    absx(k) = sqrt(sum(blk.^2));
end
end

function blk = local_get_block(xR, k)
ii = 4*(k-1)+1:4*k;
blk = xR(ii);
end

function xR = local_set_block(xR, k, blk)
ii = 4*(k-1)+1:4*k;
xR(ii) = blk(:);
end

function qc = local_qconj_pack(q)
qc = [q(1); -q(2); -q(3); -q(4)];
end

function p = local_qmul_pack(q, r)
a=q(1); b=q(2); c=q(3); d=q(4);
e=r(1); f=r(2); g=r(3); h=r(4);
p = [ a*e - b*f - c*g - d*h;
    a*f + b*e + c*h - d*g;
    a*g - b*h + c*e + d*f;
    a*h + b*g - c*f + d*e ];
end

function xR = local_right_scale_pack(xR, s)
n = numel(xR)/4;
Rs = local_Rmat(s(1),s(2),s(3),s(4));
for k=1:n
    blk = local_get_block(xR,k);
    xR = local_set_block(xR,k, Rs*blk);
end
end

% =====================================================================
%        QUATERNION RANK / NULL VIA COMPLEX EMBEDDING (PRE-PASS)
% =====================================================================
function rankA = local_rank_quat(Aq, opt)
PhiA = local_cembed(Aq);
rtol = local_ranktol(PhiA, opt);
rankA = 0.5 * rank(PhiA, rtol);
rankA = max(0, min(size(Aq,1), round(rankA)));
end

function Xq = local_null_quat(Aq, df0, opt)
n = size(Aq,1);
PhiA = local_cembed(Aq);

rtol = local_ranktol(PhiA, opt);
Zc = null(PhiA, rtol);      % complex 2n-by-m
m  = size(Zc,2);
if m == 0
    Xq = quaternion.empty(n,0);
    return;
end

cand = quaternion.empty(n,0);
for k=1:m
    u = Zc(1:n,k);
    v = Zc(n+1:end,k);
    cand(:,end+1) = quaternion(real(u), imag(u), real(v), imag(v)); %#ok<AGROW>
end

J = [zeros(n), eye(n); -eye(n), zeros(n)];
for k=1:m
    z2 = J * conj(Zc(:,k));
    u = z2(1:n);
    v = z2(n+1:end);
    cand(:,end+1) = quaternion(real(u), imag(u), real(v), imag(v)); %#ok<AGROW>
end

Xq = quaternion.empty(n,0);
for k=1:size(cand,2)
    if size(Xq,2) >= df0, break; end
    xk = cand(:,k);
    if local_rank_quat_cols([Xq, xk], opt) > local_rank_quat_cols(Xq, opt)
        Xq(:,end+1) = xk; %#ok<AGROW>
    end
end
end

function r = local_rank_quat_cols(Xq, opt)
if isempty(Xq)
    r = 0; return;
end
PhiX = local_cembed_rect(Xq);
rtol = local_ranktol(PhiX, opt);
r = 0.5*rank(PhiX, rtol);
r = max(0, min(size(Xq,2), round(r)));
end

function tol = local_ranktol(M, opt)
% Rank/null tolerance for complex embedding matrices.
% Larger RankTolFactor => more aggressive detection of singularity.
if isempty(M)
    tol = 0; return;
end
s = svd(M);
if isempty(s) || s(1)==0
    tol = 0;
else
    tol = opt.RankTolFactor * max(size(M)) * eps(s(1));
end
end

function PhiA = local_cembed(Aq)
[W,X,Y,Z] = parts(Aq);
A1 = W + 1i*X;
A2 = Y + 1i*Z;
PhiA = [A1, A2; -conj(A2), conj(A1)];
end

function PhiX = local_cembed_rect(Xq)
[W,X,Y,Z] = parts(Xq);
X1 = W + 1i*X;
X2 = Y + 1i*Z;
PhiX = [X1, X2; -conj(X2), conj(X1)];
end

% =====================================================================
% Distinct representatives (post-processing)
% =====================================================================
function [lambdaU, VU, resU, meta] = local_distinct_representatives(lambda, V, res, opt)
% Group approximately-equal lambdas and pick one representative per group.
%
% Heuristic representative choice:
%   1) Prefer substantially smaller residual (factor opt.DistinctResFactor)
%   2) Otherwise, prefer a “prettier” quaternion with more exact zeros/integers
%      (i.e., components exactly equal to 0 or exactly integer in floating point).

lambdaU = quaternion.empty(0,1);
VU      = quaternion.empty(size(V,1),0);
resU    = zeros(0,1);

meta = struct();
meta.groupOfHit = [];
meta.groups = {};
meta.representativeIndex = [];
meta.tolAbs = opt.DistinctTolAbs;
meta.tolRel = opt.DistinctTolRel;
meta.resFactor = opt.DistinctResFactor;

K = numel(lambda);
if K==0
    meta.nGroups = 0;
    return;
end

% Precompute R^4 parts
[w,x,y,z] = parts(lambda(:));
P = [w,x,y,z]; % Kx4

% Helper: distance between two 4-vectors
norm4 = @(a) sqrt(sum(a.^2,2));

repIdx = zeros(0,1);
repP   = zeros(0,4);
repRes = zeros(0,1);

groupOf = zeros(K,1);
groups = {};

for k = 1:K
    pk = P(k,:);

    % Find a matching group
    gMatch = 0;
    for g = 1:numel(repIdx)
        pref = repP(g,:);
        tol = opt.DistinctTolAbs + opt.DistinctTolRel * norm(pref);
        if norm(pk - pref) <= tol
            gMatch = g;
            break;
        end
    end

    if gMatch==0
        % New group
        repIdx(end+1,1) = k; %#ok<AGROW>
        repP(end+1,:)   = pk; %#ok<AGROW>
        repRes(end+1,1) = res(k); %#ok<AGROW>
        gMatch = numel(repIdx);
        groups{gMatch,1} = k; %#ok<AGROW>
    else
        groups{gMatch,1}(end+1,1) = k; %#ok<AGROW>

        % Decide whether to replace representative
        kBest = repIdx(gMatch);
        rBest = repRes(gMatch);
        rNew  = res(k);

        replace = false;

        if rNew < rBest / opt.DistinctResFactor
            replace = true;
        else
            % If residuals are comparable, use “prettiness” tie-break
            if rNew <= rBest * opt.DistinctResFactor
                replace = local_is_prettier(pk, P(kBest,:), rNew, rBest);
            end
        end

        if replace
            repIdx(gMatch) = k;
            repP(gMatch,:) = pk;
            repRes(gMatch) = rNew;
        end
    end

    groupOf(k) = gMatch;
end

% Build outputs in order of groups (stable)
G = numel(repIdx);
lambdaU = quaternion(w(repIdx), x(repIdx), y(repIdx), z(repIdx));
VU      = V(:, repIdx);
resU    = res(repIdx);

meta.nGroups = G;
meta.groupOfHit = groupOf;
meta.groups = groups;
meta.representativeIndex = repIdx;

end
function tf = local_is_prettier(pNew, pOld, rNew, rOld)
% Higher scores mean “nicer” representation.
% Exact zeros first, then exact integers, then smaller fractional tails.

[zNew, iNew, tailNew] = local_pretty_features(pNew);
[zOld, iOld, tailOld] = local_pretty_features(pOld);

if zNew ~= zOld
    tf = (zNew > zOld);
    return;
end
if iNew ~= iOld
    tf = (iNew > iOld);
    return;
end
if abs(tailNew - tailOld) > 0
    tf = (tailNew < tailOld);
    return;
end
% Final tie-break: smaller residual
tf = (rNew < rOld);
end

function [nZero, nInt, tail] = local_pretty_features(p)
nZero = sum(p==0);
nInt  = sum(p==round(p));
tail  = sum(abs(p - round(p)));
end
