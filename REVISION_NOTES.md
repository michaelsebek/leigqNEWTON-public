# leigqNEWTON v1.0.5 — LAA-R1-relative-2026-10-03

This revision corrects the public `leigqNEWTON` numerical core for LAA-D-26-00124.
It was published as GitHub release **v1.0.5** on **2026-10-03**, with version DOI
**[10.5281/zenodo.23124184](https://doi.org/10.5281/zenodo.23124184)**.
The numerical-core identifier remains `LAA-R1-relative-2026-10-03`.
Earlier package labels (GitHub tag 1.0.4, Contents 1.2, README 1.0) were not uniform;
this does not imply different numerical cores in the five input archives compared below.
The old DOI `10.5281/zenodo.18410141` identifies the first archived release, not v1.0.5.
After archiving, README and citation metadata on `main` were updated to record the assigned DOI.
No MATLAB source or test changes are part of that metadata-only update, and the release tag remains unchanged.

## Status and scope

The supplied GitHub main, GitHub 1.0.4 tag, audited archive, local public archive and local integrated-core archive all have the same 137 files after normalizing line endings. There are no additional local numerical fixes to merge. The main/tag files are byte-identical; the audited and local archives differ from them only by CRLF/LF in six text files.

The local integration wrappers were intentionally excluded by the user. This patch contains no wrappers and must not overwrite them. It updates the numerical core at its actual location, not a wrapper entry point with a different body.

All changes below underwent source review and independent real-array checks. The targeted MATLAB regression driver `test_leigqNEWTON_revision` was then run natively by the author in the stand-alone public environment on 2026-10-03, with the personal quaternion toolbox removed from the MATLAB path: **OVERALL: OK (29/29), 0.98 s**. This regression validation is separate from the earlier successful native `verify_revision` run validating manuscript matrix identities; neither run is a rerun of the historical benchmark campaign.

No historical numerical table, figure, timing, seed collection or benchmark has been rerun or replaced.

## Common residual and stopping convention

For a finite nonzero quaternion vector v:

    eta = norm(A*v - lambda*v,2) / ((norm(A,2)+abs(lambda))*norm(v,2)).

The matrix norm is the quaternion operator 2-norm, calculated from the standard real embedding (its spectral norm is the same); it is not the Frobenius norm of that embedding. The norm is computed once per solver call. There is no additive `+1` or `max(1,...)` floor. If A=lambda=0, eta=0 for a nonzero v. A zero or nonfinite vector never provides an acceptable eigenpair.

- `leigqNEWTON`: `Tol` (alias `ResTol`) now tests eta by default, both inside Newton and on final acceptance, including zero-eigenvalue candidates.
- `ToleranceMode='absolute'` explicitly requests a raw residual stopping criterion. This is a compatibility criterion, **not a bit-for-bit reproduction of old releases**.
- `ResidualNormalized` controls only the returned scalar residual. Setting it to false does not change the default relative acceptance test.
- `leigqNEWTON_cert_resPair`: default output remains raw; `ResidualNormalized=true` returns eta. All five output positions are unchanged.
- `leigqNEWTON_cert_resMin`: default output remains raw; the relative certificate is `sigma_min(A-lambda*I)/(norm(A,2)+abs(lambda))` and uses a unit minimizer.
- `checkNEWTON`: uses the same definitions. The minimal certificate no longer depends on the norm of a separately supplied vector, and vectors returned by the solver are no longer silently discarded and replaced by SVD vectors.
- `ZeroNullTol`: an optional extra raw-unit cap; it may tighten, never bypass, the common acceptance test.

The new private helper `private/leigqNEWTON_relres.m` evaluates the ratio with exponent splitting, avoiding overflow/underflow solely from forming the denominator. This does not promise immunity to overflow/underflow in arbitrary matrix-vector operations near floating-point limits.

## Correctness fixes

### Main Newton solver

1. Initialize trial-budget variables before an early successful return from the zero-eigenvalue pre-pass.
2. Invert the complex vector embedding correctly: if phi(x)=[u;-conj(v)] then x=u+v*j, so v=-conj(z_bottom). The previous conversion used z_bottom directly.
3. Real nullspace fallback now uses an SVD-based orthonormal basis, not `null(LA,'r')`. Quaternion array preallocation uses zero components rather than the invalid `quaternion.empty(n,m)` with both dimensions nonzero.
4. Select independent fallback columns using the full complex adjoint and right-H rank. Independence over R or over a single stacked complex vector is not the same condition.
5. Zero-prepass acceptance always checks the common certificate, even with `VerifyZeroNull=false` or a loose `ZeroNullTol`.
6. An approximate diagonal detected with a positive TriTol may guide seeding but cannot use the exact zero-residual diagonal shortcut. Only an exactly diagonal input uses that shortcut.
7. Exact diagonal entries are classified as zero only when zero, not merely when smaller than CleanTol.
8. Explicit trial budgets are not silently increased to the requested number of hits.
9. Failed line search no longer applies an untested step below MinAlpha. Invalid input and nonfinite candidates are rejected.
10. Balance the defect rows and lambda-increment variable by a matrix/eigenvalue scale, while keeping the same Newton equations in exact arithmetic. Scale the normal-equation vector-refinement matrix as well. This prevents a trivially small/large A from unbalancing defect and gauge rows.

### Local polisher

`leigqNEWTON_refine_polish` had a separate derivative error: its left-eigenvalue branch used the matrix for v*deltaLambda where the derivative requires deltaLambda*v. The corrected real coupling is right multiplication by the fixed vector entries and is shared with a direct noncommuting regression test.

The polisher now uses the common relative two-norm stopping criterion, checks an already converged initial pair before solving a potentially singular system, balances its equations, and reports stagnation separately from convergence. `res` **remains a structure**, preserving `.resInf` and `.res2` used by `refine_auto` and `refine_batch`; it adds `.relative`, `.converged`, `.returned`, and `.toleranceMode`. `TolRes`/`Tol` is relative by default; `TolStep` only detects stagnation. The old undocumented right-side gauge was not valid for a general quaternion matrix: this left-eigenpair polisher now rejects `Side='right'` rather than silently using it.

`leigqNEWTON_init_vec` remains SVD-based. Fixed-right-eigenvalue initialization requires `Gauge=false`; left-gauging was not invariant. Its help now lists the actual supported options instead of unimplemented Method/Tol choices.

### Certificates and diagnostics

Pair certificates reject zero vectors and missing vectors for a nonempty lambda list. Empty lambda/vector lists remain valid empty results. Iterative `svds` output is checked for convergence/finite values and falls back to full SVD when needed. Numeric complex inputs use the conventional (1,i) quaternion slice consistently.

## What did not change

Function names and output ordering remain stable. This is still a multistart **hit collector**, not an exhaustive-spectrum algorithm. `lambdaU` only clusters collected hits; it does not prove completeness or count infinite components. The absolute default clustering tolerance is not a scale-invariant spectral counting theorem. The zero pre-pass may return several independent vectors for one eigenvalue zero. None of these current conventions resolves the historical runner's counting conventions.

The derivative-free search stages in `refine_lambda`, `refine_auto`, and `refine_batch` retain their optimization tolerances. Those optimizer stopping criteria are not eigenpair acceptance certificates; final candidates must be assessed with the corrected certificate functions. Their calls to the corrected polisher inherit the new TolRes meaning. No broad tuning or new large benchmark is included.

## Targeted regression validation (completed)

The author has already completed this test successfully: **OVERALL: OK (29/29), 0.98 s**.
No repeat is needed solely for the post-release documentation/citation update.
Other users may optionally verify a new installation with the commands below,
starting a fresh MATLAB session after replacing MATLAB functions. Do not initialize unrelated toolboxes or clear classes/functions for this patch.

    report = test_leigqNEWTON_revision;
    save('leigqNEWTON_revision_report.mat','report');
    assert(strcmp(report.status,'OK'),'Revision regression failed.');

The driver contains 29 targeted groups: scalar/zero cases, invalid inputs, the reviewer's tiny-matrix issue, scale factors from 1e-100 to 1e100 for small examples, reporting independence, nullspace conversion and H-independence, diagonal shortcut safety, budgets, RNG restoration, polisher derivative and convergence, and diagnostics. These are small regression problems, not the historical benchmark campaign. The report records MATLAB version, platform, paths and source hashes.

For API details of changed functions, use current MATLAB help and this note. Existing rendered HTML/PDF documents and old example outputs remain historical and have not been regenerated for this revision.
