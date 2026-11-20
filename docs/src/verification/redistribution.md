# [Verification: spline-projection redistribution](@id rd-verification)

The note *Symplectic redistribution of weighted particles in phase space* implements every
analytical identity, every algorithmic step and every negative control with `SimpleSplines.jl` in
the script `verify_symplectic_redistribution.jl` (testsets R1–R5, 34 assertions) on a 1D-1V
domain ``[0, 2\pi) \times [-3.5, 3.5]``, periodic in space and clamped in velocity. The table
lists each statement with the control that must fail when a structural hypothesis is dropped.

| Testset | Statement | Positive check | Negative control |
|:--|:--|:--|:--|
| R1 | [pull-back](@ref presymplectic-form), [fold](@ref velocity-momentum-fold) | ``\|\Pi^*\omega - \omega_{\mathrm{analytic}}\| < 10^{-13}`` | Maxwellian ``D_{\bar{v}} P_{\bar{x}}`` vanishes at ``\|\bar{v}\| = v_{\mathrm{th}}`` and changes sign for ``\|\bar{v}\| > v_{\mathrm{th}}`` |
| R2 | [``L^2`` lemma](@ref rd-l2-lemma), [KKT](@ref rd-kkt-projection) | KKT spline gives ``\bar{w}_a > 0``, mass and momentum errors ``< 10^{-14}``, lower density CV | Control D: ``\min f_h^{L^2} < 0``; Control C: raw evaluation leaves ``\|\sum \bar{p}^{\mathrm{raw}}_a - \sum p_a\| > 10^{-4}`` |
| R3 | [obstructions](@ref rd-obstructions), [symmetric Hamiltonians](@ref rd-symmetric-hamiltonians) | pairwise quartic spline ``H_\theta``: ``\|\sum_a \nabla_{x_a} H_\theta\| < 10^{-15}``, ``S_{N}``-equivariance ``< 10^{-14}`` | Control A: uncoupled ``h(x, p)`` breaks ``\sum p_a``; Control B: mean-force subtraction breaks ``\Psi^*\omega = \omega``; Control C: raw ``\bar{p}^{\mathrm{raw}}`` makes ``G\theta = \Delta z`` inconsistent |
| R4 | [parameter solver](@ref rd-linear-parameter-solver) | forward Newton midpoint and symplectic Euler hit ``(\bar{x}, \bar{p})`` to ``< 10^{-13}``, IFT versus FD Jacobian ``< 10^{-9}``, ``\|\sum \bar{p}_a - \sum p_a\| < 10^{-15}`` | unregularised large displacements that violate the contraction condition cross caustics onto spurious roots |
| R5 | [closed-loop identity](@ref closed-loop-identity), [variational closure](@ref variational-closure) | variational cotangent-lift closure: ``\|D\Phi^{\mathsf{T}} \mathbb{J} D\Phi - \mathbb{J}\| < 10^{-12}`` with ``\|D\Theta\| > 10^{-2}`` | Control E: non-variational target matching gives ``\|D\Phi^{\mathsf{T}} \mathbb{J} D\Phi - \mathbb{J}\| > 10^{-3}``, matching the identity to ``< 10^{-8}`` |

## The note's script in a cold process

The research log records that the consolidated script never ran in a cold process: its numbers
for 16 and 36 markers came from an interactive session, and the four-marker check of R5 never
ran. Run unchanged in a cold process (Julia 1.13.4, SimpleSplines 0.4, the test environment of the
package), it gives **27 passes, 2 failures and 1 error** of 30 reached assertions:

- **R1 fails.** ``\|\Pi^*\omega - \omega_{\mathrm{analytic}}\| = 5.14 \cdot 10^{-11}`` against the
  bound ``10^{-13}``. The script differentiates ``\Pi`` by fourth-order differences with step
  ``10^{-5}``; ``\Pi`` is quadratic, so the truncation error vanishes, and the round-off
  ``\sim \varepsilon / h \approx 2 \cdot 10^{-11}`` of the step sets the error. The bound is
  below what the method can reach. The package differentiates ``\Pi`` exactly (error ``0``) and
  bounds the difference quotient by ``100 \varepsilon / h``.
- **R4 fails for symplectic Euler.** The forward Newton solve converges (residual
  ``1.3 \cdot 10^{-16}``) to a root ``0.184`` away from ``(\bar{x}, \bar{p})`` (scaled norm). At
  the target, ``\| \nabla_x \nabla_{\bar{p}} H_{\theta^*} \| = 1.39 > 1``: the invertibility
  condition of the proposition fails, and the identity-connected branch ends on another root.
  This is the caustic of the [remark](@ref rd-coarse-basis) (a), at the script's own relaxation
  step. The research log reports the same miss (``0.150``) as the caustic control "at the
  larger unregularised step ``0.02 \times 8``"; the script reproduces ``0.150`` for that step,
  and the target is missed at the script's step ``0.006 \times 10`` as well.
- **R5 errors.** `PairwiseSplineHamiltonian(Lx, 0.25; nx_cells = 4, np_cells = 4)` throws
  `ArgumentError: a periodic B-spline basis of degree 4 needs more than 4 cells, got n = 4`
  in SimpleSplines 0.4, before any assertion of R5. With ``5`` cells, the only change, the
  note's script runs R5: the target
  matching control gives a defect ``0.396`` and a relative error ``5.96 \cdot 10^{-10}`` of the
  identity (both pass), and the variational closure gives
  ``\|D\Phi^{\mathsf{T}} \mathbb{J} D\Phi - \mathbb{J}\| = 5.03 \cdot 10^{-10}``, which **fails**
  the bound ``10^{-12}``: the outer Jacobian of fourth-order differences with step ``10^{-5}``,
  through an inner Newton solve that itself uses difference Jacobians, has a round-off floor near
  ``10^{-10}``.

## How the package passes for the right reason

- R1 compares the pull-back by automatic differentiation with the analytic form to round-off, and
  the difference quotient with the bound of its round-off.
- R4 keeps the single symplectic Euler step as the caustic control: it asserts that the
  invertibility condition fails (``> 1``), that the step misses the target by more than ``0.1``,
  and that the target solves the step equations to round-off. The positive check runs two
  [substeps](@ref rd-solver), each with the condition below one (the larger is ``0.78``), and
  reaches the target to ``2.2 \cdot 10^{-15}``.
- R5 uses ``5`` cells, and computes the Jacobians of the closed-loop maps by implicit
  differentiation of the step equations. The variational defect is then
  ``6.5 \cdot 10^{-16}`` and the identity holds to ``3.6 \cdot 10^{-15}`` relative. The
  difference Jacobians stay as a second check, with the bound ``10^{-8}`` of their round-off.
- A self pair ``a = b`` of the pairwise Hamiltonian has the constant relative position ``0``;
  SimpleSplines 0.4 returns `NaN` partials for a dual number at the knot ``0`` of a periodic
  basis, so the package evaluates the self pair at the plain number ``0``. The values are
  unchanged.

## The numbers of the note and of the research log

`scripts/redistribution/verify_symplectic_redistribution.jl` prints every number below and checks
every statement and control (71 assertions, 52 s in a cold process). It runs the relaxation of the
script (ten steps of size ``0.006``) and the relaxation that produced the numbers of the research
log (eight steps of size ``0.02``). In the table, "log" is the research log, "script" and "log
set-up" are the two relaxations of the package; **bold** marks a quoted number the package does
not reproduce.

| Quantity | Quoted | Script set-up | Log set-up |
|:--|:--|:--|:--|
| ``\|\Pi^*\omega - \omega\|``, differences | ``< 10^{-13}`` (note) | **``3.9 \cdot 10^{-11}``** | |
| kernel-density CV, before → after | ``0.220 \to 0.158`` (log) | ``0.2199 \to 0.1939`` | ``0.2199 \to 0.1580`` |
| ``\min \bar{w}_a`` | ``0.0227`` (log) | **``0.01975``** | **``0.01680``** |
| mass error of ``\bar{w}`` | ``1.11 \cdot 10^{-16}`` (log) | ``1.1 \cdot 10^{-16}`` | ``0`` |
| momentum error of ``\bar{p}`` | ``0.0`` (log) | ``5.6 \cdot 10^{-17}`` | ``2.8 \cdot 10^{-16}`` |
| Control D, ``\min f_h^{L^2}``, 16 markers | ``< 0`` (note) | ``-0.1327`` | ``-0.1327`` |
| Control D, ``\min f_h^{L^2}``, 36 markers | ``-0.140`` (log) | **``-0.0551``** | |
| Control C, raw momentum defect | ``0.215`` (log) | **``0.01417``** | **``0.01117``** |
| Control C, least-squares residual of ``G\theta = \Delta z`` | ``> 10^{-4}`` (note) | ``3.54 \cdot 10^{-3}`` | ``2.79 \cdot 10^{-3}`` |
| rank of ``G`` | ``31`` (log) | ``31`` | ``31`` |
| ``\|\sum_a \nabla_{x_a} H_\theta\|`` | ``4.86 \cdot 10^{-17}`` (log) | ``7.6 \cdot 10^{-17}`` | ``4.86 \cdot 10^{-17}`` |
| permutation error | ``2.10 \cdot 10^{-16}`` (log) | ``1.8 \cdot 10^{-16}`` | ``2.10 \cdot 10^{-16}`` |
| Control A, momentum change | ``> 10^{-2}`` (note) | ``0.157`` | ``0.157`` |
| Control B, symplecticity defect | ``> 10^{-2}`` (note) | ``0.0914`` | ``0.0914`` |
| midpoint target error | ``8.04 \cdot 10^{-15}`` (log) | ``9.9 \cdot 10^{-15}`` | ``8.04 \cdot 10^{-15}`` |
| midpoint momentum error | ``1.67 \cdot 10^{-16}`` (log) | ``5.6 \cdot 10^{-17}`` | ``1.67 \cdot 10^{-16}`` |
| midpoint contraction number | ``< 1`` (note, sufficient) | ``2.27`` | ``2.01`` |
| midpoint IFT versus FD Jacobian | ``6.48 \cdot 10^{-10}`` (log) | ``5.8 \cdot 10^{-10}`` | ``6.48 \cdot 10^{-10}`` |
| midpoint FD symplecticity | ``2.73 \cdot 10^{-8}`` (log) | ``1.5 \cdot 10^{-8}`` | ``2.73 \cdot 10^{-8}`` |
| symplectic Euler target error, one step | ``< 10^{-13}`` (note) | **``0.184``** | ``0.150`` (log: caustic) |
| symplectic Euler, ``\|\nabla_x \nabla_{\bar{p}} H\|`` | ``< 1`` (note) | ``1.39`` | ``1.47`` |
| symplectic Euler target error, two substeps | | ``2.2 \cdot 10^{-15}`` | ``6.3 \cdot 10^{-15}`` |
| symplectic Euler IFT versus FD Jacobian | ``2.47 \cdot 10^{-10}`` (log) | ``3.5 \cdot 10^{-10}`` | **``2.53 \cdot 10^{-10}``** |
| symplectic Euler FD symplecticity | ``2.90 \cdot 10^{-8}`` (log) | ``1.3 \cdot 10^{-8}`` | **``2.72 \cdot 10^{-8}``** |
| Control E defect | ``> 10^{-3}`` (note) | ``0.396`` | |
| Control E identity, FD / exact | ``< 10^{-8}`` (note) | ``6.0 \cdot 10^{-10}`` / ``3.6 \cdot 10^{-15}`` | |
| variational closure defect, FD / exact | ``< 10^{-12}`` (note) | **``5.3 \cdot 10^{-10}``** / ``6.5 \cdot 10^{-16}`` | |
| variational closure ``\|D\Theta\|`` | ``> 10^{-2}`` (note) | ``1.04`` | |

The research log's 16-marker numbers for the relaxation, the symmetry, the midpoint step and its
Jacobian reproduce to the quoted digits with eight relaxation steps of size ``0.02``, the step that
the log names for the caustic; the script carries ten steps of size ``0.006``. The minimum weight
``0.0227``, the raw momentum defect ``0.215``, the 36-marker overshoot ``-0.140`` and the two
symplectic Euler Jacobian numbers do not reproduce with either relaxation; the set-up of the
interactive session that produced them is not recorded. The round-off numbers (errors of order
``10^{-16}``) depend on the order of the floating-point operations and are compared with the
bounds, not with their digits.

## The test files

| File | Checks | Controls that must fail |
|:--|:--|:--|
| `test/redistribution/canonical.jl` | R1: pull-back of ``\omega`` (exact and by differences), its kernel of Euler scaling fields, the factorisation freedom of a canonical map followed by a state-dependent positive split, the fold determinant in 1D and 2D, ``p_{\max}``, the two branches | a split that keeps the old weights in the velocity changes ``\bar{p}`` and breaks the pull-back (residual ``> 10^{-2}``) |
| `test/redistribution/projection.jl` | R2: ``L^2`` moments, KKT weights and moments, rank of ``A``, KKT stationarity, branch recovery | Control D (``\min f_h^{L^2} < -0.1``), Control C (raw momentum defect) |
| `test/redistribution/relaxation.jl` | R2: lower clustering, decreasing energy at every step, the log's CV ``0.220 \to 0.158`` with eight steps of size ``0.02``, the soft wall | |
| `test/redistribution/hamiltonians.jl` | R3: symmetry, derivatives against automatic differentiation, rank ``2Nd - d = 31`` and ``\mathcal{M}_0``, Fourier-moment and centroid-frame forms | Controls A, B, C |
| `test/redistribution/solver.jl` | R4: midpoint hit, symplecticity, IFT versus FD, substeps, pseudo-time path | the single symplectic Euler step: contraction number ``> 1``, target missed by more than ``0.1`` (caustic) |
| `test/redistribution/closedloop.jl` | R5: Control E (midpoint and symplectic Euler) against the identity, variational closure (exact and by differences), the stationarity of ``\theta`` | Control E (defect ``> 10^{-3}``), fixed-grid no-go (``(D\Phi)^{\mathsf{T}} \mathbb{J} (D\Phi) = 0``) |
| `test/redistribution/method.jl` | the whole redistribution through `remap`: the set-up, one midpoint step to the targets with exact weight and momentum, two symplectic Euler substeps, the canonical step at fixed parameters, the rejection of markers in more than one dimension | the single symplectic Euler step through `remap` (caustic) |

## The scripts

| Script | What it does |
|:--|:--|
| `scripts/redistribution/verify_symplectic_redistribution.jl` | runs R1–R5 for both relaxations, prints every number of the table above next to the quoted value, and checks the statements and the controls (71 assertions) |
| `scripts/redistribution/make_figures.jl` | writes the figures of the velocity–momentum fold, the spline projections, the target relaxation, the redistributed markers, the Newton convergence and the closed-loop defect (PDF), and the animations of the redistribution step and the target relaxation (MP4), to `results/redistribution/` |
