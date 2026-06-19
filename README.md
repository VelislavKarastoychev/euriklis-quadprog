# @euriklis/quadprog

A small, dependency‑free **convex quadratic‑programming** solver for JavaScript
(Node and Bun). It implements the **Goldfarb–Idnani dual active‑set method** and
ships two entry points:

- `solveQP`  — the pure‑scalar solver. Zero dependencies, deterministic, fast.
- `solveQPFast` — the same algorithm, but the cubic matrix factorisation in the
  start‑up phase runs on a `SharedArrayBuffer` worker pool. Pays off for large
  dense problems (`n ≳ 512`); below that it transparently calls `solveQP`.

Both return **bit‑for‑bit identical** results.

It solves

```
minimise    ½ xᵀ D x − dᵀ x          x ∈ ℝⁿ
subject to  A₁ᵀ x  =  b₁             (the first  meq  constraints)
            A₂ᵀ x  ≥  b₂             (the remaining ones)
```

with `D` symmetric positive‑definite. This covers portfolio optimisation with
constraints, constrained least squares, the dual problem of **support‑vector
machines**, and RBF networks, among many others.

---

## Installation

```sh
npm install @euriklis/quadprog@latest --save
```

## Usage

`A` is an `n × q` matrix whose **column `i` is the normal of constraint `i`**
(so `A[variable][constraint]`). `b` has length `q`, `d` length `n`, `D` is
`n × n`. `meq` (default `0`) is the number of leading **equality** constraints.

```js
import { solveQP } from "@euriklis/quadprog";

// minimise ½‖x‖² − dᵀx  s.t.  the three inequalities Aᵀx ≥ b
const D = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const d = [0, 5, 0];
const A = [[-4, 2, 0], [-3, 1, -2], [0, 0, 1]];
const b = [-8, 2, 0];

const r = solveQP(D, d, A, b);
// r.solution                  → [0.476190, 1.047619, 2.095238]
// r.value                     → −2.380952   (½xᵀDx − dᵀx at the optimum)
// r["Lagrangian multipliers"] → [0, 0.238095, 2.095238]
// r["active constraints"]     → [2, 1, 0]   (first two are active)
// r["count of active constraints"] → 2
```

For large dense problems, `solveQPFast` is a drop‑in async replacement:

```js
import { solveQPFast } from "@euriklis/quadprog";

const r = await solveQPFast(D, d, A, b);   // n < 512 → delegates to solveQP
```

### Return value

| field | meaning |
|---|---|
| `solution` | the minimiser `x*` (length `n`) |
| `value` | the objective `½x*ᵀDx* − dᵀx*` |
| `unconstrained_solution` | `D⁻¹d`, the minimiser ignoring all constraints |
| `Lagrangian multipliers` | one per constraint, `0` if inactive |
| `active constraints` / `count of active constraints` | which constraints bind at `x*` |
| `iterations` | `[main iterations, constraints dropped]` |
| `message` | `"No problems"`, or why it stopped |

## Examples

Runnable scripts live in [`examples/`](./examples) (run with `node` or `bun`):

| file | what it shows |
|---|---|
| [`01-basic.js`](./examples/01-basic.js) | a plain inequality‑constrained QP |
| [`02-equality.js`](./examples/02-equality.js) | an equality constraint via `meq` (projection onto a line) |
| [`03-portfolio.js`](./examples/03-portfolio.js) | long‑only minimum‑variance portfolio (budget equality + no‑short bounds) |
| [`04-fast-large.js`](./examples/04-fast-large.js) | `solveQPFast` on a large dense problem, vs `solveQP` |

```sh
node examples/03-portfolio.js
# weights  : [ '0.4530', '0.1693', '0.2673', '0.1105' ]
# Σ weights: 1.000000 (= 1, budget)
# all ≥ 0  : true
# variance : 0.020424
```

---

## The algorithm

This section explains *why* the method works. It assumes you are comfortable
with positive‑definite matrices, the Cholesky factorisation, QR / Givens
rotations, and Lagrange multipliers.

> New to the **[active set](./docs/concepts.md#the-active-set)** or the **[dual
> problem](./docs/concepts.md#the-dual-problem)**? Both get a plain‑language,
> picture‑first primer in [`docs/concepts.md`](./docs/concepts.md) — read that
> first if either term is unfamiliar.

### 1. The problem is a single point

Write the objective `f(x) = ½ xᵀ D x − dᵀ x`. Because `D` is symmetric
positive‑definite, `f` is **strictly convex**, the feasible region

```
Ω = { x : aᵢᵀx ≥ bᵢ  (i ≥ meq),  aᵢᵀx = bᵢ  (i < meq) }
```

is a convex polyhedron, and the constrained minimum — if `Ω ≠ ∅` — is the
**unique** global one. There is no question of local minima.

Its gradient is `∇f(x) = D x − d`. Setting it to zero gives the **unconstrained
minimum**

```
x⁰ = D⁻¹ d.
```

If `x⁰ ∈ Ω` we are already done. The interesting case is when `x⁰` violates some
constraints and the optimum sits on the boundary of `Ω`.

### 2. What the optimum looks like — the KKT conditions

A feasible `x*` is the optimum **iff** there exist multipliers `λᵢ ≥ 0` such that

```
   D x* − d = Σ λᵢ aᵢ          (stationarity: the gradient is a combination of
                                active constraint normals)
   λᵢ ( aᵢᵀx* − bᵢ ) = 0       (complementary slackness: λᵢ = 0 unless
                                constraint i is tight)
   λᵢ ≥ 0  for inequalities    (dual feasibility)
```

Because the problem is convex, these conditions are not just necessary but
**sufficient**. So the whole job is: find the set of constraints that are tight
at the optimum — the **[active set](./docs/concepts.md#the-active-set)** `𝒜` — together with multipliers `λ ≥ 0`
that balance the gradient. Once `𝒜` is known the solution is pure linear
algebra: minimise `f` subject to `aᵢᵀx = bᵢ` for `i ∈ 𝒜`.

### 3. Two ways to search, and why the *dual* one is used

A **primal** active‑set method keeps `x` feasible and walks along the boundary,
swapping constraints in and out until the multipliers come out non‑negative. It
needs a feasible starting point, which itself costs a solve.

Goldfarb and Idnani turn this around. The **[dual](./docs/concepts.md#the-dual-problem)** method keeps the *multiplier*
side healthy (`λ ≥ 0`) and works towards feasibility:

> Start at the unconstrained minimum `x⁰ = D⁻¹d` with an empty active set
> (`λ = 0`, trivially dual‑feasible) and, one constraint at a time, repair the
> worst violation while never letting any `λᵢ` go negative.

Its great advantage: **no feasible starting point is needed** — `x⁰` is free —
and every iterate is already the exact optimum of the QP restricted to the
current active set.

### 4. The right coordinate system: `J = R⁻¹`

Factor the (constant) matrix once,

```
D = Rᵀ R         (Cholesky, R upper‑triangular),     J := R⁻¹.
```

Then `D⁻¹ = J Jᵀ` and, crucially, `Jᵀ D J = I`: the columns of `J` form an
**orthonormal basis in the inner product defined by `D`**. Working in this basis
turns the awkward `D`‑weighted geometry into ordinary Euclidean geometry.

While the algorithm runs it maintains, for the current active set of size `k`,
a QR‑type factorisation of the **active constraint normals**
`N = [ a_{𝒜(1)} , … , a_{𝒜(k)} ]` expressed in the `J`‑basis. Split the
transformed basis as `J = (Q₁ | Q₂)`, where the first `k` columns `Q₁` span the
active normals and the remaining `Q₂` span the directions still free to move.
Adding or dropping a constraint is then **one sweep of Givens rotations** that
re‑triangularises this factor in `O(n²)`, instead of refactorising from scratch
in `O(n³)`.

### 5. One iteration

Let `x` be the current point with active set `𝒜` and multipliers `u ≥ 0`.

1. **Pick the most‑violated constraint.** Among all constraints compute the
   (norm‑scaled) violation `aₚᵀx − bₚ`; let `p` be the most negative. If none is
   negative, every KKT condition holds — **stop**, `x` is the optimum.

2. **Step direction.** With `n⁺ = aₚ` the violated normal, split it in the
   `J`‑basis into the part lying in the free subspace and the part in the active
   subspace:

   ```
   z = Q₂ Q₂ᵀ J n⁺      (primal direction: how x should move)
   r = R⁻¹ Q₁ᵀ n⁺       (dual direction:  how the active λ react)
   ```

   `z` is the direction in which moving `x` reduces the violation of `p` without
   disturbing the already‑active constraints; `r` says how the current
   multipliers change as we move.

3. **Step length** `t = min(t₁, t₂)`:

   - `t₂` (**primal**) — the distance along `z` until constraint `p` becomes
     exactly tight, `t₂ = −(aₚᵀx − bₚ) / (zᵀ aₚ)`.
   - `t₁` (**dual**) — the largest step before some active multiplier would turn
     negative, `t₁ = min_{ i : rᵢ > 0 } uᵢ / rᵢ`. The minimiser `it₁` is the
     constraint that would be driven infeasible in the dual.

4. **Take the step.**

   - **Full step** (`t₂ ≤ t₁`, and finite): move `x ← x + t₂ z`, update the
     multipliers, and **add** `p` to the active set (a Givens *update* of the
     factorisation). The active set grew by one; go to 1.
   - **Partial step** (`t₁ < t₂`): we cannot reach tightness yet — the blocking
     constraint `it₁` must leave first. Update `u ← u − t₁ r`, **drop** `it₁`
     (a Givens *down‑date*), and recompute the direction for the same `p`.
   - If `z = 0` **and** no dual step is possible, the constraints are
     inconsistent — `Ω = ∅`, report infeasible.

5. Repeat. Each main iteration either makes a constraint active for good or
   removes one that should not have been; the dual objective increases
   monotonically, so the process is **finite**.

When the loop ends, the active multipliers are read off from `u`, every inactive
constraint gets `λ = 0`, and `x` is returned as `x*`.

### 6. Cost, and what `solveQPFast` parallelises

Per iteration the work is `O(n²)` (the Givens sweeps and matrix–vector
products); there are typically `O(n + q)` iterations. The **one‑off start‑up**,
however, is `O(n³)`: the Cholesky `D = RᵀR` and the triangular inverse `J = R⁻¹`.

The active‑set loop is inherently sequential, but that start‑up is not.
`solveQPFast` keeps the identical loop and only accelerates the factorisation:

- a **blocked Cholesky** whose trailing‑matrix update `A₂₂ ← A₂₂ − L₂₁L₂₁ᵀ` is a
  matrix multiply, and
- a **blocked triangular inverse**,
  `R = [[A,B],[0,C]] ⇒ R⁻¹ = [[A⁻¹, −A⁻¹BC⁻¹],[0, C⁻¹]]`, whose only cubic term
  `−A⁻¹BC⁻¹` is two matrix multiplies,

with every matrix multiply dispatched across a `SharedArrayBuffer` worker pool.
Below `n = 512` the dispatch overhead is not worth it and `solveQPFast` simply
calls `solveQP`.

---

## Correctness

Every result is validated, value for value, against Alberto Santini's
[`quadprog`](https://github.com/albertosantini/node-quadprog) — the
long‑standing, heavily‑used 1‑indexed translation of Berwin Turlach's reference
Fortran `qpgen2` — across **864 randomised problems** (square and rectangular,
`n = 1…12`, `q = 1…24`, with and without equality constraints), plus the
analytic edge cases (a single active bound, an equality‑constrained projection,
a fully‑determined active set). Solutions and Lagrange multipliers agree to
machine precision (`≈ 10⁻¹⁵`). `solveQPFast` is checked the same way for
`n ∈ {512, 768, 1024}`.

> Note for users of earlier versions: the previous 0‑indexed port mistranslated
> several packed‑storage offsets from the 1‑based Fortran (a sentinel collision
> that prevented constraint `0` from ever activating, two `meq` off‑by‑ones, the
> triangular‑inverse strides, a Givens drop stride and loop bound, and the
> constraint count for non‑square `A`). These are fixed; the solver now matches
> the reference everywhere.

## Performance

All numbers are on an **AMD Ryzen 9 4900HS (16 threads)**, Node v20.19 / Bun
1.2.23, median ms; reproduce with the scripts in [`benchmark/`](./benchmark).

**How `solveQP` compares to Santini depends on the problem.** Both share the same
`O(n³)` factorisation but differ in the `O(n²)`‑per‑iteration active‑set loop, so
the picture splits in two:

**(a) Constraint‑heavy problems** — many inequalities ⇒ many iterations ⇒ the
loop dominates, and the dense 0‑indexed, goto‑free layout pulls ahead
(`node benchmark/constraints.js`, `q = 2n` dense random):

| n | iterations | Santini | solveQP | solveQPFast | **solveQP vs Santini** |
|---:|---:|---:|---:|---:|---:|
| 100 | 125 | 53.4 | 13.2 | 13.3 | **4.0× faster** |
| 200 | 297 | 602 | 147 | 149 | **4.1× faster** |
| 400 | 491 | 4423 | 1585 | 1617 | **2.8× faster** |

Here `solveQPFast ≈ solveQP`: it parallelises only the factorisation, which is a
minor part of a constraint‑heavy solve — the (sequential) iteration loop is the
bottleneck. `solveQPFast` earns its keep in regime (b), below.

**(b) Factorisation‑dominated problems** — few active constraints ⇒ ≈ 1
iteration ⇒ the time is almost entirely the Cholesky and triangular inverse,
which is the *same* work in both, so the two are **on par**
(`node benchmark/benchmark.js`, `n` box constraints):

| n | Santini | solveQP |
|---:|---:|---:|
| 200 | 5.5 | 5.5 |
| 400 | 50 | 49 |
| 800 | 507 | 412 |
| 1024 | 1178 | 1116 |

Regime (b) is exactly where **`solveQPFast`** helps: it runs that factorisation
on the worker pool. Same box problems, `solveQPFast` vs `solveQP` (identical
results; below `n = 512` it just calls `solveQP`, hence ≈ 1× there):

| n | Node | Bun |
|---:|---:|---:|
| 600 | 1.4× | 2.9× |
| 800 | 1.6× | 5.2× |
| 1024 | 2.2× | 4.8× |

Run `node benchmark/compare.js` for the full Node‑vs‑Bun, three‑solver table.

---

## Dependencies

None. `solveQPFast` uses only the built‑in `worker_threads` and
`SharedArrayBuffer`, available in both Node and Bun.

## License

MIT. Provided free of charge; the author is not liable for any damages arising
from its use.

## Contact

Questions, bugs, suggestions: `euriklis@hotmail.com`.
