"use strict";

/**
 * All three solvers on CONSTRAINT‑HEAVY problems (q = 2n dense random
 * inequalities). These take many active‑set iterations, so the per‑iteration
 * cost dominates — which is exactly where this package's dense, goto‑free inner
 * loop pulls ahead of Santini's reference. (On factorisation‑dominated problems
 * — few active constraints — the two are on par; see benchmark.js.)
 *
 * solveQPFast parallelises only the O(n³) start‑up, NOT the iteration loop. In
 * this regime the loop is the bottleneck, so Fast gives essentially no benefit
 * here — and below n=512 it simply calls solveQP. Its win is on
 * factorisation‑dominated problems (benchmark.js / compare.js).
 *
 *   node benchmark/constraints.js   (or:  bun …)
 */
import os from "os";
import { solveQP, solveQPFast } from "../index.js";

let santini = null;
try { const m = await import("quadprog"); santini = m.solveQP ?? m.default?.solveQP ?? m.default ?? null; } catch {}
if (!santini) { console.error("Install the reference for comparison:  npm i -D quadprog"); process.exit(1); }

const RUNTIME = typeof globalThis.Bun !== "undefined" ? `Bun ${globalThis.Bun.version}` : `Node ${process.version}`;
const CPU = os.cpus()?.[0]?.model?.trim() ?? "unknown CPU";

const now = () => Number(process.hrtime.bigint()) / 1e6;
let s;
const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
const M1 = (r) => { const m = [[]]; r.forEach((x, i) => (m[i + 1] = [0, ...x])); return m; };
const V1 = (a) => [0, ...a];
const cM = (m) => m.map((r) => r.slice());
const med = async (fn, reps) => { const t = []; for (let i = 0; i < reps; i++) { const a = now(); await fn(); t.push(now() - a); } t.sort((x, y) => x - y); return t[reps >> 1]; };

// FEASIBLE constraint-heavy problem: box −1 ≤ xᵢ ≤ 1 (q = 2n bounds) with a large
// linear term so the unconstrained minimum lands well outside the box → many
// bounds active → many active-set iterations. Guaranteed feasible (the box is
// non-empty), unlike random `Aᵀx ≥ b`, which over-constrains and can be empty —
// on an infeasible problem every solver returns garbage, so it is useless here.
function gen(n, seed) {
  s = seed;
  const B = Array.from({ length: n }, () => Array.from({ length: n }, () => rnd() * 2 - 1));
  const D = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => { let x = 0; for (let k = 0; k < n; k++) x += B[i][k] * B[j][k] / n; return x + (i === j ? 1 : 0); }));
  const d = Array.from({ length: n }, () => (rnd() * 2 - 1) * 3);
  const q = 2 * n, A = Array.from({ length: n }, () => new Array(q).fill(0));
  for (let i = 0; i < n; i++) { A[i][i] = 1; A[i][n + i] = -1; }
  return { D, d, A, b: new Array(q).fill(-1) };
}

console.log(`\nCPU: ${CPU}   Runtime: ${RUNTIME}`);
console.log("Constraint‑heavy QP (q = 2n box bounds, feasible) — median ms, lower is better\n");
console.log("    n │   q  │ iters │  Santini    solveQP       Fast │  solveQP vs Santini");
console.log("──────┼──────┼───────┼──────────────────────────────┼─────────────────────");
for (const n of [60, 100, 200, 400]) {
  const p = gen(n, 7);
  const D1 = M1(p.D), d1 = V1(p.d), A1 = M1(p.A), b1 = V1(p.b);
  const iters = solveQP(p.D, p.d, p.A, p.b).iterations[0];
  const reps = n <= 100 ? 9 : n <= 200 ? 6 : 4;
  for (let w = 0; w < 2; w++) { solveQP(p.D, p.d, p.A, p.b); await solveQPFast(p.D, p.d, p.A, p.b); santini(cM(D1), d1.slice(), cM(A1), b1.slice()); }
  const mine = await med(() => { solveQP(p.D, p.d, p.A, p.b); }, reps);
  const fast = await med(() => solveQPFast(p.D, p.d, p.A, p.b), reps);
  const sant = await med(() => { santini(cM(D1), d1.slice(), cM(A1), b1.slice()); }, reps);
  const f = (v) => v.toFixed(1).padStart(10);
  console.log(`${String(n).padStart(5)} │${String(2 * n).padStart(5)} │${String(iters).padStart(6)} │${f(sant)}${f(mine)}${f(fast)} │  ${((sant / mine).toFixed(2) + "× faster").padStart(13)}`);
}
console.log("\nFast ≈ solveQP here: it parallelises the factorisation, but a constraint‑heavy");
console.log("solve is dominated by the (sequential) iteration loop, so there is little to gain.");
console.log("Fast's win is on factorisation‑dominated problems — see `npm run benchmark`.");
process.exit(0);
