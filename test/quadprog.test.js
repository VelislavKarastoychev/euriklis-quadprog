"use strict";

/**
 * Correctness tests — run with `npm test` (works on Node and Bun).
 *
 * The core safety invariant is `solveQPFast === solveQP` to machine precision on
 * FEASIBLE problems (including large constraint‑heavy ones that exercise the
 * parallel factorisation path); `solveQP` itself is cross‑checked against the
 * reference `quadprog` (Santini) where that dev‑dependency is present.
 *
 * NOTE: tests use only FEASIBLE problems on purpose. Random `Aᵀx ≥ b` with many
 * dense constraints can be infeasible; on an infeasible problem every solver
 * returns a meaningless (constraint‑violating) point, and the chaotic active‑set
 * path then amplifies the tiny difference between the parallel and scalar
 * factorisations — that is garbage‑in/garbage‑out, not a solver disagreement.
 */
import assert from "node:assert/strict";
import { solveQP, solveQPFast } from "../index.js";

let santini = null;
try { const m = await import("quadprog"); santini = m.solveQP ?? m.default?.solveQP ?? m.default ?? null; } catch { /* optional */ }

// ── tiny test harness (no framework dependency) ──────────────────────────────
let passed = 0, failed = 0;
const results = [];
async function test(name, fn) {
  try { await fn(); passed++; results.push(`  ok   ${name}`); }
  catch (e) { failed++; results.push(`  FAIL ${name}\n         ${e.message}`); }
}
const maxDiff = (a, b, n) => { let m = 0; for (let i = 0; i < n; i++) m = Math.max(m, Math.abs(a[i] - b[i])); return m; };

// ── deterministic problem generators (all FEASIBLE) ──────────────────────────
let seed;
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
function spd(n, sd) {
  seed = sd;
  const B = Array.from({ length: n }, () => Array.from({ length: n }, () => rnd() * 2 - 1));
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => { let x = 0; for (let k = 0; k < n; k++) x += B[i][k] * B[j][k] / n; return x + (i === j ? 1 : 0); }));
}
/** Lower bounds xᵢ ≥ −1 (q = n). Few constraints bind → factorisation‑dominated. */
function boxLower(n, sd) {
  const D = spd(n, sd), d = Array.from({ length: n }, () => rnd() * 2 - 1);
  const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  return { D, d, A, b: new Array(n).fill(-1), q: n };
}
/** Box −1 ≤ xᵢ ≤ 1 (q = 2n) with a large linear term so many bounds activate —
 *  feasible AND iteration‑heavy (exercises the active‑set loop hundreds of times). */
function boxBoth(n, sd) {
  const D = spd(n, sd), d = Array.from({ length: n }, () => (rnd() * 2 - 1) * 3);
  const q = 2 * n, A = Array.from({ length: n }, () => new Array(q).fill(0));
  for (let i = 0; i < n; i++) { A[i][i] = 1; A[i][n + i] = -1; }
  return { D, d, A, b: new Array(q).fill(-1), q };
}
const M1 = (r) => { const m = [[]]; r.forEach((x, i) => (m[i + 1] = [0, ...x])); return m; };
const V1 = (a) => [0, ...a];
const minFeas = (A, b, x, n, q) => { let m = 0; for (let i = 0; i < q; i++) { let ax = 0; for (let j = 0; j < n; j++) ax += A[j][i] * x[j]; m = Math.min(m, ax - b[i]); } return m; };

// ── 1. analytic solveQP correctness ──────────────────────────────────────────
await test("solveQP — README example", () => {
  const r = solveQP([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 5, 0], [[-4, 2, 0], [-3, 1, -2], [0, 0, 1]], [-8, 2, 0]);
  assert.ok(maxDiff(r.solution, [0.476190476, 1.047619047, 2.095238095], 3) < 1e-6);
  assert.ok(Math.abs(r.value - -2.380952) < 1e-5);
});
await test("solveQP — single active bound (x ≥ 5)", () => {
  const r = solveQP([[1]], [0], [[1]], [5]);
  assert.ok(Math.abs(r.solution[0] - 5) < 1e-9 && r["count of active constraints"] === 1);
});
await test("solveQP — equality constraint (x+y = 1)", () => {
  const r = solveQP([[1, 0], [0, 1]], [0, 0], [[1], [1]], [1], 1);
  assert.ok(maxDiff(r.solution, [0.5, 0.5], 2) < 1e-9);
});

// ── 2. solveQPFast === solveQP on FEASIBLE problems (the safety invariant) ────
for (const [name, mk, sizes] of [
  ["box lower-bounds (factorisation-dominated)", boxLower, [50, 200, 512]],
  ["box both-sides (constraint-heavy, feasible)", boxBoth, [100, 300, 512]],
]) {
  for (const n of sizes) {
    await test(`solveQPFast === solveQP — ${name}, n=${n}`, async () => {
      const p = mk(n, 7);
      const sq = solveQP(p.D, p.d, p.A, p.b);
      const fa = await solveQPFast(p.D, p.d, p.A, p.b);
      assert.ok(minFeas(p.A, p.b, fa.solution, n, p.q) > -1e-9, "Fast solution must be feasible");
      assert.equal(fa.iterations[0], sq.iterations[0], `iteration count must match (${fa.iterations[0]} vs ${sq.iterations[0]})`);
      const d = maxDiff(fa.solution, sq.solution, n);
      assert.ok(d < 1e-9, `max |Δx| = ${d.toExponential(2)} exceeds 1e-9`);
    });
  }
}

// ── 3. cross-check against the reference (Santini), where available ──────────
if (santini) {
  for (const [n, mk] of [[60, boxLower], [120, boxBoth], [256, boxBoth]]) {
    await test(`solveQP === Santini — n=${n}`, () => {
      const p = mk(n, 11);
      const mine = solveQP(p.D, p.d, p.A, p.b).solution;
      const ref = santini(M1(p.D), V1(p.d), M1(p.A), V1(p.b)).solution.slice(1);
      const d = maxDiff(mine, ref, n);
      assert.ok(d < 1e-6, `max |Δx| vs Santini = ${d.toExponential(2)}`);
    });
  }
} else {
  results.push("  skip Santini cross-check (dev-dependency `quadprog` not installed)");
}

// ── report ───────────────────────────────────────────────────────────────────
console.log(results.join("\n"));
console.log(`\n${passed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
