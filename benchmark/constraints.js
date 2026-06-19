"use strict";

/**
 * solveQP vs. Santini on CONSTRAINT‑HEAVY problems (q = 2n dense random
 * inequalities). These take many active‑set iterations, so the per‑iteration
 * cost dominates — which is exactly where this package's dense, goto‑free inner
 * loop pulls ahead. (On factorisation‑dominated problems — few active
 * constraints, the time spent almost entirely in the Cholesky / inverse — the
 * two are on par; see benchmark.js.)
 *
 *   node benchmark/constraints.js   (or:  bun …)
 */
import os from "os";
import { solveQP } from "../index.js";

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
const med = (fn, reps) => { const t = []; for (let i = 0; i < reps; i++) { const a = now(); fn(); t.push(now() - a); } t.sort((x, y) => x - y); return t[reps >> 1]; };

function gen(n, seed) {
  s = seed;
  const B = Array.from({ length: n }, () => Array.from({ length: n }, () => rnd() * 2 - 1));
  const D = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => { let x = 0; for (let k = 0; k < n; k++) x += B[i][k] * B[j][k]; return x + (i === j ? n : 0); }));
  const q = 2 * n;
  const d = Array.from({ length: n }, () => rnd() * 2 - 1);
  const A = Array.from({ length: n }, () => Array.from({ length: q }, () => rnd() * 2 - 1));
  const b = Array.from({ length: q }, () => rnd() * 2 - 1);
  return { D, d, A, b };
}

console.log(`\nCPU: ${CPU}   Runtime: ${RUNTIME}`);
console.log("Constraint‑heavy QP (q = 2n dense random) — median ms\n");
console.log("    n │   q  │ iters │  Santini    solveQP │  solveQP is");
console.log("──────┼──────┼───────┼─────────────────────┼────────────");
for (const n of [30, 60, 100, 150, 200]) {
  const p = gen(n, 7);
  const D1 = M1(p.D), d1 = V1(p.d), A1 = M1(p.A), b1 = V1(p.b);
  const iters = solveQP(p.D, p.d, p.A, p.b).iterations[0];
  for (let w = 0; w < 3; w++) { solveQP(p.D, p.d, p.A, p.b); santini(cM(D1), d1.slice(), cM(A1), b1.slice()); }
  const reps = 9;
  const mine = med(() => solveQP(p.D, p.d, p.A, p.b), reps);
  const sant = med(() => santini(cM(D1), d1.slice(), cM(A1), b1.slice()), reps);
  const verdict = sant > mine ? `${(sant / mine).toFixed(2)}× faster` : `${(mine / sant).toFixed(2)}× slower`;
  console.log(`${String(n).padStart(5)} │${String(2 * n).padStart(5)} │${String(iters).padStart(6)} │${sant.toFixed(2).padStart(9)}${mine.toFixed(2).padStart(11)} │  ${verdict}`);
}
console.log("");
process.exit(0);
