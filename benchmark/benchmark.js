"use strict";

/**
 * Benchmark for the current runtime (Node or Bun). Times three solvers across a
 * range of problem sizes:
 *
 *   • Santini   — the reference npm `quadprog` (scalar, 1‑indexed)
 *   • solveQP   — this package's scalar solver
 *   • solveQPFast — this package's worker‑pool solver
 *
 * Each solver is given a fresh copy of the inputs (Santini and solveQP mutate /
 * copy internally — apples to apples). Reported time is the median over a few
 * repetitions, in milliseconds.
 *
 *   node benchmark/benchmark.js          # human table
 *   node benchmark/benchmark.js --json   # machine‑readable (used by compare.js)
 */
import os from "os";
import { solveQP, solveQPFast } from "../index.js";

// Santini's `quadprog` is an optional devDependency — skip gracefully if absent.
let santini = null;
try {
  const m = await import("quadprog");
  santini = m.solveQP ?? m.default?.solveQP ?? m.default ?? null;
} catch { /* not installed */ }

const RUNTIME = typeof globalThis.Bun !== "undefined" ? `Bun ${globalThis.Bun.version}` : `Node ${process.version}`;
const CPU = os.cpus()?.[0]?.model?.trim() ?? "unknown CPU";
const THREADS = os.cpus()?.length ?? 0;

const now = () => Number(process.hrtime.bigint()) / 1e6;

let s;
const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
function gen(n, seed) {
  s = seed;
  const B = Array.from({ length: n }, () => Array.from({ length: n }, () => rnd() * 2 - 1));
  const D = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => { let x = 0; for (let k = 0; k < n; k++) x += B[i][k] * B[j][k]; return x + (i === j ? n : 0); }));
  const d = Array.from({ length: n }, () => rnd() * 2 - 1);
  const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  const b = new Array(n).fill(-1);
  return { D, d, A, b };
}
const M1 = (r) => { const m = [[]]; r.forEach((x, i) => (m[i + 1] = [0, ...x])); return m; };
const V1 = (a) => [0, ...a];
const cloneM = (m) => m.map((r) => r.slice());

async function median(fn, reps) {
  const ts = [];
  for (let i = 0; i < reps; i++) { const t = now(); await fn(); ts.push(now() - t); }
  ts.sort((a, b) => a - b);
  return ts[Math.floor(reps / 2)];
}

const SIZES = [50, 100, 200, 400, 600, 800, 1024];
const results = [];
for (const n of SIZES) {
  const reps = n <= 200 ? 7 : n <= 600 ? 4 : 2;
  const p = gen(n, 1234);
  // warm up each path (JIT + worker pool)
  solveQP(p.D, p.d, p.A, p.b);
  await solveQPFast(p.D, p.d, p.A, p.b);
  if (santini) santini(M1(p.D), V1(p.d), M1(p.A), V1(p.b));

  const mine = await median(() => { solveQP(p.D, p.d, p.A, p.b); }, reps);
  const fast = await median(() => solveQPFast(p.D, p.d, p.A, p.b), reps);
  let sant = null;
  if (santini) {
    const D1 = M1(p.D), d1 = V1(p.d), A1 = M1(p.A), b1 = V1(p.b);
    sant = await median(() => { santini(cloneM(D1), d1.slice(), cloneM(A1), b1.slice()); }, reps);
  }
  results.push({ n, santini: sant, mine, fast });
}

if (process.argv.includes("--json")) {
  console.log(JSON.stringify({ runtime: RUNTIME, cpu: CPU, threads: THREADS, results }));
} else {
  const f = (v) => (v == null ? "    n/a" : v.toFixed(2).padStart(9));
  console.log(`\nCPU: ${CPU}  (${THREADS} threads)`);
  console.log(`Runtime: ${RUNTIME}\n`);
  console.log("    n │   Santini    solveQP       Fast │  Fast vs solveQP");
  console.log("──────┼──────────────────────────────────┼─────────────────");
  for (const r of results) {
    const ratio = (r.mine / r.fast).toFixed(2) + "×";
    console.log(`${String(r.n).padStart(5)} │${f(r.santini)} ${f(r.mine)} ${f(r.fast)} │  ${ratio.padStart(8)}`);
  }
  console.log("\n(times are median ms; lower is better)");
}

process.exit(0);
