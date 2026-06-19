"use strict";

/**
 * solveQPFast on a large dense problem. Same API and same result as solveQP,
 * but the O(n³) start‑up (Cholesky of D, then R⁻¹) runs on a worker pool.
 * For n < 512 it transparently calls solveQP, so it is always safe to use.
 *
 * Run:  node examples/04-fast-large.js   (or:  bun examples/04-fast-large.js)
 */
import { solveQP, solveQPFast } from "../index.js";

const n = 700;

// a random SPD D = BBᵀ + n·I, random linear term, simple lower bounds xᵢ ≥ −1
let s = 12345;
const rnd = () => { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; };
const B = Array.from({ length: n }, () => Array.from({ length: n }, () => rnd() * 2 - 1));
const D = Array.from({ length: n }, (_, i) =>
  Array.from({ length: n }, (_, j) => {
    let x = 0; for (let k = 0; k < n; k++) x += B[i][k] * B[j][k];
    return x + (i === j ? n : 0);
  }));
const d = Array.from({ length: n }, () => rnd() * 2 - 1);
const A = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
const b = new Array(n).fill(-1);

// Warm up: the first call spins up the worker pool (tens of ms). Time the
// steady state, which is what repeated / batch solving actually sees.
await solveQPFast(D, d, A, b);
solveQP(D, d, A, b);

const t0 = performance.now();
const fast = await solveQPFast(D, d, A, b);
const t1 = performance.now();
const slow = solveQP(D, d, A, b);
const t2 = performance.now();

let maxDiff = 0;
for (let i = 0; i < n; i++) maxDiff = Math.max(maxDiff, Math.abs(fast.solution[i] - slow.solution[i]));

console.log(`n = ${n}`);
console.log(`solveQPFast : ${(t1 - t0).toFixed(1)} ms`);
console.log(`solveQP     : ${(t2 - t1).toFixed(1)} ms`);
console.log(`speed‑up    : ${((t2 - t1) / (t1 - t0)).toFixed(2)}×`);
console.log(`max |Δx|    : ${maxDiff.toExponential(2)} (identical up to rounding)`);

process.exit(0); // release the worker pool
