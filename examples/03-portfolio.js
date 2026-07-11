"use strict";

/**
 * Long‑only minimum‑variance portfolio (a classic constrained QP).
 *
 *   minimise  ½ wᵀ Σ w        (portfolio variance)
 *   subject to  Σ wᵢ = 1      (budget — fully invested)        ← equality
 *               wᵢ ≥ 0        (no short selling)               ← inequalities
 *
 * In quadprog form: D = Σ, d = 0, and A is n × (n+1):
 *   column 0     = ones        (the budget equality, so meq = 1)
 *   columns 1..n = identity    (the non‑negativity bounds)
 * with b = [1, 0, …, 0].
 *
 * Run:  node examples/03-portfolio.js
 */
import { solveQP } from "../dist/index.js";

// covariance matrix Σ of 4 assets (symmetric positive‑definite)
const Sigma = [
  [0.040, 0.006, 0.004, 0.002],
  [0.006, 0.090, 0.008, 0.003],
  [0.004, 0.008, 0.0625, 0.005],
  [0.002, 0.003, 0.005, 0.160],
];
const n = Sigma.length;

const d = new Array(n).fill(0);
const A = Sigma.map((_, i) => {
  const row = [1];                         // budget column
  for (let j = 0; j < n; j++) row.push(i === j ? 1 : 0); // identity (wᵢ ≥ 0)
  return row;
});
const b = [1, ...new Array(n).fill(0)];

const r = solveQP(Sigma, d, A, b, /* meq = */ 1);

const w = r.solution;
console.log("weights  :", w.map((v) => v.toFixed(4)));
console.log("Σ weights:", w.reduce((s, v) => s + v, 0).toFixed(6), "(= 1, budget)");
console.log("all ≥ 0  :", w.every((v) => v >= -1e-9));
console.log("variance :", (2 * r.value + d.reduce((s, di, i) => s + di * w[i], 0)).toFixed(6));
