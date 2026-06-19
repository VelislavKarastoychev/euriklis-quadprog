"use strict";

/**
 * solveQPFast — the same Goldfarb-Idnani solver as `solveQP`, but the O(n³) init
 * (Cholesky of D, then R⁻¹) runs on a SharedArrayBuffer worker pool instead of
 * the scalar LINPACK routines. The active-set iteration itself is inherently
 * sequential and shared verbatim with `solveQP` (so both give identical results).
 *
 * The parallel factorisation only pays off once D is large enough to amortise the
 * worker dispatch (≈ n ≥ 512). Below FAST_MIN_N this transparently calls the
 * scalar `solveQP` — there is never a reason to prefer the workers for small D.
 *
 * Same signature and return shape as `solveQP(D, d, A, b, meq)`; async.
 */
import { quadprog as solveQP, qpgen2 } from "./quadProg.js";
import { sabF64, choleskyLower, triInvUpper } from "./parallel/linalg.js";

const FAST_MIN_N = 512;

export const solveQPFast = async (D, d, A, b, meq = 0) => {
  const n = D.length;
  if (n < FAST_MIN_N) return solveQP(D, d, A, b, meq);     // scalar wins; no dispatch
  const q = A[0].length;
  // qpgen2 negates equality columns of A and entries of b in place → deep-copy.
  const amat = A.map((r) => r.slice());
  const bvec = b.slice();

  // ── parallel init: D = L·Lᵀ → R = Lᵀ (upper) → J = R⁻¹ ──────────────────────
  const Df = sabF64(n * n);
  for (let i = 0; i < n; i++) { const Di = D[i]; for (let j = 0; j < n; j++) Df[i * n + j] = Di[j]; }
  await choleskyLower(Df, n);                               // Df lower ← L
  const R = sabF64(n * n);
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) R[i * n + j] = Df[j * n + i]; // R = Lᵀ
  await triInvUpper(R, n);                                  // R ← J = R⁻¹

  // J as a 2-D array (upper triangle) for qpgen2's pre-factorised path (ierr=1):
  // it derives x₀ = D⁻¹d from J and d itself (two O(n²) triangular products).
  const J = new Array(n);
  for (let i = 0; i < n; i++) {
    const row = new Array(n).fill(0);
    for (let j = i; j < n; j++) row[j] = R[i * n + j];
    J[i] = row;
  }
  return qpgen2(J, d.slice(), n, amat, bvec, q, meq, 1);
};
