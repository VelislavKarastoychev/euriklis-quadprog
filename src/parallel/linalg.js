"use strict";

/**
 * Parallel dense linear algebra for `solveQPFast` — the O(n³) work (Cholesky
 * trailing update, triangular inverse) is routed through the SharedArrayBuffer
 * matmul in `pool.js`; the cubic-but-small panel/diagonal steps stay scalar.
 *
 * All matrices are flat row-major Float64Array views over SharedArrayBuffers so
 * the workers can read/write them in place.
 */
import { gemm } from "./pool.js";

const NB = 128; // block width — also the scalar base case for the triangular inverse

export const sabF64 = (len) => new Float64Array(new SharedArrayBuffer(len * 8));

/** Scalar Cholesky of the b×b lower diagonal block at (jb,jb) in a stride-n matrix. */
const cholDiag = (A, jb, b, n) => {
  for (let j = 0; j < b; j++) {
    const dj = jb + j;
    let d = A[dj * n + dj];
    for (let k = 0; k < j; k++) { const v = A[dj * n + jb + k]; d -= v * v; }
    if (d <= 0) throw new Error("not positive definite");
    d = Math.sqrt(d);
    A[dj * n + dj] = d;
    for (let i = j + 1; i < b; i++) {
      const ri = jb + i;
      let s = A[ri * n + dj];
      for (let k = 0; k < j; k++) s -= A[ri * n + jb + k] * A[dj * n + jb + k];
      A[ri * n + dj] = s / d;
    }
  }
};

/** Scalar TRSM of the m×b panel below the diagonal block: L21 = A21·L11⁻ᵀ. */
const trsmPanel = (A, jb, b, n, m) => {
  for (let ii = 0; ii < m; ii++) {
    const i = jb + b + ii;
    for (let j = 0; j < b; j++) {
      const dj = jb + j;
      let s = A[i * n + dj];
      for (let k = 0; k < j; k++) s -= A[i * n + jb + k] * A[dj * n + jb + k];
      A[i * n + dj] = s / A[dj * n + dj];
    }
  }
};

/**
 * Blocked right-looking Cholesky, IN PLACE: the lower triangle of A (flat n×n,
 * SAB) becomes L with A = L·Lᵀ. The trailing update A22 −= L21·L21ᵀ is the only
 * Θ(n³) term and runs on the parallel matmul (syrk via transB, α=−1, β=1).
 */
export const choleskyLower = async (A, n) => {
  for (let jb = 0; jb < n; jb += NB) {
    const b = Math.min(NB, n - jb);
    cholDiag(A, jb, b, n);
    const m = n - jb - b;
    if (m <= 0) continue;
    trsmPanel(A, jb, b, n, m);
    // P = L21 (m×b), contiguous; G = P·Pᵀ (m×m) on the pool; A22 −= G.
    const P = sabF64(m * b);
    for (let i = 0; i < m; i++) for (let k = 0; k < b; k++) P[i * b + k] = A[(jb + b + i) * n + jb + k];
    const G = sabF64(m * m);
    await gemm(P.buffer, P.buffer, G.buffer, m, m, b, { transB: true, alpha: 1, beta: 0 });
    for (let i = 0; i < m; i++) {
      const base = (jb + b + i) * n + jb + b;
      const gi = i * m;
      for (let j = 0; j < m; j++) A[base + j] -= G[gi + j];
    }
  }
};

/** Scalar upper-triangular inverse (LINPACK dpori), flat stride-n, IN PLACE. */
const dporiFlat = (a, n) => {
  for (let k = 0; k < n; k++) {
    a[k * n + k] = 1.0 / a[k * n + k];
    const t = -a[k * n + k];
    for (let i = 0; i < k; i++) a[i * n + k] *= t;
    const kp1 = k + 1;
    if (kp1 >= n) break;
    for (let j = kp1; j < n; j++) {
      const tj = a[k * n + j];
      a[k * n + j] = 0.0;
      for (let i = 0; i <= k; i++) a[i * n + j] += tj * a[i * n + k];
    }
  }
};

/**
 * Inverse of a contiguous upper-triangular R (n×n, SAB), IN PLACE → R⁻¹. Blocked
 * divide-and-conquer: R = [[A,B],[0,C]] ⇒ R⁻¹ = [[A⁻¹, −A⁻¹·B·C⁻¹],[0,C⁻¹]]. The
 * diagonal blocks invert recursively (scalar at the base); the coupling term is
 * two parallel matmuls.
 */
export const triInvUpper = async (R, n) => {
  if (n <= NB) { dporiFlat(R, n); return; }
  const k = n >> 1, m = n - k;
  const A = sabF64(k * k), B = sabF64(k * m), C = sabF64(m * m);
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) A[i * k + j] = R[i * n + j];
    for (let j = 0; j < m; j++) B[i * m + j] = R[i * n + k + j];
  }
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) C[i * m + j] = R[(k + i) * n + k + j];
  await triInvUpper(A, k);
  await triInvUpper(C, m);
  const T = sabF64(k * m);
  await gemm(A.buffer, B.buffer, T.buffer, k, m, k);        // T = A⁻¹·B
  const Bp = sabF64(k * m);
  await gemm(T.buffer, C.buffer, Bp.buffer, k, m, m);       // B' = T·C⁻¹
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) R[i * n + j] = A[i * k + j];
    for (let j = 0; j < m; j++) R[i * n + k + j] = -Bp[i * m + j];
  }
  for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) R[(k + i) * n + k + j] = C[i * m + j];
};

/** Solve D·x = d for SPD D given its lower Cholesky L (D = L·Lᵀ), flat stride-n. */
export const cholSolve = (L, d, n) => {
  const y = new Float64Array(n);
  for (let i = 0; i < n; i++) {            // forward: L·y = d
    let s = d[i];
    for (let k = 0; k < i; k++) s -= L[i * n + k] * y[k];
    y[i] = s / L[i * n + i];
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {       // back: Lᵀ·x = y
    let s = y[i];
    for (let k = i + 1; k < n; k++) s -= L[k * n + i] * x[k];
    x[i] = s / L[i * n + i];
  }
  return x;
};
