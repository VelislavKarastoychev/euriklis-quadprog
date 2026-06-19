"use strict";

/**
 * Minimal worker pool for SharedArrayBuffer matmul — the only parallel primitive
 * `solveQPFast` needs (the blocked Cholesky trailing update and the triangular
 * inverse both reduce to it). Workers are spawned lazily on first use, reused
 * across calls, and `unref`'d so they never keep the process alive.
 *
 * Runtime: `worker_threads` works identically on Node and Bun.
 */
import { Worker } from "worker_threads";
import os from "os";

const POOL = Math.max(1, (os.cpus()?.length || 4) - 2);
// Below this FLOP count the dispatch round-trip costs more than the matmul, so
// run inline on the main thread instead. Mirrors the Tensor library's threshold.
const PARALLEL_FLOPS = 1 << 21; // ~2.1M

let workers = null;

const ensure = () => {
  if (workers) return workers;
  const url = new URL("./worker.js", import.meta.url);
  workers = [];
  for (let i = 0; i < POOL; i++) {
    const w = new Worker(url);
    w.unref();
    workers.push(w);
  }
  return workers;
};

/** Inline fallback: C = beta·C + alpha · A·(B or Bᵀ), rows [r0,r1). */
const inlineGemm = (A, B, C, N, K, r0, r1, transB, alpha, beta) => {
  if (transB) {
    for (let i = r0; i < r1; i++) {
      const ai = i * K, ci = i * N;
      for (let j = 0; j < N; j++) {
        let s = 0; const bj = j * K;
        for (let k = 0; k < K; k++) s += A[ai + k] * B[bj + k];
        C[ci + j] = beta * C[ci + j] + alpha * s;
      }
    }
  } else {
    for (let i = r0; i < r1; i++) {
      const ai = i * K, ci = i * N;
      for (let j = 0; j < N; j++) {
        let s = 0;
        for (let k = 0; k < K; k++) s += A[ai + k] * B[k * N + j];
        C[ci + j] = beta * C[ci + j] + alpha * s;
      }
    }
  }
};

/**
 * C[M×N] = beta·C + alpha · A[M×K] · (B or Bᵀ). All three are SharedArrayBuffers
 * (raw buffers, not views). Splits the M rows across the pool; small products run
 * inline. Resolves when C is fully written.
 */
export const gemm = async (aSab, bSab, cSab, M, N, K, opts = {}) => {
  const { transB = false, alpha = 1, beta = 0 } = opts;
  if (2 * M * N * K < PARALLEL_FLOPS) {
    inlineGemm(new Float64Array(aSab), new Float64Array(bSab), new Float64Array(cSab), N, K, 0, M, transB, alpha, beta);
    return;
  }
  const ws = ensure();
  const chunk = Math.ceil(M / ws.length);
  const tasks = [];
  for (let w = 0; w < ws.length; w++) {
    const r0 = w * chunk, r1 = Math.min(M, r0 + chunk);
    if (r0 >= r1) break;
    const worker = ws[w];
    tasks.push(new Promise((resolve) => {
      worker.once("message", resolve);
      worker.postMessage({ aSab, bSab, cSab, N, K, r0, r1, transB, alpha, beta });
    }));
  }
  await Promise.all(tasks);
};

/** Terminate the pool (lets a short-lived script exit immediately if desired). */
export const shutdown = () => {
  if (workers) { for (const w of workers) w.terminate(); workers = null; }
};
