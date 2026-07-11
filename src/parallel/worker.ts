"use strict";

/**
 * Matmul worker. Computes a row band [r0, r1) of
 *
 *     C = beta·C + alpha · A · (B or Bᵀ)
 *
 * over SharedArrayBuffer-backed Float64Array views (so the band is written in
 * place, visible to the main thread). `transB` selects A·Bᵀ (B is N×K) — used by
 * the Cholesky trailing update (syrk: C −= P·Pᵀ via alpha=−1, beta=1, transB).
 */
import { parentPort } from "worker_threads";

type GemmMessage = {
  aSab: ArrayBufferLike;
  bSab: ArrayBufferLike;
  cSab: ArrayBufferLike;
  N: number;
  K: number;
  r0: number;
  r1: number;
  transB: boolean;
  alpha: number;
  beta: number;
};

if (!parentPort) throw new Error("worker.ts must be run as a worker thread");
const port = parentPort;

port.on("message", (m: GemmMessage) => {
  const { aSab, bSab, cSab, N, K, r0, r1, transB, alpha, beta } = m;
  const A = new Float64Array(aSab), B = new Float64Array(bSab), C = new Float64Array(cSab);
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
  port.postMessage(1);
});
