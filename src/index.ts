"use strict";
export { quadprog as solveQP } from "./quadProg.js";
export { solveQPFast } from "./quadProgFast.js";
// Optional: eagerly terminate the solveQPFast worker pool. Never required —
// the workers are `unref`'d, so the process exits on its own once idle — but a
// long-running host can call this to release the threads deterministically.
export { shutdown } from "./parallel/pool.js";
export type { Vector, Matrix, QPResult } from "./types.js";
