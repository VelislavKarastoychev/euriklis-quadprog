"use strict";

/**
 * Equality‑constrained QP — the leading `meq` constraints are equalities.
 *
 *   minimise  ½‖x‖²     s.t.  x₁ + x₂ = 1
 *
 * The closed‑form answer is the projection of the origin onto the line, [½, ½].
 *
 * Run:  node examples/02-equality.js
 */
import { solveQP } from "../dist/index.js";

const D = [[1, 0], [0, 1]];
const d = [0, 0];
const A = [[1], [1]];   // single constraint, normal (1,1)
const b = [1];
const meq = 1;          // ← that constraint is an equality

const r = solveQP(D, d, A, b, meq);

console.log("solution :", r.solution.map((v) => v.toFixed(6))); // [0.500000, 0.500000]
console.log("value    :", r.value.toFixed(6));                  // 0.250000
console.log("λ        :", r["Lagrangian multipliers"].map((v) => v.toFixed(6)));
