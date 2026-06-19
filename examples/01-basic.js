"use strict";

/**
 * Basic inequality‑constrained QP.
 *
 *   minimise  ½‖x‖² − dᵀx     s.t.  Aᵀx ≥ b
 *
 * Run:  node examples/01-basic.js   (or:  bun examples/01-basic.js)
 */
import { solveQP } from "../index.js";

const D = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];    // ½‖x‖²
const d = [0, 5, 0];
const A = [[-4, 2, 0], [-3, 1, -2], [0, 0, 1]]; // column i = normal of constraint i
const b = [-8, 2, 0];

const r = solveQP(D, d, A, b);

console.log("solution :", r.solution.map((v) => v.toFixed(6)));        // [0.476190, 1.047619, 2.095238]
console.log("value    :", r.value.toFixed(6));                         // −2.380952
console.log("λ        :", r["Lagrangian multipliers"].map((v) => v.toFixed(6)));
console.log("active   :", r["active constraints"].slice(0, r["count of active constraints"]));
