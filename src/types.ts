"use strict";

/**
 * Shared types for the quadprog solver.
 *
 * The public solvers take dense row-major matrices (`number[][]`) and vectors
 * (`number[]`) and return a {@link QPResult}. These aliases document intent at
 * call sites without changing the runtime representation.
 */

/** A dense vector. */
export type Vector = number[];

/** A dense, row-major matrix: `M[i][j]` is row `i`, column `j`. */
export type Matrix = number[][];

/**
 * Result of {@link solveQP} / {@link solveQPFast}.
 *
 * The Goldfarb–Idnani dual active-set method minimises
 * `½·xᵀD·x − dᵀx` subject to `A₁ᵀx = b₁` (equalities) and `A₂ᵀx ≥ b₂`.
 */
export type QPResult = {
  /** The optimiser `x` (`n×1`). */
  solution: Vector;
  /** The objective value at the optimum. */
  value: number;
  /** The unconstrained optimiser `D⁻¹d` (`n×1`). */
  unconstrained_solution: Vector;
  /**
   * Iteration counters: `[main iterations, constraints dropped after becoming
   * active]`.
   */
  iterations: [number, number];
  /**
   * Exit code — the machine-readable companion to {@link message}:
   * - `0` — success.
   * - `1` — the constraints are inconsistent (the problem is infeasible);
   *   `solution` does not satisfy them.
   * - `2` — `D` is not positive definite; `solution` is meaningless.
   *
   * Prefer branching on this over string-matching `message`.
   */
  ierr: 0 | 1 | 2;
  /** Human-readable status. `"No problems"` on success. */
  message: string;
  /** Indices of the constraints active at the optimum. */
  "active constraints": Vector;
  /** How many constraints are active at the optimum. */
  "count of active constraints": number;
  /** The Lagrange multipliers (`q×1`). */
  "Lagrangian multipliers": Vector;
};
