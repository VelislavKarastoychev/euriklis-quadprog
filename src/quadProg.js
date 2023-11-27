"use strict";
/**
 * @param {number[][]} a - matrix
 * with dimensions lda x n
 * @param {number} lda - the leading dimension of
 * the array a.
 * @param {number} n  - the order of the matrix a.
 * @description this is a translation of the
 * subroutine dopri of fortran of the quadprog
 * algorithm. The routine "computes the inverse
 * of the factor a double precision symmetric positive
 * definite matrix using the factors computed in dopfa.
 * Modification of dpodi by BaT 05/11/95."
 * @returns {number[][]} a
 * @comment - error condition: a division by zero will
 * occur if the input factor contains a zero on the diagonal
 * and the inverse is requested. It will not occur if the
 * subroutines are called correctly and if dpoco or dpofa
 * has set variable  info to be equals to 0.
 * The original routine contains the daxpy and dscal subroutines.
 * We avoid these routines by implementing them in the
 * javascript.
 */
function dpori(a, n) {
  let i, j, k, kp1, t;
  // compute inverse(r)
  for (k = 0; k < n; k++) {
    a[k][k] = 1.0 / a[k][k];
    t = -a[k][k];
    // ~ call dscal(k  1, t, a(1, k),1)
    for (i = 0; i < k; i++) {
      a[i][k] *= t;
    }
    // end of dscal...
    kp1 = k + 1;
    if (kp1 >= n) break;
    for (j = kp1; j < n; j++) {
      t = a[k][j];
      a[k][j] = 0.0;
      // ~ call daxpy(k, t, a(1, k), 1, a(1, j), 1)
      for (i = 0; i <= k; i++) a[i][j] += t * a[i][k];
    }
  }
  return a;
}
function dpofa(a) {
  let i, j, jm1, k, t, s, n, info = [];
  n = a.length;
  for (j = 0; j < n; j++) {
    info = j;
    s = 0;
    jm1 = j - 1;
    if (jm1 < 0) {
      s = a[j][j] - s;
      if (s <= 0) {
        break;
      }
      a[j][j] = Math.sqrt(s);
    } else {
      for (k = 0; k <= jm1; k++) {
        //~ t = a[k][j] - ddot(k - 1, a[1][k], 1, a[1][j], 1);
        t = a[k][j];
        for (i = 0; i < k; i++) {
          t -= a[i][j] * a[i][k];
        }
        t /= a[k][k];
        a[k][j] = t;
        s += t * t;
      }
      s = a[j][j] - s;
      if (s <= 0) {
        break;
      }
      a[j][j] = Math.sqrt(s);
    }
    info = 0;
  }
  return [a, info];
}
/**
 * @param {number [][]} a - matrix with
 * dimensions lda, n.
 * @param {number} lda - lead dimension of the matrix a.
 * @param {number} n - order of the matrix a.
 * @param {number []} b - an array with dimension n.
 * @description dopsl solves the double precision symmetric
 * positive definite system a * x = b, using the factors computed by
 * dpoco or  dopfa.
 * @return {number [][]} [a, b].
 * @commentary error condition - a division by zero will occur
 * if the input factor contains a zero on the diagonal. Technically this
 * indicates singularity but is is usually caused by improper subroutine
 * arguments. It will not occur if the subroutines are called
 * correctly and the variable info is equals to zero (0).
 * We avoid the ddot and daxpy methods of the fortran and change them
 * with the corresponded javascript implementations.
 */
function dposl(a, n, b) {
  let i, k, kb, t;
  for (k = 0; k < n; k++) {
    //~ t = ddot(k - 1, a[1][k], 1, b[1], 1);
    t = 0;
    for (i = 0; i < k; i++) {
      t += a[i][k] * b[i];
    }
    b[k] = (b[k] - t) / a[k][k];
  }
  // solve r * x = y.
  for (kb = 0; kb < n; kb++) {
    k = n - kb - 1;
    b[k] /= a[k][k];
    t = -b[k];
    // ~daxpy(k - 1, t, a(1, k), 1, b(1), 1)
    for (i = 0; i < k; i++) {
      b[i] += t * a[i][k];
    }
  }
  return [a, b];
}
/**
 * @param {0 | 1 | 2} error_index
 * @param {number []} solution the sol vector.
 * @param {number} value the crval output of the function.
 * @param {number []} unconstrained_solution the dvec output.
 * @param {[number, number]} iterations - the iter value.
 * @param {number[][]} iact - the iact values.
 * @param {number} nact
 * @param {number[]} lagr - the lagrangian multipliers of the
 * quadratic form.
 */
function solutionQP(
  error_index,
  solution,
  value,
  unconstrained_solution,
  iterations,
  iact,
  nact,
  lagr,
) {
  let message = "";
  if (error_index === 1) message = "Constraints are inconsistent. No solution.";
  if (error_index === 2) {
    message =
      "The matrix A of the quadratic form 0.5 * x^t * A * x + b * x is not positive definite.";
  } else message = "No problems";
  return {
    solution,
    value,
    unconstrained_solution,
    iterations,
    message,
    "active constraints": iact,
    "count of active constraints": nact,
    "Lagrangian multipliers": lagr,
  };
}
/**
 * @param {Array.<Array.<number>>} dmat - an n x n symmetric
 * and positive defined matrix The matrix dmat WILL BE
 * DESTROYED ON EXIT. The user has two possibilities:
 * a) Give D (ierr === 0), in this case we use routines
 * from LINPACK to decompose D.
 * b) To get the algorithm started we need R^-1, where
 * D = R^T R. So if is cheaper to calculate R^-1 ,
 * in another way (D may be a band matrix) then with the
 * general routine, the user may pass R^-1.
 * INDICATED BY IERR NOT EQUAL TO ZERO.
 *
 * @param {number []} dvec - an n x 1 vector d from above
 * Will be DESTROYED on exit. Contains on exit the solution the
 * solution to the initial, i.e., unconstrained problem.
 * @param {number} n - the dimension of the dmat.
 * @param {number [][]} amat - a n x q matrix,
 * the matrix A from above, i.e., [ A = (A1, A2)^T ].
 * ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
 * CHANGED SIGNES ON EXIT.
 * @param {number[]} bvec - a q x 1 vector, the vector
 * of constants b in the constraints [b = (b1, b2)^T].
 * ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
 * CHANGED SIGNES ON EXIT.
 * @param {number} q - an integer that describes the number
 * of the constraints.
 * @param {number} meq - an integer that describes the number
 * of equality constraints. By definition 0 <= meq <= q.
 * @return {{
 * sol : number [],
 * lagr : number [],
 * crval : number,
 * iact : number [],
 * nact : number,
 * iter : number,
 * ierr : number
 * }}
 * @commentaries on the return:
 * 1. sol is  a n x 1 vector, that is the
 * final solution of the mathematical program, i.e. the x in
 * the notation above.
 * 2. lagr is a q x 1 vector and represents the final Lagrange
 * multipliers.
 * 3. crval is a scalar, the value of the criterion at the minimum.
 * 4. iact is a q x 1 vector, the constraints which
 * are active in the final fit (int).
 * 5. nact - a scalar, the number of constraints active in
 * the final fit (int).
 * 6. iter - 2 x 1 vector, first component gives the
 * number of "main" iterations, the second one says how
 * many constraints were deleted after they became active.
 * 7. ierr - an integer, error code on exit, if
 * ierr = 0, no problems
 * ierr = 1 , the minimization problem has no solution.
 * ierr = 2, problems with decomposing of D, in this case,
 * sol contains garbage!!!
 *
 * @description
 * Description of the qpgen2 from Berwin A. Turlach.
 * Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
 * This program is free software; you can redistribute it and / or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 * This routine uses the Goldfarb/Idnani algorithm to solve the
 * following minimization problem:
 *
 * minimize:
 *            -d^T x + 1/2 * x^T D x,
 *       s.t.:
 *             A1^T x = b1
 *             A2^T x >= b2
 * where the matrix D is assumed to be positive definite.
 * Especially, without losing of generality D is assumed
 * to be symmetric.
 * Note from the translator of the javascript edition:
 * - If the matrix D is not symmetric, then set it to the
 * matrix D := 1/2 (D + D^T) and then inside it in the subroutine.
 */
function qpgen2(dmat, dvec, n, amat, bvec, q, meq, ierr = 0) {
  /**
   * @param {number} fddmat - an integer,
   * the leading dimension of the matrix dmat.
   * @param {number} ierr - an integer, code for the status
   * of the matrix D:
   * ierr === 0: we have to decompose D,
   * ierr !== 0, D is already decomposed into R^T R and we
   * were given  R^-1.
   * @param {number [][]} work - vector with length
   * at least 2 * n _ r * (r + 5) / 2 + 2 * q + 1,
   * where r = min(n, q).
   * @param {number} i - counter
   * @param {number} j - counter
   * @param {number} l - counter
   * @param {number} l1 - counter
   * @param {number} info - an integer.
   */
  // integer parameters:
  let fddmat = n,
    fdamat = n,
    i,
    info,
    it1,
    iwnbv,
    iwrm,
    iwrv,
    iwsv,
    iwuv,
    iwzv,
    j,
    l,
    l1,
    nact,
    nvl,
    r,
    state = 1;
  // integer arrays:
  let iact = [], iter = [];
  // floating point parameters:
  let crval, gc, gs, nu, sum, t1, temp, tmpa, tmpb, tt, vsmall;
  // floating point arrays:
  let lagr = [], sol = [], work = [];
  // logical variables:
  let t1inf, t2min;
  // utility functions:
  const abs = Math.abs,
    max = Math.max,
    min = Math.min,
    sign = (a, b) => {
      return b >= 0 ? a : -a;
    },
    sqrt = Math.sqrt;
  r = min(n, q);
  l = (2 * n + 0.5 * (r * (r + 5)) + 2 * q) >> 0;
  /**
   * code gleaned from Powell's ZQPCVK routine to
   * determine a small number that can be assumed
   * to be an upper bound on the relative precision
   * of the computer arithmetic.
   */
  vsmall = 1e-60;
  do {
    vsmall += vsmall;
    tmpa = 1 + 0.1 * vsmall;
    tmpb = 1 + 0.2 * vsmall;
    if (tmpa <= 1.0) continue;
    if (tmpb <= 1.0) continue;
    break;
  } while (1);
  /**
   * store the initial dvec to calculate below
   * the unconstrained minima of the critical value.
   * Note from the author: For more time efficiency
   * we cut the loop with two bitwise shiftings.
   */
  for (i = 0; i < n >> 2; i++) {
    work[i << 2] = dvec[i << 2];
    work[(i << 2) + 1] = dvec[(i << 2) + 1];
    work[(i << 2) + 2] = dvec[(i << 2) + 2];
    work[(i << 2) + 3] = dvec[(i << 2) + 3];
  }
  if (n % 4 >= 1) work[n - 1] = dvec[n - 1];
  if (n % 4 >= 2) work[n - 2] = dvec[n - 2];
  if (n % 4 >= 3) work[n - 3] = dvec[n - 3];
  // put zeros in the other places of the work
  // vector.
  for (i = 0; i < ((l - n + 1) >> 2); i++) {
    work[n + (i << 2)] = 0.0;
    work[(i << 2) + n + 1] = 0.0;
    work[(i << 2) + n + 2] = 0.0;
    work[(i << 2) + n + 3] = 0.0;
  }
  if ((l - n + 1) % 4 >= 1) work[l] = 0.0;
  if ((l - n + 1) % 4 >= 2) work[l - 1] = 0.0;
  if ((l - n + 1) % 4 >= 3) work[l - 2] = 0.0;
  for (i = 0; i < q >> 2; i++) {
    iact[i << 2] = 0.0;
    lagr[i << 2] = 0.0;
    iact[(i << 2) + 1] = 0.0;
    lagr[(i << 2) + 1] = 0.0;
    iact[(i << 2) + 2] = 0.0;
    lagr[(i << 2) + 2] = 0.0;
    iact[(i << 2) + 3] = 0.0;
    lagr[(i << 2) + 3] = 0.0;
  }
  if (q % 4 >= 1) {
    iact[q - 1] = 0.0;
    lagr[q - 1] = 0.0;
  }
  if (q % 4 >= 2) {
    iact[q - 2] = 0.0;
    lagr[q - 2] = 0.0;
  }
  if (q % 4 >= 3) {
    iact[q - 3] = 0.0;
    lagr[q - 3] = 0.0;
  }
  /**
   * get the initial solution
   */
  if (ierr === 0) {
    [dmat, info] = dpofa(dmat);
    if (info !== 0) {
      ierr = 2;
      return solutionQP(ierr, sol, crval, dvec, iter, iact, nact, lagr);
    }
    [dmat, dvec] = dposl(dmat, n, dvec);
    dmat = dpori(dmat, n);
  } else {
    /**
     * Matrix D is already factorized, so we have to
     * multiply d first with R^T and then with the
     * R^-1. Note that the inverse matrix , i.e., R^-1
     * is stored in the upper half of the array dmat.
     */
    for (j = 0; j < n; j++) {
      sol[j] = 0.0;
      for (i = 0; i <= j; i++) {
        sol[j] += dmat[i][j] * dvec[i];
      }
    }
    for (j = 0; j < n; j++) {
      dvec[j] = 0.0;
      for (i = j; i < n; i++) {
        dvec[j] += dmat[j][i] * sol[i];
      }
    }
  }
  /**
   * set lower triangular of dmat to zero, store dvec
   * in sol and calculate value of the criterion at
   * unconstrained minima.
   */
  crval = 0.0;
  for (j = 0; j < n; j++) {
    sol[j] = dvec[j];
    crval += work[j] * sol[j];
    work[j] = 0;
    for (i = j + 1; i < n; i++) {
      dmat[i][j] = 0;
    }
  }
  crval *= -0.5;
  ierr = 0;
  /**
   * calculate some constraints, i.e.,
   * from which index on the different
   * quantities are stored in the work matrix.
   */
  iwzv = n;
  iwrv = iwzv + n;
  iwuv = iwrv + r;
  iwrm = iwuv + r + 1;
  iwsv = iwrm + 0.5 * (r * (r + 1));
  iwnbv = iwsv + q;
  /**
   * calculate the norm of each
   * column of the A matrix.
   */
  for (i = 0; i < q; i++) {
    sum = 0.0;
    for (j = 0; j < n; j++) {
      sum += amat[j][i] * amat[j][i];
    }
    work[iwnbv + i] = sqrt(sum);
  }
  nact = 0;
  iter[0] = 0;
  iter[1] = 0;
  /**
   * Note from the author:
   * To implement the goto statements of the
   * code we use the technical variable state,
   * which may get the values 1, 2 or 3, that
   * correstpond to the goto 50, goto 55 and
   * goto 797 jumps.
   */
  do {
    if (state === 1) {
      /**
       * Start a new iteration:
       */
      iter[0] += 1;
      /**
       * calculate all constraints and check which
       * are still violated for the equality constraints
       * we have to check whether the normal vector has to be
       * negated (as well as bvec in that case).
       */
      l = iwsv - 1;
      for (i = 0; i < q; i++) {
        l += 1;
        sum = -bvec[i];
        for (j = 0; j < n; j++) {
          sum += amat[j][i] * sol[j];
        }
        if (abs(sum) < vsmall) sum = 0.0;
        if (i >= meq) work[l] = sum;
        else {
          work[l] = -abs(sum);
          if (sum > 0) {
            for (j = 0; j < n; j++) {
              amat[j][i] = -amat[j][i];
            }
            bvec[i] = -bvec[i];
          }
        }
      }
      /**
       * as safeguard against rounding errors set already
       * active constraints explicitly to zero.
       */
      for (i = 0; i < nact; i++) {
        work[iwsv + iact[i]] = 0.0;
      }
      /**
       * we weight each violation by the number of non - zero
       * elements in the corresponding row of A. Then we choose
       * the violated constraint which has maximal absolute value,
       * i.e., the minimum. By obvious connecting and unconnectiong
       * take always the first constraint which is violated.
       */
      nvl = 0;
      temp = 0.0;
      for (i = 0; i < q; i++) {
        if (work[iwsv + i] < temp * work[iwnbv + i]) {
          nvl = i;
          temp = work[iwsv + i] / work[iwnbv + i];
        }
      }
      if (nvl === 0) {
        for (i = 0; i < nact; i++) {
          lagr[iact[i]] = work[iwuv + i];
        }
        break;
      }
      state = 2;
    }
    if (state === 2) {
      /**
       * calculate d = J^T n^+, where n^+ is the normal vector
       * of the violated constraint. J is stored in dmat in this
       * implementation! If we drop a constraint, we have
       * to jump back here.
       */
      for (i = 0; i < n; i++) {
        sum = 0;
        for (j = 0; j < n; j++) {
          sum += dmat[j][i] * amat[j][nvl];
        }
        work[i] = sum;
      }
      // Now calculate the z = J_2 d_2
      l1 = iwzv;
      for (i = 0; i < n; i++) {
        work[l1 + i] = 0.0;
      }
      for (j = nact; j < n; j++) {
        for (i = 0; i < n; i++) work[l1 + i] += dmat[i][j] * work[j];
      }
      /**
       * and r = R^-1 d_1, check also if r has
       * positive elements (among the entries
       * corresponding to inequalities constraints).
       */
      t1inf = true;
      for (i = nact - 1; i >= 0; i--) {
        sum = work[i];
        l = (iwrm + 0.5 * (i * (i + 3))) >> 0;
        l1 = l - i;
        for (j = i + 1; j < nact; j++) {
          sum -= work[l] * work[iwrv + j];
          l += j;
        }
        sum /= work[l1];
        work[iwrv + i] = sum;
        if (iact[i] <= meq) continue;
        if (sum <= 0) continue;
        t1inf = false;
        it1 = i;
      }
      /**
       * If r has positive elements, find the partial step
       * length t1, which is the maximum step in dual space
       * without violating dual feasibility. it1 stores in
       * which component the min of u/r, occurs.
       */
      if (!t1inf) {
        t1 = work[iwuv + it1] / work[iwrv + it1];
        for (i = 0; i < nact; i++) {
          if (iact[i] <= meq) continue;
          if (work[iwrv + i] <= 0.0) continue;
          temp = work[iwuv + i] / work[iwrv + i];
          if (temp < t1) {
            t1 = temp;
            it1 = i;
          }
        }
      }
      /**
       * test if z vector is equals to zero.
       */
      sum = 0.0;
      for (i = iwzv; i < iwzv + n; i++) {
        sum += work[i] * work[i];
      }
      if (abs(sum) <= vsmall) {
        /**
         * No step in primal space such that the new
         * constraint becomes feasible. Take step in
         * dual space and drop a constraint.
         */
        if (t1inf) {
          /**
           * No step in dual space possible, either
           * problem is not solvable.
           */
          ierr = 1;
          break;
        } else {
          /**
           * We take a partial step in dual space
           * and drop constraint it1, that is, we drop
           * the it1-th active constraint. Then we
           * continue at step 2(a) (marked by label 55).
           */
          for (i = 0; i < nact; i++) {
            work[iwuv + i] -= t1 * work[iwrv + i];
          }
          work[iwuv + nact] += t1;
          state = 3;
          continue;
        }
      } else {
        /**
         * compute full step length t2, minimum step
         * in primal space such that the constraint
         * becomes feasible. Keep sum (which is z^Tn^+)
         * to update crval below.
         */
        sum = 0.0;
        for (i = 0; i < n; i++) {
          sum += work[iwzv + i] * amat[i][nvl];
        }
        tt = -work[iwsv + nvl] / sum;
        t2min = true;
        if (!t1inf) {
          if (t1 < tt) {
            tt = t1;
            t2min = false;
          }
        }
        /**
         * take step in primal and dual space.
         */
        for (i = 0; i < n; i++) {
          sol[i] += tt * work[iwzv + i];
        }
        crval += tt * sum * (tt / 2 + work[iwuv + nact]);
        for (i = 0; i < nact; i++) {
          work[iwuv + i] -= tt * work[iwrv + i];
        }
        work[iwuv + nact] += tt;
        /**
         * if it was a full step, then we check whether
         * further constraints are violated otherwise we
         * can drop the current constraint and iterate
         * once more.
         */
        if (t2min) {
          /**
           * we took a full step. Thus add constraint nvl
           * to the list of active constraints add update
           * J and R.
           */
          ++nact;
          iact[nact - 1] = nvl;
          /**
           * to update Rwe have to put the first nact - 1
           * components of the d vector into column (nact)
           * of R.
           */
          l = iwrm + 0.5 * ((nact - 1) * nact);
          for (i = 0; i < nact - 1; i++) {
            work[l] = work[i];
            ++l;
          }
          /**
           * if now nact = n, then we just have to add
           * the last element to the new row of R. Otherwise
           * we use Givens transformations to turn the vector
           * d(nact:n) into a multiple of the first unit vector.
           * That multiple goes into the last element of the
           * new row of R and J is accordingly updated by the Givens
           * transformations.
           */
          if (nact === n) {
            work[l] = work[n - 1];
          } else {
            for (i = n - 1; i >= nact; i--) {
              /**
               * we have to find the Givens rotation
               * which will reduce the element (l1) of
               * d to zero. If it is already zero we don't
               * have to do anything, except of decreasing l1.
               */
              if (work[i] === 0.0) continue;
              gc = max(abs(work[i - 1]), abs(work[i]));
              gs = min(abs(work[i - 1]), abs(work[i]));
              temp = sign(gc * sqrt(1 + (gs * gs) / (gc * gc)), work[i - 1]);
              gc = work[i - 1] / temp;
              gs = work[i] / temp;
              /**
               * The Givens rotation is done with the matrix
               * (gc gs, gs - gc). If gc is one, then element (i)
               * of d is zero compaed with element (l1 - 1). Hence
               * we don't have to do anything. If gc is zro, then
               * we just have to switch column (i) and column (i - 1)
               * of J. Since we only switch columns in J,  we have
               * to be careful how we update d depending on the sign of
               * gs. Otherwise we have to apply the Givens ratations to
               * these columns. The i - 1 element of d has to be updated to temp.
               */
              if (gc === 1.0) continue;
              if (gc === 0.0) {
                work[i - 1] = gs * temp;
                for (j = 0; j < n; j++) {
                  temp = dmat[j][i - 1];
                  dmat[j][i - 1] = dmat[j][i];
                  dmat[j][i] = temp;
                }
              } else {
                work[i - 1] = temp;
                nu = gs / (1 + gc);
                for (j = 0; j < n; j++) {
                  temp = gc * dmat[j][i - 1] + gs * dmat[j][i];
                  dmat[j][i] = nu * (dmat[j][i - 1] + temp) - dmat[j][i];
                  dmat[j][i - 1] = temp;
                }
              }
            }
            /**
             * l is still pointin to element (nact, nact)
             * of the matrix R. So store d(nact) in R(nact, nact).
             */
            work[l] = work[nact - 1]; // or work[l] = work[nact - 1];
          }
        } else {
          /**
           * we took a partial step in dual space.
           * Thus drop constraint it1, that is we
           * drop the it1-th active constraint. Then
           * we continue at step 2(a) (marked by label 55)
           * but since the fit changed, we have to recalculate
           * now "how much" the fit violates the
           * chosen constraint naw.
           */
          sum = -bvec[nvl];
          for (j = 0; j < n; j++) {
            sum += sol[j] * amat[j][nvl];
          }
          if (nvl >= meq - 1) {
            work[iwsv + nvl] = sum;
          } else {
            work[iwsv + nvl] = -abs(sum);
            if (sum > 0) {
              for (j = 0; j < n; j++) {
                amat[j][nvl] = -amat[j][nvl];
              }
              bvec[nvl] = -bvec[nvl];
            }
          }
          state = 3;
          continue;
        }
      }
      state = 1;
    }
    /**
     * Drop constraint it1.
     */
    if (state === 3) {
      /**
       * if it1 = nact it is only neccessary to
       * update the vector u and nact.
       */
      while (it1 < nact) {
        /**
         * After updation one row of R (column of J)
         * we will also come back here.
         * We have to find the Givens rotations which
         * will reduce the element (it1 + 1, it1 + 1)
         * of R to zero. If it si already zero, we  don't
         * have to do anything except of updation u, iact,
         * and shifting column (it1 + 1) of R to column (it1)
         * l will point to element (1, it1 + 1) of R l1 will
         * point to element (it1 + 1, it1 + 1) of R.
         */
        l = iwrm + 0.5 * ((it1 + 1) * (it1 + 2));
        l1 = l + it1 + 1;
        if (work[l1] !== 0.0) {
          gc = max(abs(work[l1 - 1]), abs(work[l1]));
          gs = min(abs(work[l1 - 1]), abs(work[l1]));
          temp = sign(gc * sqrt(1 + gs * gs / (gc * gc)), work[l1 - 1]);
          gc = work[l1 - 1] / temp;
          gs = work[l1] / temp;
          /**
           * The Givens rotation is done with the matrix
           * (gc gs, gs -gc) . If gc is one, then element (it1 + 1, it1 + 1)
           * of R is zero compared with element (it1, it1 + 1).
           * Hence we don't have to do anything. If gc is zero, then we just
           * have to switch row (it1) and row (it1 + 1) of R and column (it1)
           * and column (it1 + 1) of J. Since we switch reows in R and
           * columns in J, we can ignore the sign of gs. Otherwise we have
           * to apply the Givens rotation to these rows/columns.
           */
          if (gc !== 1.0) {
            if (gc === 0) {
              for (i = it1 + 1; i < nact; i++) {
                temp = work[l1 - 1];
                work[l1 - 1] = work[l1];
                work[l1] = temp;
                l1 += i;
              }
              for (i = 0; i < n; i++) {
                temp = dmat[i][it1];
                dmat[i][it1] = dmat[i][it1 + 1];
                dmat[i][it1 + 1] = temp;
              }
            } else {
              nu = gs / (1 + gc);
              for (i = it1 + 1; i < nact; i++) {
                temp = gc * work[l1 - 1] + gs * work[l1];
                work[l1] = nu * (work[l1 - 1] + temp) - work[l1];
                work[l1 - 1] = temp;
                l1 += i;
              }
              for (i = 0; i < n; i++) {
                temp = gc * dmat[i][it1] + gs * dmat[i][it1 + 1];
                dmat[i][it1 + 1] = nu * (dmat[i][it1] + temp) -
                  dmat[i][it1 + 1];
                dmat[i][it1] = temp;
              }
            }
          }
        }
        /**
         * Shift column (it1 + 1) of R to column (it1)
         * (that is, the first it1 elements). THe position
         * of element (1, it1 + 1) of R was calculated above
         * and stored in l.
         */
        l1 = l - it1 - 1;
        for (i = 0; i <= it1; i++) {
          work[l1] = work[l];
          ++l;
          ++l1;
        }
        /**
         * update vector u and iact as necessary
         * Contunue with updation the matrices J and R.
         */
        work[iwuv + it1] = work[iwuv + it1 + 1];
        iact[it1] = iact[it1 + 1];
        ++it1;
      }
      work[iwuv + nact - 1] = work[iwuv + nact];
      work[iwuv + nact] = 0;
      iact[nact - 1] = 0;
      --nact;
      ++iter[1];
      state = 2;
    }
  } while (1);
  return solutionQP(ierr, sol, crval, dvec, iter, iact, nact, lagr);
}
/**
 *
 * @param {number [][]} D - a squared matrix with dimensions n x n.
 * @param {number[]} d - a vector with dimensions n x 1.
 * @param {number[][]} A - a rectangular matrix with dimensions q x n.
 * @param {number[]} b - a vector with dimensions q x 1.
 * @param {number?} meq - number of the equality constraints.
 * By default this parameter is set to zero.
 * @param {0 | 1 | 2?} ierr - an integer, when is zero comutes the
 * inverse of the D matrix and make factorizations, otherwise assumes that
 * the inverse matrix is alredy computed in D. By default this parameter is set to 0.
 * @returns {{
 * solution : Array.<number>,
 value : number,
 unconstrained_solution : Array.<number>,
 iterations : [number, number],
 message : string,
 'active constraints': Array.<number>,
 'count of active constraints': number,
 'Lagrangian multipliers': Array.<number>
 * }}
 *@description
 *This program is free software; you can redistribute it and / or modify it
 * under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 * This routine uses the Goldfarb/Idnani algorithm to solve the
 * following minimization problem:
 *
 * minimize:
 *            -d^T x + 1/2 * x^T D x,
 *       s.t.:
 *             A1^T x = b1
 *             A2^T x >= b2
 * where the matrix D is assumed to be positive definite.
 * Especially, without losing of generality D is assumed
 * to be symmetric.
 *
 * If the matrix D is not symmetric, then set it to the
 * matrix D := 1/2 (D + D^T) and then inside it in the subroutine.
 */
export function quadprog(D, d, A, b, meq = 0, ierr = 0) {
  // fast deep copy of the elements
  let _A = [],
    _b = [],
    _d = [],
    _D = [],
    _ierr = ierr,
    i,
    j,
    n = D.length,
    q = A.length;
  for (i = 0; i < n >> 1; i++) {
    _d[i << 1] = d[i << 1];
    _d[(i << 1) + 1] = d[(i << 1) + 1];
    _D[i << 1] = [];
    _D[(i << 1) + 1] = [];
    for (j = 0; j < n >> 1; j++) {
      _D[i << 1][j << 1] = D[i << 1][j << 1];
      _D[i << 1][(j << 1) + 1] = D[i << 1][(j << 1) + 1];
      _D[(i << 1) + 1][j << 1] = D[(i << 1) + 1][j << 1];
      _D[(i << 1) + 1][(j << 1) + 1] = D[(i << 1) + 1][(j << 1) + 1];
    }
    if (n & 1) {
      _D[i << 1][n - 1] = D[i << 1][n - 1];
      _D[(i << 1) + 1][n - 1] = D[(i << 1) + 1][n - 1];
    }
  }
  if (n & 1) {
    _d[n - 1] = d[n - 1];
    _D[n - 1] = [];
    for (i = 0; i < n >> 1; i++) {
      _D[n - 1][i << 1] = D[n - 1][i << 1];
      _D[n - 1][(i << 1) + 1] = D[n - 1][(i << 1) + 1];
    }
    _D[n - 1][n - 1] = D[n - 1][n - 1];
  }
  for (i = 0; i < q >> 1; i++) {
    _b[i << 1] = b[i << 1];
    _b[(i << 1) + 1] = b[(i << 1) + 1];
    _A[i << 1] = [];
    _A[(i << 1) + 1] = [];
    for (j = 0; j < n >> 1; j++) {
      _A[i << 1][j << 1] = A[i << 1][j << 1];
      _A[i << 1][(j << 1) + 1] = A[i << 1][(j << 1) + 1];
      _A[(i << 1) + 1][j << 1] = A[(i << 1) + 1][j << 1];
      _A[(i << 1) + 1][(j << 1) + 1] = A[(i << 1) + 1][(j << 1) + 1];
    }
    if (n & 1) {
      _A[i << 1][n - 1] = A[i << 1][n - 1];
      _A[(i << 1) + 1][n - 1] = A[(i << 1) + 1][n - 1];
    }
  }
  if (q & 1) {
    _b[q - 1] = b[q - 1];
    _A[q - 1] = [];
    for (i = 0; i < n >> 1; i++) {
      _A[q - 1][i << 1] = A[q - 1][i << 1];
      _A[q - 1][(i << 1) + 1] = A[q - 1][(i << 1) + 1];
    }
    if (n & 1) _A[q - 1][n - 1] = A[q - 1][n - 1];
  }
  return qpgen2(_D, _d, n, _A, _b, q, meq, _ierr);
}
