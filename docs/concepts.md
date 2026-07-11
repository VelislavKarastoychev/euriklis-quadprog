# Two ideas behind the solver, in plain language

The [README's algorithm section](../README.md#the-algorithm) leans on two terms
that are standard in optimization but easy to meet for the first time here: the
**active set** and the **dual problem**. This page explains both with pictures
and intuition rather than theorems. If you can picture a bowl and a fence, you
already have everything you need.

---

## The active set

Think of the objective `f(x) = ½xᵀDx − dᵀx` as a **bowl**. Because `D` is
positive‑definite the bowl curves upward in every direction, so it has exactly
one lowest point — the unconstrained minimum `x⁰ = D⁻¹d`.

Now drop the bowl inside a **fenced yard**: each constraint `aᵢᵀx ≥ bᵢ` is one
straight **wall**, and the feasible region is the area enclosed by all the walls.

Two things can happen.

```
   bowl's bottom INSIDE the yard          bowl's bottom OUTSIDE the yard
   ─────────────────────────────          ──────────────────────────────
        ┌───────────────┐                      ┌───────────────┐
        │      · x⁰      │   ← done, the       │               │
        │   (the min)   │     answer is x⁰     │            ·  │·x⁰  ← min is
        │               │                      │           x* ╲│     outside;
        └───────────────┘                      └────────────────╲    the true
                                                       wall it leans on ↑    minimum x*
                                                                            sits ON a wall
```

- If the bottom of the bowl is already inside the yard, that point **is** the
  answer. No wall matters.
- If the bottom is outside, the ball can't reach it; it rolls until it is
  **pressed against one or more walls** and can roll no further. The minimum sits
  on those walls.

A wall the solution **leans against** — one where `aᵢᵀx* = bᵢ` exactly — is
called **active**. A wall the solution **doesn't touch** — `aᵢᵀx* > bᵢ`, with
room to spare — is **inactive**: you could delete it and the answer would not
move. The collection of walls that are touched at the optimum is the

> **active set** `𝒜` — the constraints holding with equality at the solution.

### Why it is the whole game

Here is the punchline that the algorithm is built on:

> **If you already knew `𝒜`, the problem would be easy.** Minimizing the bowl
> subject to the touched walls *as equalities* (`aᵢᵀx = bᵢ` for `i ∈ 𝒜`) is a
> plain linear‑algebra solve — no inequalities left.

The entire difficulty is *discovering which walls end up active*. An
**active‑set method** does exactly that: it keeps a guess of `𝒜` and refines it
one wall at a time — adding a wall when the current point pushes through it,
dropping a wall when it turns out not to be needed — until the guess is correct.

### A tiny example

```
minimize  ½(x² + y²)      (bowl centered at the origin)
subject to  x + y ≥ 2     (one wall)
```

The unconstrained minimum is `(0, 0)`, but `0 + 0 = 0 < 2`, so it is infeasible —
the bottom of the bowl is outside the yard. The ball rolls down to the nearest
point of the line `x + y = 2`, which is `(1, 1)`. There the single constraint is
**active**, so `𝒜 = {that constraint}`. With more walls the only question is
*which* of them the ball finally rests against.

---

## The dual problem

Every minimization problem — call it the **primal** — has a shadow twin called
the **dual**, and the two meet at the same optimal value. The bridge between them
is the idea of a **price** on each constraint.

### Prices on constraints (Lagrange multipliers)

Instead of treating a wall as an absolute barrier, attach a **price**
`λᵢ ≥ 0` to it — think of it as a fine you pay per unit by which you'd like to
push past that wall. Bundle the objective and the fines into one function, the
**Lagrangian**:

```
L(x, λ) = ½xᵀDx − dᵀx  −  Σ λᵢ (aᵢᵀx − bᵢ).
```

- The **primal** view: *pick `x`* to minimize the cost while obeying every wall.
- The **dual** view: *pick the prices `λ ≥ 0`* and ask how cheap the (now
  unconstrained) Lagrangian can be made — then choose the prices that make that
  cheapest‑case as **expensive** as possible.

So the primal minimizes over points; the dual maximizes over prices.

### Why the two are the same answer here

For any prices, the dual's value never exceeds the primal's (you can't fine your
way to a better optimum than actually exists) — that is **weak duality**. For a
convex problem like this one (`D` positive‑definite), the gap closes completely —
**strong duality** — and the optimal prices `λ*` are precisely the **Lagrange
multipliers** the solver returns.

Those prices have a concrete meaning: `λᵢ` is the **shadow price** of wall `i` —
how much the optimal objective would improve if you relaxed that wall by one
unit. This is why

- an **inactive** wall has `λᵢ = 0` (a wall you don't lean on is worth nothing to
  move), and
- an **active** wall has `λᵢ > 0` (the harder it pushes back, the more you'd gain
  by easing it).

### Primal methods vs. dual methods

This is the choice that defines Goldfarb–Idnani.

| | keeps healthy… | works toward… | needs to start… |
|---|---|---|---|
| **primal** active‑set | a *feasible* point `x ∈ Ω` | optimality (`λ ≥ 0`) | already inside the yard |
| **dual** active‑set | non‑negative *prices* `λ ≥ 0` | feasibility (`x ∈ Ω`) | anywhere — `x⁰` is free |

A **primal** method stays inside the fence and slides along walls until the
prices come out non‑negative. A **dual** method does the opposite: it keeps the
prices valid (`λ ≥ 0`) and repairs feasibility.

**Goldfarb–Idnani is a dual method.** It starts at the unconstrained minimum
`x⁰ = D⁻¹d` with no active walls and all prices zero — trivially a valid set of
prices — and then, one violated wall at a time, raises that wall's price and
lets `x` move just enough, never letting any price go negative. When no wall is
violated any more, primal and dual agree and the point is optimal.

Its big practical win: **you don't need a feasible starting point.** Finding one
would itself cost a solve in a primal method; here `x⁰` is handed to you for
free.

---

*Back to the [algorithm walk‑through](../README.md#the-algorithm).*
