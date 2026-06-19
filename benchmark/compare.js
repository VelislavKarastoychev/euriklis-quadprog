"use strict";

/**
 * Runs benchmark.js under BOTH Node and Bun and prints one combined table:
 * Santini vs solveQP vs solveQPFast, side by side for the two runtimes, with the
 * CPU it ran on. Whichever runtime is missing is simply skipped.
 *
 *   node benchmark/compare.js     (or:  bun benchmark/compare.js)
 */
import { execFileSync } from "child_process";
import { fileURLToPath } from "url";
import path from "path";
import os from "os";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const bench = path.join(__dirname, "benchmark.js");

const run = (cmd) => {
  try {
    const out = execFileSync(cmd, [bench, "--json"], { encoding: "utf8", stdio: ["ignore", "pipe", "ignore"], timeout: 600000 });
    return JSON.parse(out.trim().split("\n").pop());
  } catch { return null; }
};

console.log("Running Node…"); const node = run("node");
console.log("Running Bun…");  const bun = run("bun");

if (!node && !bun) { console.error("Neither node nor bun could run the benchmark."); process.exit(1); }
const ref = node ?? bun;

console.log(`\nCPU: ${ref.cpu}  (${ref.threads} threads)`);
console.log(`Node: ${node ? node.runtime : "n/a"}   Bun: ${bun ? bun.runtime : "n/a"}`);
console.log("\nMedian solve time (ms) — lower is better\n");

const fmt = (v) => (v == null ? "   n/a" : v.toFixed(1).padStart(8));
const byN = (src) => Object.fromEntries((src?.results ?? []).map((r) => [r.n, r]));
const N = byN(node), B = byN(bun);
const sizes = (ref.results ?? []).map((r) => r.n);

console.log("        │            Node (ms)         │             Bun (ms)");
console.log("    n   │  Santini   solveQP      Fast │   Santini   solveQP      Fast");
console.log("────────┼──────────────────────────────┼──────────────────────────────");
for (const n of sizes) {
  const a = N[n] ?? {}, c = B[n] ?? {};
  console.log(
    `${String(n).padStart(6)}  │${fmt(a.santini)}${fmt(a.mine)}${fmt(a.fast)} │${fmt(c.santini)}${fmt(c.mine)}${fmt(c.fast)}`,
  );
}
console.log("\nFast vs solveQP speed‑up:");
console.log("    n   │   Node   │    Bun");
console.log("────────┼──────────┼──────────");
for (const n of sizes) {
  const a = N[n] ?? {}, c = B[n] ?? {};
  const sp = (r) => (r.mine && r.fast ? (r.mine / r.fast).toFixed(2) + "×" : "n/a");
  console.log(`${String(n).padStart(6)}  │ ${sp(a).padStart(7)}  │ ${sp(c).padStart(7)}`);
}
void os;
