import type { Attitude, Vec, Grid } from "./types";
export const add = (a: Vec, b: Vec): Vec => a.map((v, i) => v + b[i]) as Vec;
export const scale = (a: Vec, s: number): Vec => a.map((v) => v * s) as Vec;
export const dot = (a: Vec, b: Vec) => a.reduce((s, v, i) => s + v * b[i], 0);
const sub = (a: Vec, b: Vec) => add(a, scale(b, -1));
const cross = (a: Vec, b: Vec): Vec => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
export function rotate(v: Vec, roll: number, pitch: number): Vec {
  const r = (roll * Math.PI) / 180,
    p = (pitch * Math.PI) / 180;
  const x = Math.cos(r) * v[0] - Math.sin(r) * v[1],
    y = Math.sin(r) * v[0] + Math.cos(r) * v[1];
  return [
    x,
    Math.cos(p) * y - Math.sin(p) * v[2],
    Math.sin(p) * y + Math.cos(p) * v[2],
  ];
}
export function preview(a: Attitude, roll: number, pitch: number): Attitude {
  const axes = a.reference_axes.map((v) => rotate(v, roll, pitch));
  // Recover installation coordinates from old basis, then rotate using the new basis.
  const fov = a.fov.map((v) => {
    const body = a.axes.map((axis) => dot(v, axis));
    return [0, 1, 2].map((k) =>
      body.reduce((s, b, i) => s + b * axes[i][k], 0),
    ) as Vec;
  });
  return { ...a, axes, fov };
}
export function rayHit(origin: Vec, direction: Vec, grid: Grid): number | null {
  let closest = Infinity;
  const vertex = (i: number, j: number): Vec | null => {
    const v = grid.map((c) => c[i][j]);
    return v.every((x) => x !== null && Number.isFinite(x)) ? (v as Vec) : null;
  };
  const triangle = (a: Vec, b: Vec, c: Vec) => {
    const e1 = sub(b, a),
      e2 = sub(c, a),
      h = cross(direction, e2),
      det = dot(e1, h);
    if (Math.abs(det) < 1e-9) return;
    const s = sub(origin, a),
      u = dot(s, h) / det;
    if (u < 0 || u > 1) return;
    const q = cross(s, e1),
      v = dot(direction, q) / det;
    if (v < 0 || u + v > 1) return;
    const t = dot(e2, q) / det;
    if (t >= 0) closest = Math.min(closest, t);
  };
  for (let i = 0; i < grid[0].length - 1; i++)
    for (let j = 0; j < grid[0][i].length - 1; j++) {
      const a = vertex(i, j),
        b = vertex(i + 1, j),
        c = vertex(i, j + 1),
        d = vertex(i + 1, j + 1);
      if (a && b && c) triangle(a, b, c);
      if (b && c && d) triangle(b, d, c);
    }
  return Number.isFinite(closest) ? closest : null;
}
