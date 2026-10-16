import { describe, it, expect } from "vitest";
import { rotate, rayHit } from "./geometry";
import { Latest, localTime, utc } from "./api";
describe("scene math and request ordering", () => {
  it("applies fixed Z then X", () => {
    const v = rotate([1, 0, 0], 90, 90);
    expect(v[0]).toBeCloseTo(0);
    expect(v[1]).toBeCloseTo(0);
    expect(v[2]).toBeCloseTo(1);
  });
  it("rejects old successful or failed requests", () => {
    const gate = new Latest();
    const old = gate.next(),
      latest = gate.next();
    expect(gate.current(old)).toBe(false);
    expect(gate.current(latest)).toBe(true);
  });
  it("uses Beijing regardless of browser timezone", () => {
    expect(localTime("2026-10-20T16:00:00Z")).toBe("2026-10-21T00:00");
    expect(utc("2026-10-21T00:00")).toBe("2026-10-20T16:00:00.000Z");
  });
  it("finds forward mesh crossings and rejects backwards rays", () => {
    const grid = [
      [
        [1, 1],
        [1, 1],
      ],
      [
        [-1, 1],
        [-1, 1],
      ],
      [
        [-1, -1],
        [1, 1],
      ],
    ];
    expect(rayHit([0, 0, 0], [1, 0, 0], grid)).toBeCloseTo(1);
    expect(rayHit([0, 0, 0], [-1, 0, 0], grid)).toBe(null);
  });
});
