export type Vec = [number, number, number];
export interface Weather {
  pressure: number;
  by: number;
  bz: number;
  dst: number;
  beta: number;
  mach: number;
  source: string;
  observed_at: string | null;
}
export interface Orbit {
  source: "moon" | "imported";
  frame: "ECI2000" | "GSM" | "GSE";
  unit: "km" | "m" | "RE";
  samples: { time: string; position: Vec; quaternion?: number[] | null }[];
}
export interface Scene {
  version: 1;
  time: string;
  reference_time: string;
  weather: Weather;
  orbit: Orbit;
  roll: number;
  pitch: number;
  target: Vec;
  target_mode: "point" | "far" | "sheet";
  hold: "earth" | "inertial" | "tracking";
  half_thickness: number;
  magnetic: boolean;
  quality: "preview" | "fine";
}
export interface Attitude {
  time: string;
  position: Vec;
  moon: Vec;
  axes: Vec[];
  reference_axes: Vec[];
  measured_axes: Vec[] | null;
  fov: Vec[];
  sun_direction: Vec;
  target_direction: Vec;
  target: Vec | null;
  roll: number;
  pitch: number;
  solution: { roll: number; pitch: number; error: number; reachable: boolean };
  pointing_error: number;
  sun_angle: number;
  sun_ok: boolean;
  ray_miss: number | null;
  target_forward: boolean;
  quaternion: number[];
  track: Vec[];
  source: string;
  measured_error: number | null;
}
export type Grid = (number | null)[][][];
export interface Environment {
  field_time?: string;
  sheet_time?: string;
  time: string;
  weather: Weather;
  magnetopause: Grid;
  bowshock: Grid;
  sheet: Grid | null;
  field_lines: Vec[][];
  errors: string[];
}
export const initial: Scene = {
  version: 1,
  time: "2026-10-20T16:00:00Z",
  reference_time: "2026-10-20T16:00:00Z",
  weather: {
    pressure: 2,
    by: 0,
    bz: -3,
    dst: -20,
    beta: 1,
    mach: 8,
    source: "manual",
    observed_at: null,
  },
  orbit: { source: "moon", frame: "ECI2000", unit: "km", samples: [] },
  roll: -30.932074,
  pitch: -1.643111,
  target: [-30, 0, 0],
  target_mode: "point",
  hold: "earth",
  half_thickness: 3,
  magnetic: true,
  quality: "preview",
};
