# Fortran Code Analysis Summary

## 1. Code Purpose

The Fortran program `Simulation.for` processes ENA (Energetic Neutral Atom) detector data from a space-based instrument. It:

- Reads raw TDC (Time-to-Digital Converter) time data from input files
- Converts TDC time measurements to spatial coordinates (X, Y) using calibration parameters
- Calculates particle velocity and energy
- Classifies particles as Hydrogen (H) or Oxygen (O) based on time-of-flight and amplitude
- Computes azimuth and elevation angles
- Generates binned distribution data (ANGLE array)
- Outputs processed data to multiple files (VALID, INITIAL, VOID, CONTOUR)

## 2. Code Structure

### Files:
| File | Purpose |
|------|---------|
| `Simulation.for` | Main program with data processing logic |
| `Dimens.for` | Parameter definitions (array sizes) |
| `ITERANCE.NML` | Calibration parameters (boundary time values) |

### Key Parameters from ITERANCE.NML:
```
AXL=-54.917  (Probe A X left boundary, ns)
AXR=50.447   (Probe A X right boundary, ns)
AYD=-16.686  (Probe A Y down boundary, ns)
AYU=41.627   (Probe A Y up boundary, ns)
BXL=-53.405  (Probe B X left boundary, ns)
BXR=51.184   (Probe B X right boundary, ns)
```

## 3. Signal Selection and Classification Logic

The code classifies measurements into **three categories** based on data quality:

### Category 1: VOID (Noise/Invalid) → `VOID_T20250224-p1.DAT`
**Counter:** `NVO`

**Criteria for rejection (noise/invalid measurements):**

| Check | Condition | Description |
|-------|-----------|-------------|
| Zero Signal | `(X1e.EQ.0.0.AND.X2E.EQ.0.0)` OR `(Y1e.EQ.0.0.AND.Y2E.EQ.0.0)` | Both X or both Y TDC channels read zero |
| Out of Bounds (Start) | `ABS(Xs).GT.45` OR `ABS(Ys).GT.15` | Start position outside valid detector area |
| Missing Signal | `KX.EQ.1.AND.KY.EQ.1` | Both X and Y time differences out of boundary range |

**Logic flow:** These events go to `WRITE(6,111)` (VOID file) via `GOTO 11`

---

### Category 2: INITIAL (Incomplete/Edge Cases) → `INITIAL_T20250224-p1.DAT`
**Counter:** `NIN`

**Criteria for partial data (edge cases):**

| Check | Condition | Description |
|-------|-----------|-------------|
| End Position OOB | `ABS(Xe).GT.75` OR `ABS(Ye).GT.27.5` | End position outside detector boundaries |
| Invalid TOF | `DTof.LT.6.OR.DTof.GT.187` | Time-of-flight outside physical range (6-187 ns) |
| Excessive X Travel | `ABS(DX).GT.45` | X displacement exceeds detector half-width |
| Excessive Y Travel | `ABS(DY).GT.15` | Y displacement exceeds detector half-height |
| Missing End Signal | `KKX.EQ.1.AND.KKY.EQ.1` | Both X and Y end signals missing |

**Logic flow:** These events go to `WRITE(5,112)` (INITIAL file) via `GOTO 22`

---

### Category 3: VALID (Physical Signals) → `VALID_T20250224-p1.DAT`
**Counter:** `NVI`

**Criteria for valid physical signals:**

| Check | Condition | Description |
|-------|-----------|-------------|
| TOF Consistency | `Tof0.GT.2.0.AND.Tof1.LT.Tof0` | Basic time ordering check |
| Rollover Correction | Lines 154-157 handle TDC rollover | Adds 663.47 ns when end < start |
| Boundary Check | Time differences within `(XL, XR)` and `(YD, YU)` | Within calibrated detector boundaries |
| Position Reconstruction | Valid Xs, Ys calculated from time differences | Using linear transformation |
| End Position Valid | Valid Xe, Ye within detector | Using similar logic |
| TOF Range | `6.LE.DTof.LE.187` | Physically reasonable time-of-flight |
| Displacement Check | `ABS(DX).LE.45` and `ABS(DY).LE.15` | Reasonable particle trajectory |

**Logic flow:** These events go to `WRITE(4,113)` (VALID file)

---

### Particle Identification (Physical Signal Sub-classification)

For valid events, the code classifies particles based on time-of-flight:

| TOF Range | Amplitude | Classification | Energy Formula |
|-----------|-----------|----------------|----------------|
| `DTof < 40.0 ns` | Any | Hydrogen (NPar=1) | `Eng = HE * (DL/DTof)^2` |
| `DTof > 66.3 ns` | Any | Oxygen (NPar=0) | `Eng = OE * (DL/DTof)^2` |
| `40.0 ≤ DTof ≤ 66.3 ns` | `AMP > 1580` | Oxygen (NPar=0) | `Eng = OE * (DL/DTof)^2` |
| `40.0 ≤ DTof ≤ 66.3 ns` | `AMP ≤ 1580` | Hydrogen (NPar=1) | `Eng = HE * (DL/DTof)^2` |

Where:
- `HE = 5.2247` (Hydrogen energy conversion factor)
- `OE = 82.927` (Oxygen energy conversion factor)
- `AMP > 1580` indicates higher amplitude → Oxygen hypothesis

---

## 4. Coordinate Calculation Method

### Physics Model:
The code implements the linear coordinate transformation from TDC time differences to spatial coordinates:

```
CofX = X_MCP / (XR - XL)      [mm/ns]
CofY = Y_MCP / (YU - YD)      [mm/ns]
Z0X = (XR + XL) / 2           [ns]
Z0Y = (YU + YD) / 2           [ns]

X = (Δt_X - Z0X) * CofX       [mm]
Y = (Δt_Y - Z0Y) * CofY       [mm]
```

Where:
- `Δt_X = X2S - X1S` (time difference between right and left TDC channels)
- `X_MCP = 145 mm`, `Y_MCP = 55 mm` (detector dimensions)

## 4. Verification Against TDC_Coordinate_System.md

### Correct Implementations:
| Aspect | Document | Code | Status |
|--------|----------|------|--------|
| Slope formula | k = 150/(XR-XL) | CofX = 145/(XR-XL) | ✓ Correct |
| Linear form | x = k·Δt + b | X = CofX·(Δt-Z0X) | ✓ Equivalent |
| Probe A X channels | 5 (TX1), 7 (TX2) | Used in data | ✓ Consistent |
| Probe A Y channels | 1 (TY1), 3 (TY2) | Used in data | ✓ Consistent |

### Issues Found:

#### Issue 1: **CRITICAL BUG - Wrong Coefficient Usage**
- **Location:** Lines 207, 258
- **Problem:** Code uses `Cofy` instead of `CofX` for X-coordinate calculation
```fortran
Xs=(2*(Tofs-X1s)+CX-Z0X)*Cofy   ! BUG: should be CofX
Xe=(2*(Tofe-X1e)+CX-Z0X)*Cofy   ! BUG: should be CofX
```
- **Impact:** X coordinates will be incorrectly scaled when this branch is executed

#### Issue 2: **Sign Convention Mismatch**
- **Document:** `x = k·(TX1 - TX2) + b`
- **Code:** `Xs = (X2s-X1s-Z0X)*CofX` (assuming X2s=TX2, X1s=TX1)
- **Problem:** Code uses (TX2 - TX1) instead of (TX1 - TX2)
- **Impact:** X coordinate sign may be inverted depending on channel assignment

#### Issue 3: **Hardcoded Dimensions**
- **Problem:** Board dimensions (145mm, 55mm) are hardcoded instead of using parameters
- **Location:** Lines 94-95, 102-103
- **Impact:** Maintenance difficulty if detector dimensions change

#### Issue 4: **Inconsistent Variable Naming**
- **Document:** Uses A_X_LEFT_NS, A_X_RIGHT_NS, etc.
- **Code:** Uses AXL, AXR, AYD, AYU
- **Note:** Comments suggest AYU/AYD but code logic uses YU/YD inconsistently

#### Issue 5: **Time Conversion Bug**
- **Location:** Lines 86-90
- **Problem:** Hardcoded Beijing to UTC conversion (hour-8)
- **Impact:** Timezone handling is fragile and location-dependent

#### Issue 6: **Boundary Check Logic**
- **Location:** Lines 172-173, 178-179
- **Problem:** Uses strict inequalities `.GT.XL.AND..LT.XR`
- **Impact:** Events exactly at boundaries may be excluded

#### Issue 7: **Magic Numbers**
- **Problem:** Many undocumented constants throughout code
- Examples:
  - Line 154-157: `663.47` (rollover correction?)
  - Line 271: `2500` (50^2 for DL calculation)
  - Line 308: `1580` (amplitude threshold)
  - Line 40, 66.3: TOF boundaries for particle ID

#### Issue 8: **Particle ID Logic**
- **Location:** Lines 299-318
- **Problem:** Discrete TOF thresholds (40ns, 66.3ns) for H/O separation
- **Impact:** Sharp boundaries may misclassify particles near thresholds

#### Issue 9: **Missing Error Handling**
- Input file opening without existence checks
- No validation of calibration parameters
- Array bounds not checked before access

#### Issue 10: **Hardcoded File Names**
- Input/output filenames are hardcoded with specific dates
- Impact: Requires recompilation to process different data files

## 5. Recommendations

### High Priority:
1. **Fix the CofY/CofX bug** at lines 207 and 258
2. **Verify channel mapping** (X1/X2 vs TX1/TX2) to ensure correct sign convention
3. **Add parameter definitions** for board dimensions instead of hardcoding

### Medium Priority:
4. Replace magic numbers with named constants
5. Add input validation and error handling
6. Make file names configurable via namelist

### Low Priority:
7. Add comments explaining the physics (rollover correction, time calculations)
8. Consider using configuration file for particle ID thresholds

## 6. Calibration Parameter Comparison

| Parameter | Document Value | Code Value (ITERANCE.NML) | Match |
|-----------|----------------|---------------------------|-------|
| A_X_LEFT_NS | -50.44742 | -54.917 | ❌ Different |
| A_X_RIGHT_NS | 54.91726 | 50.447 | ❌ Different |
| A_Y_LEFT_NS | -34.71570 | -16.686 | ❌ Different |
| A_Y_RIGHT_NS | 23.59747 | 41.627 | ❌ Different |
| B_X_LEFT_NS | -51.18360 | -53.405 | ❌ Different |
| B_X_RIGHT_NS | 53.40479 | 51.184 | ❌ Different |

**Note:** The calibration values in the documentation differ from those in ITERANCE.NML. The code uses the values from the namelist file, which may indicate updated calibration data not reflected in the documentation.

## 7. Complete Program Flow

### 7.1 Text Flow Description

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              INITIALIZATION                                  │
│  • Read calibration parameters from ITERANCE.NML                            │
│  • Initialize counters (NVO=0, NVI=0, NIN=0)                                │
│  • Open 5 files:                                                            │
│    - Input:  ENA_202511121033_202.TXT (raw TDC data)                        │
│    - Output: VALID_T20250224-p1.DAT   (valid signals)                       │
│    - Output: INITIAL_T20250224-p1.DAT (incomplete data)                     │
│    - Output: VOID_T20250224-p1.DAT    (noise/invalid)                       │
│    - Output: CONTOUR_T20250224-p1.DAT (binned distribution)                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           MAIN LOOP (K=1 to KKK)                             │
│  For each event record in input file:                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          STEP 1: READ & PARSE DATA                           │
│  • Read: date_str, time_str, TOF01, TOF11, X1S1, X1E1, X2S1, X2E1,         │
│          Y1S1, Y1E1, Y2S1, Y2E1, AMP1                                       │
│  • Parse date/time strings                                                  │
│  • Apply timezone correction: hour = hour - 8 (Beijing → UTC)               │
│  • Select calibration based on NPM (0=Probe A, 1=Probe B)                   │
│    → Calculate CofX, CofY, Z0X, Z0Y, CX, CY, ETof                           │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       STEP 2: PREPROCESS TDC TIMES                           │
│  • IF (Tof0 > 2.0 AND Tof1 < Tof0):  (找到合理的tof采样值放到1上)                                       │
│      → Swap Tof0 and Tof1 (ensure proper ordering)                          │
│  • Rollover correction (TDC counter overflow):                              │
│      → IF (X1e < X1s): X1e = X1s + 663.47                                   │
│      → IF (X2e < X2s): X2e = X2s + 663.47                                   │
│      → IF (Y1e < Y1s): Y1e = Y1s + 663.47                                   │
│      → IF (Y2e < Y2s): Y2e = Y1s + 663.47   [BUG: Should be Y2s]            │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────┐
                    │     ZERO SIGNAL CHECK?          │
                    │  (X1e=0 AND X2E=0) OR          │
                    │  (Y1e=0 AND Y2E=0) ?            │
                    └─────────────────────────────────┘
                         │                 │
                    YES  │                 │  NO
                         ▼                 ▼
               ┌──────────────┐    ┌──────────────┐
               │ WRITE TO VOID│    │  CONTINUE TO │
               │    NVO++     │    │ START POS    │
               │  GOTO 10     │    │   CHECK      │
               └──────────────┘    └──────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     STEP 3: CALCULATE START POSITION (Xs, Ys)                │
│                                                                              │
│  X Calculation:                                                             │
│    IF (XL < (X2s-X1s) < XR):                                                │
│      → Tofx0 = (X2s + X1s - CX) / 2                                         │
│      → Xs = (X2s - X1s - Z0X) * CofX                                        │
│      → KX = 0 (valid)                                                       │
│    ELSE:                                                                    │
│      → KX = 1 (missing/edge case)                                           │
│                                                                              │
│  Y Calculation:                                                             │
│    IF (YD < (Y2s-Y1s) < YU):                                                │
│      → Tofy0 = (Y2s + Y1s - CY) / 2                                         │
│      → Ys = (Y2s - Y1s - Z0Y) * CofY                                        │
│      → KY = 0 (valid)                                                       │
│    ELSE:                                                                    │
│      → KY = 1 (missing/edge case)                                           │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────┐
                    │      BOTH X AND Y VALID?        │
                    │       (KX=0 AND KY=0)?          │
                    └─────────────────────────────────┘
                         │                 │
                    YES  │                 │  NO
                         ▼                 ▼
               ┌──────────────┐    ┌─────────────────────────┐
               │ Tofs = MAX   │    │ KX=1 AND KY=1?          │
               │  (Tofx0,Tofy0)│   └─────────────────────────┘
               │ GOTO 20      │         │              │
               │ (End calc)   │    YES   │              │  NO
               └──────────────┘         ▼              ▼
                              ┌──────────────┐  ┌──────────────────┐
                              │ WRITE TO VOID│  │ Missing X or Y:  │
                              │    NVO++     │  │ Reconstruct from │
                              │  GOTO 10     │  │ available signal │
                              └──────────────┘  └──────────────────┘
                                                         │
                                                         ▼
                              ┌──────────────────────────────────────────┐
                              │  SINGLE SIGNAL RECONSTRUCTION:           │
                              │  IF (KY=1): Use Tofs = Tofx0             │
                              │    → Reconstruct Ys from Tofs using      │
                              │      nearest Y channel                   │
                              │    → Check if other Y channel rollover   │
                              │                                          │
                              │  IF (KX=1): Use Tofs = Tofy0             │
                              │    → Reconstruct Xs from Tofs using      │
                              │      nearest X channel                   │
                              │    → Check if other X channel rollover   │
                              │    → [BUG: Uses Cofy instead of CofX]    │
                              └──────────────────────────────────────────┘
                                                         │
                                                         ▼
                    ┌─────────────────────────────────┐
                    │      VALID START POSITION?      │
                    │    (|Xs| ≤ 45 AND |Ys| ≤ 15)?   │
                    └─────────────────────────────────┘
                         │                 │
                    YES  │                 │  NO
                         ▼                 ▼
               ┌──────────────┐    ┌──────────────┐
               │  CONTINUE    │    │ WRITE TO VOID│
               │  TO STEP 4   │    │    NVO++     │
               └──────────────┘    │  GOTO 10     │
                                   └──────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: CALCULATE END POSITION (Xe, Ye)                                     │
│  [Same logic as Step 3 but for end times X1e, X2e, Y1e, Y2e]                │
│  → Produces KKX, KKY, Tofx1, Tofy1, Xe, Ye                                  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
                    ┌─────────────────────────────────┐
                    │    BOTH END SIGNALS MISSING?    │
                    │      (KKX=1 AND KKY=1)?         │
                    └─────────────────────────────────┘
                         │                 │
                    YES  │                 │  NO
                         ▼                 ▼
               ┌─────────────────┐  ┌──────────────────┐
               │ WRITE TO INITIAL│  │ Tofe = MIN       │
               │     NIN++       │  │ (Tofy1,Tofx1)    │
               │   GOTO 10       │  │ IF Tofe ≤ Tofs:  │
               └─────────────────┘  │   Tofe=MAX       │
                                    │ (Tof1,Tofx1)     │
                                    └──────────────────┘
                                                         │
                                                         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 5: VALIDATE END POSITION                                               │
│  • Reconstruct Xe/Ye if missing (similar to step 3)                         │
│  • Check: |Xe| ≤ 75 AND |Ye| ≤ 27.5                                         │
│    → IF FAIL: WRITE TO INITIAL (NIN++, GOTO 10)                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STEP 6: FINAL VALIDATION CHECKS                           │
│                                                                              │
│  1. Calculate DTof = Tofe - Tofs + ETof                                     │
│     → Check: 6 ≤ DTof ≤ 187 ns                                              │
│     → IF FAIL: WRITE TO INITIAL (NIN++, GOTO 10)                            │
│                                                                              │
│  2. Calculate DX = Xe - Xs, DY = Ye - Ys                                    │
│     → Check: |DX| ≤ 45 AND |DY| ≤ 15                                       │
│     → IF FAIL: WRITE TO INITIAL (NIN++, GOTO 10)                            │
│                                                                              │
│  3. Calculate DL = SQRT(DX² + DY² + 2500)   [Path length]                   │
│     → V = DL / DTof * 1000   [Velocity in mm/ns]                            │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     STEP 7: PARTICLE IDENTIFICATION                          │
│                                                                              │
│  IF (DTof < 40.0 ns):                                                       │
│    → NPar = 1 (Hydrogen)                                                    │
│    → Eng = HE * V²                                                          │
│                                                                              │
│  ELSE IF (DTof > 66.3 ns):                                                  │
│    → NPar = 0 (Oxygen)                                                      │
│    → Eng = OE * V²                                                          │
│                                                                              │
│  ELSE (40.0 ≤ DTof ≤ 66.3 ns):                                              │
│    → IF (AMP > 1580): NPar = 0 (Oxygen, higher amplitude)                   │
│    → IF (AMP ≤ 1580): NPar = 1 (Hydrogen, lower amplitude)                  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    STEP 8: ANGLE CALCULATION & BINNING                       │
│                                                                              │
│  Calculate Azimuth and Elevation:                                           │
│    → Azi = ATAN(DX / 50.0) * 180/PI                                         │
│    → Elv = ATAN(DY / 50.0) * 180/PI ± 5.0  (+5 for Probe A, -5 for Probe B) │
│                                                                              │
│  Bin into ANGLE(NI, NJ) array:                                              │
│    → MI = INT(Azi/3) + 16  (Azimuth bin, 1-31)                              │
│    → MJ = INT(Elv/3) + 9   (Elevation bin, 1-17)                            │
│    → ANGLE(MI, MJ) = ANGLE(MI, MJ) + 1                                      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      STEP 9: WRITE VALID EVENT                               │
│  • NVI++                                                                    │
│  • WRITE(4,113): year, month, day, hour, minute, second,                    │
│                  NPar, Tofs, Tofe, DTof, Xs, Ys, Xe, Ye,                    │
│                  Azi, Elv, Eng                                              │
│  • GOTO 10 (next event)                                                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         POST-PROCESSING (After Loop)                         │
│                                                                              │
│  1. Write binned distribution to CONTOUR file:                              │
│     → For I=1 to NI, J=1 to NJ:                                             │
│        WRITE(7,114): Azimuth, Elevation, ANGLE(I,J)                         │
│                                                                              │
│  2. Print statistics: WRITE(*,*) NVO, NVI, NIN                              │
│     → NVO: Count of VOID events (noise/invalid)                             │
│     → NVI: Count of VALID events (physical signals)                         │
│     → NIN: Count of INITIAL events (incomplete/edge cases)                  │
│                                                                              │
│  3. Close all files and STOP                                                │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 7.2 Text Flow Chart with GOTO Targets

**GOTO Reference Map:**
| Label | Line | Description | Target Action |
|-------|------|-------------|---------------|
| `GOTO 10` | Line 324 | Main loop continue | `DO 10 K=1,KKK` - next event |
| `GOTO 11` | Line 160 | Write VOID | NVO++, write to file 6, GOTO 10 |
| `GOTO 20` | Line 216 | Skip to end calc | Calculate Xe, Ye |
| `GOTO 22` | Line 241 | Write INITIAL | NIN++, write to file 5, GOTO 10 |

```
START
  │
  ▼
[Initialize] ─────────────────────────────────────────────────────────────┐
  • Read ITERANCE.NML                                                      │
  • Open files (3=input, 4=VALID, 5=INITIAL, 6=VOID, 7=CONTOUR)            │
  • NVO=0, NVI=0, NIN=0                                                    │
  │                                                                        │
  ▼                                                                        │
[DO 10 K=1,KKK] ◄──────────────────────────────────────────────────────────┤
  │                          ▲        ▲        ▲        ▲                  │
  │                          │        │        │        │                  │
  ▼                          │        │        │        │                  │
[Read TDC data]              │        │        │        │                  │
  │                          │        │        │        │                  │
  ▼                          │        │        │        │                  │
[Preprocess times]           │        │        │        │                  │
  • Swap TOF if Tof0>Tof1    │        │        │        │                  │
  • Rollover correction      │        │        │        │                  │
  │                          │        │        │        │                  │
  ▼                          │        │        │        │                  │
(Zero signal check?) ──YES──►[GOTO 11]───────┐   │        │                  │
  │ (X1e=0 & X2e=0) OR    (Line 160)         │   │        │                  │
  │ (Y1e=0 & Y2e=0)                          │   │        │                  │
  │                                          │   │        │                  │
  NO                                         │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
[Calculate Xs, Ys]                           │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
(KX=0 AND KY=0?) ──YES──►[GOTO 20]           │   │        │                  │
  │                       (Line 216)         │   │        │                  │
  │                       Tofs=MAX()         │   │        │                  │
  │                                          │   │        │                  │
  NO                                         │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
(KX=1 AND KY=1?) ──YES──►[GOTO 11]           │   │        │                  │
  │                       (Line 160)         │   │        │                  │
  │                       WRITE VOID         │   │        │                  │
  │                       NVO++              │   │        │                  │
  │                                          │   │        │                  │
  NO                                         │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
[Reconstruct from single signal]             │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
(|Xs|≤45 AND |Ys|≤15?) ──NO──►[GOTO 11]      │   │        │                  │
  │                           (Line 160)     │   │        │                  │
  │                                          │   │        │                  │
  YES                                        │   │        │                  │
  │                                          │   │        │                  │
  ▼                                          │   │        │                  │
  └─────────────────────┬────────────────────┘   │        │                  │
                        │                        │        │                  │
                        ▼                        │        │                  │
              [Calculate Xe, Ye] ◄───────────────┘        │                  │
               (Label 20, Line 216)                       │                  │
                • Calculate Tofx1, Tofy1                  │                  │
                • KKX, KKY flags                          │                  │
                        │                                 │                  │
                        ▼                                 │                  │
              (KKX=1 AND KKY=1?) ──YES──►[GOTO 22]        │                  │
                │                        (Line 241)       │                  │
                │                        WRITE INITIAL    │                  │
                │                        NIN++            │                  │
                │                                         │                  │
                NO                                        │                  │
                │                                         │                  │
                ▼                                         │                  │
              [Calculate Tofe]                            │                  │
                • Tofe = MIN(Tofy1,Tofx1)                 │                  │
                │                                         │                  │
                ▼                                         │                  │
              (|Xe|≤75 AND |Ye|≤27.5?) ──NO──►[GOTO 22]   │                  │
                │                              (Line 241) │                  │
                │                                         │                  │
                YES                                       │                  │
                │                                         │                  │
                ▼                                         │                  │
              [Final validation]                          │                  │
                • 6≤DTof≤187?                             │                  │
                • |DX|≤45?                                │                  │
                • |DY|≤15?                                │                  │
                │                                         │                  │
                ▼                                         │                  │
              (Any check FAIL?) ──YES──►[GOTO 22]         │                  │
                │                        (Line 241)       │                  │
                │                                         │                  │
                NO                                        │                  │
                │                                         │                  │
                ▼                                         │                  │
              [Particle Identification]                   │                  │
                • DTof<40ns → H                           │                  │
                • DTof>66.3ns → O                         │                  │
                • 40-66.3ns → check AMP                   │                  │
                │                                         │                  │
                ▼                                         │                  │
              [Calculate Azi, Elv]                        │                  │
              [Bin to ANGLE(MI,MJ)]                       │                  │
                │                                         │                  │
                ▼                                         │                  │
              [WRITE VALID]                               │                  │
                NVI++                                     │                  │
                WRITE(4,113) ...                          │                  │
                │                                         │                  │
                └─────────────────┬───────────────────────┘                  │
                                  │                                          │
                                  ▼                                          │
              [GOTO 10] ─────────────────────────────────────────────────────┘
               (Line 324)
                CONTINUE to next K

[After Loop]
  │
  ▼
[Write CONTOUR file]
[Print NVO, NVI, NIN]
[Close files]
  │
  ▼
STOP
```

**GOTO Target Code Details:**

| GOTO | Target Line | Code Snippet |
|------|-------------|--------------|
| `GOTO 10` | Line 324 | `10 CONTINUE` - loop end, next iteration |
| `GOTO 11` | Line 160 | `11 WRITE(6,111) ...` → `NVO=NVO+1` → `GOTO 10` |
| `GOTO 20` | Line 216 | `20 KKX=0` - start of end position calc block |
| `GOTO 22` | Line 241 | `22 WRITE(5,112) ...` → `NIN=NIN+1` → `GOTO 10` |

---

### 7.3 Decision Tree Summary

| Stage | Check | Pass Action | Fail Action |
|-------|-------|-------------|-------------|
| 1 | Zero signal? | Continue | VOID |
| 2 | Valid start signals? | Calculate Xs, Ys | Reconstruct or VOID |
| 3 | Start in bounds? | Continue to end calc | VOID |
| 4 | Valid end signals? | Calculate Xe, Ye | INITIAL |
| 5 | End in bounds? | Continue | INITIAL |
| 6 | TOF valid (6-187ns)? | Continue | INITIAL |
| 7 | Displacement valid? | Continue | INITIAL |
| 8 | All checks passed | VALID output | - |

### 7.4 Output File Summary

| File | Counter | Description | Records |
|------|---------|-------------|---------|
| VOID | NVO | Zero/missing signals, out of bounds | 15 fields: date, time, TOF, X/Y TDC, AMP |
| INITIAL | NIN | Incomplete data, edge cases, failed checks | 17 fields: + Xs, Ys reconstructed |
| VALID | NVI | Physical signals with full reconstruction | 18 fields: + NPar, Azi, Elv, Eng |
| CONTOUR | - | Binned angular distribution (31×17 grid) | 3 fields: Azimuth, Elevation, Count |

## 8. Summary

The Fortran code implements the coordinate transformation physics described in the TDC_Coordinate_System.md document correctly in principle. However, there is a **critical bug** where `Cofy` is used instead of `CofX` in two locations, which will cause incorrect X coordinate calculations. Additionally, there are concerns about sign convention consistency and hardcoded values that should be addressed for maintainability and correctness.





---

3. **Rollover:** Verify `TX_1e < TX_1s` cases get +663.47 correction
