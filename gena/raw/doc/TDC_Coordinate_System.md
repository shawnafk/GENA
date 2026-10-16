# TDC to Coordinate System Documentation

## Overview

This document describes the physics of converting TDC (Time-to-Digital Converter) time measurements to spatial coordinates (X, Y) for a position-sensitive detector using snake-line propagation delay.

---

## 0. Channel Mapping

### 0.1 Probe A Configuration

| Direction   | TDC Channel | Signal | Notes                            |
| ----------- | ----------- | ------ | -------------------------------- |
| Left (-X)   | 5           | TX1    | Time recorded at left receiver   |
| Right (+X)  | 7           | TX2    | Time recorded at right receiver  |
| Top (-Y)    | 1           | TY1    | Time recorded at top receiver    |
| Bottom (+Y) | 3           | TY2    | Time recorded at bottom receiver |

### 0.2 Probe B Configuration

| Direction   | TDC Channel | Signal | Notes                            |
| ----------- | ----------- | ------ | -------------------------------- |
| Left (-X)   | 1           | TX1    | Time recorded at left receiver   |
| Right (+X)  | 3           | TX2    | Time recorded at right receiver  |
| Top (-Y)    | 5           | TY1    | Time recorded at top receiver    |
| Bottom (+Y) | 7           | TY2    | Time recorded at bottom receiver |

---

## 1. Physical Model

### 1.1 Signal Propagation

When a particle/signal hits the board:
1. The hit triggers the TDC to **start** counting
2. Signal pulses propagate along X and Y snake lines (蛇形线) in both directions
3. When pulses reach receivers TX1 (left) and TX2 (right), TDC **stops** and records times

```
   (—)      e1 (extension)        L           e2 (extension) (+)
    ←───────────────────←─────────────────────→───────────────────→
    │                   │                     │                   │
  TX1                 x=-75                 x=+75                TX2
  (left)            (boundary)            (boundary)           (right)
                              ↑
                         Hit point at position x
```

Where:

- $e_1$ = extension length beyond left boundary
- $e_2$ = extension length beyond right boundary
- $L$ = length of the transmission line (snake line) inside the board, different from board length, $d_1+d_2 = 150 mm \ll L$

### 1.2 Time Equations

For a hit at distance **d₁** from left boundary and **d₂** from right boundary, the signal transmitted through line $L_1 + e_1$ and the rightward $L_2 + e_2$, assuming $v$ = signal propagation velocity (constant), we have

$$
TX_1 = \frac{L_1 + e_1}{v}
$$

$$
TX_2 = \frac{L_2 + e_2}{v}
$$

---



## 2. Coordinate Calculation

The coordinate on the board can be derived from time difference $\Delta t = TX_1 - TX_2$ .

### 2.1 Position in Center-Origin Coordinates

Inside the board, the width of the mcp and the length of the snake line that the signal went through has a linear relation:
$$
D = k L
$$
The coordinate (relative to center) is:

$$
x = \frac{d_1 - d_2}{2}  = \frac{k}{2}(L_1 - L_2)
$$

using equations 1 and 2

$$
x = \frac{ k }{2} \left( v (TX_1 - TX_2) - (e_1 - e_2)\right)
$$

### 2.2 Linear Form

This is a linear relationship:

$$
\boxed{x = K \cdot \Delta t + B}
$$

Where:
- $K = k v /2$ , the slope inside the board
- $B = -k(e_1 - e_2)/2$ (intercept/offset)

---

## 3. Calibration Procedure

### 3.1 Boundary Measurements

To determine the coefficiences, we have measured time difference $\Delta t = TX_1 - TX_2$ (time difference, always TX1 - TX2) at left and right boundary of the board:

| Position       | Coordinate     | Time Difference   |
| -------------- | -------------- | ----------------- |
| Left boundary  | $x_L = -75$ mm | $\Delta t_L = XL$ |
| Right boundary | $x_R = +75$ mm | $\Delta t_R = XR$ |

From the physics model:

$$
XL = \frac{e_1 - (L + e_2)}{v} = \frac{e_1 - e_2 - L}{v} = \frac{-b - L}{v}
$$

$$
XR = \frac{(L + e_1) - e_2}{v} = \frac{L + e_1 - e_2}{v} = \frac{L - b}{v}
$$

where $b =  2B/k$

---



### 3.2 Solving for Calibration Constants

#### from physics model above directly:

$$
v = \frac{2L}{XR - XL}
$$

$$
b = -\frac{XL + XR}{2} \cdot v
$$



thus
$$
K = \frac{kv}{2} = \frac{1}{2} \frac{D}{L}\cdot\frac{2L}{XR - XL} = \boxed{\frac{D}{XR - XL}}
$$

$$
B = \frac{k b}{2} = -\frac{kv}{4}(XL+XR)=\boxed{-\frac{D}{2}\frac{XL+XR}{XR-XL}}
$$



#### or from linear relation

From the two boundary equations:
$$
K \cdot XL + b = -D/2
$$

$$
K \cdot XR + b = +D/2
$$

Subtracting:

$$
K(XR - XL) = D \quad \Rightarrow \quad \boxed{K = \frac{D}{XR - XL}}
$$

Solving for intercept:

$$
b = -\frac{XL + XR}{2} \cdot K = \boxed{-\frac{D}{2}\frac{XL + XR}{XR-XL}}
$$

 

## 5. Key Parameters

| Parameters | Value | Description |
|----------|-------------|----------|
| A_X_LEFT_NS | -54.91726 | TX1-TX2 of the signal hit on the minus X boundary of the MCP; XL |
| A_X_RIGHT_NS | 50.44742 | TX1-TX2 of the signal hit on the plus X boundary of the MCP; XR |
| A_Y_LEFT_NS | -23.59747（地）\-16.68598（天） | TY1-TY2 of the signal hit on the minus Y boundary of the MCP; YL |
| A_Y_RIGHT_NS | 34.71570（地）41.62719（天） | TY1-TY2 of the signal hit on the plus Y boundary of the MCP; YR |
| B_X_LEFT_NS | \-51.18360 | same but for probe B |
| B_X_RIGHT_NS | 53.40479 | same but for probe B |
| B_Y_LEFT_NS  | \-16.68598 | same but for probe B                                        |
| B_Y_RIGHT_NS | 41.62719   | same but for probe B                                        |
| X_MCP        | 145        | MCP X Length; D in above equation                    |
| Y_MCP        | 55         | MCP Y Length; D in above equation                            |

## 6. Derived Parameters

| Parameters | Value                                               | Description                          |
| ---------- | --------------------------------------------------- | ------------------------------------ |
| A/B_X_COF  | X_MCP/(A/B_X_RIGHT_NS - A/B_X_RIGHT_NS)             | K in above equation, for X direction |
| A/B_Y_COF  | The slope k: Y_MCP/(A/B_Y_RIGHT_NS - A/B_Y_LEFT_NS) | K in above equation, for Y direction |
| A/B_X_INT  | -A/B_X_COF \* (A/B_X_LEFT_NS + A/B_X_RIGHT_NS)/2    | B in above equation, for X direction |
| A/B_Y_INT  | -A/B_Y_COF \* (A/B_Y_LEFT_NS + A/B_Y_RIGHT_NS)/2    | B in above equation, for Y direction |
| X          | A/B_X_COF \* (TX1 - TX2) + A/B_X_INT                | X coordinate (Origin at center MCP)  |
| Y          | A/B_Y_COF \* (TY1 - TY2) + A/B_Y_INT                | Y coordinate (Origin at center MCP)  |
