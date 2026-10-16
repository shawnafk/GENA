# --- deprecated statistics functions for legacy GENA particle rows --------------------------------

def channel_nonzero_count(df, c):
    if c == "X1":
        return ((df["TX_1s"] != 0) | (df["TX_1e"] != 0)).sum()
    if c == "X2":
        return ((df["TX_2s"] != 0) | (df["TX_2e"] != 0)).sum()
    if c == "Y1":
        return ((df["TY_1s"] != 0) | (df["TY_1e"] != 0)).sum()
    if c == "Y2":
        return ((df["TY_2s"] != 0) | (df["TY_2e"] != 0)).sum()
    raise ValueError("c must be one of 'X1', 'X2', 'Y1', 'Y2'.")


def hist0d_xy(l1, name, tmin=5, tmax=633):
    eff = (l1[name] > tmin) & (l1[name] < tmax)
    return window_hist0d(l1[eff], window="1min")


def pspe(a, b):
    pe = np.roots([1, -(2 * a + b), a])
    pe = pe[pe < 1]
    ps = a / pe
    return ps, pe

