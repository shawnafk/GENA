import numpy as np
import pandas as pd
import gena.filter as flt
from .decode import _normalize_probe

__all__ = ['repair_l1a', 'df1b_type']

#after repair and drop those can not be repaired, or still bad after repair.
df1b_type ={
        'UUID':          np.dtype('S36'),
        'PID':           np.uint8,
        'BM':            np.uint8,
        'BX':            np.uint16,

        'HV_C': np.float32,
        'HV_G': np.float32,
        'MJ2000_q1': np.float32,
        'MJ2000_q2': np.float32,
        'MJ2000_q3': np.float32,
        'MJ2000_Star_X': np.int32,
        'MJ2000_Star_Y': np.int32,
        'MJ2000_Star_Z': np.int32,


        'Timestamp':     np.float64,
        # Match the core L1 internal-time precision.
        'TimeInternal':  np.float64,
        'EvtLevel':      np.uint8,
        'EP':            np.float32,
        'Tof_0':         np.float32,
        'Tof_1':         np.float32,
        'TX_1s':         np.float32,
        'TX_1e':         np.float32,
        'TX_2s':         np.float32,
        'TX_2e':         np.float32,
        'TY_1s':         np.float32,
        'TY_1e':         np.float32,
        'TY_2s':         np.float32,
        'TY_2e':         np.float32,



        #mid result for repair
        'Tofs': np.float32,
        'Tofe': np.float32,
        'Tofx0': np.float32,
        'Tofy0': np.float32,
        'Tofx1': np.float32,
        'Tofy1': np.float32,
        'tmpX1': np.float32,
        'tmpX2': np.float32,
        'tmpY1': np.float32,
        'tmpY2': np.float32,
        'pop_TX_1e':      np.float32,
        'pop_TX_2e':      np.float32,
        'pop_TY_1e':      np.float32,
        'pop_TY_2e':      np.float32,
        'Type': np.uint8,
        'Process_Mask': np.uint32,

        #quality flags for repair
        'Xs_valid': np.uint8,
        'Ys_valid': np.uint8,
        'Xe_valid': np.uint8,
        'Ye_valid': np.uint8,
        'start_valid': np.uint8,
        'end_valid': np.uint8,
}




#!================================      repair_lu         ====================================
def _swap_tof01(df):
    tof0 = df['Tof_0'].copy()
    tof1 = df['Tof_1'].copy()
    idx = (tof0 > 2.5) & (tof1 < tof0)
    df.loc[idx, 'Tof_0'] = tof1
    df.loc[idx, 'Tof_1'] = tof0
    return 0

# 强制将默认参数转换为 float32
def _apply_rollover_correction(df, rollover_period=np.float32(663.47)):
    cond_x1e = (df['TX_1e'] != 0) & (df['TX_1e'] < df['TX_1s'])
    df.loc[cond_x1e, 'Process_Mask'] |= (1 << 1)
    df.loc[cond_x1e, 'TX_1e'] = df.loc[cond_x1e, 'TX_1e'] + rollover_period

    cond_x2e = (df['TX_2e'] != 0) & (df['TX_2e'] < df['TX_2s'])
    df.loc[cond_x2e, 'Process_Mask'] |= (1 << 2)
    df.loc[cond_x2e, 'TX_2e'] = df.loc[cond_x2e, 'TX_2e'] + rollover_period

    cond_y1e = (df['TY_1e'] != 0) & (df['TY_1e'] < df['TY_1s'])
    df.loc[cond_y1e, 'Process_Mask'] |= (1 << 3)
    df.loc[cond_y1e, 'TY_1e'] = df.loc[cond_y1e, 'TY_1e'] + rollover_period

    cond_y2e = (df['TY_2e'] != 0) & (df['TY_2e'] < df['TY_2s'])
    df.loc[cond_y2e, 'Process_Mask'] |= (1 << 4)
    df.loc[cond_y2e, 'TY_2e'] = df.loc[cond_y2e, 'TY_2e'] + rollover_period
    return 0

def _check_deltaT(df, params, prob):
    XL, XR = params[f'{prob.upper()}_XL'], params[f'{prob.upper()}_XR']
    YL, YR = params[f'{prob.upper()}_YL'], params[f'{prob.upper()}_YR']

    dt_x = df['TX_1s'] - df['TX_2s']
    dt_y = df['TY_1s'] - df['TY_2s']
    txs_vld = ((dt_x >= XL) & (dt_x <= XR))
    tys_vld = ((dt_y >= YL) & (dt_y <= YR))

    dt_x = df['TX_1e'] - df['TX_2e']
    dt_y = df['TY_1e'] - df['TY_2e']
    txe_vld = ((dt_x >= XL) & (dt_x <= XR))
    tye_vld = ((dt_y >= YL) & (dt_y <= YR))

    return txs_vld, tys_vld, txe_vld, tye_vld


# --- Start repair helpers ---
def _repair_y_from_x(df, cond, CX, CY, KY_slope, BY_offset):
    df.loc[cond, 'Type'] = 1

    df.loc[cond, 'Tofs'] = df.loc[cond, 'Tofx0']
    closer_to_2 = abs(df['TY_2s'] - df['Tofs']) < abs(df['TY_1s'] - df['Tofs'])

    # Y2 closer
    idx = cond & closer_to_2
    if idx.any():
        raw_TY_1s = df['TY_1s'].copy() # 备份原始数据，避免污染
        df.loc[idx, 'TY_2s'] = df.loc[idx, 'TY_2s'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TY_1s'] = CY - df.loc[idx, 'TY_2s']
        df.loc[idx, 'tmpY1'] = KY_slope * (df.loc[idx, 'TY_1s'] - df.loc[idx, 'TY_2s']) + BY_offset
        df.loc[idx,'Process_Mask'] |= (1 << 5)

        df.loc[idx, 'pop_TY_1e'] = np.float32(0.0) # 使用 float32 的 0.0
        y1e_ok = idx & ((raw_TY_1s - df['Tofs']) > CY)
        df.loc[y1e_ok, 'pop_TY_1e'] = raw_TY_1s[y1e_ok] # 使用原始数据判断和恢复
        df.loc[y1e_ok,'Process_Mask'] |= (1 << 6)

    # Y1 closer
    idx = cond & (~ closer_to_2)
    if idx.any():
        raw_TY_2s = df['TY_2s'].copy() # 备份原始数据
        df.loc[idx, 'TY_1s'] = df.loc[idx, 'TY_1s'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TY_2s'] = CY - df.loc[idx, 'TY_1s']
        df.loc[idx, 'tmpY1'] = KY_slope * (df.loc[idx, 'TY_1s'] - df.loc[idx, 'TY_2s']) + BY_offset
        df.loc[idx,'Process_Mask'] |= (1 << 7)

        df.loc[idx, 'pop_TY_2e'] = np.float32(0.0)
        y2e_ok = idx & ((raw_TY_2s - df['Tofs']) > CY)
        df.loc[y2e_ok, 'pop_TY_2e'] = raw_TY_2s[y2e_ok] # 使用原始数据判断和恢复
        df.loc[y2e_ok,'Process_Mask'] |= (1 << 8)


def _repair_x_from_y(df, cond, CX, CY, KX_slope, BX_offset):
    df.loc[cond, 'Type'] = 2
    df.loc[cond, 'Tofs'] = df.loc[cond, 'Tofy0']
    closer_to_2 = abs(df['TX_2s'] - df['Tofs']) < abs(df['TX_1s'] - df['Tofs'])

    # X2 closer
    idx = cond & closer_to_2
    if idx.any():
        raw_TX_1s = df['TX_1s'].copy()
        df.loc[idx, 'TX_2s'] = df.loc[idx, 'TX_2s'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TX_1s'] = CX - df.loc[idx, 'TX_2s']
        df.loc[idx, 'tmpX1'] = KX_slope * (df.loc[idx, 'TX_1s'] - df.loc[idx, 'TX_2s']) + BX_offset
        df.loc[idx,'Process_Mask'] |= (1 << 9)

        df.loc[idx, 'pop_TX_1e'] = np.float32(0.0)
        x1e_ok = idx & ((raw_TX_1s - df['Tofs']) > CX)
        df.loc[x1e_ok, 'pop_TX_1e'] = raw_TX_1s[x1e_ok]
        df.loc[x1e_ok,'Process_Mask'] |= (1 << 10)

    # X1 closer
    idx = cond & (~ closer_to_2)
    if idx.any():
        raw_TX_2s = df['TX_2s'].copy()
        df.loc[idx, 'TX_1s'] = df.loc[idx, 'TX_1s'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TX_2s'] = CX - df.loc[idx, 'TX_1s']
        df.loc[idx, 'tmpX1'] = KX_slope * (df.loc[idx, 'TX_1s'] - df.loc[idx, 'TX_2s']) + BX_offset
        df.loc[idx,'Process_Mask'] |= (1 << 11)

        df.loc[idx, 'pop_TX_2e'] = np.float32(0.0)
        x2e_ok = idx & ((raw_TX_2s - df['Tofs']) > CX)
        df.loc[x2e_ok, 'pop_TX_2e'] = raw_TX_2s[x2e_ok]

# --- End repair helpers ---

def _repair_y_from_x_end(df, cond, CX, CY, KY_slope, BY_offset):
    df.loc[cond, 'Type'] = 3
    df.loc[cond, 'Tofe'] = df.loc[cond, 'Tofx1']
    closer_to_2 = abs(df['TY_2e'] - df['Tofe']) < abs(df['TY_1e'] - df['Tofe'])

    # Y2 closer
    idx = cond & closer_to_2
    if idx.any():
        pure_board_ty2e = df.loc[idx, 'TY_2e'] - df.loc[idx, 'Tofe']
        pure_board_ty1e = CY - pure_board_ty2e
        Dtof = df.loc[idx, 'Tofe'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TY_1e'] =  Dtof + pure_board_ty1e
        df.loc[idx, 'tmpY2'] = KY_slope * (df.loc[idx, 'TY_1e'] - df.loc[idx, 'TY_2e']) + BY_offset
        df.loc[idx,'Process_Mask'] |= (1 << 12)

    # Y1 closer
    idx = cond & (~ closer_to_2)
    if idx.any():
        pure_board_ty1e = df.loc[idx, 'TY_1e'] - df.loc[idx, 'Tofe']
        pure_board_ty2e = CY - pure_board_ty1e
        Dtof = df.loc[idx, 'Tofe'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TY_2e'] = Dtof + pure_board_ty2e
        df.loc[idx, 'tmpY2'] = KY_slope * (df.loc[idx, 'TY_1e'] - df.loc[idx, 'TY_2e']) + BY_offset
        df.loc[idx,'Process_Mask'] |= (1 << 13)


def _repair_x_from_y_end(df, cond, CX, CY, KX_slope, BX_offset):
    df.loc[cond, 'Type'] = 4
    df.loc[cond, 'Tofe'] = df.loc[cond, 'Tofy1']
    closer_to_2 = abs(df['TX_2e'] - df['Tofe']) < abs(df['TX_1e'] - df['Tofe'])

    # X2 closer
    idx = cond & closer_to_2
    if idx.any():
        pure_board_tx2e = df.loc[idx, 'TX_2e'] - df.loc[idx, 'Tofe']
        pure_board_tx1e = CX - pure_board_tx2e
        Dtof = df.loc[idx, 'Tofe'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TX_1e'] = Dtof + pure_board_tx1e
        df.loc[idx, 'tmpX2'] = KX_slope * (df.loc[idx, 'TX_1e'] - df.loc[idx, 'TX_2e']) + BX_offset
        df.loc[idx,'Process_Mask'] |= (1 << 14)

    # X1 closer
    idx = cond & (~ closer_to_2)
    if idx.any():
        pure_board_tx1e = df.loc[idx, 'TX_1e'] - df.loc[idx, 'Tofe']
        pure_board_tx2e = CX - pure_board_tx1e
        Dtof = df.loc[idx, 'Tofe'] - df.loc[idx, 'Tofs']
        df.loc[idx, 'TX_1e'] = Dtof + pure_board_tx2e
        df.loc[idx, 'tmpX2'] = KX_slope * (df.loc[idx, 'TX_1e'] - df.loc[idx, 'TX_2e']) + BX_offset
        df.loc[idx,'Process_Mask'] |= (1 << 15)


def _repair_cascading_end(df, cond, CX, CY, XL, XR, YL, YR,
                          KX_slope, BX_offset, KY_slope, BY_offset):
    df.loc[cond, 'Type'] = 5
    remaining = cond.copy()

    # Fallback 1: recovered X1 end signal
    fb = remaining & (df['pop_TX_1e'] != 0)
    if fb.any():
        new_dt = df['pop_TX_1e'] - df['TX_2e']
        valid = fb & (new_dt >= XL) & (new_dt <= XR)
        if valid.any():
            #df.loc[valid, 'Tofx1'] = (df.loc[valid, 'pop_TX_1e'] + df.loc[valid, 'TX_2e'] - CX) / fp32_2
            #df.loc[valid, 'TX_1e'] = df.loc[valid, 'pop_TX_1e'] - df.loc[valid, 'Tofx1']
            #df.loc[valid, 'TX_2e'] = df.loc[valid, 'TX_2e'] - df.loc[valid, 'Tofx1']
            df.loc[valid, 'TX_1e'] = df.loc[valid, 'pop_TX_1e']
            df.loc[valid, 'Tofx1'] = (df.loc[valid, 'TX_1e'] + df.loc[valid, 'TX_2e'] - CX) / fp32_2
            _repair_y_from_x_end(df, valid, CX, CY, KY_slope, BY_offset)
            df.loc[valid,'Process_Mask'] |= (1 << 16)
            remaining &= ~valid

    # Fallback 2: recovered X2 end signal
    fb = remaining & (df['pop_TX_2e'] != 0)
    if fb.any():
        new_dt = df['TX_1e'] - df['pop_TX_2e']
        valid = fb & (new_dt >= XL) & (new_dt <= XR)
        if valid.any():
            df.loc[valid, 'TX_2e'] = df.loc[valid, 'pop_TX_2e']
            df.loc[valid, 'Tofx1'] = (df.loc[valid, 'TX_1e'] + df.loc[valid, 'pop_TX_2e'] - CX) / fp32_2
            df.loc[valid, 'tmpX2'] = KX_slope * (df.loc[valid, 'TX_1e'] - df.loc[valid, 'pop_TX_2e']) + BX_offset
            _repair_y_from_x_end(df, valid, CX, CY, KY_slope, BY_offset)
            df.loc[valid,'Process_Mask'] |= (1 << 17)
            remaining &= ~valid

    # Fallback 3: recovered Y1 end signal
    fb = remaining & (df['pop_TY_1e'] != 0)
    if fb.any():
        new_dt = df['pop_TY_1e'] - df['TY_2e']
        valid = fb & (new_dt >= YL) & (new_dt <= YR)
        if valid.any():
            df.loc[valid, 'TY_1e'] = df.loc[valid, 'pop_TY_1e']
            df.loc[valid, 'Tofy1'] = (df.loc[valid, 'TY_1e'] + df.loc[valid, 'pop_TY_2e'] - CX) / fp32_2
            df.loc[valid, 'tmpY2'] = KY_slope * (df.loc[valid, 'TY_1e'] - df.loc[valid, 'TY_2e']) + BY_offset
            _repair_x_from_y_end(df, valid, CX, CY, KX_slope, BX_offset)
            df.loc[valid,'Process_Mask'] |= (1 << 18)
            remaining &= ~valid

    # Fallback 4: recovered Y2 end signal
    fb = remaining & (df['pop_TY_2e'] != 0)
    if fb.any():
        new_dt = df['TY_1e'] - df['pop_TY_2e']
        valid = fb & (new_dt >= YL) & (new_dt <= YR)
        if valid.any():
            df.loc[valid, 'TY_2e'] = df.loc[valid, 'pop_TY_2e']
            df.loc[valid, 'Tofy1'] = (df.loc[valid, 'TY_1e'] + df.loc[valid, 'pop_TY_2e'] - CX) / fp32_2
            df.loc[valid, 'tmpY2'] = KY_slope * (df.loc[valid, 'TY_1e'] - df.loc[valid, 'TY_2e']) + BY_offset
            _repair_x_from_y_end(df, valid, CX, CY, KX_slope, BX_offset)
            df.loc[valid,'Process_Mask'] |= (1 << 19)
            remaining &= ~valid

def _repair_start(df, params, prob):
    p = prob.upper()
    XL, XR = params[f'{p}_XL'], params[f'{p}_XR']
    YL, YR = params[f'{p}_YL'], params[f'{p}_YR']
    KX_slope = params[f'{p}_KX']
    KY_slope = params[f'{p}_KY']
    BX_offset = params[f'{p}_BX']
    BY_offset = params[f'{p}_BY']

    CX = max(abs(XR), abs(XL))
    CY = max(abs(YR), abs(YL))

    x_vld = (df['Xs_valid'] == 1)
    df.loc[x_vld, 'Tofx0'] = (df.loc[x_vld, 'TX_2s'] + df.loc[x_vld, 'TX_1s']  - CX) / fp32_2
    df.loc[x_vld, 'TX_1s'] = df.loc[x_vld, 'TX_1s'] - df.loc[x_vld, 'Tofx0']
    df.loc[x_vld, 'TX_2s'] = df.loc[x_vld, 'TX_2s'] - df.loc[x_vld, 'Tofx0']
    df.loc[x_vld, 'tmpX1'] = KX_slope * (df.loc[x_vld, 'TX_1s'] - df.loc[x_vld, 'TX_2s']) + BX_offset

    y_vld = (df['Ys_valid'] == 1)
    df.loc[y_vld, 'Tofy0'] = (df.loc[y_vld, 'TY_2s'] + df.loc[y_vld, 'TY_1s']  - CY) / fp32_2
    df.loc[y_vld, 'TY_1s'] = df.loc[y_vld, 'TY_1s'] - df.loc[y_vld, 'Tofy0']
    df.loc[y_vld, 'TY_2s'] = df.loc[y_vld, 'TY_2s'] - df.loc[y_vld, 'Tofy0']
    df.loc[y_vld, 'tmpY1'] = KY_slope * (df.loc[y_vld, 'TY_1s'] - df.loc[y_vld, 'TY_2s']) + BY_offset

    df.loc[x_vld & y_vld, 'Tofs'] = np.maximum(df.loc[x_vld & y_vld, 'Tofy0'], df.loc[x_vld & y_vld, 'Tofx0'],dtype=np.float32)

    _repair_y_from_x(df, x_vld & (~ y_vld), CX, CY, KY_slope, BY_offset)
    _repair_x_from_y(df, (~ x_vld) & y_vld, CX, CY, KX_slope, BX_offset)

    return 0


def _repair_end(df, params, prob):
    p = prob.upper()
    XL, XR = params[f'{p}_XL'], params[f'{p}_XR']
    YL, YR = params[f'{p}_YL'], params[f'{p}_YR']
    KX_slope = params[f'{p}_KX']
    KY_slope = params[f'{p}_KY']
    BX_offset = params[f'{p}_BX']
    BY_offset = params[f'{p}_BY']

    CX = max(abs(XR), abs(XL))
    CY = max(abs(YR), abs(YL))

    x_vld = (df['Xe_valid']==1) & (df['start_valid'] ==1)
    df.loc[x_vld, 'Tofx1'] = (df.loc[x_vld, 'TX_2e'] + df.loc[x_vld, 'TX_1e']  - CX) / fp32_2
    df.loc[x_vld, 'TX_1e'] = df.loc[x_vld, 'TX_1e'] - df.loc[x_vld, 'Tofx1']
    df.loc[x_vld, 'TX_2e'] = df.loc[x_vld, 'TX_2e'] - df.loc[x_vld, 'Tofx1']
    df.loc[x_vld, 'tmpX2'] = KX_slope * (df.loc[x_vld, 'TX_1e'] - df.loc[x_vld, 'TX_2e']) + BX_offset

    y_vld = (df['Ye_valid']==1) & (df['start_valid'] ==1)
    df.loc[y_vld, 'Tofy1'] = (df.loc[y_vld, 'TY_2e'] + df.loc[y_vld, 'TY_1e']  - CY) / fp32_2
    df.loc[y_vld, 'TY_1e'] = df.loc[y_vld, 'TY_1e'] - df.loc[y_vld, 'Tofy1']
    df.loc[y_vld, 'TY_2e'] = df.loc[y_vld, 'TY_2e'] - df.loc[y_vld, 'Tofy1']
    df.loc[y_vld, 'tmpY2'] = KY_slope * (df.loc[y_vld, 'TY_1e'] - df.loc[y_vld, 'TY_2e']) + BY_offset

    df.loc[x_vld & y_vld, 'Tofe'] = np.minimum(df.loc[x_vld & y_vld, 'Tofy1'], df.loc[x_vld & y_vld, 'Tofx1'],dtype=np.float32)

    _repair_y_from_x_end(df, x_vld & (~ y_vld), CX, CY, KY_slope, BY_offset)
    _repair_x_from_y_end(df, (~ x_vld) & y_vld, CX, CY, KX_slope, BX_offset)

    _repair_cascading_end(df, (~ x_vld) & (~ y_vld), CX, CY, XL, XR, YL, YR,
                          KX_slope, BX_offset, KY_slope, BY_offset)

    tof_reverse = df['Tofe'] < df['Tofs']
    if tof_reverse.any():
        tofs_tmp = df.loc[tof_reverse, 'Tofs'].copy()
        df.loc[tof_reverse, 'Tofs'] = df.loc[tof_reverse, 'Tofe']
        df.loc[tof_reverse, 'Tofe'] = tofs_tmp
    return 0


#! ===========================================  to do repair simple -================================
#for start over lim or sum not match
def _repair_start_2(df, params, prob, dir = 'X'):
    p = prob.upper()
    C = params[f'{p}_{dir}_SUM']

    ck1s = flt.check('length',dir=dir,edge='1',channel='s',params=params,prob=prob)
    ck1e = flt.check('length',dir=dir,edge='1',channel='e',params=params,prob=prob)
    ck2s = flt.check('length',dir=dir,edge='2',channel='s',params=params,prob=prob)
    ck2e = flt.check('length',dir=dir,edge='2',channel='e',params=params,prob=prob)
    cksum_ss = flt.check('sum', dir=dir, channel=('s','s'), params=params,prob=prob)
    cksum_se = flt.check('sum', dir=dir, channel=('s','e'), params=params,prob=prob)
    cksum_es = flt.check('sum', dir=dir, channel=('e','s'), params=params,prob=prob)
    cksum_ee = flt.check('sum', dir=dir, channel=('e','e'), params=params,prob=prob)

    # 2e 才是真 2s。此时原本的 2e 绝对是垃圾（丢弃）。
    cond_func = ~cksum_ss & ck1s & ~ck2s & ck2e & cksum_se
    cond1 = cond_func(df)
    df.loc[cond1,f'pop_T{dir}_2e'] = df.loc[cond1,f'T{dir}_2s']
    df.loc[cond1,f'T{dir}_2s'] = df.loc[cond1,f'T{dir}_2e']
    df.loc[cond1,'Type'] = 1
       # 1e 才是真 1s。此时原本的 1e 绝对是垃圾（丢弃）。
    cond_func = ~cksum_ss & ck2s & ~ck1s & ck1e & cksum_es
    cond2 = cond_func(df)
    df.loc[cond2,f'pop_T{dir}_1e'] = df.loc[cond2,f'T{dir}_1s']
    df.loc[cond2,f'T{dir}_1s'] = df.loc[cond2,f'T{dir}_1e']
    df.loc[cond1,'Type'] = 2

    cond_func = (~cksum_ss) & ck1s & (~ck2s) & (~cksum_se)
    cond_recon_2s = cond_func(df)
    df.loc[cond_recon_2s, f'pop_T{dir}_2s'] = df.loc[cond_recon_2s, f'T{dir}_2s']
    df.loc[cond_recon_2s, f'T{dir}_2s'] = C - df.loc[cond_recon_2s, f'T{dir}_1s']
    df.loc[cond_recon_2s, 'Type'] = 11

    cond_func = ~cksum_ss & ck2s & (~ck1s) & (~cksum_es)
    cond_recon_1s = cond_func(df)
    df.loc[cond_recon_1s, f'pop_T{dir}_1s'] = df.loc[cond_recon_1s, f'T{dir}_1s']
    df.loc[cond_recon_1s, f'T{dir}_1s'] = C - df.loc[cond_recon_1s, f'T{dir}_2s']
    df.loc[cond_recon_1s, 'Type'] = 12


    # all e are real s。
    cond_func = ~cksum_ss & ~ck1s & ~ck2s & ck2e & ck1e & cksum_ee
    cond4 = cond_func(df)
    df.loc[cond4,f'pop_T{dir}_1e'] = df.loc[cond4,f'T{dir}_1s']
    df.loc[cond4,f'pop_T{dir}_2e'] = df.loc[cond4,f'T{dir}_2s']
    df.loc[cond4,f'T{dir}_1s'] = df.loc[cond4,f'T{dir}_1e']
    df.loc[cond4,f'T{dir}_2s'] = df.loc[cond4,f'T{dir}_2e']
    df.loc[cond_recon_1s, 'Type'] = 21

    other = ~(cond1|cond2|cond_recon_1s|cond_recon_2s|cond4)
    df.loc[other,'Type'] = 255


def _repair_method_1(df1a, params, prob):
    if len(df1a) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in df1b_type.items()})

    df1b = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in df1b_type.items()})

    base_cols = ['UUID', 'PID', 'BM', 'BX', 'Timestamp', 'TimeInternal', 'EvtLevel',
                 'HV_C', 'HV_G',
                 'MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3',
                 'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z',
                 'EP', 'Tof_0', 'Tof_1',
                 'TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
                 'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']
    for k in base_cols:
        df1b[k] = df1a[k]

    df1b['Type'] = np.uint8(255)
    df1b['Process_Mask'] = np.uint32(0)
    df1b.loc[flt.time_tof_eff()(df1b),'Type'] = np.uint8(0)

    quality_cols = ['Xs_valid', 'Ys_valid', 'Xe_valid', 'Ye_valid', 'start_valid','end_valid']
    for k in quality_cols:
        df1b[k] = 1

    mid = ['Tofs', 'Tofe', 'Tofx0', 'Tofy0', 'Tofx1', 'Tofy1', 'tmpX1', 'tmpX2', 'tmpY1', 'tmpY2']
    for k in mid:
        df1b[k] = np.float32(np.nan) # 使用 float32 包装 nan

    for k in ['pop_TX_1e', 'pop_TX_2e', 'pop_TY_1e', 'pop_TY_2e']:
        df1b[k] = np.float32(0.0)

    _swap_tof01(df1b)
    _apply_rollover_correction(df1b)

    df1b['Xs_valid'],df1b['Ys_valid'],df1b['Xe_valid'],df1b['Ye_valid'] = _check_deltaT(df1b, params, prob)

    _repair_start(df1b, params, prob)

    X1, Y1 = ns_2_xy(df1b, 0, params=params,prob=prob)
    df1b['start_valid'] = (abs(X1) < 45) & (abs(Y1) < 15)

    _repair_end(df1b,params, prob)

    X2, Y2 = ns_2_xy(df1b, 1, params=params,prob=prob)
    df1b['end_valid'] = (abs(X2) < 75) & (abs(Y2) < 30)

    df1b.loc[~(df1b['start_valid'] & df1b['start_valid']),'Type'] = 255
    return df1b


def _repair_method_2(df1a, params, prob):
    if len(df1a) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in df1b_type.items()})

    df1b = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in df1b_type.items()})

    base_cols = ['UUID', 'PID', 'BM', 'BX', 'Timestamp', 'TimeInternal', 'EvtLevel',
                 'HV_C', 'HV_G',
                 'MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3',
                 'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z',
                 'EP', 'Tof_0', 'Tof_1',
                 'TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
                 'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']
    for k in base_cols:
        df1b[k] = df1a[k]

    df1b['Type'] = np.uint8(255)
    df1b['Process_Mask'] = np.uint32(0)
    df1b.loc[flt.time_tof_eff()(df1b),'Type'] = np.uint8(0)

    quality_cols = ['Xs_valid', 'Ys_valid', 'Xe_valid', 'Ye_valid', 'start_valid','end_valid']
    for k in quality_cols:
        df1b[k] = 1

    mid = ['Tofs', 'Tofe', 'Tofx0', 'Tofy0', 'Tofx1', 'Tofy1', 'tmpX1', 'tmpX2', 'tmpY1', 'tmpY2']
    for k in mid:
        df1b[k] = np.float32(np.nan) # 使用 float32 包装 nan

    for k in ['pop_TX_1e', 'pop_TX_2e', 'pop_TY_1e', 'pop_TY_2e']:
        df1b[k] = np.float32(0.0)
    '''
    below check for 1s + 2s do not match, eigher only 1s is wrong, or only 2s is wrong, or both are wrong
    #if condtion 1
    #the lim condition is not sufficient, < lim could also a wrong sigal, e.g. a small noise, > lim must be wrong;
    if 1s > lim and 2s < lim:
        #assuming 1s is not good, 2s is good
        #check 1e can be used as 1s, and pop 1s as potential 1e
        pop_1e = 1s
        #signal match condition is much more stronger
        #if 1e matched with 2s, use 1e as 1s
        if 1e < lim and 1e + 2s match C:
            1s = 1e
        else:
            if 1e can not match 2s, it could be:
                A. original 1e might not a start signal (1e>lim) and  ~B
                B. 2s is broken too, and ~A
                C. both.
            if A, means no valid start signal for edge 1, so we need find correct 2s to reconstructe 1s, also we need test 1e is true stop sigal or noise (how), also we have poped 1s as pop_1e to be candicate stop signal, need test pop_1e as well.
            if B, need test if 2e can match with origal 1e; if 2e can, use 2e as 2s, pop 2s as 2e (but 2s is small, proberply not a stop sigal, just do it and check later) and test pop_2e is correct or noise.
            if C. need test and reconstructure all, 1s and 2s and 1e and 2e
    if 1s < lim and 2s >lim:
    ...
    if 1s > lim and 2s > lim
        #test 1s - 2s with in L and R #test 1e + 2e match C and 1e - 2e within L, R # test (1s+2s) - (1e+2e) is positive
        then revert
        if not, try rollover repair
    '''
    #so you can recalculate start for edge s, using correct "1s" (which is 1e orignally, just because 1s is broken), but you are not sure the current "1s" is correct? or, you can do check on edge 2 from the beginning, and fix 1 with 2s, will the two fix with 1/2 be the same result, how do you know

    df1b['Xs_valid'],df1b['Ys_valid'],df1b['Xe_valid'],df1b['Ye_valid'] = _check_deltaT(df1b, params, prob)


    _repair_start_2(df1b, params, prob, 'X')
    _repair_start_2(df1b, params, prob, 'Y')

    X1, Y1 = ns_2_xy(df1b, 0, params=params,prob=prob)
    df1b['start_valid'] = (abs(X1) < 45) & (abs(Y1) < 15)

    _repair_end(df1b,params, prob)

    X2, Y2 = ns_2_xy(df1b, 1, params=params,prob=prob)
    df1b['end_valid'] = (abs(X2) < 75) & (abs(Y2) < 30)

    df1b.loc[~(df1b['start_valid'] & df1b['start_valid']),'Type'] = 255
    return df1b

#get repaired level 1a data
def repair_l1a(df1a, params, prob,method=2):
    """Reject legacy repair algorithms until their timing model is validated.

    The historical implementations remain below their helper definitions
    for audit. Method 1 destroys the common start/stop time reference;
    method 2 does not initialize the start reference and overwrites quality
    labels. Neither is a supported source of scientific L1B/L2B products.
    """
    _normalize_probe(prob)
    if method not in (1, 2):
        raise ValueError("method must be 1 or 2")
    raise NotImplementedError(
        f"Legacy repair method {method} is disabled: its start/stop timing "
        "and quality flags are not validated. Use unmodified L1A/L2 products; "
        "do not promote reconstructed events to golden data."
    )



'''

cond1 = (df['X1_start (ns)'] > df['X1_end (ns)']) & (df['X1_end (ns)'] !=0)
cond2 = (df['X2_start (ns)'] > df['X2_end (ns)']) & (df['X2_end (ns)'] !=0)
subdf = df[cond1 & cond2]
a = (subdf['X1_start (ns)'] +  subdf['X2_start (ns)'])
b = (subdf['X2_end (ns)'] +  subdf['X1_end (ns)'] + 8192 * 0.081)


def revise_df(df):
    # 创建副本避免修改原数据
    df = df.copy()

    # 修正X1时间
    idx1 = (df['X1_end (ns)'] - df['X1_start (ns)'])<0
    df.loc[idx1, 'X1_end (ns)'] = df.loc[idx1, 'X1_end (ns)']  + 8192 * 0.081

    # 修正X2时间
    idx2 = (df['X2_end (ns)'] - df['X2_start (ns)'])<0
    df.loc[idx2, 'X2_end (ns)'] = df.loc[idx2, 'X2_end (ns)']  + 8192 * 0.081

    # 修正Y1时间
    idx3 = (df['Y1_end (ns)'] - df['Y1_start (ns)'])<0
    df.loc[idx3, 'Y1_end (ns)'] = df.loc[idx3, 'Y1_end (ns)']  + 8192 * 0.081  # 修正idx1->idx3

    # 修正Y2时间
    idx4 = (df['Y2_end (ns)'] - df['Y2_start (ns)'])<0
    df.loc[idx4, 'Y2_end (ns)'] = df.loc[idx4, 'Y2_end (ns)']  + 8192 * 0.081  # 修正idx1->idx4

    return df  # 添加返回值

'''
