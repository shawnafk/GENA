import sys
import numpy as np
import pandas as pd
from . import constants as const
import gena.filter as flt

from . import structure as st
from . import utils

fp32_2 = np.float32(2.0)

#df0 to df1 filter1 df1_flt to df1a filter2 to df1a_flt to df2 
def tdc_by_index(tdc):
    index = (tdc <= 94)
    # 显式指定 dtype 为 float32，避免默认的 float64
    converted = np.empty(len(tdc), dtype=np.float32)
    converted[index] = 0
    # 常数强制转为 float32
    converted[~ index] = (tdc[~ index] - 94) * np.float32(0.081)
    return converted

def ep_by_index(ep):
    index = (ep == 4095) | (ep ==0)
    # 显式指定 dtype 为 float32
    converted = np.empty(len(ep), dtype=np.float32)
    converted[index] = 0
    # 常数强制转为 float32
    converted[~ index] = (2044 - ep[~ index]) * np.float32(100.0 / 123.0)
    return converted
    
#snid 1 orbit
#snid 2 ground
#this fuction need return, between different levels, make df1 first
def l0_2_l1(df0, prob = 'a', snid = 1):
    '''
    Convert L0 raw to L1  
    '''
    #get num A and B
    if prob == 'a':
        index = df0['prob_flag'] == 0
        comb_xy_channel = list(zip( ['TX_1','TX_2','TY_1','TY_2'], ['ch_5','ch_7','ch_1','ch_3']))
    elif prob == 'b':
        index = df0['prob_flag'] == 1
        if snid == 1:
            comb_xy_channel = list(zip(['TX_1','TX_2','TY_1','TY_2'], ['ch_1','ch_3','ch_5','ch_7']))
        elif snid == 2:
            comb_xy_channel = list(zip(['TX_1','TX_2','TY_1','TY_2'], ['ch_7','ch_5','ch_3','ch_1']))
        else:
            print('wrong satellite id, should be 1 (orbit) or 2 (ground)')
            sys.exit(1)
    else:
        print('wrong probe flag, should be a or b')
        sys.exit(1)
    #init L1a and L1b dataframe with that size
    #init dict instead of dataframe to avoid index problem
    df1 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1_type.items()})
    for k in ['UUID','PID','BM','BX','Timestamp','TimeInternal','EvtLevel']:
        df1[k] = df0.loc[index,k]
    
    for tof,channel in zip(
            ['Tof_0','Tof_1','Tof_2','Tof_3'],
            ['ch_tof_0','ch_tof_1','ch_tof_2','ch_tof_3']
    ):
        df1[tof] = tdc_by_index(df0.loc[index,channel])
    
    #convert fly time to corresponting channel in A B probes
    #use list avoid drain iterative zip
    for se in ['s','e']: 
        for xy,channel in comb_xy_channel:
            df1[xy+se] = tdc_by_index(df0.loc[index,channel+se])
    
    df1['EP'] = ep_by_index(df0.loc[index,'EP'])

    return df1

#also need return
def l1_2_l1a(df1, pdf):
    '''
    Convert L1 to L1A: Select 8-channel events and interpolate PDF parameters.
    L1A = "full 8-channel subset of L1" - only rows with all 8 TDC channels non-zero.
    Quaternion interpolation uses SLERP (not linear).
    '''
    if len(df1) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1a_type.items()})

    # Initialize output DataFrame with L1A schema
    df1a_out = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1a_type.items()})

    # Copy base columns from L1
    base_cols = ['UUID', 'PID', 'BM', 'BX', 'Timestamp', 'TimeInternal', 'EvtLevel', 'EP',
                 'Tof_0', 'Tof_1', 'Tof_2', 'Tof_3',
                 'TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
                 'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']
 
    # 需要检查不为 0 的列
    check_cols = ['TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
              'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']

    # 一键生成 mask：所有列 !=0
    mask = (df1[check_cols] != 0).all(axis=1)

    # 直接筛选 + 重置索引
    df1a_out = df1.loc[mask, base_cols].reset_index(drop=True)

    if len(df1a_out) == 0:
        return df1a_out

    # Interpolate HV parameters linearly
    df1a_out['HV_C'] = utils.get_para(pdf, 'HV_C', df1.loc[mask,'Timestamp'].values)
    df1a_out['HV_G'] = utils.get_para(pdf, 'HV_G', df1.loc[mask,'Timestamp'].values)

    # Quaternion interpolation via SLERP
    event_rots = utils.get_rot(pdf, df1.loc[mask,'Timestamp'].values)
    quats = event_rots.as_quat()  # shape (N, 4) in xyzw format
    df1a_out['MJ2000_q1'] = quats[:, 0]
    df1a_out['MJ2000_q2'] = quats[:, 1]
    df1a_out['MJ2000_q3'] = quats[:, 2]

    # Satellite position: cubic spline interpolation for trajectory consistency
    pos = utils.interp_trajectory(pdf, df1.loc[mask,'Timestamp'].values)
    df1a_out['MJ2000_Star_X'] = pos[:, 0]
    df1a_out['MJ2000_Star_Y'] = pos[:, 1]
    df1a_out['MJ2000_Star_Z'] = pos[:, 2]

    return df1a_out

#zhengjiangshan: be carefull, this only works for tdc being triggered normally 
def f_tof_mxy(df):
    sum_start_x1_x2 = df['TX_1s'] + df['TX_2s']
    sum_end_x1_x2 = df['TX_1e'] + df['TX_2e']
    tof_x = (sum_end_x1_x2 - sum_start_x1_x2)/2
    
    sum_start_y1_y2 = df['TY_1s'] + df['TY_2s']
    sum_end_y1_y2 = df['TY_1e'] + df['TY_2e']
    tof_y = (sum_end_y1_y2 - sum_start_y1_y2)/2
    tof_m = df['Tof_1']
    return tof_m, tof_x, tof_y

def f_tof_ep(df):
    '''Calculate electron time-of-flight through acceleration and flight regions.'''
    l_acc = 3.1 * 1e-3     # acceleration region length (m)
    l_fly = 47 * 1e-3      # flight region length (m)
    energy = abs(df['HV_C'] - df['HV_G'])
    ve = (2 * energy / const.MASS_E) ** 0.5 * const.SPEED_C
    time_acc = 2 * l_acc / ve * 1e9
    time_fly = l_fly / ve * 1e9
    # 强制拍扁回 float32
    return (time_acc + time_fly).astype('float32')

#1 is H, 2 is O, 3 is no tof, default 0 is Other
def ENA_spec(df):
    '''Classify particle species from TOF and EP columns in DataFrame.'''
    tof = df['Tof']
    ep = df['EP']
    fastest_O = 28.79
    slowest_H = 161.55
    spec = np.zeros(len(df), dtype=np.float32) # 指定类型
    spec[np.isnan(tof)] = 3
    faster_than_FO = (tof < fastest_O)
    slower_than_SH = (tof > slowest_H)
    between_SH_FO  = (~ faster_than_FO) & (~ slower_than_SH)
    spec[faster_than_FO] = 1
    spec[slower_than_SH] = 2
    spec[between_SH_FO & (ep <= 1580)] = 1   # H (low amplitude or no EP measurement)
    spec[between_SH_FO & (ep > 1580)] = 2    # O (high amplitude)
    return spec

#zhengjiangshan: when consider lu's method, this function may also be modified
def ENA_flytime(df):
    tof_m = df['Tof_M']
    tof_x = df['Tof_X']
    tof_y = df['Tof_Y']
    tof = np.zeros(len(tof_m), dtype=np.float32) # 指定类型
    
    cond_m = (tof_m > 0) 
    cond_x = (abs(tof_m - tof_x) < 5)
    cond_y = (abs(tof_m - tof_y) < 5)
    cond_xy = (abs(tof_x - tof_y) < 5) | ( (abs(tof_x/tof_y) > 0.9) & (abs(tof_x/tof_y) < 1.1))
    
    cond1 = cond_m & cond_x & cond_y
    tof[cond1] = tof_m[cond1]
    
    cond2 =( ~ cond1)  & cond_xy
    tof[cond2] = np.minimum(tof_x[cond2], tof_y[cond2])
    
    cond3 = ( ~ cond1) | ( ~ cond2)
    tof[cond3] = np.minimum(tof_x[cond3], tof_y[cond3])
    # 确保返回的是 float32
    return tof.astype('float32')

#EV
def ENA_energy(df):
    x1 = df['X1']
    x2 = df['X2']
    y1 = df['Y1']
    y2 = df['Y2']
    tof = df['Tof']
    spec = df['Spec']
    fly_distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + 50**2)
    fly_speed = fly_distance * 1e-3 / (tof*1e-9)/const.SPEED_C
    mass = np.ones(len(tof), dtype=np.float32)*const.MASS_H # 指定类型
    mass[spec == 2] = const.MASS_O
    mass[spec == 3] = 0
    energy = 1/2 * mass * fly_speed**2 /1000
    # 拍扁回 float32
    return energy.astype('float32')


# ---------------------------------------------------------------------------
# Position conversion
# ---------------------------------------------------------------------------


#level 2
def l1_2_l2(df1,params=const.init_param(1),prob='a'):
    if len(df1) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df2_type.items()})
    tof_m, tof_x, tof_y = f_tof_mxy(df1)
    df2 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df2_type.items()})
    for k in ['UUID','PID','BM','BX','Timestamp','TimeInternal','EvtLevel','EP','HV_G','HV_C','MJ2000_q1','MJ2000_q2','MJ2000_q3','MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z']:
        df2[k] = df1[k]

    tof_ep = f_tof_ep(df2)
    df2['Tof_X']  = tof_x
    df2['Tof_Y']  = tof_y
    df2['Tof_M'] = tof_m

    df2['Tof'] = ENA_flytime(df2) + tof_ep
    df2['Spec'] = ENA_spec(df2)
    
    x1,y1 = params['ns_2_xy'](df1,0,prob=prob)
    x2,y2 = params['ns_2_xy'](df1,1,prob=prob)

    df2['X1'] = x1
    df2['Y1'] = y1
    df2['X2'] = x2
    df2['Y2'] = y2

    df2['Energy'] = ENA_energy(df2)
    
    df2['PX'] = (np.arctan2(x1-x2,50)*180/np.pi).astype('float32')
    df2['PY'] = (np.arctan2(y1-y2,50)*180/np.pi).astype('float32')

    return df2

#!=====================================      repair_lu                     ====================================
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
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1b_type.items()})

    df1b = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1b_type.items()})

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
    df1b.loc[flt.golden_std()(df1b),'Type'] = np.uint8(0)
    
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
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1b_type.items()})

    df1b = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1b_type.items()})

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
    df1b.loc[flt.golden_std()(df1b),'Type'] = np.uint8(0)
    
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


def l1a_2_l1b(df1a, params, prob,method=2):
    if len(df1a) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1b_type.items()})
    if method == 1:
        return _repair_method_1(df1a,params,prob)
    if method == 2:
        return _repair_method_2(df1a,params,prob)
    else:
        sys.exit(1)
        
        
        
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
