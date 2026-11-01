import sys
import numpy as np
import pandas as pd
from . import constants as const
import gena.filter as flt
from . import structure as st

fp32_2 = np.float32(2.0)


def _normalize_probe(prob):
    # Imported lazily: decode is the public owner, but it imports conversion.
    from .decode import _normalize_probe as normalize
    return normalize(prob)

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
    """Convert ADC codes without unsigned subtraction wraparound.

    Codes 0 and 4095 retain the telemetry convention for no EP measurement.
    Out-of-range or non-finite codes are invalid, represented by NaN.
    """
    values = np.asarray(ep, dtype=np.float64)
    missing = (values == 4095) | (values == 0)
    valid = np.isfinite(values) & (values > 0) & (values < 4095)
    converted = np.full(values.shape, np.nan, dtype=np.float32)
    converted[missing] = 0
    converted[valid] = (2044.0 - values[valid]) * (100.0 / 123.0)
    return converted



    
#snid 1 orbit
#snid 2 ground
#this fuction need return, between different levels, make df1 first
def l0_2_l1(df0, prob = 'a', snid = 1):
    '''
    Convert L0 raw to L1  
    '''
    prob = _normalize_probe(prob)
    if snid not in (1, 2):
        raise ValueError("snid must be 1 (orbit) or 2 (ground)")
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
    #init L1a and L1b dataframe with that size
    #init dict instead of dataframe to avoid index problem
    df1 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1_type.items()})
    for k in ['UUID','PacketID','EventID','PID','BM','BX','Timestamp','TimeInternal','EvtLevel']:
        if k in df0.columns:
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
class EngineeringIndex:
    """Packet-level lookup for PDF engineering values.

    Both HV and attitude/position are recorded per packet in the PDF. They are
    joined by PacketID; there is no event-time interpolation. That keeps all
    events in one packet on exactly the packet's recorded engineering values.
    """

    def __init__(self, pdf):
        if 'PacketID' not in pdf.columns:
            self.hv = None
            self.pose = None
            return
        deduped = pdf.drop_duplicates('PacketID').set_index('PacketID')
        self.hv = deduped[['HV_C', 'HV_G']]
        self.pose = deduped[[
            'MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3',
            'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z',
            'PosePresent',
        ]]

    def hv_at(self, packet_ids):
        if self.hv is None:
            return None
        return self.hv.reindex(np.asarray(packet_ids))

    def pose_at(self, packet_ids):
        """Return packet-level quaternion/position for the given PacketID rows."""
        if self.pose is None:
            return None
        pose = self.pose.reindex(np.asarray(packet_ids))
        present = pose['PosePresent'].fillna(False).to_numpy(dtype=bool)
        q = pose[['MJ2000_q1', 'MJ2000_q2', 'MJ2000_q3']].to_numpy(dtype=float)
        pos = pose[['MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z']].to_numpy(dtype=float)
        norm2 = np.sum(q * q, axis=1)
        valid = (present & np.isfinite(q).all(axis=1)
                 & np.isfinite(pos).all(axis=1) & (norm2 <= 1.0 + 1e-6))
        quat = np.full((len(packet_ids), 4), np.nan, dtype=float)
        result_pos = np.full((len(packet_ids), 3), np.nan, dtype=float)
        if valid.any():
            quat[valid, :3] = q[valid]
            quat[valid, 3] = np.sqrt(np.maximum(0.0, 1.0 - norm2[valid]))
            result_pos[valid] = pos[valid]
        return quat, result_pos



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
    """Electron flight correction in ns; invalid acceleration gives NaN."""
    l_acc = 3.1 * 1e-3     # acceleration region length (m)
    l_fly = 47 * 1e-3      # flight region length (m)
    energy = np.abs(np.asarray(df['HV_C'], dtype=np.float64)
                    - np.asarray(df['HV_G'], dtype=np.float64))
    valid = np.isfinite(energy) & (energy > 0)
    correction = np.full(len(df), np.nan, dtype=np.float32)
    ve = np.sqrt(2 * energy[valid] / const.MASS_E) * const.SPEED_C
    correction[valid] = (2 * l_acc + l_fly) / ve * 1e9
    return correction

#1 is H, 2 is O, 3 is no tof, default 0 is Other
def ENA_spec(df):
    """Classify valid TOF: 1=H, 2=O, 3=invalid TOF, 0=undetermined.

    The existing EP=0 hydrogen assumption is retained in the ambiguous TOF
    band. Non-finite EP in that band leaves the species undetermined.
    """
    tof = np.asarray(df['Tof'], dtype=np.float64)
    ep = np.asarray(df['EP'], dtype=np.float64)
    fastest_O = 28.79
    slowest_H = 161.55
    valid = np.isfinite(tof) & (tof > 0)
    spec = np.zeros(len(df), dtype=np.uint8)
    spec[~valid] = 3
    faster_than_FO = valid & (tof < fastest_O)
    slower_than_SH = valid & (tof > slowest_H)
    between_SH_FO = valid & ~faster_than_FO & ~slower_than_SH & np.isfinite(ep)
    spec[faster_than_FO] = 1
    spec[slower_than_SH] = 2
    spec[between_SH_FO & (ep <= 1580)] = 1   # H (low amplitude or no EP measurement)
    spec[between_SH_FO & (ep > 1580)] = 2    # O (high amplitude)
    return spec

#zhengjiangshan: when consider lu's method, this function may also be modified
def ENA_flytime(df):
    """Select a supported positive, finite flight time in ns.

    Prefer the main TOF when both axes agree with it within 5 ns. Otherwise
    use the smaller axis TOF only when X/Y agree within 5 ns or their ratio
    lies strictly between 0.9 and 1.1. Unsupported estimates remain NaN;
    this function does not drop event rows or repair recorded channels.
    """
    tof_m = np.asarray(df['Tof_M'], dtype=np.float64)
    tof_x = np.asarray(df['Tof_X'], dtype=np.float64)
    tof_y = np.asarray(df['Tof_Y'], dtype=np.float64)
    valid_m = np.isfinite(tof_m) & (tof_m > 0)
    valid_xy = np.isfinite(tof_x) & (tof_x > 0) & np.isfinite(tof_y) & (tof_y > 0)
    with np.errstate(invalid='ignore', divide='ignore'):
        main = valid_m & valid_xy & (np.abs(tof_m - tof_x) < 5) & (np.abs(tof_m - tof_y) < 5)
        ratio = tof_x / tof_y
        consistent_xy = (np.abs(tof_x - tof_y) < 5) | ((ratio > 0.9) & (ratio < 1.1))
    axes = ~main & valid_xy & consistent_xy
    tof = np.full(len(df), np.nan, dtype=np.float32)
    tof[main] = tof_m[main]
    tof[axes] = np.minimum(tof_x[axes], tof_y[axes])
    return tof

#EV
def ENA_energy(df):
    """Kinetic energy in keV; unsupported species or geometry gives NaN."""
    xyz = np.asarray(df[['X1', 'X2', 'Y1', 'Y2']], dtype=np.float64)
    tof = np.asarray(df['Tof'], dtype=np.float64)
    spec = np.asarray(df['Spec'])
    valid = np.isfinite(xyz).all(axis=1) & np.isfinite(tof) & (tof > 0) & np.isin(spec, [1, 2])
    energy = np.full(len(df), np.nan, dtype=np.float32)
    x1, x2, y1, y2 = xyz[valid].T
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + 50**2)
        speed = distance * 1e-3 / (tof[valid]*1e-9) / const.SPEED_C
        mass = np.where(spec[valid] == 2, const.MASS_O, const.MASS_H)
        values = 0.5 * mass * speed**2 / 1000
    good = np.isfinite(values) & (values > 0) & (values <= np.finfo(np.float32).max)
    energy[np.flatnonzero(valid)[good]] = values[good]
    return energy


# ---------------------------------------------------------------------------
# Position conversion
# ---------------------------------------------------------------------------



def l1_2_l1a(df1, engineering):
    """Convert L1 to L1A with global engineering lookup.

    Only rows with all eight finite position channels are retained. HV is
    looked up exactly by PacketID; attitude and position are interpolated
    globally from the PDF table.
    """
    if len(df1) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1a_type.items()})

    # Initialize output DataFrame with L1A schema
    df1a_out = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df1a_type.items()})

    # Copy the L1 event columns that L1A and L2 consume.
    base_cols = ['UUID', 'PacketID', 'EventID', 'PID', 'BM', 'BX',
                 'Timestamp', 'TimeInternal', 'EvtLevel', 'EP',
                 'Tof_0', 'Tof_1', 'Tof_2', 'Tof_3',
                 'TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
                 'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']
    base_cols = [c for c in base_cols if c in df1.columns]

    check_cols = ['TX_1s', 'TX_1e', 'TX_2s', 'TX_2e',
                  'TY_1s', 'TY_1e', 'TY_2s', 'TY_2e']

    # All eight position channels must be present and finite. A NaN must not
    # be treated as a non-zero measurement.
    mask = ((df1[check_cols] != 0).all(axis=1)
            & np.isfinite(df1[check_cols]).all(axis=1))

    df1a_out = df1.loc[mask, base_cols].reset_index(drop=True)
    if len(df1a_out) == 0:
        return df1a_out

    if not isinstance(engineering, EngineeringIndex):
        raise TypeError('engineering must be an EngineeringIndex')
    if 'PacketID' not in df1.columns:
        raise ValueError('L1A conversion requires PacketID')

    packet_ids = df1.loc[mask, 'PacketID'].to_numpy()
    hv = engineering.hv_at(packet_ids)
    if hv is None:
        raise ValueError('EngineeringIndex has no PacketID/HV table')
    df1a_out['HV_C'] = hv['HV_C'].to_numpy(dtype=np.float32)
    df1a_out['HV_G'] = hv['HV_G'].to_numpy(dtype=np.float32)

    # Attitude and position are packet-level values, looked up by PacketID.
    quats, pos = engineering.pose_at(packet_ids)
    if quats is None:
        raise ValueError('EngineeringIndex has no PacketID/pose table')
    df1a_out['MJ2000_q1'] = quats[:, 0].astype(np.float32)
    df1a_out['MJ2000_q2'] = quats[:, 1].astype(np.float32)
    df1a_out['MJ2000_q3'] = quats[:, 2].astype(np.float32)
    df1a_out['MJ2000_Star_X'] = pos[:, 0].astype(np.float32)
    df1a_out['MJ2000_Star_Y'] = pos[:, 1].astype(np.float32)
    df1a_out['MJ2000_Star_Z'] = pos[:, 2].astype(np.float32)

    return df1a_out


# L1A -> L2
def l1a_2_l2(l1a, params, prob='a'):
    """Convert L1A event rows without selection or channel repair.

    Supply the same calibration mapping used to decode the selected probe
    (including ground/orbit selection). Unsupported TOF, species and energy
    remain explicitly invalid rather than generating plausible energies.
    """
    prob = _normalize_probe(prob)
    if params is None:
        raise ValueError('params is required')
    if len(l1a) == 0:
        return pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df2_type.items()})
    tof_m, tof_x, tof_y = f_tof_mxy(l1a)
    df2 = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in st.df2_type.items()})
    for k in ['UUID','PacketID','EventID','PID','BM','BX','Timestamp','TimeInternal','EvtLevel','EP','HV_G','HV_C','MJ2000_q1','MJ2000_q2','MJ2000_q3','MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z']:
        if k in l1a.columns:
            df2[k] = l1a[k]

    tof_ep = f_tof_ep(df2)
    df2['Tof_X']  = tof_x
    df2['Tof_Y']  = tof_y
    df2['Tof_M'] = tof_m

    df2['Tof'] = ENA_flytime(df2) + tof_ep
    df2['Spec'] = ENA_spec(df2)
    
    x1,y1 = params['ns_2_xy'](l1a,0,prob=prob)
    x2,y2 = params['ns_2_xy'](l1a,1,prob=prob)

    df2['X1'] = x1
    df2['Y1'] = y1
    df2['X2'] = x2
    df2['Y2'] = y2

    df2['Energy'] = ENA_energy(df2)
    
    df2['PX'] = (np.arctan2(x1-x2,50)*180/np.pi).astype('float32')
    df2['PY'] = (np.arctan2(y1-y2,50)*180/np.pi).astype('float32')

    return df2
