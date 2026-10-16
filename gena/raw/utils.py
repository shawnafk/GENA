import pandas as pd
import numpy as np

from scipy import interpolate
i1d = interpolate.interp1d
from scipy.spatial.transform import Slerp
from scipy.spatial.transform import Rotation as Rot

def init_pqt(df_type,filename):
    df = pd.DataFrame({col: pd.Series(dtype=dtype) for col, dtype in df_type.items()})
    df.to_parquet(filename, engine='fastparquet',compression='gzip')
    return 0


def save_pqt(df,filename):
    #df = dict_2_df(df_dict,df_type)
    df.to_parquet(filename, engine='fastparquet',append=True)
    #pd.DataFrame(dict).astype(df_type).to_parquet(filename, engine='fastparquet')
    return 0

from scipy import interpolate
#!!! this interpolation function need a lot of pdf data, if only few pdf data, interpolation may fail
def func_interp(pdf,name):
    x = pdf['Timestamp']
    y = pdf[name]
    idx = y!=0
    return i1d(x[idx],y[idx],bounds_error=False, fill_value="extrapolate")


def get_para(pdf,name,trange=None):
    if trange is None:
        return func_interp(pdf,name)(pdf['Timestamp'])
    else:
        return func_interp(pdf,name)(trange)

def load_loc(pdf,gap=1):
    t = pdf['Timestamp']
    sx = pdf['MJ2000_Star_X']
    sy = pdf['MJ2000_Star_Y']
    sz = pdf['MJ2000_Star_Z']
    idx = (sx!=0) & (sy!=0) &  (sz!=0)
    #return np.column_stack([sx,sy,sz])/const.RE/1000
    #don't change unit here, keep the origin
    return t[idx][::gap] ,np.column_stack([sx[idx],sy[idx],sz[idx]])[::gap,:]

'''
File ~/syncthing/GENA-CODE/merge/load.py:420, in load_rot(pdf, gap)
    418 w_valid = np.sqrt(1-q1_valid**2-q2_valid**2-q3_valid**2)
    419 quat_ary = np.stack([ q1_valid, q2_valid, q3_valid,w_valid], axis=1)  # (N, 4)
--> 420 return t_valid[::gap],Rot.from_quat(quat_ary[::gap,:])

File _rotation.pyx:999, in scipy.spatial.transform._rotation.Rotation.from_quat()

File _rotation.pyx:857, in scipy.spatial.transform._rotation.Rotation.__init__()

ValueError: Found zero norm quaternions in `quat`.
'''

def load_rot(pdf,gap=1,norm_tol=0.1):
    t = pdf['Timestamp']
    q1 = pdf['MJ2000_q1']
    q2 = pdf['MJ2000_q2']
    q3 = pdf['MJ2000_q3']
    idx = (q1!=0) & (q2!=0) &  (q3!=0)
    #t_valid  = t[idx]
    #q1_valid = q1[idx]
    #q2_valid = q2[idx]
    #q3_valid = q3[idx]
    #w_valid = np.sqrt(1-q1_valid**2-q2_valid**2-q3_valid**2)
    #quat_ary = np.stack([ q1_valid, q2_valid, q3_valid,w_valid], axis=1)  # (N, 4)
    #return t_valid[::gap],Rot.from_quat(quat_ary[::gap,:])
    # 2. 计算平方和，避免根号下为负数（数值裁剪）
    q_sq_sum = q1**2 + q2**2 + q3**2
    # 裁剪平方和到 [0, 1]（防止因噪声导致 q_sq_sum > 1 或 <0）
    q_sq_sum_clipped = np.clip(q_sq_sum, 0.0, 1.0)
    
    # 3. 计算 w 分量（四元数实部）
    w = np.sqrt(1 - q_sq_sum_clipped)
    
    # 4. 堆叠四元数数组（注意：scipy 支持 (x,y,z,w) 或 (w,x,y,z)，需确认数据格式！）
    # 这里假设你的数据是 (x=q1, y=q2, z=q3, w)，若实际是 (w,q1,q2,q3)，需调整顺序为 [w, q1, q2, q3]
    quat_ary = np.stack([q1, q2, q3, w], axis=1)  # shape=(N,4)
    
    # 5. 计算每个四元数的范数，筛选“范数接近1”的有效数据（核心修复！）
    quat_norm = np.linalg.norm(quat_ary, axis=1)  # 计算每个四元数的范数（shape=(N,)）
    # 有效条件：范数在 [1-norm_tol, 1+norm_tol] 内，且非 nan/inf
    valid_mask = (quat_norm >= (1 - norm_tol)) & (quat_norm <= (1 + norm_tol)) & np.isfinite(quat_norm) & idx
    
    # 6. 应用有效掩码，筛选四元数和时间戳
    t_valid = t[valid_mask]
    quat_valid = quat_ary[valid_mask]
    
    # 7. 四元数归一化（确保范数严格=1，彻底解决零范数问题）
    quat_normalized = quat_valid / np.linalg.norm(quat_valid, axis=1, keepdims=True)
    
    t_sampled = t_valid[::gap]
    quat_sampled = quat_normalized[::gap]
    
    rot = Rot.from_quat(quat_sampled)
    
    # ===================== 新增：按时间戳严格排序 =====================
    # 转换为 numpy 数组方便操作
    t_array = np.array(t_sampled)
    quat_array = rot.as_quat()
    
    # 按时间戳排序
    sorted_indices = np.argsort(t_array)
    t_sorted = t_array[sorted_indices]
    quat_sorted = quat_array[sorted_indices]
    
    # 去重：只保留唯一时间戳（取第一个出现的）
    _, unique_indices = np.unique(t_sorted, return_index=True)
    t_unique = t_sorted[unique_indices]
    quat_unique = quat_sorted[unique_indices]
    
    # 重新生成 Rotation 对象
    rot_sorted = Rot.from_quat(quat_unique)
    
    return t_unique.tolist(), rot_sorted
    #return t_sampled.to_list(), rot

#Slerp exterpolation for rotation
#take list/array as parameter
def i_rot_1(t,rot,t_eval):
    #if t_eval is out of t range, extend rot with the edge value
    if t_eval[0] < t[0]:
        t = np.insert(t,0,t_eval[0])
        rot = Rot.concatenate([Rot.from_quat([rot.as_quat()[0]]), rot])
    slerp = Slerp(t, rot)
    new_rot = slerp(t_eval)
    return new_rot


def i_rot(t, rot, t_eval):
    # ===================== 修复：确保时间严格递增 =====================
    t_array = np.array(t)
    quat_array = rot.as_quat()
    
    # 1. 按时间排序
    sorted_idx = np.argsort(t_array)
    t_sorted = t_array[sorted_idx]
    quat_sorted = quat_array[sorted_idx]
    
    # 2. 去重（强制唯一时间戳）
    _, unique_idx = np.unique(t_sorted, return_index=True)
    t_clean = t_sorted[unique_idx]
    quat_clean = quat_sorted[unique_idx]
    
    # 3. 检查是否严格递增（防止极端情况）
    if len(t_clean) < 2:
        # 数据太少，直接返回边界值
        return Rot.concatenate([Rot.from_quat([quat_clean[0]])] * len(t_eval))
    
    # 4. 生成干净的 Rotation 对象
    rot_clean = Rot.from_quat(quat_clean)
    
    # 5. 边界外推（保留你原来的逻辑）
    if t_eval[0] < t_clean[0]:
        t_clean = np.insert(t_clean, 0, t_eval[0])
        rot_clean = Rot.concatenate([Rot.from_quat([quat_clean[0]]), rot_clean])
    if t_eval[-1] > t_clean[-1]:
        t_clean = np.append(t_clean, t_eval[-1])
        rot_clean = Rot.concatenate([rot_clean, Rot.from_quat([quat_clean[-1]])])
    
    # 6. 安全初始化 Slerp
    slerp = Slerp(t_clean, rot_clean)
    new_rot = slerp(t_eval)
    return new_rot

def get_rot(pdf,t_eval):
    ts,rot = load_rot(pdf)
    q4 = i_rot(ts,rot,t_eval)
    return q4

def interp_trajectory_1(pdf, t_eval):
    '''
    Interpolate satellite position (MJ2000_Star_X/Y/Z) in spherical coordinates
    and convert back to Cartesian. This respects orbital trajectory structure
    better than independent linear interpolation of each Cartesian component.
    '''
    t = pdf['Timestamp']
    x = pdf['MJ2000_Star_X'].values.astype(np.float64)
    y = pdf['MJ2000_Star_Y'].values.astype(np.float64)
    z = pdf['MJ2000_Star_Z'].values.astype(np.float64)
    # Filter out zero/invalid positions
    idx = (x != 0) & (y != 0) & (z != 0)
    if idx.sum() < 2:
        # Not enough valid data, fall back to simple linear interp
        return (func_interp(pdf, 'MJ2000_Star_X')(t_eval),
                func_interp(pdf, 'MJ2000_Star_Y')(t_eval),
                func_interp(pdf, 'MJ2000_Star_Z')(t_eval))
    t_v, x_v, y_v, z_v = t[idx], x[idx], y[idx], z[idx]
    # Convert to spherical for smoother interpolation
    r = np.sqrt(x_v**2 + y_v**2 + z_v**2)
    lat = np.arcsin(np.clip(z_v / r, -1, 1))
    lon = np.arctan2(y_v, x_v)
    # Interpolate in spherical
    r_i = i1d(t_v, r, bounds_error=False, fill_value='extrapolate')(t_eval)
    lat_i = i1d(t_v, lat, bounds_error=False, fill_value='extrapolate')(t_eval)
    lon_i = i1d(t_v, lon, bounds_error=False, fill_value='extrapolate')(t_eval)
    # Convert back to Cartesian
    xi = r_i * np.cos(lat_i) * np.cos(lon_i)
    yi = r_i * np.cos(lat_i) * np.sin(lon_i)
    zi = r_i * np.sin(lat_i)
    return xi, yi, zi
    
def interp_trajectory_2(pdf, t_eval):
    """
    Interpolate satellite trajectory positions under gravity.
    
    X, Y, Z coordinates constitute a trajectory under gravity.
    They should be interpolated considering physical constraints rather than
    simple linear interpolation.
    
    This uses a simplified orbital mechanics approach:
    1. Fit a smooth curve through known points
    2. Ensure the interpolation respects gravitational continuity
    
    Args:
        pdf: DataFrame with 'Timestamp', 'MJ2000_Star_X', 'MJ2000_Star_Y', 'MJ2000_Star_Z'
        t_eval: timestamps to evaluate at
        
    Returns:
        np.array of shape (N, 3) with interpolated X, Y, Z positions
    """
    from scipy.interpolate import CubicSpline
    
    t = pdf['Timestamp'].values
    x = pdf['MJ2000_Star_X'].values
    y = pdf['MJ2000_Star_Y'].values
    z = pdf['MJ2000_Star_Z'].values
    
    # Filter valid points
    valid = (x != 0) & (y != 0) & (z != 0) & np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    t_valid = t[valid]
    x_valid = x[valid]
    y_valid = y[valid]
    z_valid = z[valid]
    
    if len(t_valid) < 4:
        # Not enough points for cubic, fall back to linear
        return np.column_stack([
            func_interp(pdf, 'MJ2000_Star_X')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Y')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Z')(t_eval)
        ])
    
    # Use cubic spline for smoother interpolation that respects continuity
    # This better approximates orbital trajectories under gravity
    cs_x = CubicSpline(t_valid, x_valid)
    cs_y = CubicSpline(t_valid, y_valid)
    cs_z = CubicSpline(t_valid, z_valid)
    
    return np.column_stack([cs_x(t_eval), cs_y(t_eval), cs_z(t_eval)])

def interp_trajectory(pdf, t_eval):
    """
    Interpolate satellite trajectory positions under gravity.
    """
    from scipy.interpolate import CubicSpline
    
    t = pdf['Timestamp'].values
    x = pdf['MJ2000_Star_X'].values
    y = pdf['MJ2000_Star_Y'].values
    z = pdf['MJ2000_Star_Z'].values
    
    # Filter valid points
    valid = (x != 0) & (y != 0) & (z != 0) & np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    t_valid = t[valid]
    x_valid = x[valid]
    y_valid = y[valid]
    z_valid = z[valid]
    
    if len(t_valid) < 4:
        # Not enough points for cubic, fall back to linear
        return np.column_stack([
            func_interp(pdf, 'MJ2000_Star_X')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Y')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Z')(t_eval)
        ])
    
    # ===================== 修复：时间排序 + 去重 =====================
    # 1. 把时间和坐标绑在一起
    data = np.column_stack([t_valid, x_valid, y_valid, z_valid])
    
    # 2. 按时间戳排序
    sorted_indices = np.argsort(data[:, 0])
    data_sorted = data[sorted_indices]
    
    # 3. 去重：只保留唯一时间戳
    _, unique_indices = np.unique(data_sorted[:, 0], return_index=True)
    data_clean = data_sorted[unique_indices]
    
    # 4. 拆分回时间和坐标
    t_clean = data_clean[:, 0]
    x_clean = data_clean[:, 1]
    y_clean = data_clean[:, 2]
    z_clean = data_clean[:, 3]
    
    # 5. 再次检查数据量
    if len(t_clean) < 4:
        return np.column_stack([
            func_interp(pdf, 'MJ2000_Star_X')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Y')(t_eval),
            func_interp(pdf, 'MJ2000_Star_Z')(t_eval)
        ])
    
    # Use cubic spline for smoother interpolation
    cs_x = CubicSpline(t_clean, x_clean)
    cs_y = CubicSpline(t_clean, y_clean)
    cs_z = CubicSpline(t_clean, z_clean)
    
    return np.column_stack([cs_x(t_eval), cs_y(t_eval), cs_z(t_eval)])