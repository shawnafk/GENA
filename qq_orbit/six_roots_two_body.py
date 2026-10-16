import numpy as np
import math
from datetime import datetime, timedelta

# 月球引力参数 (km³/s²)
LUNAR_MU = 4902.8001

def solve_kepler(M, e, tol=1e-12, max_iter=200):
    """求解开普勒方程，得到偏近点角E"""
    E = M if e < 0.8 else np.pi * np.sign(M)
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        f_prime = 1 - e * np.cos(E)
        delta = f / f_prime
        E -= delta
        if abs(delta) < tol:
            return E
    raise RuntimeError(f"开普勒方程未收敛（M={M}, e={e}）")

def orbital_to_lunar_equatorial(elements):
    """将轨道六根数转换为月心赤道J2000坐标系的位置和速度"""
    a, e, i, Omega, omega, M = elements
    
    # 计算偏近点角和真近点角
    E = solve_kepler(M, e)
    cos_f = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_f = (np.sqrt(1 - e**2) * np.sin(E)) / (1 - e * np.cos(E))
    f = np.arctan2(sin_f, cos_f)
    
    # 轨道半径和角动量
    r = a * (1 - e**2) / (1 + e * cos_f)
    h = np.sqrt(LUNAR_MU * a * (1 - e**2))
    
    # 近月点坐标系状态
    r_peri = np.array([r * cos_f, r * sin_f, 0.0])
    v_peri = np.array([
        -LUNAR_MU / h * sin_f,
        LUNAR_MU / h * (e + cos_f),
        0.0
    ])
    
    # 旋转矩阵：近月点坐标系 → 月心赤道J2000坐标系
    def rot_z(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    def rot_x(theta):
        c, s = np.cos(theta), np.sin(theta)
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    
    R = rot_z(-Omega) @ rot_x(-i) @ rot_z(-omega)
    r_lunar = R @ r_peri
    v_lunar = R @ v_peri
    
    return r_lunar, v_lunar

def predict_state_at_time(known_elements, known_time, target_time):
    """
    预测目标时刻的航天器状态
    
    参数:
        known_elements: 已知时刻的轨道六根数 [a, e, i, Omega, omega, M_known]（角度为弧度）
        known_time: 已知六根数对应的时刻（datetime对象）
        target_time: 目标时刻（datetime对象）
    
    返回:
        r_target: 目标时刻的位置矢量 (km)
        v_target: 目标时刻的速度矢量 (km/s)
    """
    a, e, i, Omega, omega, M_known = known_elements
    
    # 1. 计算平均角速度 n (rad/s)
    n = np.sqrt(LUNAR_MU / (a**3))
    
    # 2. 计算时间差（转换为秒）
    delta_t_known = (known_time - datetime(1970, 1, 1)).total_seconds()  # 已知时刻的时间戳（秒）
    delta_t_target = (target_time - datetime(1970, 1, 1)).total_seconds()  # 目标时刻的时间戳（秒）
    
    # 3. 反推近月点时刻 t0（秒级时间戳）
    t0 = delta_t_known - M_known / n  # 由 M_known = n*(delta_t_known - t0) 推导
    
    # 4. 计算目标时刻的平近点角 M_target
    delta_t = delta_t_target - t0  # 目标时刻与近月点的时间差
    M_target = (n * delta_t) % (2 * np.pi)  # 取模确保在 [0, 2π) 范围内
    
    # 5. 用新的 M_target 更新六根数，计算目标状态
    target_elements = [a, e, i, Omega, omega, M_target]
    r_target, v_target = orbital_to_lunar_equatorial(target_elements)
    
    return r_target, v_target



r = []
v = []

dt = np.linspace(0,86400,360)
t0 = pd.to_datetime(df.iloc[7,1])  # 已知时刻的时间
for dt in dt:
    t = t0 + timedelta(seconds=dt)
    r_target, v_target = predict_state_at_time(df.iloc[7,2:], t0,  t)
    r.append(r_target)
    v.append(v_target)