from concurrent.futures import ProcessPoolExecutor, as_completed
#3d view of tail event in paraview

#for qq, you will need xyz axes/ fov cone with four edge arrow of prob A an prob B.
#only input should be time
import os
import numpy as np
import pyvista as pv
from common import RE, SPEED_C
from omni import load
from common import coord as gcd
from astros import ephemeris as emhs 
from scipy.interpolate import interp1d
import pandas as pd
from geopack import geopack, t96

out_dir = '/mnt/c/Users/Jiangshan/ena_paraview_data'
trange = ['20251204 00:00', '20251209 23:59']



# ==========================================
# 0. 加载 OMNI 太阳风参数
# ==========================================
t_window = '10min'

# 【修复】必须加载 SYM_H (或 Dst)，因为 T96 需要它
omni_vars = ['Vx', 'Vy', 'Vz', 'Pressure', 'BY_GSM', 'BZ_GSM', 'SYM_H']
omni_dict = load.load_spedas('omni', trange, download=False, omni_vars=omni_vars)

# 这里我用模拟数据代替你的 load_omni，请在运行时替换为你自己的 DataFrame
date = pd.date_range(trange[0], trange[1], freq=t_window)

def interp_omni(dict,t):
    ts_ori = np.array([_t.timestamp() for _t in dict['x']])
    y_ori = np.array(dict['y'])
    notnan = ~np.isnan(y_ori)
    f = interp1d(ts_ori[notnan], y_ori[notnan], kind='linear', fill_value='extrapolate')
    return np.array([f(_t.timestamp()) for _t in t])

pressure = interp_omni(omni_dict['Pressure'], date)
pressure[pressure <= 0] = 0.001
vx = interp_omni(omni_dict['Vx'], date)
vy = interp_omni(omni_dict['Vy'], date)
vz = interp_omni(omni_dict['Vz'], date)
imfbygsm = interp_omni(omni_dict['BY_GSM'], date)
imfbzgsm = interp_omni(omni_dict['BZ_GSM'], date)
#symh = interp_omni(omni_dict['SYM_H'], date)
dst_kyoto = load.load_dst(trange)  # T96 需要 SYM_H 或 Dst
dst = interp_omni(dst_kyoto, date)


# ==========================================
# 1. 初始化 3D 空间网格 (GSM坐标系)
# ==========================================
Nx = 80; Ny = 60; Nz = 60
xgsm = np.linspace(-70, 10, Nx)
ygsm = np.linspace(-30, 30, Ny)
zgsm = np.linspace(-30, 30, Nz)
X, Y, Z = np.meshgrid(xgsm, ygsm, zgsm, indexing='ij')
X_flat = X.flatten(order='F')
Y_flat = Y.flatten(order='F')
Z_flat = Z.flatten(order='F')

ave_ut = np.average(np.array([_t.timestamp() for _t in date]))

## average T96
geopack.recalc(ave_ut)
ps = geopack.psi 
parmod = np.array([
    np.average(pressure), 
    np.average(dst), 
    np.average(imfbygsm), 
    np.average(imfbzgsm), 0.0])
b0xgsm, b0ygsm, b0zgsm = geopack.dip(X_flat, Y_flat, Z_flat)

dbxgsm = np.zeros_like(X_flat)
dbygsm = np.zeros_like(Y_flat)
dbzgsm = np.zeros_like(Z_flat)

# Python 级循环 (在此进程内独享 CPU 核心)
for k in range(len(X_flat)):
    dbxgsm[k], dbygsm[k], dbzgsm[k] = t96.t96(parmod, ps, X_flat[k], Y_flat[k], Z_flat[k])

bxgsm = b0xgsm + dbxgsm
bygsm = b0ygsm + dbygsm
bzgsm = b0zgsm + dbzgsm
    


# ==========================================
# 1. 定义独立的工作进程函数 (Worker)
# ==========================================
def process_single_frame(args):
    """
    独立计算某一帧的 3D 磁场并保存。
    为了极速的进程间通信，我们只传入 1D 的坐标轴 xgsm, ygsm, zgsm，在内部重新网格化。
    """
    i, ut, pdyn, _dst, byimf, bzimf, xgsm, ygsm, zgsm, out_dir = args
    
    # 1. 进程内独立构建网格，避免传输巨大数组，极大地节省内存和进程通信时间
    X, Y, Z = np.meshgrid(xgsm, ygsm, zgsm, indexing='ij')
    X_flat = X.flatten(order='F')
    Y_flat = Y.flatten(order='F')
    Z_flat = Z.flatten(order='F')

    # 2. 更新当前【独立进程】的 geopack 偶极倾角
    geopack.recalc(ut)
    ps = geopack.psi 
    
    parmod = np.array([pdyn, _dst, byimf, bzimf, 0.0])
    
    # 3. 计算场
    b0xgsm, b0ygsm, b0zgsm = geopack.dip(X_flat, Y_flat, Z_flat)
    
    dbxgsm = np.zeros_like(X_flat)
    dbygsm = np.zeros_like(Y_flat)
    dbzgsm = np.zeros_like(Z_flat)
    
    # Python 级循环 (在此进程内独享 CPU 核心)
    for k in range(len(X_flat)):
        dbxgsm[k], dbygsm[k], dbzgsm[k] = t96.t96(parmod, ps, X_flat[k], Y_flat[k], Z_flat[k])
    
    bxgsm = b0xgsm + dbxgsm
    bygsm = b0ygsm + dbygsm
    bzgsm = b0zgsm + dbzgsm
    
    # 4. 进程内独立构建 PyVista 对象
    # (VTK 对象跨进程传递会导致段错误，必须在 Worker 内部生成和销毁)
    grid = pv.StructuredGrid(X, Y, Z)
    b_vectors = np.column_stack((bxgsm, bygsm, bzgsm))
    grid.point_data["B_field"] = b_vectors
    
    # 5. 独立保存文件
    frame_file = f"t96_frame_{i:04d}.vts"
    grid.save(os.path.join(out_dir, frame_file))
    
    # 打印进度 (包含进程 PID，你可以看到多个核心在狂奔)
    print(f"[PID {os.getpid()}] 帧 {i:04d} 完成: Pdyn={pdyn:.2f}, Dst={_dst:.0f}")
    
    # 返回这一帧的 PVD 字符串和索引 (因为并行完成顺序是乱的，需要 i 来排序)
    pvd_line = f'  <DataSet timestep="{ut}" group="" part="0" file="{frame_file}"/>'
    return i, pvd_line

# ==========================================
# 2. 主程序入口 (多进程强制要求)
# ==========================================
if __name__ == '__main__':
    # 假设你上面的坐标定义、OMNI数据加载等代码都在这里...
    # xgsm, ygsm, zgsm = ...
    # date, pressure, dst, imfbygsm, imfbzgsm = ...
    
    out_dir = '/mnt/c/Users/Jiangshan/ena_paraview_data'
    os.makedirs(out_dir, exist_ok=True) 

    pvd_lines_dict = {} # 用于存储乱序返回的 PVD 字符串

    # 1. 组装所有需要并行计算的任务包裹 (Arguments)
    tasks = []
    for i, current_time in enumerate(date):
        ut = current_time.timestamp()
        pdyn = float(pressure[i])
        _dst = float(dst[i])
        byimf = float(imfbygsm[i])
        bzimf = float(imfbzgsm[i])
        
        # 将任务打包
        task_args = (i, ut, pdyn, _dst, byimf, bzimf, xgsm, ygsm, zgsm, out_dir)
        tasks.append(task_args)

    print(f"🚀 开始并行计算 3D 磁场，共 {len(tasks)} 帧...")
    
    # 2. 启动多进程池 (默认使用机器所有的物理+逻辑核心)
    # 如果你不想把 CPU 跑满导致电脑卡顿，可以指定 max_workers=4 或 8
    with ProcessPoolExecutor() as executor:
        # 提交所有任务
        futures = {executor.submit(process_single_frame, args): args for args in tasks}
        
        # 收集结果
        for future in as_completed(futures):
            try:
                frame_idx, pvd_str = future.result()
                pvd_lines_dict[frame_idx] = pvd_str
            except Exception as exc:
                print(f"❌ 某一帧生成失败: {exc}")

    # 3. 按时间顺序重组 PVD 文件 (极其关键，因为完成顺序是乱的)
    print("正在生成最终的 PVD 时间序列文件...")
    final_pvd_lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '<Collection>'
    ]
    
    # 按照帧号 i (0 到 863) 从小到大依次写入
    for i in range(len(tasks)):
        if i in pvd_lines_dict:
            final_pvd_lines.append(pvd_lines_dict[i])
            
    final_pvd_lines.extend(['</Collection>', '</VTKFile>'])
    
    with open(os.path.join(out_dir, 't96_magnetic_field.pvd'), 'w') as f:
        f.write('\n'.join(final_pvd_lines))

    print(f"✅ 全部计算完毕！完美收工！")
