'''
This script will take a time range and
1. automatically determines all observations intverals by checking the time gap in the parameter dataframe, Large grain level, gap = 3600s.
2.  within the Large grain level, it will list the space environment parameters, such as the solar wind (density, velocity, magnetic field), geomagnetic indices (Kp, Dst, AE, AL, AU), and the satellite position and pointing direction (intersection point with GSM XZ plane at Y = 0) variabtions. Also we will give ENA counts rate here, per 10mins base on our experience, but note that since the pointing may change, and also the satellite may turn on and off, and not to mention the parameter settings, the counts rate may not be continuous.
3. then, within Large grain level, it will check again the small grain level, gap = 300s (that would means the instrument has continuous observations without turn on and off), and give the average parameters about the MCP and Anode parameters. 
4. finally, with in each small observation interval, it will generate the quick summary plot of the observed ENA result (strict level) including: 
    - ENA localtions passed through the Carbon film and reached MCP/Anode
    - All range fov/energy counts rate v.s. time
    - Energy distribution
    - Satellite quick 3D look in GSM coordinates
    - FOV distribution plus projected distribution on the GSM XZ/XY/YZ plane cross 0.0.0 (Earth mode), GSM YZ plane at -65 RE (Tail mode)
    all above plots will be given in total Species, H/O and their ratio
'''
import numpy as np
import pandas as pd
from merge import statistic as sts
from merge import load
from merge.astros import ephemeris as epms
from merge import coord as gcd
from merge import constants as const
from merge import project

from matplotlib import pyplot as plt
from merge.visulize import quicklook as ql
from merge.visulize import aux_plot
#contain all date
trange = ['20170607', '20260301']

#grain level, gap = 3600s
L_Level = 3600
S_Level = 300
ENA_DATA_LEVEL = 'l2s_a'
ENA_DATA_LEVEL2 = 'l1a_a'
HIST_TWINDOW = '10min'
POINTING_THRESHOLD = 10
SUMMARY_COLUMNS = [
    'observation_index',
    'parent_l_level_index',
    'sub_observation_index',
    'start_time',
    'end_time',
    'duration_seconds',
    'pointing_x_min',
    'pointing_x_max',
    'pointing_x_mean',
    'pointing_z_min',
    'pointing_z_max',
    'pointing_z_mean',
    'earth_centered',
    'ena_strict_total_counts',
]

COORD = 'GSM'
l2a = load.load_gena_from(ENA_DATA_LEVEL,trange)
pdf = load.load_gena_from('pdf',trange)
num_obsvations, observation_time_interval,observation_pdf = pdf.split(L_Level)
print(f'total obvervation times {num_obsvations}')
import os

def _pointing_projection(pdf_frame, coord=COORD):
    if len(pdf_frame) == 0:
        return pd.DataFrame(columns=['Date', 'PX', 'PZ'])

    ts_list, eci2000 = load.load_loc(pdf_frame,gap=1)
    try:
        ts_rot, st_rot = load.load_rot(pdf_frame,gap=1)
    except ValueError:
        return pd.DataFrame(columns=['Date', 'PX', 'PZ'])
    n = min(len(ts_list), len(ts_rot), len(eci2000))
    if n == 0:
        return pd.DataFrame(columns=['Date', 'PX', 'PZ'])

    ts_list = list(ts_list)[:n]
    eci2000 = eci2000[:n]
    st_rot = st_rot[:n]
    pd_date = pd.to_datetime(ts_list,unit='s').to_list()
    qq_loc =  eci2000/const.RE
    moon_loc = epms.moon(pd_date)/const.RE*1000
    start  = gcd.gen(pd_date,qq_loc+moon_loc,'ECI2000',0)
    start.convert(coord)
    qq_mz_axis = gcd.gen(pd_date, np.tile([0,0,-1], (n, 1)), 'QQ', st_rot)
    qq_mz_axis.convert(coord)
    plane_normal = np.array([0,1,0])
    plane_point = np.array([0,0,0])
    start_coord = np.atleast_2d(start[coord])
    direction_coord = np.atleast_2d(qq_mz_axis[coord])
    denom = np.dot(direction_coord, plane_normal)
    numer = np.dot(plane_point - start_coord, plane_normal)
    t = np.full(denom.shape, np.nan, dtype=np.float64)
    valid = np.abs(denom) > 1e-9
    t[valid] = numer[valid] / denom[valid]
    proj = start_coord + t[:, np.newaxis] * direction_coord
    return pd.DataFrame({
        'Date': pd_date,
        'PX': proj[:,0],
        'PZ': proj[:,2],
    })

#how can we orgnized each "functions" in this for loop, and control whether we active this function or not.
#for example, functions are plot env and ENA counts, plot 2d distributions like fov, loc, layer, plot start and end signal
def show_observations(observation_time_interval):
    root_folder = "ALL"
    if not os.path.exists(root_folder):
        os.makedirs(root_folder)
        print(f"✅ 已创建根文件夹: {root_folder}")

    # 遍历每个时间区间，生成子文件夹
    for idx, (start_ts, end_ts) in enumerate(observation_time_interval, start=1):
        # 计算时间差（秒数）
        duration_seconds = int((end_ts - start_ts).total_seconds())
        
        # 格式化开始时间：YYYYMMDDTHHMMSS
        start_time_str = start_ts.strftime("%Y%m%dT%H%M%S")
        end_time_str = end_ts.strftime("%Y%m%dT%H%M%S")
        
        # 拼接文件夹名称：索引_开始时间_持续秒数
        folder_name = f"{idx}_{start_time_str}_{end_time_str}_{duration_seconds}"
        
        # 完整路径
        folder_path = os.path.join(root_folder, folder_name)
        
        # 创建子文件夹
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"📁 已创建: {folder_path}")
        #process solar wind and geomagnetic parameters, and save to the folder, all in csv list
        
        trange = [start_ts.strftime("%Y-%m-%dT%H:%M:%S"), end_ts.strftime("%Y-%m-%dT%H:%M:%S")]
        
        omni_vars = [
        'BY_GSM',
        'BZ_GSM',
        'flow_speed',
        'proton_density',
        ]
        omni = load.load_spedas('omni',trange,omni_vars=omni_vars,download=False)
        dst = load.load_dst(trange)
        ae = load.load_ae(trange,key='ae')
        al = load.load_ae(trange,key='al')
        au = load.load_ae(trange,key='au')
        
        #return empty if no data
        tot_ena = load.load_gena_from(ENA_DATA_LEVEL,trange)
        tot_ena_l1a = load.load_gena_from(ENA_DATA_LEVEL2,trange)
        if len(tot_ena !=0):
            t0, h0 = sts.window_hist0d(tot_ena,window=HIST_TWINDOW)
        else:
            t0=[]
            h0=[]
        if len(tot_ena_l1a !=0):
            t1, h1 = sts.window_hist0d(tot_ena_l1a,window=HIST_TWINDOW)
        else:
            t1=[]
            h1=[]
            
        pdf = load.load_gena_from('pdf',trange)
        ts_list,eci2000 = load.load_loc(pdf,gap=1)
        ts_list,st_rot = load.load_rot(pdf,gap=1)
        pd_date = pd.to_datetime(ts_list,unit='s').to_list()
        qq_loc =  eci2000/const.RE
        moon_loc = epms.moon(pd_date)/const.RE*1000
        start  = gcd.gen(pd_date,qq_loc+moon_loc,'ECI2000',0)
        start.convert(COORD)
        qq_mz_axis = gcd.gen(pd_date, np.tile([0,0,-1], (len(ts_list), 1)), 'QQ', st_rot)
        qq_mz_axis.convert(COORD) 
        proj = project.project(start[COORD],qq_mz_axis[COORD],np.array([0,1,0]),np.array([0,0,0]))
        #plot figure 1
        #sw velocity
        #density
        #IMF BZ
        #IMF BY
        #Dst
        #AE, AL, AU
        #ENA counts rate
        f,axs = plt.subplots(13,1,figsize=(12,30),sharex=True)
        aux_plot.qk_plot(omni['flow_speed'],f=f,ax=axs[0],color='k')
        aux_plot.qk_plot(omni['proton_density'],f=f,ax=axs[1],color='k')
        aux_plot.qk_plot(omni['BZ_GSM'],f=f,ax=axs[2],color='k')
        aux_plot.qk_plot(omni['BY_GSM'],f=f,ax=axs[3],color='k')
        aux_plot.qk_plot(dst,f=f,ax=axs[4],color='k')
        aux_plot.qk_plot(ae,f=f,ax=axs[5],color='r',label='AE')
        aux_plot.qk_plot(al,f=f,ax=axs[5],color='g',label='AL')
        aux_plot.qk_plot(au,f=f,ax=axs[5],color='b',label='AU')
        axs[5].legend(loc='upper right')
        
        axs[6].plot(t1,h1)
        axs[7].plot(t0,h0)
        
        axs[8].plot(pd_date,start[COORD][:,0])
        axs[9].plot(pd_date,start[COORD][:,1])
        axs[10].plot(pd_date,start[COORD][:,2])
        
        axs[11].plot(pd_date,proj[:,0])
        axs[12].plot(pd_date,proj[:,2])
        
        for ax,yname in zip(axs,[r'$v_{sw} (km/s)$',r'$n_{sw} (cm^{-3})$',r'$B^{GSM}_Z (nT)$',r'$B^{GSM}_Y (nT)$','Dst','AE/AL/AU',f'VLD8 /{HIST_TWINDOW}',f'Strict /{HIST_TWINDOW}',r'$X^{QQ} (RE)$',r'$Y^{QQ} (RE)$',r'$Z^{QQ} (RE)$',r'$PX^{qq} (RE)$',r'$PZ^{qq} (RE)$']):
            ax.set_ylabel(yname,rotation=0,labelpad=30)
            ax.grid('on')
        
        
        axs[-1].set_xlim(start_ts-pd.Timedelta('10min'),end_ts+pd.Timedelta('10min'))
        
        f.autofmt_xdate()
        f.savefig(os.path.join(folder_path,'environment_parameters.pdf'),dpi=300,bbox_inches='tight')
        
        plt.close(f)

        if(len(pdf)==0):
            continue
        _,_,subpdf = pdf.split(S_Level)
        
        print(trange)       
        print(len(subpdf))       
        for i,_subpdf in enumerate(subpdf):
            if(len(_subpdf)==0 or len(tot_ena) == 0):
                continue
            tstart = _subpdf.iloc[0]['Date']
            tend = _subpdf.iloc[len(_subpdf)-1]['Date']
            duration =  int((tend - tstart).total_seconds())
            txttitle = f'{tstart.strftime("%Y-%m-%dT%H:%M:%S")}_{tend.strftime("%Y-%m-%dT%H:%M:%S")}_{duration}'
            
            _sub = tot_ena.between([tstart.strftime("%Y-%m-%dT%H:%M:%S"),tend.strftime("%Y-%m-%dT%H:%M:%S")])[0]
            _sub_l1a = tot_ena_l1a.between([tstart.strftime("%Y-%m-%dT%H:%M:%S"),tend.strftime("%Y-%m-%dT%H:%M:%S")])[0]
            f,ax = ql.show_E_lin(_sub)
            ax.set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_E_lin.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)

            f,ax = ql.show_fov(_sub)
            ax.set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_FOV.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)

            f,ax = plt.subplots(2,1,sharex=True,sharey=True)
            ql.show_layer(_sub,1,f=f,ax=ax[0])
            ql.show_layer(_sub,2,f=f,ax=ax[1])
            ax[0].set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_LOC.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)

            f,ax = ql.show_start(_sub_l1a)
            ax.set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_start.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)

            f,ax = ql.show_end(_sub_l1a)
            ax.set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_end.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)
            
            f,ax = ql.show_joint(_sub_l1a)
            ax.set_title(txttitle)
            f.savefig(os.path.join(folder_path,str(i)+'_joint.jpg'),dpi=300,bbox_inches='tight')
            plt.close(f)

            #energy spectrum


def generate_observation_report(observation_time_interval):
    """
    生成所有小观测汇总的 Excel 报告 (.xlsx)，每行一个小观测。
    列: NL, ns, num, start_time, end_time, duration_seconds,
        ena_total_counts, ena_avg_counts_per_sec, PX_mean, PZ_mean, is_earth_mode
    同一大观测 (NL) 的行使用相同背景色标记。过滤 ena_total_counts == 0 的行。
    """
    from openpyxl import Workbook
    from openpyxl.styles import PatternFill

    color_list = plt.cm.tab20(np.linspace(0, 1, 20))[:, :3]

    rows = []
    nl_color_map = {}

    for NL, (start_ts, end_ts) in enumerate(observation_time_interval, start=1):
        trange = [start_ts.strftime("%Y-%m-%dT%H:%M:%S"),
                  end_ts.strftime("%Y-%m-%dT%H:%M:%S")]

        tot_ena = load.load_gena_from(ENA_DATA_LEVEL, trange)
        pdf_large = load.load_gena_from('pdf', trange)

        if len(pdf_large) == 0 or len(tot_ena) == 0:
            continue

        _, _, subpdfs = pdf_large.split(S_Level)

        valid_rows_nl = []

        for ns_raw, _subpdf in enumerate(subpdfs, start=1):
            if len(_subpdf) == 0:
                continue

            tstart = _subpdf.iloc[0]['Date']
            tend = _subpdf.iloc[-1]['Date']
            duration = int((tend - tstart).total_seconds())

            _sub_ena = tot_ena.between(
                [tstart.strftime("%Y-%m-%dT%H:%M:%S"),
                 tend.strftime("%Y-%m-%dT%H:%M:%S")])[0]
            ena_total = len(_sub_ena)
            if ena_total == 0:
                continue
            ena_avg = ena_total / duration

            proj_df = _pointing_projection(_subpdf)
            px_vals = proj_df['PX'].dropna()
            pz_vals = proj_df['PZ'].dropna()
            if len(px_vals) > 0:
                px_mean = float(px_vals.mean())
                pz_mean = float(pz_vals.mean())
                dist = np.sqrt(px_mean ** 2 + pz_mean ** 2)
                is_earth = bool(dist <= POINTING_THRESHOLD)
            else:
                px_mean = None
                pz_mean = None
                is_earth = None

            valid_rows_nl.append({
                'NL': NL,
                'start_time': tstart.strftime("%Y-%m-%dT%H:%M:%S"),
                'end_time': tend.strftime("%Y-%m-%dT%H:%M:%S"),
                'duration_seconds': duration,
                'ena_total_counts': ena_total,
                'ena_avg_counts_per_sec': round(ena_avg, 6),
                'PX_mean': round(px_mean, 3) if px_mean is not None else None,
                'PZ_mean': round(pz_mean, 3) if pz_mean is not None else None,
                'is_earth_mode': is_earth,
            })

        num = len(valid_rows_nl)
        if num == 0:
            continue

        color_idx = (NL - 1) % len(color_list)
        r, g, b = color_list[color_idx]
        hex_color = f"{int(r * 255):02X}{int(g * 255):02X}{int(b * 255):02X}"
        nl_color_map[NL] = hex_color

        for ns, row_data in enumerate(valid_rows_nl, start=1):
            row_data['ns'] = ns
            row_data['num'] = num
            rows.append(row_data)

        print(f"NL={NL}: {num} small observations, color=#{hex_color}")

    if not rows:
        print("⚠ 没有符合条件的观测数据")
        return None

    wb = Workbook()
    ws = wb.active
    ws.title = "Observation Report"

    columns = ['NL', 'ns', 'num', 'start_time', 'end_time',
               'duration_seconds', 'ena_total_counts', 'ena_avg_counts_per_sec',
               'PX_mean', 'PZ_mean', 'is_earth_mode']
    ws.append(columns)

    for row_data in rows:
        ws.append([row_data[c] for c in columns])

    for row_idx in range(2, ws.max_row + 1):
        nl_val = ws.cell(row=row_idx, column=1).value
        hex_color = nl_color_map.get(nl_val)
        if hex_color is None:
            continue
        fill = PatternFill(start_color=hex_color, end_color=hex_color, fill_type="solid")
        for col in range(1, len(columns) + 1):
            ws.cell(row=row_idx, column=col).fill = fill

    wb.save('observation_report.xlsx')
    print(f"报告已保存: observation_report.xlsx, 共 {len(rows)} 行")
    return pd.DataFrame(rows)


report_df = generate_observation_report(observation_time_interval)
print(report_df.head(20))
