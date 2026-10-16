#datetime_str = [ f"{int(t[0])}-{int(t[1])}-{int(t[2])} {int(t[3])}:{int(t[4])}:{int(t[5])}" for t in data[:,:6]]
#
#ax[0].plot(df['time'],df['a'],'k-',marker='x')
#ax[1].plot(df['time'],df['e'],'k-',marker='x')
#ax[2].plot(df['time'],df['i'],'k-',marker='x')
#ax[3].plot(df['time'],df['long_of_node'],'k-',marker='x')
#ax[4].plot(df['time'],df['arg_of_perigee'],'k-',marker='x')
#ax[5].plot(df['time'],df['M'],'k-',marker='x')
from astropy import units as u
from poliastro.bodies import Moon, Earth, Sun
from poliastro.twobody import Orbit
from astropy.time import Time
def mean_to_true_anomaly(M, e, tol=1e-12, max_iter=100):
    """
    将平近点角转换为真近点角

    参数:
        M: 平近点角 (弧度)，范围应在 [0, 2π)
        e: 轨道偏心率，e ∈ [0, 1)（椭圆轨道）
        tol: 数值迭代的收敛精度
        max_iter: 最大迭代次数

    返回:
        f: 真近点角 (弧度)，范围在 [0, 2π)
    """
    # 确保平近点角在 [0, 2π) 范围内
    M = M % (2 * np.pi)

    # 步骤1: 求解开普勒方程 M = E - e*sin(E)，得到偏近点角E
    # 选择初始值（对大偏心率更稳定）
    if e < 0.8:
        E = M
    else:
        E = np.pi if M > np.pi else 0.0

    # 牛顿迭代法求解开普勒方程
    for _ in range(max_iter):
        f_kepler = E - e * np.sin(E) - M  # 开普勒方程残差
        f_prime = 1 - e * np.cos(E)       # 残差的导数
        delta = f_kepler / f_prime        # 迭代修正量
        E -= delta

        if abs(delta) < tol:
            break
    else:
        raise RuntimeError(f"开普勒方程未收敛（M={M:.4f}, e={e:.4f}）")

    # 步骤2: 将偏近点角E转换为真近点角f
    # 使用三角函数关系避免反正切可能的象限错误
    cos_f = (np.cos(E) - e) / (1 - e * np.cos(E))
    sin_f = (np.sqrt(1 - e**2) * np.sin(E)) / (1 - e * np.cos(E))

    # 计算真近点角并确保在 [0, 2π) 范围内
    f = np.arctan2(sin_f, cos_f)

    return f

import numpy as np
import pandas as pd 
df = pd.read_csv('QQ_24H.csv')
#7-9
oid=7
a = df['a'][oid] <<u.km
ecc = df['e'][oid] <<u.one
inc = df['i'][oid] <<u.deg
raan = df['Omega'][oid] <<u.deg
argp = df['omega'][oid] <<u.deg
nu = np.degrees(mean_to_true_anomaly(np.radians(df['M'][oid]),df['e'][oid])) << u.deg

# ---- predict orbit -----
#given six roots at a moment
initial = Orbit.from_classical(Moon, a, ecc, inc, raan, argp, nu, epoch=Time(df['time'][oid]))

qq_oct = initial.to_ephem()

r,v = qq_oct.rv()
f,ax = subplots(3,1,sharex=True)
for i in range(3):
    ax[i].plot(pd.to_datetime(qq_oct.epochs.value), r[:,i])
    
    
from astropy.time import Time, TimeDelta
from poliastro.twobody.sampling import EpochsArray

tofs = TimeDelta(np.linspace(0, 2 * u.day, num=1000))
L15_qq = initial.to_ephem(EpochsArray(initial.epoch + tofs))
loc_qq,_ = L15_qq.rv()
#t_qq = L15_qq.epochs.unix
from merge import astros, coord as gcod

start_time = "2025-09-06 00:00:00"
end_time = "2025-09-08 00:00:00"
full_datetime = pd.date_range(start=start_time, end=end_time, freq='min')
full_moon = gcod.gen(full_datetime, astros.moon(full_datetime), 'ECI2000', 0)
full_moon.convert('GSM')

loc_moon = astros.moon(L15_qq.epochs)
global_qq = gcod.gen(pd.to_datetime(L15_qq.epochs.value), loc_moon + loc_qq.value, 'ECI2000', 0)
global_qq.convert('GSM')

f,ax=subplots(3,2,sharex=True)
for i in range(3): ax[i,0].plot(full_datetime, full_moon['ECI2000'][:,i]/6371,label='Moon')
for i in range(3): ax[i,1].plot(full_datetime, full_moon['GSM'][:,i]/6371,label='Moon')

for i in range(3): ax[i,0].plot(pd.to_datetime(L15_qq.epochs.value), global_qq['ECI2000'][:,i]/6371,label='QQ')
for i in range(3): ax[i,1].plot(pd.to_datetime(L15_qq.epochs.value), global_qq['GSM'][:,i]/6371,label='QQ')

for _ in ax.flatten():
    _.legend()

for i in range(3): ax[i,0].set_ylabel([r'$X_{J2000}$',r'$Y_{J2000}$',r'$Z_{J2000}$'][i])
for i in range(3): ax[i,1].set_ylabel([r'$X_{GSM}$',r'$Y_{GSM}$',r'$Z_{GSM}$'][i])



from lunar_python import Solar, Lunar
from datetime import datetime, timedelta

def get_lunar_fifteenth_in_year():
    """获取未来一年内所有农历十五的公历日期"""
    lunar_fifteenths = []
    today = datetime.now()
    
    # 检查未来365天
    for i in range(366):  # 加1天确保覆盖完整一年
        # 计算未来第i天的公历日期
        current_date = today + timedelta(days=i)
        
        # 转换为农历
        solar = Solar.fromYmd(current_date.year, current_date.month, current_date.day)
        lunar = solar.getLunar()
        
        # 检查是否是农历十五
        if lunar.getDay() == 15:
            # 存储公历日期
            lunar_fifteenths.append(current_date.strftime('%Y-%m-%d'))
    return lunar_fifteenths

# 获取结果
fifteenth_dates = get_lunar_fifteenth_in_year()

#orbs=[]
#for i in range(12):
#    orbs.append(initial.propagate(28.*i<<u.day, method=CowellPropagator(f=f)))
pre_orb = initial.propagate(Time(fifteenth_dates[0])-3*u.day, method=CowellPropagator(f=f))


tofs = TimeDelta(np.linspace(0, 6 * u.day, num=1000))
qq_ephem_L15 = pre_orb.to_ephem(
    EpochsArray(pre_orb.epoch + tofs, method=CowellPropagator(rtol=1e-6, f=f)),
)

f,ax=subplots(5,1,sharex=True)
t_orb = [pd.to_datetime(orbs[i].epoch.value) for i in range(30)]
ax[0].plot(t_orb, [orbs[i].a.value for i in range(30)],marker='o')
ax[1].plot(t_orb, [orbs[i].ecc.value for i in range(30)],marker='o')
ax[2].plot(t_orb, [np.degrees(orbs[i].inc.value ) for i in range(30)],marker='o')
ax[3].plot(t_orb, [np.degrees(orbs[i].raan.value) for i in range(30)],marker='o')
ax[4].plot(t_orb, [np.degrees(orbs[i].argp.value) for i in range(30)],marker='o')


t_obs = pd.to_datetime(df['time'])
ax[0].plot(t_obs,df['a'],'kx')
ax[1].plot(t_obs,df['e'],'kx')
ax[2].plot(t_obs,df['i'],'kx')
ax[3].plot(t_obs,df['Omega'],'kx')
ax[4].plot(t_obs,df['omega'],'kx')
#ax[5].plot([pd.to_datetime(orbs[i].epoch.value) for i in range(30)], [orbs[i].a.value for i in range(30)])



