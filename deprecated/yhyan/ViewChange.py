from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

import astropy.units as units
from sunpy.coordinates import frames

from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

# timestamp = ['2025-04-14T01:00:00', '2025-04-14T04:00:00', '2025-04-14T07:00:00', '2025-04-14T10:00:00',
#              '2025-04-15T01:00:00', '2025-04-15T04:00:00', '2025-04-15T07:00:00', '2025-04-15T10:00:00',
#              '2025-04-16T01:00:00', '2025-04-16T04:00:00', '2025-04-16T07:00:00', '2025-04-16T10:00:00']
timestamp = '2025-04-16T10:00:00'
solar_system_ephemeris.set('input/de432s.bsp')

moon_pos_gei = get_body('moon', Time(timestamp)).transform_to(frames.GeocentricEarthEquatorial).cartesian
earth_pos_gei = get_body('earth', Time(timestamp)).transform_to(frames.GeocentricEarthEquatorial).cartesian
sun_pos_gei = get_body('sun', Time(timestamp)).transform_to(frames.GeocentricEarthEquatorial).cartesian

sun = sun_pos_gei.xyz.value
earth = earth_pos_gei.xyz.value
moon = moon_pos_gei.xyz.value

# 初始视线方向（从月球指向地球）
v0 = earth - moon
v0 /= np.linalg.norm(v0)

# 定义日地连线方向
L = earth - sun
L /= np.linalg.norm(L)

# 计算日-月-地平面的法向量
n = np.cross(earth - moon, sun - moon)
n /= np.linalg.norm(n)

# 生成在平面内且垂直于日地连线的目标方向
vt = np.cross(n, L)
vt /= np.linalg.norm(vt)

# 计算旋转轴和旋转角
u = np.cross(v0, vt)
u /= np.linalg.norm(u)
theta = np.arccos(np.clip(np.dot(v0, vt), -1.0, 1.0))

# 生成四元数
q = np.array([
    np.cos(theta / 2),
    u[0] * np.sin(theta / 2),
    u[1] * np.sin(theta / 2),
    u[2] * np.sin(theta / 2)
])

negZt = vt + moon
negZt_gei = SkyCoord(x=negZt[0], y=negZt[1], z=negZt[2], representation_type='cartesian', frame = frames.GeocentricEarthEquatorial, obstime=Time(timestamp), unit=[units.km, units.km, units.km])
negZt_icrs = negZt_gei.transform_to('icrs').cartesian
moon_pos_icrs = get_body('moon', Time(timestamp)).transform_to('icrs').cartesian

negZ0 = v0 + moon
negZ0_gei = SkyCoord(x=negZ0[0], y=negZ0[1], z=negZ0[2], representation_type='cartesian', frame = frames.GeocentricEarthEquatorial, obstime=Time(timestamp), unit=[units.km, units.km, units.km])
negZ0_icrs = negZ0_gei.transform_to('icrs').cartesian
# moon_pos_icrs = get_body('moon', Time(timestamp)).transform_to('icrs').cartesian

print('Timestamp:', timestamp)  
print("negZ0:", negZ0_icrs-moon_pos_icrs)
print("NegZt:", negZt_icrs-moon_pos_icrs)
# calculate the quaternion
def get_quaternion(v0, vt):
    u = np.cross(v0, vt)
    u /= np.linalg.norm(u)
    theta = np.arccos(np.clip(np.dot(v0, vt), -1.0, 1.0))
    q = np.array([
        np.cos(theta / 2),
        u[0] * np.sin(theta / 2),
        u[1] * np.sin(theta / 2),
        u[2] * np.sin(theta / 2)
    ])
    return q

q = get_quaternion((negZ0_icrs-moon_pos_icrs).xyz.value, (negZt_icrs-moon_pos_icrs).xyz.value)
print("Quaternion:", q)


# plot the result in ICRS frame
scale = 1e5
negZ0_scaled = (negZ0_icrs-moon_pos_icrs).xyz.value
negZt_scaled = (negZt_icrs-moon_pos_icrs).xyz.value
moon_scaled = moon_pos_icrs.xyz.value / scale
sun_scaled = get_body('sun', Time(timestamp)).transform_to('icrs').cartesian.xyz.value / scale
earth_scaled = get_body('earth', Time(timestamp)).transform_to('icrs').cartesian.xyz.value / scale

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# plot the moon, negZ0 and negZt
ax.scatter(*sun_scaled, color='orange', s=100, label='Sun')
ax.scatter(*earth_scaled, color='blue', s=100, label='Earth')
ax.scatter(*moon_scaled, color='gray', s=100, label='Moon')

ax.quiver(*moon_scaled, *negZ0_scaled, color='red', length=1, arrow_length_ratio=0.2, label='Initial view')
ax.quiver(*moon_scaled, *negZt_scaled, color='green', length=1, arrow_length_ratio=0.2, label='Adjusted view')
ax.plot([sun_scaled[0], earth_scaled[0]], [sun_scaled[1], earth_scaled[1]], [sun_scaled[2], earth_scaled[2]],
        color='orange', linestyle='--', label='Sun-Earth line')

# set the plot parameters
ax.set_xlim([moon_scaled[0]-5, moon_scaled[0]+5])
ax.set_ylim([moon_scaled[1]-5, moon_scaled[1]+5])
ax.set_zlim([moon_scaled[2]-5, moon_scaled[2]+5])
# ax.set_xlim([-5, 5])
# ax.set_ylim([-5, 5])
# ax.set_zlim([-5, 5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# add the legend
ax.legend()

plt.title('Satellite View Adjustment in ICRS Frame')
plt.show()

v0 = (negZ0_icrs-moon_pos_icrs).xyz.value
vt = R.from_quat(q, scalar_first=True).apply(v0)



# 缩放因子（单位转换成更易读的数量级）
scale = 1e5

# 缩放后的位置
sun_scaled = sun / scale
earth_scaled = earth / scale
moon_scaled = moon / scale

# 缩放后的方向
v0_scaled = v0
vt_scaled = vt

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 绘制太阳、地球和月球
ax.scatter(*sun_scaled, color='orange', s=100, label='Sun')
ax.scatter(*earth_scaled, color='blue', s=100, label='Earth')
ax.scatter(*moon_scaled, color='gray', s=100, label='Moon')

# 绘制初始视线方向
ax.quiver(*moon_scaled, *v0_scaled, color='red', length=1, arrow_length_ratio=0.2, label='Initial view')

# 绘制调整后的视线方向
ax.quiver(*moon_scaled, *vt_scaled, color='green', length=1, arrow_length_ratio=0.2, label='Adjusted view')

# 绘制日地连线
ax.plot([sun_scaled[0], earth_scaled[0]], [sun_scaled[1], earth_scaled[1]], [sun_scaled[2], earth_scaled[2]],
        color='orange', linestyle='--', label='Sun-Earth line')

# 设置图形参数
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# 添加图例
ax.legend()

plt.title('Satellite View Adjustment in Geocentric Frame')
plt.show()


