import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# 经纬度范围
lng_bins = np.arange(0-180, 180 + 5, 5)
lat_bins = np.arange(-90, 90 + 5, 5)

# 创建摩尔维特投影图
figure, ax = plt.subplots(figsize=(10, 5), dpi=300, subplot_kw={'projection': ccrs.Mollweide()})

# 绘制白色带状图
ax.plot(ribbon[0], ribbon[1] + 20, transform=ccrs.PlateCarree(), c='white', linestyle='--', linewidth=1)
ax.plot(ribbon[0], ribbon[1], transform=ccrs.PlateCarree(), c='white', linewidth=1)
ax.plot(ribbon[0], ribbon[1] - 20,transform=ccrs.PlateCarree(), c='white', linestyle='--', linewidth=1)

# 绘制二维直方图（将经纬度转换为弧度）
hist = ax.hist2d(mod(ec_lbd/pi*180,360)-180, ec_beta/pi*180,  bins=[lng_bins, lat_bins], cmin=0,  transform=ccrs.PlateCarree())

ax.set_global()
ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
plt.title("Mollweide")
plt.xlabel("Ecliptic Longitude (°)")
plt.ylabel("Ecliptic Latitude (°)")

plt.show()


# 颜色条
plt.colorbar(hist[3], ax=ax, label='Counts')

# 绘制观测点
ax.scatter(view_direction_list[:, 0]-180, view_direction_list[:, 1],
           c='red', s=1, transform=ccrs.PlateCarree())

# 绘制四个角点连线框
# ax.plot(points_x-180, points_y, marker='.', linestyle='-', color='blue', linewidth=1, transform=ccrs.PlateCarree())

# 设置图像范围和标签
