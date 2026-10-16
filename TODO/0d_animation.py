import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel('Lat')
ax.set_ylabel('Phi')
ax.set_title('ARTEMIS IN FOV')
ax.grid(True)

x = px
y = py
t = pd_time
# 初始化点和时间文本
point, = ax.scatter([], [], 'ro', markersize=10)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

# 初始化函数
def init():
    point.set_data([], [])
    time_text.set_text('')
    return point, time_text

# 更新函数
def update(frame):
    point.set_data(x[frame], y[frame])
    time_text.set_text(f'time {str(t[frame])}')
    return point, time_text

# 创建动画
ani = FuncAnimation(fig, update, frames=len(t), init_func=init,
                    interval=100, blit=True)

# 显示动画
plt.tight_layout()
plt.show()

# 如果需要保存动画，取消下面一行的注释
# ani.save('point_movement.gif', writer='pillow', fps=10)    
