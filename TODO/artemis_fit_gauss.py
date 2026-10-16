import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def mid(x):
    return (x[1:]+x[:-1])/2
# 定义二维高斯函数
def gaussian_2d(xy, A, x0, y0, sigma_x, sigma_y):
    x, y = xy
    return A * np.exp(-((x-x0)**2/(2*sigma_x**2) + (y-y0)** 2/(2*sigma_y**2))) / (np.pi * sigma_x * sigma_y)

# 数据准备：x, y为坐标数组，counts为计数数组
def gauss_fit(subdf):
    h,x,y = np.histogram2d(subdf['PX'],subdf['PY'],[np.arange(-30,30,2.5),np.arange(-15,15,2.5)])
    x= mid(x)
    y= mid(y)
    Y,X = np.meshgrid(y, x)
    counts = h

    # 初始化参数
    A_init = np.max(counts)
    x0_init = np.sum(X * counts) / np.sum(counts)
    y0_init = np.sum(Y * counts) / np.sum(counts)
    sigma_x_init = 1.0  # 手动设定初始值
    sigma_y_init = 1.0  # 手动设定初始值

    # 拟合参数
    params, covariance = curve_fit(
        gaussian_2d, 
        (X.ravel(), Y.ravel()), 
        counts.ravel(),
        p0=[A_init, x0_init, y0_init, sigma_x_init, sigma_y_init]
    )
    A, x0, y0, sigma_x, sigma_y = params

    # 计算FWHM
    fwhm_x = 2.3548 * sigma_x
    fwhm_y = 2.3548 * sigma_y
    return A, x0, y0, fwhm_x, fwhm_y



#method2 

from pysr import PySRRegressor

model = PySRRegressor(
    maxsize=20,
    niterations=40,  # < Increase me for better results
    binary_operators=["+", "*", "-", "/"],
    unary_operators=[
        "exp",
        "inv(x) = 1/x",
        # ^ Custom operator (julia syntax)
    ],
    extra_sympy_mappings={"inv": lambda x: 1 / x},
    # ^ Define operator for SymPy as well
    elementwise_loss="loss(prediction, target) = (prediction - target)^2",
    # ^ Custom loss function (julia syntax)
)


import numpy as np
from pysr import PySRRegressor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# 重塑数据为 PySR 接受的格式：(n_samples, 2) 输入和 (n_samples,) 输出
X_flat = X.flatten()
Y_flat = Y.flatten()
Z_flat = h.flatten()
input_data = np.column_stack((X_flat, Y_flat))

# 初始化 PySR 模型
model = PySRRegressor(
    niterations=100,  # 增加迭代次数以提高发现复杂函数的概率
    binary_operators=["+", "*", "/", "-", "^"],  # 二元运算符
    unary_operators=["exp", "log", "sqrt", "sin", "cos"],  # 一元函数
    maxsize=20,  # 限制表达式最大复杂度
    populations=5,  # 并行搜索的种群数
    procs=4,  # 使用的 CPU 核心数
    verbosity=0,  # 控制输出详细程度
    elementwise_loss="loss(prediction, target) = (prediction - target)^2",  # 损失函数：均方误差
)

# 训练模型
model.fit(input_data, Z_flat)

# 输出最佳发现的表达式
print("最佳表达式：", model.get_best()["equation"])

# 评估拟合效果
Z_pred = model.predict(input_data).reshape(X.shape)

# 可视化（3D 图）
fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
ax1.set_title('真实函数')

ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, Y, Z_pred, cmap='plasma', alpha=0.8)
ax2.set_title('PySR 拟合结果')

plt.tight_layout()
plt.show()