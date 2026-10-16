from scipy.optimize import root_scalar
import numpy as np
# 1. 定义目标函数
def f(x,r_sat,degs):
    earth_sat =  - r_sat
    obj_sat = np.array([x,0,0]) - r_sat
    return obj_sat.dot(earth_sat) / (np.linalg.norm(obj_sat) * np.linalg.norm(earth_sat)) - np.cos(np.radians(degs))


r_sat = np.array([-49,26,-17])
degs = 35
result = root_scalar( f, args=(r_sat, degs), bracket=[0, -100],   method='bisect'  )


#verify
earth_sat =  - r_sat
obj_sat = np.array([result.root,0,0]) - r_sat
np.degrees(np.arccos(obj_sat.dot(earth_sat) / (np.linalg.norm(obj_sat) * np.linalg.norm(earth_sat))))

print("方程的根:", result.root)
print("函数值 f(根):", result.fun)
print("迭代次数:", result.iterations) 