#basic libraries
import numpy as np

#astropy lib for astro locations
from astropy import time as AT
from astropy.coordinates import  get_body
from scipy.spatial.transform import Rotation as Rot

#calculate current sheet location
from geopack import geopack, t96

RE=6371
AU = 149597870700/1e3

#self defined functions
#va, vb, 3*1 numpy array
#return angle in radian
def angle_bwteen_vec(va,vb):
    return np.arccos(np.dot(va,vb)/(np.linalg.norm(va)*np.linalg.norm(vb)))

#rot x,y,z, use euler angle, see new gio/coord lib
#always from earth looking mode, with long edge (mechanic X axis) parallel to Equator plane
#return new pointing, to Y=0 plane
def gen_rot(u0, angle):
    q1 = u0[0]*np.sin(angle/2)
    q2 = u0[1]*np.sin(angle/2)
    q3 = u0[2]*np.sin(angle/2)
    w = np.cos(angle/2)
    quat_ary = np.array([q1, q2, q3, w])
    return Rot.from_quat(quat_ary)

def rot_inertia(angle,minus_z_vec,axis):
    Rot = gen_rot(axis,angle)
    new_minus_z  =  Rot.as_matrix()@minus_z_vec
    return new_minus_z

def find_intersection_with_Y0(direction_vec, start_point):
    """
    计算由起点和方向矢量定义的直线与Y=0平面的交点
    
    参数:
        direction_vec: 3x1 numpy数组，表示直线的方向矢量
        start_point: 3x1 numpy数组，表示直线的起点（默认原点）
        
    返回:
        intersection: 3x1 numpy数组，表示交点坐标（x, 0, z）
        若直线与Y=0平面平行或在平面上，则返回None
    """
    # 确保输入是3x1的数组
    direction_vec = np.array(direction_vec).reshape(3, 1)
    start_point = np.array(start_point).reshape(3, 1)
    
    # 提取方向矢量的Y分量和起点的Y分量
    dir_y = direction_vec[1, 0]
    start_y = start_point[1, 0]
    
    # 检查直线是否与Y=0平面平行（方向矢量Y分量为0）
    if np.isclose(dir_y, 0):
        if np.isclose(start_y, 0):
            print("直线在Y=0平面上，有无穷多个交点")
        else:
            print("直线与Y=0平面平行，无交点")
        return None
    
    # 计算参数t：直线方程为 start_point + t*direction_vec，求Y=0时的t
    # 即 start_y + t*dir_y = 0 → t = -start_y / dir_y
    t = -start_y / dir_y
    
    # 计算交点（此时Y分量必为0）
    intersection = start_point + t * direction_vec
    return intersection[0,0],intersection[2,0]


from spacepy import coordinates as coord
from spacepy.time import Ticktock
def convert(isot,value,in_type='ECI2000',out_type='GSM'):
    temp = coord.Coords(value,in_type,'car',ticks=Ticktock(isot)).convert(out_type,'car')
    return np.array([temp.x,temp.y,temp.z]).flatten()

#Calculation Flow
#ICRS coord
#'2025-09-09 00:00:00'
#tailward is +, headward is -
#northward is +, southward is -

def print_info(isot_str,angle_tail,angle_updown,loc_sat=None):
    
    isot = AT.Time(isot_str)
    
    loc_sat=None
    #use moon location to approximate queqiao location
    if loc_sat == None:
        moon = get_body('moon',isot).represent_as('cartesian')
        loc_sat = moon.xyz.value*AU/RE
    loc_sat_gsm = convert(isot,loc_sat)
    
    print('----------------------------------------------')
    print('UTC Time is', isot_str)
    #decide the sign of rotation angle based on lunar date
    #to do
    
    
    #decide the sign of rotation angle based on GSM_Y
    moon_half = 0
    if loc_sat_gsm[1]<0:
        moon_half = 1
        print('Second half of the month, A in the North')
    else:
        print('First half of the month, A in the South')
    if moon_half == 0:
        new_angle_tail = -angle_tail
        new_angle_updown = -angle_updown
    else:
        new_angle_tail = angle_tail
        new_angle_updown = angle_updown

    
    print('QueQiao location is (GSM)', loc_sat_gsm)

    loc_sun = get_body('sun',isot).represent_as('cartesian').xyz.value*AU/RE
    sun_sat_vec = loc_sun - loc_sat

    #By default, looking at Earth
    loc_point=np.array([0,0,0])
    minus_z_vec  = loc_point - loc_sat

    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    print("Initial minus z pointing at GSM Y0 plane at", X0, Z0, "RE")
    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Initial Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')

    #---------- rotation along ICRS Z ---
    minus_z_vec = rot_inertia(np.radians(new_angle_tail),minus_z_vec,np.array([0,0,1]))
    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)

    print("After rotate tailward ", str(new_angle_tail), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")
    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')

    #---------- rotation along sun-earth line ---
    sun_earth = loc_sun - np.array([0,0,0])
    norm_sun_earth = sun_earth/np.linalg.norm(sun_earth)
    minus_z_vec = rot_inertia(np.radians(new_angle_updown),minus_z_vec,norm_sun_earth)

    minus_z_vec_gsm = convert(isot,minus_z_vec)
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    print("After rotate northward ",str(new_angle_updown), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")

    sunavoid = angle_bwteen_vec(sun_sat_vec,minus_z_vec)
    print("Sun avoid angle is", np.degrees(sunavoid), "degs")
    print('----------------------------------------------')


from scipy.spatial.transform import Rotation as Rot
def get_rotation_quaternion_B2A(z_B_A, x_sign=1):
    z_B_A = np.array(z_B_A, dtype=np.float64)
    norm_z = np.linalg.norm(z_B_A)
    if norm_z < 1e-9:
        raise ValueError("B 的 Z 轴矢量不能为零向量！")
    u_z = z_B_A / norm_z  # (u_z1, u_z2, u_z3)

    # 步骤 2：确定 B 的 X 轴 u_x（满足与 A 的 Z 轴垂直 + 与 B 的 Z 轴垂直）
    # 设 u_x = (Xx, Xy, 0)，由 u_x · u_z = 0 得 Xy = - (u_z1 / u_z2) * Xx（假设 u_z2 ≠ 0）
    if abs(u_z[1]) < 1e-9:  # 若 u_z2=0，则 u_z1 ≠ 0（否则 u_z 为 (0,0,±1)）
        Xx = 0.0
        Xy = x_sign * 1.0  # 此时 u_z1 ≠ 0，Xy 取 ±1
    else:
        Xx = x_sign * 1.0  # Xx 符号由 x_sign 确定
        Xy = - (u_z[0] / u_z[1]) * Xx
    x_B_A = np.array([Xx, Xy, 0.0])
    # 归一化 X 轴
    norm_x = np.linalg.norm(x_B_A)
    u_x = x_B_A / norm_x  # (u_x1, u_x2, u_x3=0)

    # 步骤 3：用右手定则计算 B 的 Y 轴 u_y = u_z × u_x
    u_y = np.cross(u_z, u_x)  # 叉乘，天然正交且单位

    # 步骤 4：构建旋转矩阵 R_BA（列向量：u_x, u_y, u_z）
    R_BA = np.column_stack([u_x, u_y, u_z])

    r11, r12, r13 = R_BA[0, 0], R_BA[0, 1], R_BA[0, 2]
    r21, r22, r23 = R_BA[1, 0], R_BA[1, 1], R_BA[1, 2]
    r31, r32, r33 = R_BA[2, 0], R_BA[2, 1], R_BA[2, 2]

    w = 0.5 * np.sqrt(1 + r11 + r22 + r33)
    if abs(w) < 1e-9:  
        x = 0.5 * np.sqrt(1 + r11 - r22 - r33)
        y = r12 / (4 * x)
        z = r13 / (4 * x)
        w = (r23 - r32) / (4 * x)
    else:
        x = (r23 - r32) / (4 * w)
        y = (r31 - r13) / (4 * w)
        z = (r12 - r21) / (4 * w)

    q_BA = np.array([x, y, z, w])
    q_BA = q_BA / np.linalg.norm(q_BA)

    return Rot.from_quat(q_BA)


def print_gse(isot_str,angle_tail,angle_updown,loc_sat=None):
    
    isot = AT.Time(isot_str)
    
    loc_sat=None
    #use moon location to approximate queqiao location
    if loc_sat == None:
        moon = get_body('moon',isot).represent_as('cartesian')
        loc_sat = moon.xyz.value*AU/RE
    loc_sat_gse = convert(isot,loc_sat,out_type='GSE')
    #By default, looking at Earth
    loc_point=np.array([0,0,0])
    minus_z_vec_gse  = loc_point - loc_sat_gse

    
    print('----------------------------------------------')
    print('UTC Time is', isot_str)
    print('QueQiao location (J2000) is', loc_sat)
    #decide the sign of rotation angle based on lunar date
    #to do
    
    #decide the sign of rotation angle based on GSM_Y
    moon_half = 0
    if loc_sat_gse[1]<0:
        moon_half = 1
        print('Second half of the month, A in the North')
    else:
        moon_half=0
        print('First half of the month, A in the South')
    if moon_half == 0:
        new_angle_tail = -angle_tail
        new_angle_updown = -angle_updown
    else:
        new_angle_tail = angle_tail
        new_angle_updown = angle_updown

    #---------- rotation along GSE Z ---
    minus_z_vec_gse = rot_inertia(np.radians(new_angle_tail),minus_z_vec_gse,np.array([0,0,1]))
    #---------- rotation along GSE -X line ---
    minus_z_vec_gse = rot_inertia(np.radians(new_angle_updown),minus_z_vec_gse,np.array([-1,0,0]))

    #q4 num
    q4_mec2gse = get_rotation_quaternion_B2A(minus_z_vec_gse, x_sign=1)

    #1. sunavoid
    loc_sun = get_body('sun',isot).represent_as('cartesian').xyz.value*AU/RE
    sun_sat_vec = loc_sun - loc_sat
    sun_sat_vec_gse = convert(isot,sun_sat_vec,out_type='GSE')
    sunavoid = angle_bwteen_vec(sun_sat_vec_gse,minus_z_vec_gse)
    print("Sun avoid angle is", np.degrees(sunavoid), "degs")
    

    #2. -z intersection
    minus_z_vec_gsm = convert(isot,minus_z_vec_gse,in_type='GSE',out_type='GSM')
    
    loc_sat_gsm = convert(isot,loc_sat,out_type='GSM')
    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    print("After rotate",str(new_angle_tail),'and',str(new_angle_updown), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")

    point_A = gen_rot(np.array([1,0,0]), np.radians(5)).as_matrix()@np.array([0,0,-1])
    point_A_gse = q4_mec2gse.as_matrix()@point_A
    point_A_gsm = convert(isot,point_A_gse,in_type='GSE',out_type='GSM')
    X1,Z1 = find_intersection_with_Y0(point_A_gsm,loc_sat_gsm)
    print("After rotate",str(new_angle_tail),'and',str(new_angle_updown), "degs, A pointing GSM Y0 plane at", X1, Z1, "RE")

    #print('----------------------------------------------')
    return 0


def print_gsm(isot_str,angle_tail,angle_updown,loc_sat=None):
    
    isot = AT.Time(isot_str)
    
    loc_sat=None
    #use moon location to approximate queqiao location
    if loc_sat == None:
        moon = get_body('moon',isot).represent_as('cartesian')
        loc_sat = moon.xyz.value*AU/RE
    loc_sat_gsm = convert(isot,loc_sat,out_type='GSM')
    
    #print('----------------------------------------------')
    #print('UTC Time is', isot_str)
    #decide the sign of rotation angle based on lunar date
    #to do
    
    
    #decide the sign of rotation angle based on GSM_Y
    moon_half = 0
    if loc_sat_gsm[1]<0:
        moon_half = 1
        #print('Second half of the month, A in the North')
    else:
        moon_half=0
        #print('First half of the month, A in the South')
    if moon_half == 0:
        new_angle_tail = -angle_tail
        new_angle_updown = -angle_updown
    else:
        new_angle_tail = angle_tail
        new_angle_updown = angle_updown

    
    #print('QueQiao location is (GSE)', loc_sat_gse)

    #By default, looking at Earth
    loc_point=np.array([0,0,0])
    minus_z_vec_gsm  = loc_point - loc_sat_gsm
    #print("Initial minus z pointing at GSM Y0 plane at", X0, Z0, "RE")
    #print('----------------------------------------------')

    #---------- rotation along GSM Z ---
    minus_z_vec_gsm = rot_inertia(np.radians(new_angle_tail),minus_z_vec_gsm,np.array([0,0,1]))
    #---------- rotation along GSM X line ---
    minus_z_vec_gsm = rot_inertia(np.radians(new_angle_updown),minus_z_vec_gsm,np.array([1,0,0]))

    X0,Z0 = find_intersection_with_Y0(minus_z_vec_gsm,loc_sat_gsm)
    #print("After rotate",str(new_angle_tail),'and',str(new_angle_updown), "degs, minus z pointing GSM Y0 plane at", X0, Z0, "RE")

    #print('----------------------------------------------')
    return X0,Z0

X=[]
Z=[]
i=0
for _t,_sat in zip(pd.to_datetime(t,unit='s'),coord_qq['ECI2000']):
    print(i)
    i+=1
    _x,_z = print_gse(_t,30,2.5,_sat)
    X.append(_x)
    Z.append(_z)
